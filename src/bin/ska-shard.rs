//! `ska-shard` — split a ska2 `.skf` into minimizer bins, or concatenate bins.
//!
//! Workflow this enables (bounds `ska merge` peak memory to ~1/n):
//!   1. `ska build` each sample → per-sample `.skf`
//!   2. `ska-shard split sample.skf -n N`  (run per sample)
//!   3. `ska merge` bin i across all samples (parallel, ~1/n of the key space)
//!   4. `ska-shard concat merged_bin*.skf -o merged.skf`

use std::collections::{BTreeSet, HashMap};
use std::path::{Path, PathBuf};

use anyhow::{bail, ensure, Context, Result};
use clap::{Parser, Subcommand};
use ndarray::{concatenate, Array2, Axis};

use ska_minimizer_split::dispatch_skf_width;
use ska_minimizer_split::minimizer::{decode_flank, minimizer_bin};
use ska_minimizer_split::skf::{SkaArray, SkfInt};

#[derive(Parser)]
#[command(
    name = "ska-shard",
    about = "Split a ska2 .skf into minimizer bins, or concatenate bins back together",
    version
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Split one .skf into N bins by split-kmer minimizer value.
    Split(SplitArgs),
    /// Concatenate per-bin .skf files (sharing identical samples) into one .skf.
    Concat(ConcatArgs),
}

#[derive(Parser)]
struct SplitArgs {
    /// Input .skf file.
    input: PathBuf,
    /// Number of bins to split into.
    #[arg(short = 'n', long = "bins")]
    bins: usize,
    /// Minimizer (l-mer) length; must be <= k-1.
    #[arg(short = 'l', long = "minimizer-len", default_value_t = 9)]
    minimizer_len: usize,
    /// Output prefix; bins are written as <prefix>.<i>.skf. Defaults to the
    /// input file stem.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: Option<String>,
}

#[derive(Parser)]
struct ConcatArgs {
    /// Per-bin .skf inputs to concatenate (>= 1).
    #[arg(required = true)]
    inputs: Vec<PathBuf>,
    /// Output merged .skf.
    #[arg(short = 'o', long = "output")]
    output: PathBuf,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Split(args) => {
            dispatch_skf_width!(&args.input, IntT => run_split::<IntT>(&args))
        }
        Command::Concat(args) => {
            dispatch_skf_width!(&args.inputs[0], IntT => run_concat::<IntT>(&args))
        }
    }
}

fn run_split<IntT: SkfInt>(args: &SplitArgs) -> Result<()> {
    ensure!(args.bins >= 1, "--bins must be >= 1");
    let arr: SkaArray<IntT> = SkaArray::load(&args.input)?;
    let k = arr.k;
    let flank_len = k - 1;
    ensure!(
        args.minimizer_len <= flank_len,
        "--minimizer-len {} exceeds k-1 = {} for this skf (k={})",
        args.minimizer_len,
        flank_len,
        k
    );

    let ncols = arr.names.len();
    let n = args.bins;

    // Per-bin row buffers.
    let mut keys: Vec<Vec<IntT>> = vec![Vec::new(); n];
    let mut counts: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut variants: Vec<Array2<u8>> = (0..n).map(|_| Array2::zeros((0, ncols))).collect();

    for i in 0..arr.n_kmers() {
        let flank = decode_flank(arr.split_kmers[i].into(), k);
        let b = minimizer_bin(&flank, args.minimizer_len, n)?;
        keys[b].push(arr.split_kmers[i]);
        counts[b].push(arr.variant_count[i]);
        variants[b]
            .push_row(arr.variants.row(i))
            .expect("row width matches ncols");
    }

    let prefix = args.out_prefix.clone().unwrap_or_else(|| default_prefix(&args.input));

    for b in 0..n {
        let out = PathBuf::from(format!("{prefix}.{b}.skf"));
        let bin = SkaArray::<IntT> {
            k: arr.k,
            rc: arr.rc,
            names: arr.names.clone(),
            split_kmers: std::mem::take(&mut keys[b]),
            variants: std::mem::replace(&mut variants[b], Array2::zeros((0, ncols))),
            variant_count: std::mem::take(&mut counts[b]),
            ska_version: arr.ska_version.clone(),
            k_bits: arr.k_bits,
        };
        println!("{} -> {} split k-mers", out.display(), bin.n_kmers());
        bin.save(&out)?;
    }

    println!(
        "split {} ({} split k-mers) into {} bins",
        args.input.display(),
        arr.n_kmers(),
        n
    );
    Ok(())
}

fn run_concat<IntT: SkfInt>(args: &ConcatArgs) -> Result<()> {
    // Load every bin. concat necessarily materialises the full merged matrix, so
    // this is the same memory profile as stacking the bins one at a time.
    let bins: Vec<SkaArray<IntT>> = args
        .inputs
        .iter()
        .map(|p| SkaArray::load(p))
        .collect::<Result<_>>()?;

    let merged = concat_bins(&bins)?;

    println!(
        "concatenated {} bins -> {} ({} split k-mers, {} samples)",
        args.inputs.len(),
        args.output.display(),
        merged.n_kmers(),
        merged.names.len()
    );
    merged.save(&args.output)?;
    Ok(())
}

/// Concatenate per-bin SKFs into one, aligning sample columns by NAME.
///
/// The bins partition the split-kmer space, so their rows are simply stacked.
/// Their `variants` columns, however, may be in different sample orders (each
/// per-bin `ska merge` orders columns by its own input order) — so columns are
/// realigned by sample name to a canonical (sorted-union) order, with any sample
/// absent from a bin filled with the SKA missing byte `b'-'`. Only structural
/// parameters (`k`/`rc`/`k_bits`) must match across bins.
fn concat_bins<IntT: SkfInt>(bins: &[SkaArray<IntT>]) -> Result<SkaArray<IntT>> {
    let first = &bins[0];

    for (i, bin) in bins.iter().enumerate().skip(1) {
        if bin.k != first.k || bin.rc != first.rc || bin.k_bits != first.k_bits {
            bail!(
                "bin {i} has incompatible parameters (k/rc/k_bits) with the first input"
            );
        }
    }

    // Canonical sample order = sorted union of every bin's names. Deterministic, and
    // tolerant of a bin missing a sample that had no k-mers in its minimizer range.
    let canonical: Vec<String> = bins
        .iter()
        .flat_map(|b| b.names.iter().cloned())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let ncols = canonical.len();
    let canon_idx: HashMap<&str, usize> = canonical
        .iter()
        .enumerate()
        .map(|(j, name)| (name.as_str(), j))
        .collect();

    // Remap each bin's variant columns into canonical positions by sample name,
    // filling samples absent from a bin with the SKA missing byte b'-'.
    let mut split_kmers: Vec<IntT> = Vec::new();
    let mut variant_count: Vec<usize> = Vec::new();
    let mut blocks: Vec<Array2<u8>> = Vec::with_capacity(bins.len());
    let mut warned_missing = false;
    for bin in bins {
        if bin.names.len() != ncols && !warned_missing {
            eprintln!(
                "[ska-shard concat] note: bins have differing sample sets; \
                 absent samples are filled with '-' (missing)"
            );
            warned_missing = true;
        }
        let mut block = Array2::from_elem((bin.n_kmers(), ncols), b'-');
        for (c, name) in bin.names.iter().enumerate() {
            let j = canon_idx[name.as_str()];
            block.column_mut(j).assign(&bin.variants.column(c));
        }
        blocks.push(block);
        // variant_count is per-k-mer (count of non-missing bases); reordering
        // columns or adding all-missing columns leaves it unchanged.
        split_kmers.extend(bin.split_kmers.iter().copied());
        variant_count.extend(bin.variant_count.iter().copied());
    }

    let views: Vec<_> = blocks.iter().map(|a| a.view()).collect();
    let variants = if views.is_empty() {
        Array2::zeros((0, ncols))
    } else {
        concatenate(Axis(0), &views).context("stacking variant blocks")?
    };

    Ok(SkaArray::<IntT> {
        k: first.k,
        rc: first.rc,
        names: canonical,
        split_kmers,
        variants,
        variant_count,
        ska_version: first.ska_version.clone(),
        k_bits: first.k_bits,
    })
}

/// Default split output prefix = input file stem (strips the `.skf`).
fn default_prefix(input: &Path) -> String {
    input
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "shard".to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    /// Build a minimal u64 SkaArray bin from sample names and a row-major variant
    /// matrix (rows = k-mers, cols = samples, in `names` order).
    fn bin(names: &[&str], keys: &[u64], variants: Array2<u8>) -> SkaArray<u64> {
        let variant_count = (0..variants.nrows())
            .map(|r| variants.row(r).iter().filter(|&&b| b != b'-').count())
            .collect();
        SkaArray {
            k: 31,
            rc: true,
            names: names.iter().map(|s| s.to_string()).collect(),
            split_kmers: keys.to_vec(),
            variants,
            variant_count,
            ska_version: "test".to_string(),
            k_bits: 64,
        }
    }

    #[test]
    fn concat_realigns_columns_by_name() {
        // Same samples, DIFFERENT column orders across the two bins.
        let b0 = bin(&["A", "B", "C"], &[10, 11], array![[b'A', b'C', b'G'], [b'T', b'T', b'T']]);
        let b1 = bin(&["C", "B", "A"], &[20], array![[b'G', b'C', b'A']]);

        let merged = concat_bins(&[b0, b1]).unwrap();

        // Canonical order is the sorted union: A, B, C.
        assert_eq!(merged.names, vec!["A", "B", "C"]);
        // Rows stacked; b1's columns realigned to [A, B, C].
        assert_eq!(
            merged.variants,
            array![
                [b'A', b'C', b'G'], // b0 row 0
                [b'T', b'T', b'T'], // b0 row 1
                [b'A', b'C', b'G'], // b1 row 0, realigned C,B,A -> A,B,C
            ]
        );
        assert_eq!(merged.split_kmers, vec![10, 11, 20]);
    }

    #[test]
    fn concat_fills_absent_sample_with_missing() {
        // Second bin is missing sample C entirely.
        let b0 = bin(&["A", "B", "C"], &[10], array![[b'A', b'C', b'G']]);
        let b1 = bin(&["B", "A"], &[20], array![[b'T', b'A']]);

        let merged = concat_bins(&[b0, b1]).unwrap();

        assert_eq!(merged.names, vec!["A", "B", "C"]);
        assert_eq!(
            merged.variants,
            array![
                [b'A', b'C', b'G'],  // b0
                [b'A', b'T', b'-'],  // b1: A,B realigned; C absent -> '-'
            ]
        );
        // variant_count is unchanged per row (missing fill adds no non-missing base).
        assert_eq!(merged.variant_count, vec![3, 2]);
    }

    #[test]
    fn concat_rejects_incompatible_k() {
        let mut b1 = bin(&["A"], &[20], array![[b'A']]);
        b1.k = 21;
        let b0 = bin(&["A"], &[10], array![[b'A']]);
        assert!(concat_bins(&[b0, b1]).is_err());
    }
}
