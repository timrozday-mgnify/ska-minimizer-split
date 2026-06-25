//! `ska-shard` — split a ska2 `.skf` into hash bins, subset it, or concatenate bins.
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
use ska_minimizer_split::minimizer::{decode_flank, flank_bin, flank_hash};
use ska_minimizer_split::skf::{SkaArray, SkfInt};

#[derive(Parser)]
#[command(
    name = "ska-shard",
    about = "Split a ska2 .skf into hash bins, subset it, or concatenate bins back together",
    version
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Split one .skf into N bins by split-kmer hash value.
    Split(SplitArgs),
    /// Write a sparse hash-selected subset of one .skf.
    Subset(SubsetArgs),
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
    /// Output prefix; bins are written as <prefix>.<i>.skf. Defaults to the
    /// input file stem.
    #[arg(short = 'o', long = "out-prefix")]
    out_prefix: Option<String>,
}

#[derive(Parser)]
struct SubsetArgs {
    /// Input .skf file.
    input: PathBuf,
    /// Output subset .skf.
    #[arg(short = 'o', long = "output")]
    output: PathBuf,
    /// Fraction of the ntHash domain to retain, from 0 to 1.
    #[arg(long = "sparsity")]
    sparsity: Option<f64>,
    /// Inclusive lower hash bound.
    #[arg(long = "min-hash")]
    min_hash: Option<u64>,
    /// Inclusive upper hash bound.
    #[arg(long = "max-hash")]
    max_hash: Option<u64>,
    /// Inclusive lower bound as a fraction of u64::MAX.
    #[arg(long = "min-proportion")]
    min_proportion: Option<f64>,
    /// Inclusive upper bound as a fraction of u64::MAX.
    #[arg(long = "max-proportion")]
    max_proportion: Option<f64>,
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
        Command::Subset(args) => {
            dispatch_skf_width!(&args.input, IntT => run_subset::<IntT>(&args))
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

    let ncols = arr.names.len();
    let n = args.bins;

    log_hash_context(
        "split",
        &args.input,
        None,
        k,
        arr.k_bits,
        arr.n_kmers(),
        arr.names.len(),
    );

    // Per-bin row buffers.
    let mut keys: Vec<Vec<IntT>> = vec![Vec::new(); n];
    let mut counts: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut variants: Vec<Array2<u8>> = (0..n).map(|_| Array2::zeros((0, ncols))).collect();

    for i in 0..arr.n_kmers() {
        let flank = decode_flank(arr.split_kmers[i].into(), k);
        let b = flank_bin(&flank, n)?;
        keys[b].push(arr.split_kmers[i]);
        counts[b].push(arr.variant_count[i]);
        variants[b]
            .push_row(arr.variants.row(i))
            .expect("row width matches ncols");
    }

    let prefix = args
        .out_prefix
        .clone()
        .unwrap_or_else(|| default_prefix(&args.input));

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

fn run_subset<IntT: SkfInt>(args: &SubsetArgs) -> Result<()> {
    let selection = HashSelection::from_args(args)?;
    let arr: SkaArray<IntT> = SkaArray::load(&args.input)?;
    let subset = subset_array(&arr, selection)?;

    log_hash_context(
        "subset",
        &args.input,
        Some(&args.output),
        arr.k,
        arr.k_bits,
        arr.n_kmers(),
        arr.names.len(),
    );
    println!("selection: {}", selection.describe());
    if selection.range_overrides_sparsity() {
        println!(
            "selection note: explicit hash/proportion bounds were supplied; \
             --sparsity was ignored"
        );
    }
    println!(
        "subset retained {} of {} split k-mers (dropped {}, retained fraction {:.6})",
        subset.n_kmers(),
        arr.n_kmers(),
        arr.n_kmers() - subset.n_kmers(),
        retained_fraction(subset.n_kmers(), arr.n_kmers())
    );

    subset.save(&args.output)?;
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

#[derive(Debug, Clone, Copy, PartialEq)]
enum HashSelection {
    Empty {
        ignored_sparsity: bool,
    },
    Range {
        min: u64,
        max: u64,
        explicit_range: bool,
        ignored_sparsity: bool,
    },
}

impl HashSelection {
    fn from_args(args: &SubsetArgs) -> Result<Self> {
        validate_optional_proportion(args.sparsity, "--sparsity")?;
        validate_optional_proportion(args.min_proportion, "--min-proportion")?;
        validate_optional_proportion(args.max_proportion, "--max-proportion")?;

        ensure!(
            !(args.min_hash.is_some() && args.min_proportion.is_some()),
            "use only one of --min-hash or --min-proportion"
        );
        ensure!(
            !(args.max_hash.is_some() && args.max_proportion.is_some()),
            "use only one of --max-hash or --max-proportion"
        );

        let explicit_range = args.min_hash.is_some()
            || args.max_hash.is_some()
            || args.min_proportion.is_some()
            || args.max_proportion.is_some();
        let ignored_sparsity = explicit_range && args.sparsity.is_some();

        if explicit_range {
            let min = args
                .min_hash
                .or_else(|| args.min_proportion.map(proportion_to_hash))
                .unwrap_or(0);
            let max = args
                .max_hash
                .or_else(|| args.max_proportion.map(proportion_to_hash))
                .unwrap_or(u64::MAX);
            ensure!(
                min <= max,
                "minimum hash bound must be <= maximum hash bound"
            );
            return Ok(Self::Range {
                min,
                max,
                explicit_range,
                ignored_sparsity,
            });
        }

        let sparsity = args
            .sparsity
            .context("provide --sparsity or at least one explicit hash/proportion bound")?;
        if sparsity == 0.0 {
            return Ok(Self::Empty { ignored_sparsity });
        }
        Ok(Self::Range {
            min: 0,
            max: proportion_to_hash(sparsity),
            explicit_range,
            ignored_sparsity,
        })
    }

    fn contains(self, hash: u64) -> bool {
        match self {
            Self::Empty { .. } => false,
            Self::Range { min, max, .. } => min <= hash && hash <= max,
        }
    }

    fn describe(self) -> String {
        match self {
            Self::Empty { .. } => "sparsity mode, retaining no hash values".to_string(),
            Self::Range {
                min,
                max,
                explicit_range,
                ..
            } => {
                let mode = if explicit_range {
                    "explicit range mode"
                } else {
                    "sparsity mode"
                };
                format!("{mode}, inclusive hash range {min}..={max}")
            }
        }
    }

    fn range_overrides_sparsity(self) -> bool {
        match self {
            Self::Empty { ignored_sparsity } => ignored_sparsity,
            Self::Range {
                ignored_sparsity, ..
            } => ignored_sparsity,
        }
    }
}

fn subset_array<IntT: SkfInt>(
    arr: &SkaArray<IntT>,
    selection: HashSelection,
) -> Result<SkaArray<IntT>> {
    let ncols = arr.names.len();
    let mut split_kmers = Vec::new();
    let mut variant_count = Vec::new();
    let mut variants = Array2::zeros((0, ncols));

    for i in 0..arr.n_kmers() {
        let flank = decode_flank(arr.split_kmers[i].into(), arr.k);
        let hash = flank_hash(&flank)?;
        if selection.contains(hash) {
            split_kmers.push(arr.split_kmers[i]);
            variant_count.push(arr.variant_count[i]);
            variants
                .push_row(arr.variants.row(i))
                .expect("row width matches ncols");
        }
    }

    Ok(SkaArray::<IntT> {
        k: arr.k,
        rc: arr.rc,
        names: arr.names.clone(),
        split_kmers,
        variants,
        variant_count,
        ska_version: arr.ska_version.clone(),
        k_bits: arr.k_bits,
    })
}

fn proportion_to_hash(proportion: f64) -> u64 {
    if proportion <= 0.0 {
        0
    } else if proportion >= 1.0 {
        u64::MAX
    } else {
        (proportion * u64::MAX as f64).floor() as u64
    }
}

fn validate_optional_proportion(value: Option<f64>, flag: &str) -> Result<()> {
    if let Some(value) = value {
        ensure!(
            (0.0..=1.0).contains(&value),
            "{flag} must be between 0 and 1"
        );
    }
    Ok(())
}

fn retained_fraction(retained: usize, total: usize) -> f64 {
    if total == 0 {
        0.0
    } else {
        retained as f64 / total as f64
    }
}

fn log_hash_context(
    action: &str,
    input: &Path,
    output: Option<&Path>,
    k: usize,
    k_bits: u32,
    n_kmers: usize,
    n_samples: usize,
) {
    let flank_len = k.saturating_sub(1);
    match output {
        Some(output) => println!(
            "{action} {} -> {} ({} split k-mers, {} samples)",
            input.display(),
            output.display(),
            n_kmers,
            n_samples
        ),
        None => println!(
            "{action} {} ({} split k-mers, {} samples)",
            input.display(),
            n_kmers,
            n_samples
        ),
    }
    println!(
        "hash context: k={k}, flank length={flank_len}, k_bits={k_bits}, hash domain=0..={}",
        u64::MAX
    );
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
            bail!("bin {i} has incompatible parameters (k/rc/k_bits) with the first input");
        }
    }

    // Canonical sample order = sorted union of every bin's names. Deterministic, and
    // tolerant of a bin missing a sample that had no k-mers in its hash range.
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

    fn subset_args() -> SubsetArgs {
        SubsetArgs {
            input: PathBuf::from("in.skf"),
            output: PathBuf::from("out.skf"),
            sparsity: None,
            min_hash: None,
            max_hash: None,
            min_proportion: None,
            max_proportion: None,
        }
    }

    #[test]
    fn concat_realigns_columns_by_name() {
        // Same samples, DIFFERENT column orders across the two bins.
        let b0 = bin(
            &["A", "B", "C"],
            &[10, 11],
            array![[b'A', b'C', b'G'], [b'T', b'T', b'T']],
        );
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
                [b'A', b'C', b'G'], // b0
                [b'A', b'T', b'-'], // b1: A,B realigned; C absent -> '-'
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

    #[test]
    fn sparsity_zero_selects_empty_hash_set() {
        let mut args = subset_args();
        args.sparsity = Some(0.0);
        let selection = HashSelection::from_args(&args).unwrap();
        assert!(!selection.contains(0));
        assert!(!selection.contains(u64::MAX));
    }

    #[test]
    fn sparsity_one_selects_full_hash_set() {
        let mut args = subset_args();
        args.sparsity = Some(1.0);
        let selection = HashSelection::from_args(&args).unwrap();
        assert!(selection.contains(0));
        assert!(selection.contains(u64::MAX));
    }

    #[test]
    fn explicit_range_overrides_sparsity() {
        let mut args = subset_args();
        args.sparsity = Some(1.0);
        args.min_hash = Some(10);
        args.max_hash = Some(20);
        let selection = HashSelection::from_args(&args).unwrap();
        assert!(selection.range_overrides_sparsity());
        assert!(!selection.contains(9));
        assert!(selection.contains(10));
        assert!(selection.contains(20));
        assert!(!selection.contains(21));
    }

    #[test]
    fn proportional_bounds_map_to_hash_bounds() {
        let mut args = subset_args();
        args.min_proportion = Some(0.0);
        args.max_proportion = Some(1.0);
        let selection = HashSelection::from_args(&args).unwrap();
        assert!(selection.contains(0));
        assert!(selection.contains(u64::MAX));
    }

    #[test]
    fn conflicting_bound_types_are_rejected() {
        let mut args = subset_args();
        args.min_hash = Some(1);
        args.min_proportion = Some(0.25);
        assert!(HashSelection::from_args(&args).is_err());
    }

    #[test]
    fn subset_preserves_selected_rows_and_metadata() {
        let arr = bin(
            &["A", "B"],
            &[10, 20, 30],
            array![[b'A', b'C'], [b'G', b'T'], [b'C', b'C']],
        );

        let selection = HashSelection::Range {
            min: 0,
            max: u64::MAX,
            explicit_range: true,
            ignored_sparsity: false,
        };
        let subset = subset_array(&arr, selection).unwrap();

        assert_eq!(subset.k, arr.k);
        assert_eq!(subset.rc, arr.rc);
        assert_eq!(subset.names, arr.names);
        assert_eq!(subset.split_kmers, vec![10, 20, 30]);
        assert_eq!(subset.variant_count, vec![2, 2, 2]);
        assert_eq!(
            subset.variants,
            array![[b'A', b'C'], [b'G', b'T'], [b'C', b'C']]
        );
    }

    #[test]
    fn subset_empty_selection_keeps_columns() {
        let arr = bin(&["A", "B"], &[10], array![[b'A', b'C']]);

        let subset = subset_array(
            &arr,
            HashSelection::Empty {
                ignored_sparsity: false,
            },
        )
        .unwrap();

        assert_eq!(subset.names, vec!["A", "B"]);
        assert_eq!(subset.split_kmers, Vec::<u64>::new());
        assert_eq!(subset.variant_count, Vec::<usize>::new());
        assert_eq!(subset.variants.shape(), &[0, 2]);
    }
}
