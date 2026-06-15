//! `ska-shard` — split a ska2 `.skf` into minimizer bins, or concatenate bins.
//!
//! Workflow this enables (bounds `ska merge` peak memory to ~1/n):
//!   1. `ska build` each sample → per-sample `.skf`
//!   2. `ska-shard split sample.skf -n N`  (run per sample)
//!   3. `ska merge` bin i across all samples (parallel, ~1/n of the key space)
//!   4. `ska-shard concat merged_bin*.skf -o merged.skf`

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
    let first: SkaArray<IntT> = SkaArray::load(&args.inputs[0])?;
    let ncols = first.names.len();

    let mut split_kmers: Vec<IntT> = first.split_kmers.clone();
    let mut variant_count: Vec<usize> = first.variant_count.clone();
    let mut variant_blocks: Vec<Array2<u8>> = vec![first.variants.clone()];

    for path in &args.inputs[1..] {
        let arr: SkaArray<IntT> = SkaArray::load(path)?;
        check_compatible(&first, &arr, path)?;
        split_kmers.extend(arr.split_kmers);
        variant_count.extend(arr.variant_count);
        variant_blocks.push(arr.variants);
    }

    // Vertically stack the per-bin middle-base blocks (columns already aligned
    // by the identical sample ordering check above).
    let views: Vec<_> = variant_blocks.iter().map(|a| a.view()).collect();
    let variants = if views.is_empty() {
        Array2::zeros((0, ncols))
    } else {
        concatenate(Axis(0), &views).context("stacking variant blocks")?
    };

    let merged = SkaArray::<IntT> {
        k: first.k,
        rc: first.rc,
        names: first.names.clone(),
        split_kmers,
        variants,
        variant_count,
        ska_version: first.ska_version.clone(),
        k_bits: first.k_bits,
    };

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

/// Reject inputs whose k / rc / k_bits / sample set or ordering differ, since
/// concatenation aligns `variants` columns by position.
fn check_compatible<IntT: SkfInt>(
    a: &SkaArray<IntT>,
    b: &SkaArray<IntT>,
    b_path: &Path,
) -> Result<()> {
    if a.k != b.k || a.rc != b.rc || a.k_bits != b.k_bits {
        bail!(
            "{} has incompatible parameters (k/rc/k_bits) with the first input",
            b_path.display()
        );
    }
    if a.names != b.names {
        bail!(
            "{} has a different sample set or ordering than the first input; \
             all bins must come from the same ordered `ska merge` inputs",
            b_path.display()
        );
    }
    Ok(())
}

/// Default split output prefix = input file stem (strips the `.skf`).
fn default_prefix(input: &Path) -> String {
    input
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "shard".to_string())
}
