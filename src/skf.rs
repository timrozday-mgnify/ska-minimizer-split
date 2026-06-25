//! Read/write ska2 `.skf` files.
//!
//! A `.skf` is a snappy-framed CBOR stream encoding ska.rust's `MergeSkaArray`
//! struct. We mirror that struct field-for-field so that a file written here is
//! indistinguishable from one written by ska2 itself. ciborium serialises a
//! derived struct as a CBOR map keyed by field name, so it is the field *names*
//! and *types* that must match — order is irrelevant.
//!
//! See ska.rust `src/merge_ska_array.rs` (struct + `save`/`load`).

use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

use anyhow::{Context, Result};
use ndarray::Array2;
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};

/// In-memory mirror of ska2's `MergeSkaArray<IntT>`.
///
/// `IntT` is `u64` for k ≤ 31 and `u128` for 31 < k ≤ 63 — the same split-kmer
/// integer width ska2 chose when the file was written. Keeping the loaded width
/// is what makes the round-trip byte-faithful.
#[derive(Serialize, Deserialize, Clone)]
pub struct SkaArray<IntT> {
    /// K-mer size.
    pub k: usize,
    /// Whether reverse-complement split k-mers were used.
    pub rc: bool,
    /// Sample names (columns of `variants`).
    pub names: Vec<String>,
    /// Canonical split-kmer integer keys (rows of `variants`).
    pub split_kmers: Vec<IntT>,
    /// Middle bases: rows ↔ `split_kmers`, columns ↔ `names`.
    pub variants: Array2<u8>,
    /// Count of non-missing bases per split k-mer.
    pub variant_count: Vec<usize>,
    /// ska version string embedded by the writer.
    pub ska_version: String,
    /// Bits used per split k-mer.
    pub k_bits: u32,
}

/// Trait alias for the two concrete integer widths a `.skf` can use.
///
/// `Into<u128>` lets [`crate::minimizer`] decode either width through a single
/// code path; `Ord`/`Copy` cover the row-shuffling the tools perform.
pub trait SkfInt: Serialize + DeserializeOwned + Copy + Ord + Into<u128> + 'static {}
impl SkfInt for u64 {}
impl SkfInt for u128 {}

impl<IntT: SkfInt> SkaArray<IntT> {
    /// Load a `.skf` (snappy → CBOR) into the chosen integer width.
    pub fn load(path: &Path) -> Result<Self> {
        let reader = BufReader::new(
            File::open(path).with_context(|| format!("opening {}", path.display()))?,
        );
        let decompress = snap::read::FrameDecoder::new(reader);
        let obj: Self = ciborium::de::from_reader(decompress)
            .with_context(|| format!("decoding skf {}", path.display()))?;
        Ok(obj)
    }

    /// Write a `.skf` (CBOR → snappy), matching ska2's `save` exactly.
    pub fn save(&self, path: &Path) -> Result<()> {
        let writer = BufWriter::new(
            File::create(path).with_context(|| format!("creating {}", path.display()))?,
        );
        let mut compress = snap::write::FrameEncoder::new(writer);
        ciborium::ser::into_writer(self, &mut compress)
            .with_context(|| format!("encoding skf {}", path.display()))?;
        Ok(())
    }

    /// Number of split k-mers (rows).
    pub fn n_kmers(&self) -> usize {
        self.split_kmers.len()
    }
}

/// Minimal header used only to read `k` without committing to an integer width.
///
/// serde ignores the other CBOR map entries, so this deserialises from a full
/// `.skf` regardless of whether it is a u64 or u128 file.
#[derive(Deserialize)]
struct KPeek {
    k: usize,
}

/// Read just the k-mer size from a `.skf`, to pick the integer width.
pub fn peek_k(path: &Path) -> Result<usize> {
    let reader =
        BufReader::new(File::open(path).with_context(|| format!("opening {}", path.display()))?);
    let decompress = snap::read::FrameDecoder::new(reader);
    let header: KPeek = ciborium::de::from_reader(decompress)
        .with_context(|| format!("reading k from {}", path.display()))?;
    Ok(header.k)
}

/// Does this k use 64-bit split k-mers? (ska2: k ≤ 31 → u64, else u128.)
pub fn is_u64(k: usize) -> bool {
    k <= 31
}

/// Run `$body` with `SkaArray<$T>` selected from the k of `$path`.
///
/// Rust has no generic closures, so width dispatch is a macro: it peeks k, picks
/// `u64`/`u128`, and binds the concrete type as `$T` inside `$body`.
#[macro_export]
macro_rules! dispatch_skf_width {
    ($path:expr, $T:ident => $body:expr) => {{
        let k = $crate::skf::peek_k($path)?;
        if $crate::skf::is_u64(k) {
            type $T = u64;
            $body
        } else {
            type $T = u128;
            $body
        }
    }};
}
