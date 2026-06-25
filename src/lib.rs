//! Library for `ska-shard`: read/write ska2 `.skf` files and bin their split
//! k-mers by full-flank hash value.
//!
//! - [`skf`] mirrors ska2's on-disk `MergeSkaArray` for faithful round-trips.
//! - [`minimizer`] hashes split-kmer flanks and assigns them to bins deterministically.

pub mod minimizer;
pub mod skf;
