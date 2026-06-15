//! Library for `ska-shard`: read/write ska2 `.skf` files and bin their split
//! k-mers by minimizer value.
//!
//! - [`skf`] mirrors ska2's on-disk `MergeSkaArray` for faithful round-trips.
//! - [`minimizer`] assigns each split k-mer to a bin deterministically.

pub mod minimizer;
pub mod skf;
