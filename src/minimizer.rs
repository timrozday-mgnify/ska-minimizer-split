//! Bin assignment for split k-mers by hashing their full flank.
//!
//! A SKA split-kmer key encodes the `(k-1)` flanking bases (the variable middle
//! base lives in the `variants` array, not the key). We decode those flanks back
//! to ACGT and take the canonical ntHash of the **entire** flank sequence, then
//! map the hash to a bin. Because the key is canonical and the same in every
//! sample, this assignment is deterministic per biological split k-mer — the
//! property that makes "shard → merge each bin → concatenate" equivalent to a
//! direct merge.
//!
//! The decode mirrors ska.rust `src/ska_dict/bit_encoding.rs::decode_kmer`; the
//! hash uses the canonical ntHash scheme from the `nthash` crate.

use anyhow::{ensure, Result};
use nthash::NtHashIterator;

/// 2-bit base code → ASCII, matching ska2's `LETTER_CODE`.
const LETTER_CODE: [u8; 4] = [b'A', b'C', b'T', b'G'];

/// Decode the `(k-1)` flanking bases of a split k-mer key to an ASCII sequence.
///
/// Port of ska2 `decode_kmer`: the key holds an upper and a lower half of
/// `half_k = (k-1)/2` bases each; within a half the bases are decoded LSB-first
/// then reversed. The returned `Vec` is `upper ++ lower`, length `k-1`.
pub fn decode_flank(key: u128, k: usize) -> Vec<u8> {
    let half_k = (k - 1) / 2;
    let shift = half_k * 2;
    let lower_mask: u128 = if shift >= 128 { u128::MAX } else { (1u128 << shift) - 1 };

    let decode_half = |mut bits: u128| -> Vec<u8> {
        let mut half = Vec::with_capacity(half_k);
        for _ in 0..half_k {
            half.push(LETTER_CODE[(bits & 0x3) as usize]);
            bits >>= 2;
        }
        half.reverse();
        half
    };

    let mut flank = decode_half((key >> shift) & lower_mask);
    flank.extend(decode_half(key & lower_mask));
    flank
}

/// Assign a flank sequence to one of `n` bins via the canonical ntHash of the full flank.
///
/// Hashing the entire `(k-1)` flank (rather than taking a minimizer over shorter
/// windows) eliminates a length parameter and gives better hash entropy.
pub fn flank_bin(flank: &[u8], n: usize) -> Result<usize> {
    ensure!(n >= 1, "number of bins must be >= 1");
    ensure!(!flank.is_empty(), "flank must be non-empty");
    let mut iter = NtHashIterator::new(flank, flank.len())
        .map_err(|e| anyhow::anyhow!("ntHash init failed: {e}"))?;
    let hash = iter.next().expect("full-flank window always exists");
    Ok((hash % n as u64) as usize)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Encode `(k-1)` ACGT bases into a split-kmer key the way ska2 would, so we
    /// can check `decode_flank` is its exact inverse.
    fn encode_flank(seq: &[u8], k: usize) -> u128 {
        let half_k = (k - 1) / 2;
        assert_eq!(seq.len(), 2 * half_k);
        let encode_base = |b: u8| ((b >> 1) & 0x3) as u128; // ska2 encode_base
        let mut key: u128 = 0;
        for &b in seq {
            key = (key << 2) | encode_base(b);
        }
        key
    }

    #[test]
    fn decode_is_inverse_of_encode() {
        let k = 31; // half_k = 15, flank length 30
        let seq = b"ACGTACGTACGTACGGTCAGTCAGTCAGTC"; // 30 bases, ACGT only
        let key = encode_flank(seq, k);
        assert_eq!(decode_flank(key, k), seq.to_vec());
    }

    #[test]
    fn decode_known_small_key() {
        // k=5 → half_k=2, flank length 4. Bases A=0,C=1,T=2,G=3 (ska2 codes).
        // seq "ACGT": codes 0,1,3,2 → key = 0b00_01_11_10 = 30.
        let k = 5;
        let key = encode_flank(b"ACGT", k);
        assert_eq!(key, 0b00_01_11_10);
        assert_eq!(decode_flank(key, k), b"ACGT".to_vec());
    }

    #[test]
    fn bin_is_deterministic_and_in_range() {
        let flank = b"ACGTACGTACGTACGGTCAGTCAGTCAGTC";
        let n = 7;
        let a = flank_bin(flank, n).unwrap();
        let b = flank_bin(flank, n).unwrap();
        assert_eq!(a, b);
        assert!(a < n);
    }
}
