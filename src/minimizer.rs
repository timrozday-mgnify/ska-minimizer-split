//! Minimizer-based bin assignment for split k-mers.
//!
//! A SKA split-kmer key encodes the `(k-1)` flanking bases (the variable middle
//! base lives in the `variants` array, not the key). We decode those flanks back
//! to ACGT, take the canonical ntHash minimizer of the `l`-mer windows, and map
//! it to a bin. Because the key is canonical and the same in every sample, this
//! assignment is deterministic per biological split k-mer — the property that
//! makes "shard → merge each bin → concatenate" equivalent to a direct merge.
//!
//! The decode mirrors ska.rust `src/ska_dict/bit_encoding.rs::decode_kmer`; the
//! ntHash scheme mirrors rust-mdbg's use of canonical `ntc64`.

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

/// Assign a flank sequence to one of `n` bins via its canonical ntHash minimizer.
///
/// `l` is the minimizer (l-mer) length and must satisfy `l <= flank.len()`.
pub fn minimizer_bin(flank: &[u8], l: usize, n: usize) -> Result<usize> {
    ensure!(n >= 1, "number of bins must be >= 1");
    ensure!(
        l >= 1 && l <= flank.len(),
        "minimizer length {l} must be in 1..={}",
        flank.len()
    );
    let iter = NtHashIterator::new(flank, l)
        .map_err(|e| anyhow::anyhow!("ntHash init failed: {e}"))?;
    let min_hash = iter.min().expect("at least one l-mer window");
    Ok((min_hash % n as u64) as usize)
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
        let a = minimizer_bin(flank, 9, n).unwrap();
        let b = minimizer_bin(flank, 9, n).unwrap();
        assert_eq!(a, b);
        assert!(a < n);
    }

    #[test]
    fn rejects_oversized_minimizer() {
        assert!(minimizer_bin(b"ACGT", 5, 4).is_err());
    }
}
