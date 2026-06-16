# ska-minimizer-split

`ska-shard` — an accessory to [ska2](https://github.com/bacpop/ska.rust) that
partitions a `.skf` (split k-mer file) into `n` bins by **minimizer value**, and
concatenates per-bin `.skf` files back into one.

## Why

`ska merge` deserialises every input `.skf` fully into memory, so peak RAM scales
with the whole split-k-mer space. A SKA split-kmer key is *canonical* — a
deterministic function of the biological k-mer, independent of which sample it
came from — so any pure function of the key (here, a minimizer hash) assigns the
**same** k-mer to the **same** bin in every sample. That makes the merge
partitionable:

```
1. ska build each sample                  -> sample.skf
2. ska-shard split sample.skf -n N         (run per sample) -> sample.0.skf .. sample.{N-1}.skf
3. ska merge bin i across all samples      (parallel; each merge touches ~1/N of the key space)
4. ska-shard concat merged_bin*.skf -o merged.skf
```

The resulting `merged.skf` is equivalent to a direct `ska merge` of all samples
(verified: identical `ska align` output), but each per-bin merge holds only
~1/N of the k-mers, bounding peak memory.

## Usage

```bash
# Split one .skf into N bins (written as <prefix>.<i>.skf; prefix defaults to input stem)
ska-shard split sample.skf -n 8 [-l 9] [-o prefix]

# Concatenate per-bin .skf files (must share identical, identically-ordered samples)
ska-shard concat merged_bin0.skf merged_bin1.skf ... -o merged.skf
```

- `-n/--bins` number of bins.
- `-l/--minimizer-len` minimizer (l-mer) length, must be `<= k-1` (default 9).
- Split is per-file; run it on each sample's `.skf` with the **same** `-n` and
  `-l` so corresponding bins are mergeable.

## Container

A prebuilt image is published to GHCR on each version tag:

```bash
# Docker / Podman
docker pull ghcr.io/timrozday-mgnify/ska-minimizer-split:0.1.1
docker run --rm ghcr.io/timrozday-mgnify/ska-minimizer-split:0.1.1 ska-shard --help

# Singularity / Apptainer — pull the prebuilt SIF (no on-node OCI->SIF conversion)
singularity pull oras://ghcr.io/timrozday-mgnify/ska-minimizer-split:0.1.1-sif
```

The image is built on Alpine with a fully static (musl) binary, and a ready-made
`.sif` is published via ORAS so HPC nodes download and mount it directly instead
of converting the OCI image themselves.

Build it locally with `docker build -t ghcr.io/timrozday-mgnify/ska-minimizer-split:0.1.1 .`.
The [subspecies-phylogeny](https://github.com/timrozday-mgnify/subspecies-phylogeny)
pipeline consumes this image in its `SKA2_SHARD_SPLIT` / `SKA2_SHARD_CONCAT`
Nextflow modules.

## How it works

- **SKF I/O** mirrors ska2's `MergeSkaArray` (snappy-framed CBOR via `ciborium`
  + `snap`, `ndarray` for the variant matrix), so files are byte-compatible and
  round-trip through ska2 unchanged. Integer width (u64 for k≤31, u128 for
  31<k≤63) is detected from the file's `k`.
- **Minimizers** decode each split-kmer key to its `(k-1)` flanking bases and
  take the canonical [ntHash](https://crates.io/crates/nthash) minimizer of the
  `l`-mer windows, then `bin = minimizer % n`. Same scheme as
  [rust-mdbg](https://github.com/timrozday-mgnify/rust-mdbg).
- **Concat** validates matching `k`/`rc`/`k_bits`/sample-ordering, then row-wise
  concatenates the bins. SKA imposes no ordering on split k-mers, so no sort is
  needed.

## Testing

```bash
cargo test                 # unit tests (decode/encode inverse, binning, validation)
```

End-to-end verification against real ska2 (requires Docker): build three
genomes into `.skf`s, shard each into bins, `ska merge` per bin, `ska-shard
concat`, and confirm `ska align` output matches a direct merge column-for-column.
