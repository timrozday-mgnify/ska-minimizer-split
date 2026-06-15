# Build the ska-shard binary, then ship it on a slim Debian base.
# Runtime keeps a shell (/bin/bash) so Nextflow's task wrapper can run — do not
# switch to distroless.
# Needs a Cargo new enough for edition2024 transitive deps (Rust >= 1.85).
FROM rust:1-slim AS builder
WORKDIR /build
COPY . .
RUN cargo build --release --locked

FROM debian:bookworm-slim
LABEL org.opencontainers.image.source=https://github.com/timrozday-mgnify/ska-minimizer-split
LABEL org.opencontainers.image.description="ska-shard: split/concat ska2 .skf files by minimizer"
# procps provides `ps`, which Nextflow needs to collect task metrics.
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
COPY --from=builder /build/target/release/ska-shard /usr/local/bin/ska-shard
CMD ["ska-shard", "--help"]
