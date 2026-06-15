# Build a fully static (musl) ska-shard binary, then ship it on Alpine.
#
# Why this shape for OCI-runtime compatibility (Singularity/Apptainer --oci):
#   * A static musl binary has no libc/base coupling, so it runs on any base.
#   * Alpine is a small, ubiquitous, OCI-runtime-friendly base with a clean
#     /dev + rootfs that converts and mounts reliably under the OCI runtime.
#   * bash is added because Nextflow's task wrapper uses /bin/bash; procps gives
#     a real `ps`, which Nextflow needs to collect task metrics.
FROM rust:1-alpine AS builder
RUN apk add --no-cache musl-dev
WORKDIR /build
COPY . .
RUN cargo build --release --locked

FROM alpine:3.20
LABEL org.opencontainers.image.source=https://github.com/timrozday-mgnify/ska-minimizer-split
LABEL org.opencontainers.image.description="ska-shard: split/concat ska2 .skf files by minimizer"
RUN apk add --no-cache bash procps
COPY --from=builder /build/target/release/ska-shard /usr/local/bin/ska-shard
CMD ["ska-shard", "--help"]
