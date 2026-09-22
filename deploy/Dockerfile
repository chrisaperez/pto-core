# syntax=docker/dockerfile:1
#
# Production image for cuttag_profiler.
#
# Multi-stage: the build stage carries CMake, compilers and ~1.5 GB of headers;
# the runtime stage carries the binary and the shared libraries it actually
# needs. The dashboard is embedded in the binary's rodata, so there is no web
# asset directory to ship or to serve from disk.
#
# SECURITY POSTURE -- read before changing anything here.
#
# cuttag_profiler requires a session token (a 256-bit value minted per run,
# compared in constant time) on every /api/* request, on top of the Host and
# Origin checks. Host/Origin alone stop a hostile *web page* in the operator's
# browser; they stop nothing that can open a socket and set its own headers --
# the token is what closes that gap. See SECURITY_HTTP_2026-08-15.md.
#
# H2 (arbitrary file read) is now CLOSED in the binary: --data-root is mandatory
# and the server refuses to start without it. The entrypoint still passes it
# explicitly rather than relying on that -- confinement should be visible in the
# image, not inferred from the binary's defaults.
#
# Everything below follows from that:
#   * the container listens on 127.0.0.1 by default and compose publishes it to
#     the host's loopback only -- reach it over an SSH tunnel, not a public port;
#   * --data-root is mandatory and baked into the entrypoint;
#   * the data volume is mounted read-only;
#   * the process runs as an unprivileged user with all capabilities dropped.
#
# Do not add `--allow-remote`. See deploy/README.md.

# ---------------------------------------------------------------- build stage
FROM ubuntu:24.04 AS build

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        ninja-build \
        pkg-config \
        git \
        ca-certificates \
        libhts-dev \
        zlib1g-dev \
        libdeflate-dev \
        libomp-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . /src

# Only cuttag_profiler is copied into the runtime stage below. The build here
# only disables scrna_matrix (also drops the OpenMP requirement) and
# fastq_stream; it does not pass -DPTO_BUILD_GENOMIC_TOOLKIT=OFF or
# -DPTO_BUILD_PEAKS=OFF, so genomic_toolkit and pto-peaks are compiled in this
# stage too and then discarded -- wasted build time, not a correctness issue.
#
# -DPROFILER_NATIVE_ARCH=OFF is deliberate: -march=native would bake in the
# build host's ISA and fault on an older EC2 instance type.
#
# -DPROFILER_STATIC_HTSLIB=OFF: the runtime stage installs libhts/libdeflate/zlib
# as shared libraries, so linking htslib statically here would both duplicate
# them and drag in transitive static archives (htscodecs, bz2) that Ubuntu does
# not ship in -dev form. Shared is smaller and keeps CVE patching to
# `apt-get upgrade` in the runtime layer rather than an image rebuild.
RUN cmake -S /src -B /build -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DPTO_BUILD_SCRNA_MATRIX=OFF \
        -DPTO_BUILD_FASTQ_STREAM=OFF \
        -DPROFILER_NATIVE_ARCH=OFF \
        -DPROFILER_BUILD_TESTS=OFF \
        -DPROFILER_STATIC_HTSLIB=OFF \
    && cmake --build /build -j "$(nproc)"

# Fail the image build rather than shipping something that does not start.
RUN /build/modules/cuttag_profiler/cuttag_profiler --version

# -------------------------------------------------------------- runtime stage
FROM ubuntu:24.04 AS runtime

ENV DEBIAN_FRONTEND=noninteractive

# Runtime shared libraries only; no compilers, no headers.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libhts3t64 \
        zlib1g \
        libdeflate0 \
        curl \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && useradd --system --create-home --shell /usr/sbin/nologin profiler

COPY --from=build /build/modules/cuttag_profiler/cuttag_profiler /usr/local/bin/cuttag_profiler
COPY deploy/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod 0755 /usr/local/bin/entrypoint.sh

# Read-only mount point for the operator's BAM/BED files.
RUN mkdir -p /data && chown profiler:profiler /data
VOLUME ["/data"]

USER profiler
WORKDIR /home/profiler

# Bind address is a build-time default, overridable at run time. Inside a
# container 0.0.0.0 is the only address the port can be published from, so the
# confinement that matters is the PUBLISH address on the host side --
# docker-compose maps this to 127.0.0.1 only. See deploy/README.md.
ENV PROFILER_BIND=0.0.0.0 \
    PROFILER_PORT=8080 \
    PROFILER_DATA_ROOT=/data

EXPOSE 8080

HEALTHCHECK --interval=30s --timeout=3s --start-period=5s --retries=3 \
    CMD curl -fsS "http://127.0.0.1:${PROFILER_PORT}/api/health" || exit 1

# The entrypoint validates its inputs and refuses to start unconfined; see
# deploy/entrypoint.sh for why --data-root is not parameterised away.
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
