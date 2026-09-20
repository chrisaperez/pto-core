#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
#
# Build htslib (and libdeflate) inside a manylinux container, STATIC ONLY.
#
# Run by cibuildwheel's `before-all` on Linux. manylinux has no htslib package,
# and htslib is not optional for this wheel: `modules/cuttag_profiler`'s
# CMakeLists says `pkg_check_modules(HTSLIB REQUIRED htslib)` -- BAM reading is
# that module's whole premise, not an add-on.
#
# WHY STATIC ONLY, AND WHY THAT IS ENFORCED BY DELETING THE .so
# -------------------------------------------------------------
# `auditwheel` vendors shared libraries for the extension modules it finds in a
# wheel. This wheel contains no extension module -- it contains three ordinary
# ELF executables -- and relying on auditwheel to notice and patch those is
# relying on behaviour that is incidental to its purpose. A missing vendored
# libhts.so does not fail the build; it fails at the user's first BAM, on a
# machine we cannot see.
#
# Only `cuttag_profiler` has a PROFILER_STATIC_HTSLIB option;
# `genomic_toolkit` and `fastq_stream` link whatever `pkg-config --libs htslib`
# resolves to. So rather than adding an option to two more modules -- a change
# to module internals this package has no business making -- the .so is removed
# after install. `-lhts` then has exactly one thing to resolve to, and all three
# binaries are static with respect to htslib whether or not they asked to be.
#
# ci.yml's `static-link` job already proves this configuration links and passes
# ctest, and asserts the result carries no libhts-family .so dependency. This
# script puts the wheel on that same footing.
#
# libcurl is DISABLED, deliberately and twice over. It is what drags nghttp2 and
# brotli into a static link -- the exact closure problem ci.yml records as the
# reason static htslib is unsupported on macOS -- and it is also htslib's remote
# CRAM reference fetcher, which every deployment of these tools already refuses:
# `job_specs._htslib_offline_env` pins REF_PATH locally because a VPC with no
# NAT turns that fetch into a hang rather than an error. A cluster node behind a
# proxy behaves the same way. No libcurl, no hang.
set -euo pipefail

HTSLIB_VERSION="${PTO_HTSLIB_VERSION:-1.21}"
LIBDEFLATE_VERSION="${PTO_LIBDEFLATE_VERSION:-1.22}"
PREFIX="${PTO_DEPS_PREFIX:-/usr/local}"

log() { printf '\n[htslib-build] %s\n' "$*" >&2; }

# `yum` is present on manylinux_2_28 (AlmaLinux 8) and manylinux2014.
log "installing build dependencies"
yum install -y bzip2-devel xz-devel zlib-devel gcc gcc-c++ make autoconf tar >/dev/null

work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
cd "$work"

# Optional digest pinning.
#
# No SHA256 is hardcoded here on purpose: a constant that nobody in this
# repository has verified against the upstream release is not a supply-chain
# control, it is a decoration that breaks the build the first time upstream
# re-rolls a tarball. Set PTO_HTSLIB_SHA256 in the workflow (or a release
# checklist) to a digest you actually checked, and it is enforced.
verify() {
  local file="$1" expected="${2:-}"
  [ -n "$expected" ] || { log "no digest pinned for $file (set PTO_HTSLIB_SHA256)"; return 0; }
  echo "${expected}  ${file}" | sha256sum -c -
}

log "fetching libdeflate ${LIBDEFLATE_VERSION}"
curl -fsSL --retry 3 -o libdeflate.tar.gz \
  "https://github.com/ebiggers/libdeflate/archive/refs/tags/v${LIBDEFLATE_VERSION}.tar.gz"
verify libdeflate.tar.gz "${PTO_LIBDEFLATE_SHA256:-}"
tar xzf libdeflate.tar.gz
cmake -S "libdeflate-${LIBDEFLATE_VERSION}" -B build-libdeflate \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
  -DLIBDEFLATE_BUILD_SHARED_LIB=OFF \
  -DLIBDEFLATE_BUILD_GZIP=OFF
cmake --build build-libdeflate -j"$(nproc)"
cmake --install build-libdeflate

log "fetching htslib ${HTSLIB_VERSION}"
curl -fsSL --retry 3 -o htslib.tar.bz2 \
  "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2"
verify htslib.tar.bz2 "${PTO_HTSLIB_SHA256:-}"
tar xjf htslib.tar.bz2
cd "htslib-${HTSLIB_VERSION}"

# --disable-plugins keeps the reader in-process: a plugin is a dlopen at
# runtime, which a statically linked, relocatable wheel cannot satisfy.
./configure \
  --prefix="${PREFIX}" \
  --with-libdeflate \
  --disable-libcurl \
  --disable-plugins \
  --disable-gcs \
  --disable-s3 \
  CFLAGS="-fPIC -O2"
make -j"$(nproc)"
make install

log "removing shared htslib so -lhts can only resolve statically"
rm -fv "${PREFIX}"/lib/libhts.so*

# Prove the intent held rather than assuming it. A .pc that still names a
# shared library, or a missing archive, would otherwise surface as a wheel that
# builds and then cannot open a BAM on the user's machine.
test -f "${PREFIX}/lib/libhts.a" || { log "ERROR: libhts.a was not installed"; exit 1; }
if compgen -G "${PREFIX}/lib/libhts.so*" >/dev/null; then
  log "ERROR: a shared libhts survived; the binaries would link against it"
  exit 1
fi

ldconfig 2>/dev/null || true
log "htslib ${HTSLIB_VERSION} installed static-only under ${PREFIX}"
PKG_CONFIG_PATH="${PREFIX}/lib/pkgconfig" pkg-config --modversion htslib >&2
