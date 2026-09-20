#!/usr/bin/env bash
#
# Builds a static htslib with every network transport compiled out.
#
# Why this exists: a stock htslib archive contains hfile_libcurl.o and
# hfile_s3.o, so the linked binary carries an HTTP/S3 client even though the
# profiler never constructs a remote URL. For a deployment into a regulated
# genomics environment, "the code path is unreachable" is a weaker claim than
# "the code is not in the binary". This build removes it, which also drops the
# libcurl and OpenSSL link dependencies entirely.
#
# Afterwards, point the project's configure step at the result:
#
#   scripts/build_htslib_minimal.sh --prefix "$PWD/third_party/htslib-install"
#   cmake -S . -B build \
#     -DCMAKE_PREFIX_PATH="$PWD/third_party/htslib-install" \
#     -DPKG_CONFIG_PATH="$PWD/third_party/htslib-install/lib/pkgconfig"
#
set -euo pipefail

VERSION="1.22"
PREFIX="$PWD/third_party/htslib-install"
BUILD_DIR="$PWD/third_party/htslib-build"
TARBALL=""
JOBS="$( (command -v nproc >/dev/null && nproc) || sysctl -n hw.ncpu )"

usage() {
  cat <<'EOF'
usage: scripts/build_htslib_minimal.sh [options]

  --version V     htslib release to fetch            [1.22]
  --prefix DIR    install prefix                     [./third_party/htslib-install]
  --build-dir DIR scratch build directory            [./third_party/htslib-build]
  --tarball FILE  use a local htslib-X.Y.tar.bz2 instead of downloading
                  (required on a machine with no outbound network)
  --jobs N        parallel make jobs                 [all cores]
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --version) VERSION="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --build-dir) BUILD_DIR="$2"; shift 2 ;;
    --tarball) TARBALL="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

if [[ -n "$TARBALL" ]]; then
  echo "==> unpacking $TARBALL"
  tar xf "$TARBALL"
else
  URL="https://github.com/samtools/htslib/releases/download/${VERSION}/htslib-${VERSION}.tar.bz2"
  echo "==> downloading $URL"
  echo "    (on an air-gapped host, fetch this elsewhere and pass --tarball)"
  curl -fsSL "$URL" -o "htslib-${VERSION}.tar.bz2"
  tar xf "htslib-${VERSION}.tar.bz2"
fi

cd "htslib-${VERSION}"

echo "==> configuring without any network transport"
# --disable-libcurl removes the HTTPS/S3/GCS transports outright, which is the
# point of this script. --disable-plugins keeps everything in the archive so
# there is no runtime .so lookup. --disable-lzma/--disable-bz2 could shrink it
# further but would refuse legitimate CRAM inputs, so both stay on.
./configure \
  --prefix="$PREFIX" \
  --disable-libcurl \
  --disable-s3 \
  --disable-gcs \
  --disable-plugins \
  --enable-libdeflate \
  --without-libdeflate 2>/dev/null || \
./configure \
  --prefix="$PREFIX" \
  --disable-libcurl \
  --disable-s3 \
  --disable-gcs \
  --disable-plugins

make -j"$JOBS"
make install

echo
echo "==> verifying that no network symbols survived"
ARCHIVE="$PREFIX/lib/libhts.a"
if command -v nm >/dev/null 2>&1; then
  if nm -uo "$ARCHIVE" 2>/dev/null | grep -q "curl_easy"; then
    echo "WARNING: libcurl symbols are still referenced by $ARCHIVE" >&2
    exit 1
  fi
  echo "    no curl_* symbols referenced by $ARCHIVE"
fi

echo
echo "htslib ${VERSION} installed to $PREFIX"
echo "reconfigure with:"
echo "  cmake -S . -B build -DCMAKE_PREFIX_PATH=$PREFIX \\"
echo "        -DPKG_CONFIG_PATH=$PREFIX/lib/pkgconfig"
