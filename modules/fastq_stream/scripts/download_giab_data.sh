#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
#
# Fetch the GIAB HG002 (NA24385) Illumina NovaSeq 2x150bp PCR-free WGS dataset
# used for the benchmarks in docs/PAPER_DRAFT.tex.
#
# The full pair is ~100 GB compressed. Use --subset to pull a bounded prefix
# instead, which is enough to reproduce the throughput ranking without the
# storage cost.
#
# NOTE ON URLS: public genomics mirrors reorganise their layouts periodically.
# This script verifies every URL with a HEAD request *before* downloading and
# fails loudly with the offending URL rather than silently producing a
# truncated or empty file. If a URL has moved, check the GIAB index at
#   https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/
# and update SOURCES below.

set -euo pipefail

DEST="${DEST:-data/giab}"
SUBSET_READS=0
JOBS=4
VERIFY_ONLY=0

usage() {
  cat <<'EOF'
usage: download_giab_data.sh [options]

  --dest DIR         download directory [data/giab]
  --subset N         keep only the first N reads per mate (0 = full file)
  --jobs N           parallel connections per file [4]
  --verify-only      check URL availability and exit
  -h, --help         this message

Environment:
  DEST               same as --dest
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dest) DEST="$2"; shift 2 ;;
    --subset) SUBSET_READS="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    --verify-only) VERIFY_ONLY=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

# HG002 NovaSeq 6000, PCR-free, 2x150bp. Hosted by Google Brain Genomics as
# part of the DeepVariant training data release; this is the same library GIAB
# distributes and is the one cited in the manuscript.
BASE="https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/pcr_free/30x"
SOURCES=(
  "${BASE}/HG002.novaseq.pcr-free.30x.R1.fastq.gz"
  "${BASE}/HG002.novaseq.pcr-free.30x.R2.fastq.gz"
)

have() { command -v "$1" >/dev/null 2>&1; }

if ! have curl; then
  echo "error: curl is required" >&2
  exit 1
fi

echo "==> verifying source URLs"
missing=0
for url in "${SOURCES[@]}"; do
  # -L follows redirects; -f makes HTTP errors a non-zero exit.
  if size=$(curl -sIL -f "$url" | awk 'BEGIN{IGNORECASE=1} /^content-length:/ {v=$2} END{print v+0}'); then
    if [[ "${size}" -gt 0 ]]; then
      printf '    ok   %12d bytes  %s\n' "$size" "$url"
    else
      printf '    WARN (no content-length) %s\n' "$url"
    fi
  else
    printf '    FAIL %s\n' "$url" >&2
    missing=1
  fi
done

if [[ "$missing" -ne 0 ]]; then
  cat >&2 <<EOF

error: one or more source URLs are unreachable.

  These are public mirrors and their paths do change. Do not proceed with a
  partial download -- benchmark numbers from a truncated FASTQ are worthless.
  Locate the current path for the HG002 NovaSeq PCR-free 30x FASTQs and update
  SOURCES in this script.
EOF
  exit 1
fi

[[ "$VERIFY_ONLY" -eq 1 ]] && { echo "==> verify-only, done"; exit 0; }

mkdir -p "$DEST"

fetch() {
  local url="$1" out="$2"
  if [[ -f "$out" ]]; then
    echo "    exists, skipping: $out"
    return
  fi
  echo "    fetching $(basename "$out")"
  if have aria2c; then
    aria2c -x "$JOBS" -s "$JOBS" -d "$(dirname "$out")" -o "$(basename "$out")" "$url"
  else
    # --fail so an HTML error page never lands on disk as a .gz
    curl -fL --retry 3 --retry-delay 5 -o "${out}.part" "$url"
    mv "${out}.part" "$out"
  fi
}

echo "==> downloading to ${DEST}"
for url in "${SOURCES[@]}"; do
  fetch "$url" "${DEST}/$(basename "$url")"
done

if [[ "$SUBSET_READS" -gt 0 ]]; then
  echo "==> building ${SUBSET_READS}-read subsets"
  lines=$(( SUBSET_READS * 4 ))
  for url in "${SOURCES[@]}"; do
    f="${DEST}/$(basename "$url")"
    sub="${f%.fastq.gz}.subset.fastq.gz"
    [[ -f "$sub" ]] && { echo "    exists, skipping: $sub"; continue; }
    echo "    $(basename "$sub")"
    # head closing the pipe early is expected; tolerate the resulting SIGPIPE.
    ( gzip -dc "$f" || true ) | head -n "$lines" | gzip -c > "$sub"
  done
fi

# A BGZF copy is required to exercise the block-parallel path. Plain gzip from
# the sequencer is single-member and cannot be inflated in parallel by any
# tool -- see docs/ARCHITECTURE.md section 3.
if have bgzip; then
  echo "==> building BGZF copies (enables the block-parallel path)"
  for f in "${DEST}"/*.fastq.gz; do
    [[ "$f" == *.bgzf.gz ]] && continue
    out="${f%.fastq.gz}.bgzf.gz"
    [[ -f "$out" ]] && continue
    echo "    $(basename "$out")"
    gzip -dc "$f" | bgzip -@ "$JOBS" -c > "$out"
  done
else
  echo "==> bgzip not found; skipping BGZF copies"
  echo "    install htslib to benchmark the block-parallel path"
fi

echo "==> done"
du -sh "$DEST" 2>/dev/null || true
