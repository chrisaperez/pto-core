#!/usr/bin/env bash
#
# Benchmarking protocol for the Application Note: runs cuttag_profiler and
# `deeptools computeMatrix reference-point` over the same BAM and region set,
# records wall time / peak RSS / CPU for each, and quantifies how closely the
# two agree bin-for-bin.
#
# Nothing here fabricates a number: every figure quoted in the manuscript comes
# out of results/benchmark.tsv and results/accuracy.txt produced by this script.
#
#   scripts/benchmark.sh --bam data/demo.bam --regions data/promoters.bed
#
set -euo pipefail

BAM=""
REGIONS=""
OUT_DIR="results"
THREADS="$( (command -v nproc >/dev/null && nproc) || sysctl -n hw.ncpu )"
REPLICATES=3
UPSTREAM=2000
DOWNSTREAM=2000
BIN_SIZE=50
PROFILER="./build/cuttag_profiler"
COMPUTE_MATRIX="computeMatrix"

usage() {
  cat <<'EOF'
usage: scripts/benchmark.sh --bam FILE --regions FILE [options]

  --bam FILE           indexed BAM to profile           (required)
  --regions FILE       BED of loci                      (required)
  --out-dir DIR        where results land               [results]
  --threads N          worker threads for both tools    [all cores]
  --replicates N       timed repeats per tool           [3]
  --upstream BP        window before the anchor         [2000]
  --downstream BP      window after the anchor          [2000]
  --bin-size BP        bin width                        [50]
  --profiler PATH      cuttag_profiler binary           [./build/cuttag_profiler]
  --compute-matrix P   deeptools computeMatrix binary   [computeMatrix]
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2 ;;
    --regions) REGIONS="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --replicates) REPLICATES="$2"; shift 2 ;;
    --upstream) UPSTREAM="$2"; shift 2 ;;
    --downstream) DOWNSTREAM="$2"; shift 2 ;;
    --bin-size) BIN_SIZE="$2"; shift 2 ;;
    --profiler) PROFILER="$2"; shift 2 ;;
    --compute-matrix) COMPUTE_MATRIX="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$BAM" && -n "$REGIONS" ]] || { usage; exit 1; }
[[ -x "$PROFILER" ]] || { echo "no profiler binary at $PROFILER" >&2; exit 1; }

mkdir -p "$OUT_DIR"
RESULTS="$OUT_DIR/benchmark.tsv"

# GNU time reports peak RSS; macOS /usr/bin/time -l does too, in a different
# format. Normalise both to kilobytes.
measure() {
  local label="$1"; shift
  local logfile
  logfile="$(mktemp)"

  if /usr/bin/time -l true >/dev/null 2>&1; then
    /usr/bin/time -l "$@" >/dev/null 2>"$logfile" || true
    local wall rss
    wall="$(awk '/real/ {print $1}' "$logfile" | head -1)"
    rss="$(awk '/maximum resident set size/ {printf "%.0f", $1/1024}' "$logfile")"
  else
    /usr/bin/time -v "$@" >/dev/null 2>"$logfile" || true
    local wall rss
    wall="$(awk -F': ' '/Elapsed \(wall clock\)/ {print $2}' "$logfile" |
            awk -F: '{ if (NF==3) print $1*3600+$2*60+$3; else if (NF==2) print $1*60+$2; else print $1 }')"
    rss="$(awk -F': ' '/Maximum resident set size/ {print $2}' "$logfile")"
  fi

  rm -f "$logfile"
  printf '%s\t%s\t%s\n' "$label" "${wall:-NA}" "${rss:-NA}"
}

echo "region file: $(wc -l < "$REGIONS") loci"
echo "threads:     $THREADS"
echo "replicates:  $REPLICATES"
echo

printf 'tool\treplicate\twall_seconds\tpeak_rss_kb\n' > "$RESULTS"

echo "--- cuttag_profiler ---"
for r in $(seq 1 "$REPLICATES"); do
  read -r _ wall rss < <(measure profiler \
    "$PROFILER" profile \
      --bam "$BAM" --regions "$REGIONS" \
      --upstream "$UPSTREAM" --downstream "$DOWNSTREAM" \
      --bin-size "$BIN_SIZE" --normalization CPM \
      --threads "$THREADS" --quiet \
      --out-matrix "$OUT_DIR/profiler_matrix.tsv")
  printf 'cuttag_profiler\t%s\t%s\t%s\n' "$r" "$wall" "$rss" >> "$RESULTS"
  echo "  replicate $r: ${wall}s, ${rss} kB"
done

# computeMatrix consumes a bigWig, not a BAM, so the comparable unit of work
# is the two-stage pipeline a bench scientist actually runs: bamCoverage to
# build the coverage track, then computeMatrix to window it. Both stages are
# timed, and the sum is what the speedup is quoted against; the computeMatrix
# stage is also reported alone so a reader with a cached bigWig can judge it.
BAM_COVERAGE="${COMPUTE_MATRIX%computeMatrix}bamCoverage"
if command -v "$COMPUTE_MATRIX" >/dev/null 2>&1 && command -v "$BAM_COVERAGE" >/dev/null 2>&1; then
  BIGWIG="$OUT_DIR/$(basename "${BAM%.bam}").bw"

  echo "--- deeptools bamCoverage ---"
  for r in $(seq 1 "$REPLICATES"); do
    read -r _ wall rss < <(measure bamcoverage \
      "$BAM_COVERAGE" \
        -b "$BAM" -o "$BIGWIG" \
        --binSize "$BIN_SIZE" \
        --normalizeUsing CPM --extendReads \
        -p "$THREADS")
    printf 'deeptools_bamCoverage\t%s\t%s\t%s\n' "$r" "$wall" "$rss" >> "$RESULTS"
    echo "  replicate $r: ${wall}s, ${rss} kB"
  done

  echo "--- deeptools computeMatrix ---"
  for r in $(seq 1 "$REPLICATES"); do
    read -r _ wall rss < <(measure computematrix \
      "$COMPUTE_MATRIX" reference-point \
        --referencePoint TSS \
        -S "$BIGWIG" \
        -R "$REGIONS" \
        -a "$DOWNSTREAM" -b "$UPSTREAM" \
        --binSize "$BIN_SIZE" \
        --numberOfProcessors "$THREADS" \
        --outFileName "$OUT_DIR/deeptools_matrix.gz" \
        --outFileNameMatrix "$OUT_DIR/deeptools_matrix.tsv")
    printf 'deeptools_computeMatrix\t%s\t%s\t%s\n' "$r" "$wall" "$rss" >> "$RESULTS"
    echo "  replicate $r: ${wall}s, ${rss} kB"
  done
else
  echo "deeptools not on PATH - skipping the baseline arm." >&2
  echo "Install with: pip install deeptools" >&2
fi

echo
echo "--- summary ---"
awk -F'\t' 'NR>1 && $3 != "NA" { sum[$1] += $3; n[$1]++; if ($4+0 > rss[$1]) rss[$1] = $4 }
     END {
       printf "%-26s %10s %12s\n", "stage", "mean_s", "peak_rss_MB";
       for (t in sum) printf "%-26s %10.2f %12.1f\n", t, sum[t]/n[t], rss[t]/1024;

       ours = ("cuttag_profiler" in sum) ? sum["cuttag_profiler"]/n["cuttag_profiler"] : 0;
       cov  = ("deeptools_bamCoverage" in sum) ? sum["deeptools_bamCoverage"]/n["deeptools_bamCoverage"] : 0;
       cm   = ("deeptools_computeMatrix" in sum) ? sum["deeptools_computeMatrix"]/n["deeptools_computeMatrix"] : 0;

       if (ours > 0 && cm > 0) {
         printf "\n%-26s %10.2f\n", "deeptools pipeline total", cov + cm;
         printf "%-26s %9.1fx\n", "speedup (end to end)", (cov + cm) / ours;
         printf "%-26s %9.1fx\n", "speedup (computeMatrix)", cm / ours;
       }
     }' "$RESULTS" | tee "$OUT_DIR/summary.txt"

# --- Accuracy -------------------------------------------------------------
# Compares the two meta-profiles bin by bin. The tools normalise differently,
# so the comparison is on the Pearson correlation of the curves plus the
# largest relative deviation after rescaling to a common mean.
if [[ -s "$OUT_DIR/deeptools_matrix.tsv" && -s "$OUT_DIR/profiler_matrix.tsv" ]]; then
  python3 scripts/compare_matrices.py \
    --profiler "$OUT_DIR/profiler_matrix.tsv" \
    --deeptools "$OUT_DIR/deeptools_matrix.tsv" \
    | tee "$OUT_DIR/accuracy.txt"
fi

echo
echo "results written to $RESULTS"
