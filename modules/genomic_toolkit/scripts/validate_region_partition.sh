#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
#
# Differential check: `--region` must PARTITION a file, not merely filter it.
#
# This is the property a map-reduce over regions depends on, and it is not the
# property htslib's iterator has. `sam_itr_querys` returns every record
# OVERLAPPING a region, so a fragment straddling a tile boundary is returned by
# the queries on both sides. Summing such shards double-counts, and the result
# is a plausible number that is simply too big -- no error, no warning.
#
# `stream_bam` therefore assigns a fragment to the region containing its START.
# This script asserts the consequence: for a set of tiles covering a contig,
# the per-tile counts sum to exactly the whole-contig counts. It deliberately
# builds a fixture whose fragments straddle a boundary, because a random BAM
# almost never produces one and a test that cannot fail is not a test -- the
# samtools cross-check below prints the over-count this rule avoids.
#
# usage: validate_region_partition.sh <genomic_toolkit>
set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "usage: $0 <genomic_toolkit-binary>" >&2
    exit 64
fi
BINARY="$1"

for tool in samtools python3; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "skip: $tool is not installed" >&2
        exit 0
    fi
done

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

BOUNDARY=50000
python3 - "$WORK" "$BOUNDARY" <<'PY'
import sys
work, boundary = sys.argv[1], int(sys.argv[2])
lines = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:100000"]
# Fragments deliberately straddling `boundary`, plus a spread of ordinary ones.
starts = [boundary - 100, boundary - 50, boundary - 20, boundary - 1, boundary, boundary + 50]
starts += list(range(1000, 40000, 2500)) + list(range(60000, 95000, 2500))
for n, st in enumerate(sorted(starts)):
    tlen = 300
    mate = st + tlen - 100
    lines.append(f"s{n}\t99\tchr1\t{st+1}\t60\t100M\t=\t{mate+1}\t{tlen}\t{'A'*100}\t{'I'*100}")
    lines.append(f"s{n}\t147\tchr1\t{mate+1}\t60\t100M\t=\t{st+1}\t{-tlen}\t{'A'*100}\t{'I'*100}")
open(f"{work}/in.sam", "w").write("\n".join(lines) + "\n")
open(f"{work}/peaks.bed", "w").write("chr1\t0\t100000\n")
PY

samtools sort -o "$WORK/in.bam" "$WORK/in.sam" 2>/dev/null
samtools index "$WORK/in.bam"

field() { python3 -c "import json,sys;print(json.load(sys.stdin)['$1'])"; }

LEFT="chr1:1-${BOUNDARY}"
RIGHT="chr1:$((BOUNDARY + 1))-100000"

echo "== samtools, overlap semantics (what NOT to do) =="
SA=$(samtools view -c -f 64 "$WORK/in.bam" "$LEFT")
SB=$(samtools view -c -f 64 "$WORK/in.bam" "$RIGHT")
STRUE=$(samtools view -c -f 64 "$WORK/in.bam" chr1)
echo "   left=$SA right=$SB sum=$((SA + SB))  true=$STRUE"
if [ "$((SA + SB))" -eq "$STRUE" ]; then
    echo "FAIL: fixture does not straddle the boundary, so this check proves nothing" >&2
    exit 1
fi
echo "   (overlap semantics over-counts by $((SA + SB - STRUE)) -- this is the bug being avoided)"

echo "== genomic_toolkit --region, start ownership =="
status=0
for CMD in "sizes" "frip"; do
    if [ "$CMD" = "frip" ]; then
        EXTRA=(--peaks "$WORK/peaks.bed" --keep-dups)
        KEYS="fragments in_peaks bases_in_peaks bases_total"
    else
        EXTRA=()
        KEYS="fragments records_read"
    fi
    # `${EXTRA[@]+...}`, not a bare `"${EXTRA[@]}"`: macOS ships bash 3.2,
    # where expanding an EMPTY array under `set -u` is an unbound-variable
    # error. CI has a macOS runner, so this is a real portability constraint
    # rather than a hypothetical one.
    set -- ${EXTRA[@]+"${EXTRA[@]}"}

    for KEY in $KEYS; do
        W=$("$BINARY" "$CMD" "$WORK/in.bam" "$@" --json | field "$KEY")
        L=$("$BINARY" "$CMD" "$WORK/in.bam" "$@" --region "$LEFT"  --json | field "$KEY")
        R=$("$BINARY" "$CMD" "$WORK/in.bam" "$@" --region "$RIGHT" --json | field "$KEY")
        if [ "$((L + R))" -eq "$W" ]; then
            printf '   OK   %-6s %-16s %s + %s == %s\n' "$CMD" "$KEY" "$L" "$R" "$W"
        else
            printf '   FAIL %-6s %-16s %s + %s != %s\n' "$CMD" "$KEY" "$L" "$R" "$W" >&2
            status=1
        fi
    done
done

# base_frip must be reconstructible from the two emitted terms, which is the
# entire reason they are emitted.
python3 - "$BINARY" "$WORK" "$LEFT" "$RIGHT" <<'PY' || status=1
import json, subprocess, sys
binary, work, left, right = sys.argv[1:5]
def run(*extra):
    return json.loads(subprocess.run(
        [binary, "frip", f"{work}/in.bam", "--peaks", f"{work}/peaks.bed",
         "--keep-dups", "--json", *extra],
        capture_output=True, text=True, check=True).stdout)

whole = run()
parts = [run("--region", left), run("--region", right)]
merged = sum(p["bases_in_peaks"] for p in parts) / sum(p["bases_total"] for p in parts)
if abs(merged - whole["base_frip"]) > 1e-9:
    print(f"   FAIL base_frip  merged {merged:.9f} != whole {whole['base_frip']:.9f}",
          file=sys.stderr)
    sys.exit(1)
print(f"   OK   frip   base_frip        merged {merged:.6f} == whole")
PY

exit $status
