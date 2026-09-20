#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
#
# Regression check: FRiP on a BAM must not depend on the ORDER in which the
# peak file lists its contigs.
#
# `frip` pre-seeds the stream's contig dictionary with the peak file's contigs
# before the BAM header is read. `stream_bam` used to write the raw BAM tid into
# each Fragment on the assumption that the dictionary and the header agreed on
# numbering -- true only when the dictionary started empty. On a multi-contig
# BAM whose peak file listed contigs in any other order, every fragment was
# looked up against the wrong contig's peaks: a silently wrong FRiP, exit 0.
# Found on ENCODE K562 ATAC-seq (ENCFF121ZQX, chr22 slice with its full
# 195-contig header): 0 fragments in peaks where bedtools counts 22,473.
#
# Every in-tree BAM fixture has ONE contig, where the two numberings coincide,
# which is why nothing caught it. This fixture has two, and a peak file that
# lists them backwards.
#
# usage: validate_bam_contig_order.sh <genomic_toolkit>
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

# chr1: 10 fragments at 20000..29000, none in a peak.
# chr2: 10 fragments at 1000..10000, the first 4 inside chr2's peak.
# By construction: fragments=20, in_peaks=4.
#
# The two contigs' fragments must sit at DIFFERENT positions. With identical
# layouts a swapped lookup sends chr1's fragments into chr2's peak and scores
# the same 4 -- the first version of this script did exactly that, and passed
# against the pre-fix binary.
python3 - "$WORK" <<'PY'
import sys
work = sys.argv[1]
lines = ["@HD\tVN:1.6\tSO:coordinate",
         "@SQ\tSN:chr1\tLN:100000", "@SQ\tSN:chr2\tLN:100000"]
for chrom, first in (("chr1", 20000), ("chr2", 1000)):
    for n, st in enumerate(range(first, first + 10000, 1000)):
        tlen, mate = 300, st + 200
        q = f"{chrom}_{n}"
        lines.append(f"{q}\t99\t{chrom}\t{st+1}\t60\t100M\t=\t{mate+1}\t{tlen}\t{'A'*100}\t{'I'*100}")
        lines.append(f"{q}\t147\t{chrom}\t{mate+1}\t60\t100M\t=\t{st+1}\t{-tlen}\t{'A'*100}\t{'I'*100}")
open(f"{work}/in.sam", "w").write("\n".join(lines) + "\n")
# chr2 covers fragments starting 1000..4000; chr1's peak touches nothing.
peaks = {"chr1": "chr1\t90000\t91000\n", "chr2": "chr2\t900\t4500\n"}
open(f"{work}/header_order.bed", "w").write(peaks["chr1"] + peaks["chr2"])
open(f"{work}/reversed.bed", "w").write(peaks["chr2"] + peaks["chr1"])
PY

samtools sort -o "$WORK/in.bam" "$WORK/in.sam" 2>/dev/null
samtools index "$WORK/in.bam"

status=0
for PEAKS in header_order reversed; do
    # `|| true`: a pre-fix binary exits 2 here, and the count is what we report.
    OUT=$("$BINARY" frip "$WORK/in.bam" --peaks "$WORK/$PEAKS.bed" --keep-dups --json 2>/dev/null || true)
    GOT=$(printf '%s' "$OUT" | python3 -c "import json,sys
try: d=json.load(sys.stdin); print(d['fragments'], d['in_peaks'])
except Exception: print('? ?')")
    if [ "$GOT" = "20 4" ]; then
        printf '   OK   peaks in %-12s fragments/in_peaks = %s\n' "$PEAKS" "$GOT"
    else
        printf '   FAIL peaks in %-12s fragments/in_peaks = %s (want 20 4)\n' "$PEAKS" "$GOT" >&2
        status=1
    fi
done
exit $status
