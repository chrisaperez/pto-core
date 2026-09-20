#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
#
# Differential check: genomic_toolkit's FRiP against the samtools/bedtools
# pipeline it replaces, on the same BAM and peak file.
#
# The point is the *count*, not the timing. A speedup over a pipeline that
# computes something subtly different is worthless, so this asserts exact
# equality of "fragments overlapping a peak" and fails loudly otherwise. Run it
# before believing any change to peak_set.hpp or fragment_stream.cpp.
#
# usage: validate_against_bedtools.sh <genomic_toolkit> <in.bam> <peaks.bed>
set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "usage: $0 <genomic_toolkit-binary> <in.bam> <peaks.bed>" >&2
    exit 64
fi

BINARY="$1"
BAM="$2"
PEAKS="$3"

for tool in samtools bedtools awk sort; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "skip: $tool is not installed" >&2
        exit 0
    fi
done

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

echo "== reference pipeline (samtools sort -n | bedtools bamtobed | sort | bedtools intersect) =="
REF_START=$(date +%s)
samtools sort -n -o "$WORK/ns.bam" "$BAM"
# bamtobed -bedpe needs name-sorted input; the fragment is the span from the
# leftmost mate start to the rightmost mate end, which is what the tool's
# BEDPE path reconstructs too.
bedtools bamtobed -bedpe -i "$WORK/ns.bam" 2>/dev/null \
    | awk 'BEGIN { OFS = "\t" } $1 == $4 {
          s = ($2 < $5 ? $2 : $5); e = ($3 > $6 ? $3 : $6); print $1, s, e
      }' \
    | sort -k1,1 -k2,2n > "$WORK/frags.bed"
REF_HITS=$(bedtools intersect -u -a "$WORK/frags.bed" -b "$PEAKS" | wc -l | tr -d ' ')
REF_TOTAL=$(wc -l < "$WORK/frags.bed" | tr -d ' ')
REF_END=$(date +%s)
REF_BYTES=$(($(wc -c < "$WORK/ns.bam") + $(wc -c < "$WORK/frags.bed")))

echo "== genomic_toolkit =="
NEW_START=$(date +%s)
# --keep-dups, because the reference pipeline does no duplicate marking. Both
# sides must gate fragment length identically: the reference applies none, so
# --max-length 0 disables this tool's default 1000 bp cap.
JSON=$("$BINARY" frip "$BAM" --peaks "$PEAKS" --keep-dups --max-length 0 --json)
NEW_END=$(date +%s)

NEW_HITS=$(printf '%s' "$JSON" | sed -n 's/.*"in_peaks":\([0-9]*\).*/\1/p')
NEW_TOTAL=$(printf '%s' "$JSON" | sed -n 's/.*"fragments":\([0-9]*\).*/\1/p')

printf '\n%-28s %12s %12s\n' "" "reference" "toolkit"
printf '%-28s %12s %12s\n' "fragments"        "$REF_TOTAL"  "$NEW_TOTAL"
printf '%-28s %12s %12s\n' "fragments in peaks" "$REF_HITS" "$NEW_HITS"
printf '%-28s %12s %12s\n' "seconds"          "$((REF_END - REF_START))" "$((NEW_END - NEW_START))"
printf '%-28s %12s %12s\n' "intermediate bytes" "$REF_BYTES" "0"

STATUS=0
if [ "$REF_TOTAL" != "$NEW_TOTAL" ]; then
    echo "MISMATCH: fragment counts differ ($REF_TOTAL vs $NEW_TOTAL)" >&2
    STATUS=1
fi
if [ "$REF_HITS" != "$NEW_HITS" ]; then
    echo "MISMATCH: in-peak counts differ ($REF_HITS vs $NEW_HITS)" >&2
    STATUS=1
fi
[ "$STATUS" -eq 0 ] && echo "counts agree exactly"
exit "$STATUS"
