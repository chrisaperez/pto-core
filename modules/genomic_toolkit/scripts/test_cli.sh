#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
#
# Black-box checks of the genomic_toolkit CLI.
#
# The unit suites link the library, not src/main.cpp, so everything that lives
# there -- argument parsing, the exit codes pto-cloud and Nextflow key on, and
# the refusals `frip` and `markdup` make after streaming -- is reachable only by
# running the binary. Every case below was a defect on 2026-09-11 (abort,
# silent narrowing, silent undercount, silent partial BAM); the audit record is
# docs/AUDIT_2026-09-11_genomic_toolkit.md.
#
# Text inputs always run. The BAM cases need a build with htslib, plus samtools
# and python3 to write the fixture; without them they SKIP, loudly.
#
# usage: test_cli.sh <genomic_toolkit>
set -uo pipefail  # not -e: exit codes are the thing under test

if [ "$#" -ne 1 ]; then
    echo "usage: $0 <genomic_toolkit-binary>" >&2
    exit 64
fi
BIN="$1"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

checks=0
failures=0
ok()   { checks=$((checks + 1)); printf '   OK   %s\n' "$1"; }
fail() { checks=$((checks + 1)); failures=$((failures + 1)); printf '   FAIL %s\n' "$1" >&2; }

# expect CODE DESCRIPTION command...   (stdout+stderr land in $WORK/out)
expect() {
    local want="$1" desc="$2"
    shift 2
    "$@" >"$WORK/out" 2>&1
    local got=$?
    if [ "$got" -eq "$want" ]; then
        ok "$desc (exit $got)"
    else
        fail "$desc: exit $got, wanted $want"
        sed 's/^/          /' "$WORK/out" >&2
    fi
}

# says NEEDLE DESCRIPTION   -- about the most recent `expect`
says() {
    if grep -qF -- "$1" "$WORK/out"; then
        ok "$2"
    else
        fail "$2: output lacks '$1'"
    fi
}

F="$WORK/f.bed"
P="$WORK/p.bed"
printf 'chr1\t100\t300\nchr1\t200\t260\nchr1\t500\t900\n' >"$F"
printf 'chr1\t0\t1000\n' >"$P"

echo "== numeric options: refused by name, never aborted, never narrowed"
for bad in abc 100x 4294967396 99999999999999999999 -50 0x10 ''; do
    expect 1 "--min-length '$bad' is refused" "$BIN" sizes "$F" --min-length "$bad" --json
    says "expects an integer" "--min-length '$bad' names the problem"
done
expect 1 "--max-length 4294967296 is refused (was narrowed to 0, 'off')" \
    "$BIN" sizes "$F" --max-length 4294967296 --json
expect 1 "--min-mapq 300 is refused" "$BIN" sizes "$F" --min-mapq 300 --json
expect 1 "--min-mapq -1 is refused" "$BIN" sizes "$F" --min-mapq -1 --json
expect 1 "--genome-size 0x10 is refused (was read as 0)" \
    "$BIN" frip "$F" --peaks "$P" --genome-size 0x10 --json
expect 1 "--genome-size -1 is refused" "$BIN" frip "$F" --peaks "$P" --genome-size -1 --json
expect 0 "--min-length at the int32 ceiling is accepted" \
    "$BIN" sizes "$F" --min-length 2147483647 --max-length 0 --json
expect 0 "--min-mapq 255 is accepted" "$BIN" sizes "$F" --min-mapq 255 --json
expect 0 "--genome-size 3100000000 is accepted" \
    "$BIN" frip "$F" --peaks "$P" --genome-size 3100000000 --json

echo "== frip refuses an unsorted input when it is marking duplicates"
U="$WORK/unsorted.bed"
printf 'chr1\t100\t300\nchr1\t500\t900\nchr1\t100\t300\n' >"$U"
expect 2 "frip on unsorted input" "$BIN" frip "$U" --peaks "$P" --json
says "out of coordinate order" "frip says why"
expect 0 "frip --unsorted on the same input" "$BIN" frip "$U" --peaks "$P" --unsorted --json
says '"duplicates_excluded":1' "frip --unsorted finds the duplicate"
expect 0 "frip --keep-dups on the same input" "$BIN" frip "$U" --peaks "$P" --keep-dups --json
expect 2 "markdup on unsorted input (unchanged)" "$BIN" markdup "$U" --json

echo "== an input with no parseable record is an error, not an empty library"
BAD="$WORK/bad.bed"
printf 'chr1\t300\t100\nchr1\tabc\t200\n' >"$BAD"
expect 2 "sizes over nothing but malformed records" "$BIN" sizes "$BAD" --json
says "none of the 2 records" "sizes says how many"
expect 2 "markdup over nothing but malformed records" "$BIN" markdup "$BAD" --json
expect 2 "frip over nothing but malformed records" "$BIN" frip "$BAD" --peaks "$P" --json
MIXED="$WORK/mixed.bed"
printf 'chr1\t100\t300\nchr1\tabc\t200\n' >"$MIXED"
expect 0 "a tolerated handful of malformed lines still runs" "$BIN" sizes "$MIXED" --json
: >"$WORK/empty.bed"
expect 0 "an empty file is not 'all malformed'" "$BIN" sizes "$WORK/empty.bed" --json

echo "== column 5: a score by default, a count only with --with-counts"
B6="$WORK/scored.bed"
printf 'chr1\t100\t300\tr1\t60\t+\nchr1\t400\t600\tr2\t.\t-\n' >"$B6"
expect 0 "a BED6 with a score column" "$BIN" sizes "$B6" --json
says '"fragments":2' "a BED6 is one fragment per interval, not score-many"
expect 0 "frip over the same BED6" "$BIN" frip "$B6" --peaks "$P" --keep-dups --json
says '"fragments":2' "frip does not multiply by the score either"
C="$WORK/copies.tsv"
printf 'chr1\t100\t300\tBC\t5\nchr1\t100\t300\tBC\tabc\nchr1\t100\t300\tBC\t0\n' >"$C"
expect 0 "--with-counts over good and corrupt counts" "$BIN" sizes "$C" --with-counts --json
says '"fragments":5' "--with-counts counts only the good line's 5 copies"
says '"malformed":2' "--with-counts reports the corrupt counts as malformed"
expect 0 "the same file without --with-counts" "$BIN" sizes "$C" --json
says '"fragments":3' "without --with-counts each line is one fragment"
printf 'chr1\t100\t300\tBC\t99999999999999\n' >"$WORK/copies_big.tsv"
expect 2 "--with-counts refuses an absurd count, not clamped to a million" \
    "$BIN" sizes "$WORK/copies_big.tsv" --with-counts --json
expect 0 "without --with-counts the same line is one fragment" \
    "$BIN" sizes "$WORK/copies_big.tsv" --json

echo "== peak files"
printf 'chr1\t150\t150\nchr1\t400\t400\n' >"$WORK/zero_peaks.bed"
expect 1 "a peak file of only zero-length intervals" "$BIN" frip "$F" --peaks "$WORK/zero_peaks.bed" --json
says "zero-length" "the error names zero-length intervals, not contig naming"
printf 'chr1\t0\t1000\r\n' >"$WORK/crlf_peaks.bed"
expect 0 "a CRLF peak file" "$BIN" frip "$F" --peaks "$WORK/crlf_peaks.bed" --keep-dups --json
says '"in_peaks":3' "the CRLF peak file gives the LF answer"

echo "== BAM"
if "$BIN" --help | grep -q "unavailable: built without htslib"; then
    printf '  SKIP BAM cases (this build has no htslib)\n'
elif ! command -v samtools >/dev/null 2>&1 || ! command -v python3 >/dev/null 2>&1; then
    printf '  SKIP BAM cases (samtools and python3 are needed to write the fixture)\n'
else
    # Enough pairs that the BAM spans many BGZF blocks, so there is a block
    # boundary to cut at that is not the end of the file.
    python3 - "$WORK/in.sam" <<'PY'
import sys
lines = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:50000000"]
seq, qual = "A" * 100, "I" * 100
for n in range(40000):
    st = 1000 + n * 300
    mate = st + 200
    lines.append(f"p{n}\t99\tchr1\t{st+1}\t60\t100M\t=\t{mate+1}\t300\t{seq}\t{qual}")
    lines.append(f"p{n}\t147\tchr1\t{mate+1}\t60\t100M\t=\t{st+1}\t-300\t{seq}\t{qual}")
open(sys.argv[1], "w").write("\n".join(lines) + "\n")
PY
    samtools view -b -o "$WORK/in.bam" "$WORK/in.sam" 2>/dev/null
    python3 - "$WORK/in.bam" "$WORK" <<'PY'
import struct, sys
src, work = sys.argv[1], sys.argv[2]
data = open(src, "rb").read()
offs, pos = [], 0
while pos < len(data):
    assert data[pos:pos + 2] == b"\x1f\x8b", f"not BGZF at {pos}"
    offs.append(pos)
    pos += struct.unpack_from("<H", data, pos + 16)[0] + 1
assert len(offs) > 4, f"fixture has only {len(offs)} BGZF blocks"
mid = offs[len(offs) // 2]
open(f"{work}/cut_at_block.bam", "wb").write(data[:mid])
open(f"{work}/no_eof.bam", "wb").write(data[:-28])
open(f"{work}/cut_mid_block.bam", "wb").write(data[:mid + 1000])
PY
    expect 0 "the intact BAM" "$BIN" sizes "$WORK/in.bam" --json
    says '"fragments":40000' "the intact BAM yields every pair"
    expect 1 "a BAM cut at a BGZF block boundary is refused" "$BIN" sizes "$WORK/cut_at_block.bam" --json
    says "end-of-file marker" "the error says the file is truncated"
    expect 1 "frip refuses the same truncated BAM" \
        "$BIN" frip "$WORK/cut_at_block.bam" --peaks "$P" --json
    expect 1 "a BAM missing only its EOF block is refused" "$BIN" sizes "$WORK/no_eof.bam" --json
    expect 1 "a BAM cut mid-block is refused" "$BIN" sizes "$WORK/cut_mid_block.bam" --json
fi

printf '%d checks, %d failure(s)\n' "$checks" "$failures"
[ "$failures" -eq 0 ]
