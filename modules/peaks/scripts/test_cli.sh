#!/bin/sh
# SPDX-License-Identifier: MIT
#
# Drives the pto-peaks BINARY as a black box.
#
# The unit suites cover the window and the caller as headers. This covers what
# only the executable has: argument parsing, the BED reader's fail-closed
# checks, the narrowPeak writer, and the exit codes a caller's retry logic
# depends on -- 2 means "the input's shape made the answer wrong", is
# deterministic, and is never retried.
#
# POSIX sh and coreutils only. No python, no bedtools: this must run on the
# same bare cluster the binary is built for.
set -eu

BIN="${1:?usage: test_cli.sh /path/to/pto-peaks [/path/to/test_peaks_bam_streamer]}"
# The BAM fixture writer. Optional: a build without BGZF support has no BAM to
# make, and the BAM cases below skip loudly rather than failing.
BAMGEN="${2:-}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT INT TERM

fail() {
    printf '  FAIL %s\n' "$1" >&2
    FAILURES=$((FAILURES + 1))
}
ok() { printf '  ok   %s\n' "$1"; }
FAILURES=0

expect_status() {
    want="$1"
    name="$2"
    shift 2
    set +e
    "$@" >"$WORK/out" 2>"$WORK/err"
    got=$?
    set -e
    if [ "$got" -eq "$want" ]; then
        ok "$name"
    else
        fail "$name: exit $got, wanted $want"
        sed 's/^/        /' "$WORK/err" | head -3 >&2
    fi
}

printf 'chr1\t100000\nchr2\t50000\n' >"$WORK/sizes.txt"

# ---------------------------------------------------------------------------
# A background tiling with one dense peak planted at 50000.
# ---------------------------------------------------------------------------
awk 'BEGIN {
        for (s = 0; s + 200 <= 99800; s += 200) printf "chr1\t%d\t%d\n", s, s + 200
     }' >"$WORK/bg.bed"
awk 'BEGIN { for (i = 0; i < 300; i++) printf "chr1\t50000\t50200\n" }' >"$WORK/peak.bed"
sort -k1,1 -k2,2n "$WORK/bg.bed" "$WORK/peak.bed" >"$WORK/frags.bed"

printf '== peaks/cli ==\n'

expect_status 0 "--help exits 0" "$BIN" --help
expect_status 0 "--version exits 0" "$BIN" --version
expect_status 1 "--chrom-sizes is required" "$BIN" "$WORK/frags.bed"

# ---------------------------------------------------------------------------
# --chrom-sizes itself is malformed: each defect names the field and value
# rather than one generic "malformed contig record" for all four causes.
# ---------------------------------------------------------------------------
printf '\t100000\n' >"$WORK/sizes_empty_name.txt"
expect_status 1 "an empty contig name in --chrom-sizes exits 1" \
    "$BIN" --chrom-sizes "$WORK/sizes_empty_name.txt" "$WORK/frags.bed"
if grep -q "empty contig name (field 1)" "$WORK/err"; then
    ok "  ...and names the field"
else
    fail "  ...and names the field"
    sed 's/^/        /' "$WORK/err" | head -3 >&2
fi

printf 'chr1\tXX\n' >"$WORK/sizes_nonint.txt"
expect_status 1 "a non-integer length in --chrom-sizes exits 1" \
    "$BIN" --chrom-sizes "$WORK/sizes_nonint.txt" "$WORK/frags.bed"
if grep -q "length field 'XX' is not an integer" "$WORK/err"; then
    ok "  ...and quotes the bad value"
else
    fail "  ...and quotes the bad value"
    sed 's/^/        /' "$WORK/err" | head -3 >&2
fi

printf 'chr1\t0\n' >"$WORK/sizes_zero.txt"
expect_status 1 "a zero length in --chrom-sizes exits 1" \
    "$BIN" --chrom-sizes "$WORK/sizes_zero.txt" "$WORK/frags.bed"
if grep -q "length 0 must be positive" "$WORK/err"; then
    ok "  ...and shows the value"
else
    fail "  ...and shows the value"
    sed 's/^/        /' "$WORK/err" | head -3 >&2
fi

printf 'chr1\t-500\n' >"$WORK/sizes_negative.txt"
expect_status 1 "a negative length in --chrom-sizes exits 1" \
    "$BIN" --chrom-sizes "$WORK/sizes_negative.txt" "$WORK/frags.bed"
if grep -q "length -500 must be positive" "$WORK/err"; then
    ok "  ...and shows the value"
else
    fail "  ...and shows the value"
    sed 's/^/        /' "$WORK/err" | head -3 >&2
fi

printf 'chr1\t99999999999999\n' >"$WORK/sizes_huge.txt"
expect_status 1 "a length past int32 max in --chrom-sizes exits 1" \
    "$BIN" --chrom-sizes "$WORK/sizes_huge.txt" "$WORK/frags.bed"
if grep -q "exceeds the maximum representable contig length" "$WORK/err"; then
    ok "  ...and says so rather than silently truncating"
else
    fail "  ...and says so rather than silently truncating"
    sed 's/^/        /' "$WORK/err" | head -3 >&2
fi
expect_status 1 "unknown option is rejected" "$BIN" --chrom-sizes "$WORK/sizes.txt" --nope

# ---------------------------------------------------------------------------
# The planted peak.
# ---------------------------------------------------------------------------
if "$BIN" --chrom-sizes "$WORK/sizes.txt" --cutoff 5 --min-length 50 \
        --max-gap 100 "$WORK/frags.bed" >"$WORK/peaks.narrowPeak" 2>"$WORK/err"; then
    ok "calls peaks on planted data"
else
    fail "calls peaks on planted data"
fi

n=$(wc -l <"$WORK/peaks.narrowPeak" | tr -d ' ')
if [ "$n" -eq 1 ]; then
    ok "one peak from one planted site"
else
    fail "one peak from one planted site: got $n"
    cat "$WORK/peaks.narrowPeak" >&2
fi

cols=$(head -1 "$WORK/peaks.narrowPeak" | awk -F'\t' '{print NF}')
if [ "$cols" = "10" ]; then
    ok "narrowPeak has 10 columns"
else
    fail "narrowPeak has 10 columns: got $cols"
fi

start=$(head -1 "$WORK/peaks.narrowPeak" | cut -f2)
end=$(head -1 "$WORK/peaks.narrowPeak" | cut -f3)
summit=$(head -1 "$WORK/peaks.narrowPeak" | cut -f10)
if [ "$start" = "50000" ] && [ "$end" = "50200" ]; then
    ok "peak boundaries match the planted site"
else
    fail "peak boundaries match the planted site: got $start-$end"
fi
# The pileup is flat across the plateau, so the summit is its midpoint.
if [ "$summit" -ge 80 ] && [ "$summit" -le 120 ]; then
    ok "summit is in the middle of the plateau"
else
    fail "summit is in the middle of the plateau: offset $summit"
fi

# ---------------------------------------------------------------------------
# The same data on stdin.
# ---------------------------------------------------------------------------
if "$BIN" --chrom-sizes "$WORK/sizes.txt" --cutoff 5 --min-length 50 --max-gap 100 \
        <"$WORK/frags.bed" >"$WORK/stdin.narrowPeak" 2>/dev/null &&
        cmp -s "$WORK/stdin.narrowPeak" "$WORK/peaks.narrowPeak"; then
    ok "stdin gives the same answer as a file"
else
    fail "stdin gives the same answer as a file"
fi

# ---------------------------------------------------------------------------
# Fail-closed input handling. Every one of these must be exit 2: the answer
# would be wrong, and it would be wrong the same way on a retry.
# ---------------------------------------------------------------------------
printf 'chr1\t2000\t2100\nchr1\t1000\t1100\n' >"$WORK/unsorted.bed"
expect_status 2 "unsorted records exit 2" \
    "$BIN" --chrom-sizes "$WORK/sizes.txt" "$WORK/unsorted.bed"

printf 'chrZ\t100\t200\n' >"$WORK/badchrom.bed"
expect_status 2 "a contig missing from the sizes file exits 2" \
    "$BIN" --chrom-sizes "$WORK/sizes.txt" "$WORK/badchrom.bed"

printf 'chr1\t99900\t100100\n' >"$WORK/oob.bed"
expect_status 2 "an interval past the contig end exits 2" \
    "$BIN" --chrom-sizes "$WORK/sizes.txt" "$WORK/oob.bed"
if grep -q "end coordinate 100100 exceeds contig 'chr1' length 100000" "$WORK/err"; then
    ok "  ...and names the field, the value, and the contig length"
else
    fail "  ...and names the field, the value, and the contig length"
    sed 's/^/        /' "$WORK/err" | head -3 >&2
fi

# Negative start and an inverted interval are different defects from "past the
# contig end" -- upstream tools that emit them are broken in different ways,
# so each gets its own message rather than one "interval outside the contig"
# line for all three (this is the Step 3 hardening pass: descriptive,
# field-specific diagnostics, not a generic rejection).
printf 'chr1\t-100\t200\n' >"$WORK/negative.bed"
expect_status 2 "a negative start coordinate exits 2" \
    "$BIN" --chrom-sizes "$WORK/sizes.txt" "$WORK/negative.bed"
if grep -q "negative start coordinate -100 (field 2)" "$WORK/err"; then
    ok "  ...and names the field and the value"
else
    fail "  ...and names the field and the value"
    sed 's/^/        /' "$WORK/err" | head -3 >&2
fi

printf 'chr1\t500\t200\n' >"$WORK/inverted.bed"
expect_status 2 "an inverted interval (start > end) exits 2" \
    "$BIN" --chrom-sizes "$WORK/sizes.txt" "$WORK/inverted.bed"
if grep -q "inverted or empty interval \[500, 200)" "$WORK/err"; then
    ok "  ...and shows both endpoints"
else
    fail "  ...and shows both endpoints"
    sed 's/^/        /' "$WORK/err" | head -3 >&2
fi

printf 'chr1\t100\t100\n' >"$WORK/empty_interval.bed"
expect_status 2 "a zero-length interval (start == end) exits 2" \
    "$BIN" --chrom-sizes "$WORK/sizes.txt" "$WORK/empty_interval.bed"

printf 'chr1\t100\tXX\n' >"$WORK/junk.bed"
expect_status 2 "a malformed coordinate exits 2" \
    "$BIN" --chrom-sizes "$WORK/sizes.txt" "$WORK/junk.bed"

# from_chars writes the parsed prefix before reporting trailing junk, so "100x"
# must be refused rather than silently becoming 100.
printf 'chr1\t100\t200x\n' >"$WORK/trailing.bed"
expect_status 2 "a coordinate with trailing junk exits 2" \
    "$BIN" --chrom-sizes "$WORK/sizes.txt" "$WORK/trailing.bed"

printf 'chr1\t100\t200\nchr2\t100\t200\nchr1\t300\t400\n' >"$WORK/interleaved.bed"
expect_status 2 "a contig reappearing later exits 2" \
    "$BIN" --chrom-sizes "$WORK/sizes.txt" "$WORK/interleaved.bed"

# ---------------------------------------------------------------------------
# Degenerate but legal input.
# ---------------------------------------------------------------------------
: >"$WORK/empty.bed"
if "$BIN" --chrom-sizes "$WORK/sizes.txt" "$WORK/empty.bed" >"$WORK/empty.out" 2>/dev/null; then
    if [ ! -s "$WORK/empty.out" ]; then
        ok "empty input yields no peaks and exits 0"
    else
        fail "empty input yields no peaks and exits 0: output not empty"
    fi
else
    fail "empty input yields no peaks and exits 0"
fi

# Comments and a track line are skipped, not parsed.
{ printf '# a comment\n'; printf 'track name=x\n'; cat "$WORK/frags.bed"; } >"$WORK/commented.bed"
if "$BIN" --chrom-sizes "$WORK/sizes.txt" --cutoff 5 --min-length 50 --max-gap 100 \
        "$WORK/commented.bed" >"$WORK/commented.out" 2>/dev/null &&
        cmp -s "$WORK/commented.out" "$WORK/peaks.narrowPeak"; then
    ok "comment and track lines are skipped"
else
    fail "comment and track lines are skipped"
fi

# min-length discards the call entirely.
if "$BIN" --chrom-sizes "$WORK/sizes.txt" --cutoff 5 --min-length 5000 --max-gap 100 \
        "$WORK/frags.bed" >"$WORK/long.out" 2>/dev/null && [ ! -s "$WORK/long.out" ]; then
    ok "--min-length discards a short call"
else
    fail "--min-length discards a short call"
fi

# --extend-peaks pads the REPORTED interval symmetrically and nothing else:
# the planted peak is [50000, 50200) with its summit near the middle, so a
# 40 bp pad must give [49960, 50240) with the summit offset up by exactly 40.
if "$BIN" --chrom-sizes "$WORK/sizes.txt" --cutoff 5 --min-length 50 --max-gap 100 \
        --extend-peaks 40 "$WORK/frags.bed" >"$WORK/ext.narrowPeak" 2>/dev/null; then
    en=$(wc -l <"$WORK/ext.narrowPeak" | tr -d ' ')
    ecols=$(head -1 "$WORK/ext.narrowPeak" | awk -F'\t' '{print NF}')
    estart=$(head -1 "$WORK/ext.narrowPeak" | cut -f2)
    eend=$(head -1 "$WORK/ext.narrowPeak" | cut -f3)
    esummit=$(head -1 "$WORK/ext.narrowPeak" | cut -f10)
    if [ "$en" -eq 1 ] && [ "$ecols" = "10" ] && [ "$estart" = "49960" ] &&
            [ "$eend" = "50240" ] && [ "$esummit" -ge 120 ] && [ "$esummit" -le 160 ]; then
        ok "--extend-peaks pads the reported interval symmetrically"
    else
        fail "--extend-peaks pads the reported interval: got $estart-$eend summit $esummit ($en rows, $ecols cols)"
    fi
else
    fail "--extend-peaks pads the reported interval: non-zero exit"
fi

expect_status 1 "a negative --extend-peaks is rejected" \
    "$BIN" --chrom-sizes "$WORK/sizes.txt" --extend-peaks -1 "$WORK/frags.bed"

# The 3' edge is clamped to the contig length: chr1 is 100000 here and the
# planted peak ends at 50200, so a pad wider than the tail must still stop at
# 100000 rather than emitting a coordinate past the chromosome.
awk 'BEGIN { for (i = 0; i < 300; i++) printf "chr1\t99700\t99900\n" }' >"$WORK/edge.bed"
sort -k1,1 -k2,2n "$WORK/bg.bed" "$WORK/edge.bed" >"$WORK/edgefrags.bed"
if "$BIN" --chrom-sizes "$WORK/sizes.txt" --cutoff 5 --min-length 50 --max-gap 100 \
        --extend-peaks 5000 "$WORK/edgefrags.bed" >"$WORK/edge.narrowPeak" 2>/dev/null; then
    cend=$(head -1 "$WORK/edge.narrowPeak" | cut -f3)
    if [ "$cend" -le 100000 ]; then
        ok "--extend-peaks clamps the 3' edge to the contig length"
    else
        fail "--extend-peaks clamps the 3' edge: got end $cend past 100000"
    fi
else
    fail "--extend-peaks clamps the 3' edge: non-zero exit"
fi

# A stricter FDR than the score cutoff implies is the only thing that can
# filter; at 1e-300 nothing survives.
if "$BIN" --chrom-sizes "$WORK/sizes.txt" --cutoff 5 --min-length 50 --max-gap 100 \
        --qvalue 1e-300 "$WORK/frags.bed" >"$WORK/strict.out" 2>/dev/null &&
        [ ! -s "$WORK/strict.out" ]; then
    ok "--qvalue gates the output"
else
    fail "--qvalue gates the output"
fi

expect_status 1 "a nonsense --qvalue is rejected" \
    "$BIN" --chrom-sizes "$WORK/sizes.txt" --qvalue 2 "$WORK/frags.bed"

# A contig name longer than the writer's line buffer. narrowPeak is written
# with snprintf, which reports the length the line WOULD have needed rather
# than the length it wrote -- so a name this program does not choose, and which
# a BAM header may set to 64 KiB, used to be written back out of a 512-byte
# stack buffer along with whatever followed it. The check is that the line is
# whole and has its ten columns, not merely that nothing crashed: a clamp would
# also stop the over-read and would silently emit a truncated track.
LONGNAME="chrUn_$(printf 'S%.0s' $(seq 1 600))"
printf '%s\t100000\n' "$LONGNAME" >"$WORK/long.sizes"
awk -v n="$LONGNAME" 'BEGIN { for (i = 0; i < 300; ++i) printf "%s\t%d\t%d\n", n, 1000 + i, 1120 + i }' \
    >"$WORK/long.bed"
if "$BIN" --chrom-sizes "$WORK/long.sizes" --cutoff 1 --min-length 20 \
        "$WORK/long.bed" >"$WORK/long.narrowPeak" 2>/dev/null; then
    lcols=$(head -1 "$WORK/long.narrowPeak" | awk -F'\t' '{print NF}')
    lname=$(head -1 "$WORK/long.narrowPeak" | cut -f1)
    if [ "$lcols" = "10" ] && [ "$lname" = "$LONGNAME" ]; then
        ok "a contig name past the line buffer is written whole"
    else
        fail "a contig name past the line buffer is written whole: $lcols columns, name length ${#lname}"
    fi
else
    fail "a contig name past the line buffer is written whole: non-zero exit"
fi

# ---------------------------------------------------------------------------
# Torture regressions, 2026-09-10 (docs/TORTURE_2026-09-10.md, pto-peaks).
# ---------------------------------------------------------------------------

# An output that is an input. Each of these used to exit 0 with the input
# replaced by the narrowPeak computed from it.
cp "$WORK/frags.bed" "$WORK/alias.bed"
cp "$WORK/sizes.txt" "$WORK/alias.sizes"
expect_status 1 "-o the input BED is refused" \
    "$BIN" --chrom-sizes "$WORK/sizes.txt" -o "$WORK/alias.bed" "$WORK/alias.bed"
expect_status 1 "-o the --chrom-sizes file is refused" \
    "$BIN" --chrom-sizes "$WORK/alias.sizes" -o "$WORK/alias.sizes" "$WORK/frags.bed"
# The single quotes are deliberate: $1..$3 expand in the inner sh, not here.
# shellcheck disable=SC2016
expect_status 1 "-o the file stdin reads is refused" \
    sh -c '"$1" --chrom-sizes "$2" -o "$3" <"$3"' sh "$BIN" "$WORK/sizes.txt" "$WORK/alias.bed"
if cmp -s "$WORK/frags.bed" "$WORK/alias.bed" && cmp -s "$WORK/sizes.txt" "$WORK/alias.sizes"; then
    ok "refused runs leave their inputs intact"
else
    fail "refused runs leave their inputs intact"
fi

# A narrowPeak write that fails partway: a 512-byte file size limit with
# SIGXFSZ ignored, so write() returns EFBIG. Thirty planted peaks is well past it.
#
# Under ThreadSanitizer this ulimit kills the binary before it even reaches
# main: TSan's own runtime needs somewhere between 128 KiB and 512 KiB of file
# size headroom just to initialize (confirmed by sweeping the limit), well
# above the 512 bytes this test deliberately imposes -- so it SIGBUSes inside
# TSan's own strlen interceptor during Options::Options()'s first std::string
# construction, before any argument parsing or file I/O happens. Raising the
# limit enough to clear that floor also clears the ~2 KiB this fixture writes,
# which would stop the test from ever exercising EFBIG at all; the fixture
# would need tens of thousands of peaks to still overflow a 512 KiB+ limit,
# which is a real cost to a CLI smoke test under TSan's own overhead for no
# corresponding gain. Same category as the hnswlib deadlock false positives in
# docker/tsan.supp: a tool/environment limitation, not an application defect.
awk 'BEGIN { for (s = 0; s + 200 <= 99800; s += 200) printf "chr1\t%d\t%d\n", s, s + 200
             for (p = 0; p < 30; p++) for (i = 0; i < 300; i++) printf "chr1\t%d\t%d\n", 2000 + p * 3000, 2200 + p * 3000 }' |
    sort -k1,1 -k2,2n >"$WORK/many.bed"
if [ -n "${TSAN_OPTIONS:-}" ]; then
    printf '  SKIP a narrowPeak write failure exits 1 (TSan runtime needs more file-size headroom than this test grants; see comment above)\n'
else
    expect_status 1 "a narrowPeak write failure exits 1" \
        sh -c 'trap "" XFSZ; ulimit -f 1; exec "$@"' sh \
        "$BIN" --chrom-sizes "$WORK/sizes.txt" --cutoff 5 --min-length 50 --max-gap 100 \
        -o "$WORK/big.narrowPeak" "$WORK/many.bed"
fi

# Non-finite numbers. Each used to run to exit 0 with zero peaks.
for v in nan inf -inf; do
    expect_status 1 "--lambda-bg $v is rejected" \
        "$BIN" --chrom-sizes "$WORK/sizes.txt" --lambda-bg "$v" "$WORK/frags.bed"
done
expect_status 1 "--cutoff inf is rejected" \
    "$BIN" --chrom-sizes "$WORK/sizes.txt" --cutoff inf "$WORK/frags.bed"

# CRLF. A chrom.sizes with CRLF endings failed every line; a CRLF blank line in
# the BED was a malformed record. Both must now give the LF answer exactly.
printf 'chr1\t100000\r\nchr2\t50000\r\n' >"$WORK/crlf.sizes"
if "$BIN" --chrom-sizes "$WORK/crlf.sizes" --cutoff 5 --min-length 50 --max-gap 100 \
        "$WORK/frags.bed" >"$WORK/crlf.narrowPeak" 2>/dev/null &&
        cmp -s "$WORK/crlf.narrowPeak" "$WORK/peaks.narrowPeak"; then
    ok "a CRLF chrom.sizes gives the LF answer"
else
    fail "a CRLF chrom.sizes gives the LF answer"
fi
{ head -n 5 "$WORK/frags.bed"; printf '\r\n'; tail -n +6 "$WORK/frags.bed"; } >"$WORK/crblank.bed"
if "$BIN" --chrom-sizes "$WORK/sizes.txt" --cutoff 5 --min-length 50 --max-gap 100 \
        "$WORK/crblank.bed" >"$WORK/crblank.narrowPeak" 2>/dev/null &&
        cmp -s "$WORK/crblank.narrowPeak" "$WORK/peaks.narrowPeak"; then
    ok "a CRLF blank line in the BED is skipped"
else
    fail "a CRLF blank line in the BED is skipped"
fi

# A line with no end. Read into memory whole before, and refused only after.
head -c 3000000 /dev/zero | tr '\000' 'A' >"$WORK/noline.bed"
expect_status 2 "a line longer than 1 MiB exits 2" \
    "$BIN" --chrom-sizes "$WORK/sizes.txt" "$WORK/noline.bed"

# A compressed fragment BED is refused, and the message says what it is.
gzip -c "$WORK/frags.bed" >"$WORK/frags.bed.gz"
expect_status 2 "a gzipped fragment BED exits 2" \
    "$BIN" --chrom-sizes "$WORK/sizes.txt" "$WORK/frags.bed.gz"
if grep -q "decompress it first" "$WORK/err"; then
    ok "a gzipped fragment BED is named as one"
else
    fail "a gzipped fragment BED is named as one"
fi

# ---------------------------------------------------------------------------
# BAM input. The fixture is a background tiling plus the same planted peak, so
# the two front ends can be compared on the same signal.
# ---------------------------------------------------------------------------
# The fixture must actually be BGZF: a build without a decompressor emits a
# diagnostic instead, and a non-empty file is not the same as a BAM.
if [ -z "$BAMGEN" ] || ! "$BAMGEN" --emit-fixture >"$WORK/fix.bam" 2>/dev/null ||
        [ ! -s "$WORK/fix.bam" ] ||
        [ "$(dd if="$WORK/fix.bam" bs=1 count=2 2>/dev/null | od -An -tx1 | tr -d ' ')" != "1f8b" ]; then
    printf '  SKIP BAM cases (no fixture writer, or this build has no BGZF support)\n'
else
    # A BAM carries its own contig lengths, so --chrom-sizes is not required.
    if "$BIN" --cutoff 5 --min-length 50 --max-gap 100 "$WORK/fix.bam" \
            >"$WORK/bam.narrowPeak" 2>/dev/null; then
        ok "BAM input needs no --chrom-sizes"
    else
        fail "BAM input needs no --chrom-sizes"
    fi

    bstart=$(head -1 "$WORK/bam.narrowPeak" | cut -f2)
    bend=$(head -1 "$WORK/bam.narrowPeak" | cut -f3)
    if [ "$bstart" = "50000" ] && [ "$bend" = "50200" ]; then
        ok "BAM finds the planted peak at the right coordinates"
    else
        fail "BAM finds the planted peak: got $bstart-$bend"
    fi

    # The format is detected from the first byte, so a pipe works too.
    if "$BIN" --cutoff 5 --min-length 50 --max-gap 100 <"$WORK/fix.bam" \
            >"$WORK/bamstdin.narrowPeak" 2>/dev/null &&
            cmp -s "$WORK/bamstdin.narrowPeak" "$WORK/bam.narrowPeak"; then
        ok "BAM on stdin is detected and matches"
    else
        fail "BAM on stdin is detected and matches"
    fi

    # A BAM cut at a block boundary decodes perfectly and is missing the end of
    # the chromosome. That must not be reported as success.
    head -c 20000 "$WORK/fix.bam" >"$WORK/cut.bam"
    expect_status 2 "a truncated BAM exits 2" "$BIN" "$WORK/cut.bam"

    # One flipped byte inside a compressed block: the CRC is what catches it.
    cp "$WORK/fix.bam" "$WORK/corrupt.bam"
    printf '\377' | dd of="$WORK/corrupt.bam" bs=1 seek=900 count=1 conv=notrunc 2>/dev/null
    expect_status 2 "a corrupted BAM block exits 2" "$BIN" "$WORK/corrupt.bam"

    # Not a BAM at all, but starts with the gzip magic.
    # Octal, not \xHH: POSIX printf has no hex escape and dash emits it literally.
    printf '\037\213\010\000rubbish' >"$WORK/notbam.gz"
    expect_status 2 "gzip that is not BGZF exits 2" "$BIN" "$WORK/notbam.gz"

    # Every alignment in the fixture has mapq 30.
    if "$BIN" --cutoff 5 --min-length 50 --max-gap 100 --min-mapq 60 \
            "$WORK/fix.bam" >"$WORK/mapq.narrowPeak" 2>/dev/null &&
            [ ! -s "$WORK/mapq.narrowPeak" ]; then
        ok "--min-mapq filters BAM alignments"
    else
        fail "--min-mapq filters BAM alignments"
    fi

    expect_status 1 "a nonsense --min-mapq is rejected" \
        "$BIN" --min-mapq 999 "$WORK/fix.bam"
fi

printf '%d failure(s)\n' "$FAILURES"
[ "$FAILURES" -eq 0 ]
