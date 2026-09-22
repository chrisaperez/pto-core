#!/bin/sh
# SPDX-License-Identifier: MIT
#
# fastq_stream CLI, driven as a black box. POSIX sh + coreutils + gzip.
#
# Covers the defects in docs/TORTURE_2026-09-10.md that live in main.cpp or
# only show through a real process: argument validation, report files, and
# inputs that are not regular files. test_trimmer covers the library.
#
#   sh tests/test_cli.sh path/to/fastq_stream
set -u

BIN=${1:?usage: test_cli.sh path/to/fastq_stream}
WORK=$(mktemp -d "${TMPDIR:-/tmp}/fq_cli.XXXXXX")
trap 'rm -rf "$WORK"' EXIT INT TERM
FAILS=0
CHECKS=0

pass() { CHECKS=$((CHECKS + 1)); }
fail() { CHECKS=$((CHECKS + 1)); FAILS=$((FAILS + 1)); echo "FAIL: $*" >&2; }

# expect_exit CODE DESCRIPTION -- CMD...
expect_exit() {
  want=$1; what=$2; shift 3
  "$@" >"$WORK/stdout" 2>"$WORK/stderr"
  got=$?
  if [ "$got" -eq "$want" ]; then pass; else
    fail "$what: exit $got, wanted $want"; sed 's/^/    /' "$WORK/stderr" >&2
  fi
}

# A FASTQ large enough to outrun the reader's 1 MiB first read. SEQ and QUAL
# are cut from one length so the fixture cannot drift out of step with itself.
SEQ=ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
QUAL=$(printf '%s' "$SEQ" | sed 's/./I/g')
i=0
while [ $i -lt 12000 ]; do
  printf '@r%d\n%s\n+\n%s\n' $i "$SEQ" "$QUAL"
  i=$((i + 1))
done >"$WORK/in.fq"
READS=12000
expect_exit 0 "fixture is valid" -- "$BIN" -i "$WORK/in.fq" --qc-only

# A background FIFO peer blocks in open() until the binary opens the other end.
# If the binary exits early instead, a bare `wait` hangs the suite forever, so
# every peer is reaped by pid and killed first.
reap() { kill "$1" 2>/dev/null; wait "$1" 2>/dev/null; }

reads_in() { sed -n 's/.*"reads_in": \([0-9]*\).*/\1/p' "$1"; }

# --- output aliasing input ----------------------------------------------------
cp "$WORK/in.fq" "$WORK/same.fq"
expect_exit 1 "-o same as -i is refused" -- "$BIN" -i "$WORK/same.fq" -o "$WORK/same.fq"
if cmp -s "$WORK/in.fq" "$WORK/same.fq"; then pass; else fail "-i X -o X destroyed the input"; fi
# Through the standard streams: the shell does not truncate an input redirect,
# and `--stdout >> X` would read its own output back.
# The single quotes are deliberate: $1/$2 expand in the inner sh, not here.
# shellcheck disable=SC2016
expect_exit 1 "-o the file stdin reads is refused" -- \
  sh -c '"$1" -i - -o "$2" < "$2"' sh "$BIN" "$WORK/same.fq"
if cmp -s "$WORK/in.fq" "$WORK/same.fq"; then pass; else fail "-o X < X destroyed the input"; fi
# shellcheck disable=SC2016
expect_exit 1 "--stdout appended to the input is refused" -- \
  sh -c '"$1" -i "$2" --stdout >> "$2"' sh "$BIN" "$WORK/same.fq"
if cmp -s "$WORK/in.fq" "$WORK/same.fq"; then pass; else fail "--stdout >> X changed the input"; fi

# --- non-seekable input: FIFO, the shape of `-i <(zcat ...)` --------------------
mkfifo "$WORK/in.fifo"
cat "$WORK/in.fq" >"$WORK/in.fifo" &
peer=$!
expect_exit 0 "FIFO input" -- "$BIN" -i "$WORK/in.fifo" --qc-only -j "$WORK/fifo.json"
reap $peer
if [ "$(reads_in "$WORK/fifo.json")" = "$READS" ]; then pass; else
  fail "FIFO input lost data: reads_in=$(reads_in "$WORK/fifo.json"), wanted $READS"; fi

# --- adapter options that used to be silently substituted ---------------------
expect_exit 2 "empty --adapter"            -- "$BIN" -i "$WORK/in.fq" --qc-only -a ''
expect_exit 2 "non-ACGTN --adapter"        -- "$BIN" -i "$WORK/in.fq" --qc-only -a 'AGAT-CGG'
expect_exit 2 "--adapter-mismatch > 0.5"   -- "$BIN" -i "$WORK/in.fq" --qc-only --adapter-mismatch 0.9
expect_exit 2 "--adapter-seed > overlap"   -- "$BIN" -i "$WORK/in.fq" --qc-only --adapter-seed 12
expect_exit 2 "overlap > default adapter"  -- "$BIN" -i "$WORK/in.fq" --qc-only --adapter-overlap 14
expect_exit 2 "overlap > given adapter"    -- "$BIN" -i "$WORK/in.fq" --qc-only -a AGATCG --adapter-overlap 7
expect_exit 0 "consistent adapter options" -- "$BIN" -i "$WORK/in.fq" --qc-only -a AGATCGGAAGAGC --adapter-seed 8 --adapter-overlap 10 --adapter-mismatch 0.5
expect_exit 0 "--no-adapter skips adapter validation" -- "$BIN" -i "$WORK/in.fq" --qc-only --no-adapter --adapter-overlap 100

# --- report files -------------------------------------------------------------
# A JSON write that fails partway, deterministically: a 512-byte file size
# limit with SIGXFSZ ignored (an ignored disposition survives exec), so the
# write returns EFBIG. A pipe whose reader exits early is the obvious way to
# stage this and is racy -- the report fits in the pipe buffer.
#
# Under ThreadSanitizer this ulimit kills the binary before main runs at all:
# TSan's own runtime needs somewhere between 128 KiB and 512 KiB of file size
# headroom just to initialize (confirmed by sweeping the limit against
# pto-peaks's identical construction -- same TSan runtime, same failure), well
# above the 512 bytes this test deliberately imposes, so it SIGBUSes inside
# TSan's own strlen interceptor before any argument parsing or file I/O
# happens. Raising the limit enough to clear that floor also clears the
# report this fixture actually writes, which would stop the test from ever
# exercising EFBIG at all. Same category as the hnswlib deadlock false
# positives in docker/tsan.supp: a tool/environment limitation, not an
# application defect.
if [ -n "${TSAN_OPTIONS:-}" ]; then
  echo "  SKIP JSON report write failure (TSan runtime needs more file-size headroom than this test grants; see comment above)"
else
  expect_exit 1 "JSON report write failure" -- \
    sh -c 'trap "" XFSZ; ulimit -f 1; exec "$@"' sh \
    "$BIN" -i "$WORK/in.fq" --qc-only --report /dev/null -j "$WORK/big.json"
fi
expect_exit 1 "unwritable report path" -- "$BIN" -i "$WORK/in.fq" --qc-only --report "$WORK/no/such/dir/r.txt"

CR=$(printf '\r')
cp "$WORK/in.fq" "$WORK/a${CR}b.fq"
expect_exit 0 "control character in input path" -- "$BIN" -i "$WORK/a${CR}b.fq" --qc-only -j "$WORK/cr.json"
if LC_ALL=C grep -q "$CR" "$WORK/cr.json"; then fail "raw CR in JSON report"; else pass; fi

# --- truncated gzip -----------------------------------------------------------
gzip -c "$WORK/in.fq" >"$WORK/in.fq.gz"
size=$(wc -c <"$WORK/in.fq.gz")
head -c $((size - 8)) "$WORK/in.fq.gz" >"$WORK/notrailer.fq.gz"
head -c 10 "$WORK/in.fq.gz" >"$WORK/header.fq.gz"
expect_exit 0 "whole gzip"               -- "$BIN" -i "$WORK/in.fq.gz" --qc-only
expect_exit 1 "gzip missing CRC trailer" -- "$BIN" -i "$WORK/notrailer.fq.gz" --qc-only
expect_exit 1 "gzip header only"         -- "$BIN" -i "$WORK/header.fq.gz" --qc-only

# --- last record without a final newline ----------------------------------------
printf '@r\nACGT\n+\nIIII' >"$WORK/nonl.fq"
expect_exit 0 "missing final newline" -- "$BIN" -i "$WORK/nonl.fq" -o "$WORK/nonl.out" -l 0 --no-adapter
if [ "$(cat "$WORK/nonl.out"; echo x)" = "$(printf '@r\nACGT\n+\nIIII\nx')" ]; then pass; else
  fail "missing final newline: output not a terminated record"; fi

# --- invalid quality bytes ------------------------------------------------------
printf '@r\nACGT\n+\nII\001I\n' >"$WORK/ctrl.fq"
expect_exit 1 "control byte in QUAL" -- "$BIN" -i "$WORK/ctrl.fq" --qc-only

# --- audit 2026-09-11 (docs/AUDIT_2026-09-11_fastq_stream.md) ------------------
# Phred+64: 'B' is Q2 and 'h' is Q40 at offset 64. Undeclared, the run reports
# (JSON written) and exits 2; declared, it converts to Phred+33; an explicit 33
# overrides the detection for the rare Phred+33 file that looks the same.
Q64="BBBB$(printf '%s' "$SEQ" | cut -c5- | sed 's/./h/g')"
i=0
while [ $i -lt 200 ]; do
  printf '@p%d\n%s\n+\n%s\n' $i "$SEQ" "$Q64"
  i=$((i + 1))
done >"$WORK/p64.fq"
expect_exit 2 "undeclared Phred+64 exits 2" -- "$BIN" -i "$WORK/p64.fq" --qc-only -j "$WORK/p64.json"
if [ -s "$WORK/p64.json" ]; then pass; else fail "undeclared Phred+64: report not written before exit 2"; fi
expect_exit 0 "--phred-offset 64 converts" -- \
  "$BIN" -i "$WORK/p64.fq" -o "$WORK/p64.out" --phred-offset 64 -l 0 --no-adapter --window 0
WANT33="####$(printf '%s' "$SEQ" | cut -c5- | sed 's/./I/g')"
if [ "$(sed -n 4p "$WORK/p64.out")" = "$WANT33" ]; then pass; else fail "--phred-offset 64: output not Phred+33"; fi
expect_exit 0 "--phred-offset 33 overrides detection" -- "$BIN" -i "$WORK/p64.fq" --qc-only --phred-offset 33
# High-accuracy Phred+33 with no base below Q31 must not be flagged: AVITI's
# Q44-Q50 and HiFi's Q93 both tripped the detector's first version.
AVITI="$(printf '%s' "$SEQ" | sed 's/./M/g')"   # Q44
HIFI="$(printf '%s' "$SEQ" | sed 's/./~/g; s/^..../AAAA/')"   # Q93, a few Q32
i=0
while [ $i -lt 200 ]; do
  printf '@a%d\n%s\n+\n%s\n' $i "$SEQ" "$AVITI"
  i=$((i + 1))
done >"$WORK/aviti.fq"
sed "s/^$AVITI\$/$HIFI/" "$WORK/aviti.fq" >"$WORK/hifi.fq"
expect_exit 0 "AVITI-like Phred+33 (Q44) is not flagged" -- "$BIN" -i "$WORK/aviti.fq" --qc-only
expect_exit 0 "HiFi-like Phred+33 (Q93) is not flagged" -- "$BIN" -i "$WORK/hifi.fq" --qc-only
expect_exit 2 "--phred-offset 50 is rejected" -- "$BIN" -i "$WORK/in.fq" --qc-only --phred-offset 50

# A control byte in a read name.
printf '@r\033[2Jx\nACGT\n+\nIIII\n' >"$WORK/esc.fq"
expect_exit 1 "control byte in a header" -- "$BIN" -i "$WORK/esc.fq" -o "$WORK/esc.out" -l 0 --no-adapter

# CR-only line endings are named, not reported as over-long reads.
printf '@r\rACGT\r+\rIIII\r@s\rACGT\r+\rIIII\r' >"$WORK/cr.fq"
expect_exit 1 "CR-only line endings" -- "$BIN" -i "$WORK/cr.fq" --qc-only
if grep -q "CR-only" "$WORK/stderr"; then pass; else fail "CR-only line endings: message does not say so"; fi

# An 80 kB record is refused wherever it sits, not only when it straddles a chunk.
LONG=$(head -c 40000 /dev/zero | tr '\000' 'A')
LONGQ=$(printf '%s' "$LONG" | sed 's/./I/g')
printf '@long\n%s\n+\n%s\n' "$LONG" "$LONGQ" >"$WORK/long_first.fq"
cat "$WORK/in.fq" >>"$WORK/long_first.fq"
expect_exit 1 "an over-limit record at offset 0" -- "$BIN" -i "$WORK/long_first.fq" --qc-only

echo "test_cli: $((CHECKS - FAILS))/$CHECKS checks passed"
[ "$FAILS" -eq 0 ]
