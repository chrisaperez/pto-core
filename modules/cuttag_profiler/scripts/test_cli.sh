#!/usr/bin/env bash
# SPDX-License-Identifier: MIT
#
# Black-box checks of `cuttag_profiler profile`.
#
# The unit suites link profiler_core, not src/main.cpp, so the CLI's argument
# handling and its exit codes are reachable only by running the binary. Every
# case below was a defect on 2026-09-11 -- a skipped-region-diluted profile, an
# all-skipped run exiting 0, option values narrowed or unchecked, a window
# whose bin arithmetic overflowed, a GTF name taken from the wrong attribute.
# The audit record is docs/AUDIT_2026-09-11_cuttag_profiler.md.
#
# The BAM is written here with samtools, because the module's demo data is not
# tracked; without samtools and python3 the script SKIPs, loudly.
#
# usage: test_cli.sh <cuttag_profiler>
set -uo pipefail  # not -e: exit codes are the thing under test

if [ "$#" -ne 1 ]; then
    echo "usage: $0 <cuttag_profiler-binary>" >&2
    exit 64
fi
BIN="$1"

for tool in samtools python3; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "skip: $tool is not installed (needed to write the test BAM)" >&2
        exit 0
    fi
done

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

says() {
    if grep -qF -- "$1" "$WORK/out"; then
        ok "$2"
    else
        fail "$2: output lacks '$1'"
    fi
}

# Three 20-read pile-ups at 20 kb, 40 kb and 60 kb on a 200 kb chr1.
python3 - "$WORK" <<'PY'
import sys
work = sys.argv[1]
lines = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:200000"]
seq, qual = "A" * 100, "I" * 100
n = 0
for tss in (20000, 40000, 60000):
    for _ in range(20):
        lines.append(f"r{n}\t0\tchr1\t{tss - 50 + 1}\t60\t100M\t*\t0\t0\t{seq}\t{qual}")
        n += 1
open(f"{work}/in.sam", "w").write("\n".join(lines) + "\n")
good = "".join(f"chr1\t{t}\t{t + 200}\tg{t}\t0\t+\n" for t in (20000, 40000, 60000))
absent = ("chrAbsent\t20000\t20200\tmissing1\t0\t+\n"
          "chrAbsent\t40000\t40200\t.\t0\t+\n"
          "chrAlsoAbsent\t100\t300\tmissing3\t0\t+\n")
open(f"{work}/good.bed", "w").write(good)
open(f"{work}/mixed.bed", "w").write(good + absent)
open(f"{work}/absent.bed", "w").write(absent)
open(f"{work}/names.gtf", "w").write(
    "chr1\tsrc\tgene\t20001\t20200\t.\t+\t.\t"
    "gene_id \"G1\"; havana_gene_name \"WRONG\"; gene_name \"RIGHT\";\n")
PY
samtools sort -o "$WORK/in.bam" "$WORK/in.sam" 2>/dev/null
samtools index "$WORK/in.bam"

run() {
    "$BIN" profile --bam "$WORK/in.bam" --threads 1 --quiet \
        --upstream 1000 --downstream 1000 --bin-size 100 --normalization raw "$@"
}

echo "== skipped regions"
expect 0 "the regions that exist" run --regions "$WORK/good.bed" \
    --out-matrix "$WORK/good_m.tsv" --out-profile "$WORK/good_p.tsv"
expect 0 "the same regions plus three on absent contigs" run --regions "$WORK/mixed.bed" \
    --out-matrix "$WORK/mixed_m.tsv" --out-profile "$WORK/mixed_p.tsv"
if cmp -s "$WORK/good_p.tsv" "$WORK/mixed_p.tsv"; then
    ok "skipped regions leave the meta-profile unchanged"
else
    fail "skipped regions changed the meta-profile"
    paste "$WORK/good_p.tsv" "$WORK/mixed_p.tsv" | sed -n '9,12p' | sed 's/^/          /' >&2
fi
empty_names=$(awk -F'\t' 'NR > 1 && $1 == ""' "$WORK/mixed_m.tsv" | wc -l | tr -d ' ')
if [ "$empty_names" = "0" ]; then
    ok "every matrix row, skipped or not, has a name"
else
    fail "$empty_names matrix rows have an empty name"
fi
if grep -q '^chrAbsent:40000-40200	' "$WORK/mixed_m.tsv"; then
    ok "an unnamed skipped region is named by its coordinates"
else
    fail "an unnamed skipped region is not named by its coordinates"
fi
expect 2 "every region on an absent contig" run --regions "$WORK/absent.bed" \
    --out-matrix "$WORK/absent_m.tsv" --out-profile "$WORK/absent_p.tsv"
says "chromosome-naming mismatch" "the error names the likely cause"
says "'chrAbsent'" "the error names an absent contig"
if [ -s "$WORK/absent_m.tsv" ] && [ -s "$WORK/absent_p.tsv" ]; then
    ok "the outputs are still written, as evidence"
else
    fail "the outputs were not written before exit 2"
fi

echo "== option values: refused by name, never narrowed"
for bad in 4294967296 300 -1; do
    expect 1 "--min-mapq $bad" run --regions "$WORK/good.bed" --out-profile "$WORK/p.tsv" --min-mapq "$bad"
    says "between 0 and 255" "--min-mapq $bad names the range"
done
for bad in 4294967297 -3; do
    expect 1 "--threads $bad" "$BIN" profile --bam "$WORK/in.bam" --regions "$WORK/good.bed" \
        --out-profile "$WORK/p.tsv" --quiet --threads "$bad"
done
expect 1 "--max-fragment -5 (silently disabled the gate)" run --regions "$WORK/good.bed" \
    --out-profile "$WORK/p.tsv" --max-fragment -5
expect 1 "--extend-reads -100" run --regions "$WORK/good.bed" \
    --out-profile "$WORK/p.tsv" --extend-reads -100 --no-extend-fragment
expect 1 "a window whose bin arithmetic overflowed int64" "$BIN" profile --bam "$WORK/in.bam" \
    --regions "$WORK/good.bed" --out-profile "$WORK/p.tsv" --quiet --threads 1 \
    --upstream 4000000000000000000 --downstream 1 --bin-size 1000000000000
says "bp limit" "the window error names the limit"
expect 0 "--min-mapq 255 is accepted" run --regions "$WORK/good.bed" \
    --out-profile "$WORK/p.tsv" --min-mapq 255

echo "== GTF names"
expect 0 "a GTF with havana_gene_name before gene_name" run --regions "$WORK/names.gtf" \
    --out-matrix "$WORK/gtf_m.tsv"
if [ "$(sed -n 2p "$WORK/gtf_m.tsv" | cut -f1)" = "RIGHT" ]; then
    ok "the region is named by gene_name, not havana_gene_name"
else
    fail "the region is named '$(sed -n 2p "$WORK/gtf_m.tsv" | cut -f1)', not RIGHT"
fi

printf '%d checks, %d failure(s)\n' "$checks" "$failures"
[ "$failures" -eq 0 ]
