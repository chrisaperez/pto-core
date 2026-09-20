# Audit — `genomic_toolkit`, 2026-09-11

Audit target 3: the `genomic_toolkit` module — its three subcommands (`sizes`,
`markdup`, `frip`), the text readers (fragment BED / 10x `fragments.tsv`,
BEDPE, peak BED), the htslib BAM path, CLI option parsing, and the 64-bit
accumulators. Built as shipped (Release, htslib 1.22, OpenMP) and under
AddressSanitizer + UndefinedBehaviorSanitizer, on Apple M4 Pro / Apple clang.

Every finding was reproduced against the pre-fix binary before it was patched.
`scripts/test_cli.sh` fails **37 of 56** checks on the pre-fix binary and passes
56/56 after.

**Out of scope, and why.** `cuttag_profiler` is a separate module, not part of
`genomic_toolkit`; its BED/GTF reader is still the unfuzzed parser that
`TORTURE_2026-08-17.md` lists as a residual. The "matrix serializers" named in
the audit brief are `cuttag_profiler`'s too — `genomic_toolkit` writes one
optional two-column histogram and nothing else.

## Method

1. **Sanitizer baseline first.** The 8 existing suites under
   `-fsanitize=address,undefined -fno-sanitize-recover=all -O1 -g`: **8/8 clean**.
   None of the findings below came from a sanitizer report on the existing
   suites; all came from driving hostile inputs at the binary.
2. **Black-box probes** of every numeric option, every text format (CRLF,
   zero-length, inverted, negative, non-numeric, partially numeric, past int32,
   empty, overlapping, unsorted, corrupt count column) and of BAMs that were
   truncated mid-block, truncated at a BGZF block boundary, stripped of the EOF
   block, bit-flipped, header-only, empty and fake.
3. **Reachability measured, not estimated**, for the 64-bit accumulators.
4. Fix, then pin each fix with a unit test where the behaviour lives in the
   library and with `scripts/test_cli.sh` where it lives in `src/main.cpp`
   (which no unit suite links).

## Findings

G1 is the multi-contig FRiP fault recorded in `TORTURE_2026-09-10.md`.

| # | Severity | Defect | Status |
|---|---|---|---|
| G2 | **High** | A BGZF input (BAM, `.gz`) cut at a block boundary reports statistics for the surviving prefix, exit 0 | Fixed |
| G3 | **High** | `frip` on unsorted input leaves duplicates unmarked and counted, exit 0 (`markdup` refuses the same input) | Fixed |
| G4 | Medium | Numeric options: SIGABRT on a non-number; silent int32 narrowing, prefix reads and hex-as-zero; no range checks | Fixed |
| G5 | Medium | Fragment-BED count column fails open: a corrupt count is one copy, an absurd one is clamped to a million | Fixed |
| G6 | Medium | FRiP's signed base totals overflow (UB) and the size sum wraps, reachable from a crafted file in minutes | Fixed |
| G7 | Low | BEDPE: a mate with an inverted interval is accepted; a MAPQ outside 0–255 is filed as a filter drop | Fixed |
| G8 | Low | A file in which every record is malformed is reported as an empty library, exit 0 | Fixed |
| G9 | Low | A peak file of only zero-length intervals is diagnosed as a chromosome-naming mismatch | Fixed |
| G10 | **High** | A BED5/BED6 score column was read as a read-pair count, multiplying every statistic by the score, exit 0 | Fixed (second pass) |

### G2 — a truncated BAM is half a library (High)

ENCODE `ENCFF121ZQX` chr22 slice: 423 BGZF blocks. Cut at block 211:

```
$ genomic_toolkit sizes cut_at_block.bam --json
[W::bam_hdr_read] EOF marker is absent. The input is probably truncated
{"records_read":66296,"fragments":33112, ...}          # of 132,934 / 66,406
$ echo $?
0
```

A cut *inside* a block already failed (htslib reports the short read); a cut at
a boundary decodes cleanly, and htslib's warning goes to stderr where a
pipeline does not look. An interrupted `cp`/`aws s3 cp` produces exactly this.
**Fix:** `stream_bam` calls `bgzf_check_EOF` right after the header, before any
record is read, whenever `hts_get_format(fp)->compression == bgzf` (not
`is_bgzf`, which htslib documents as not implying BGZF); a missing marker is
exit 1 with no partial report. Same rule `fastq_stream` and `pto-peaks` already
enforce. **Test:** `test_cli.sh` writes a 40,000-pair BAM and checks the intact
file, a block-boundary cut (`sizes` and `frip`), a missing EOF block alone, and
a mid-block cut.

### G3 — `frip` did not check what `markdup` checks (High)

```
chr1 100 300 / chr1 500 900 / chr1 100 300        # a duplicate, out of order
frip             -> fragments 3, duplicates_excluded 0, exit 0
frip --unsorted  -> fragments 2, duplicates_excluded 1
markdup          -> exit 2, "records arrived out of coordinate order"
```

Both subcommands run the same `DuplicateMarker` in sorted-stream mode; only
`markdup` read its `out_of_order` counter. The FRiP printed includes fragments
it claims to have excluded. This is the pattern `REVIEW_2026-08-15.md`'s
retrospective names: a guard present in one entry point and missing in its
sibling. **Fix:** `frip` exits 2 on out-of-order input whenever it is marking
duplicates, after printing its report, with `--unsorted` / `--keep-dups` named
in the message. **Test:** `test_cli.sh`.

### G4 — numeric options (Medium)

`std::stol`/`std::stoi`/`std::stoll` followed by `static_cast`:

| input | before | after |
|---|---|---|
| `--min-length abc` | uncaught `std::invalid_argument` from argument parsing, outside `main`'s `try`: **SIGABRT, exit 134** | exit 1, named |
| `--min-length 99999999999999999999` | `std::out_of_range`, **exit 134** | exit 1 |
| `--min-length 100x` | read as **100**, exit 0 | exit 1 |
| `--min-length 4294967396` | narrowed to **100**, exit 0 | exit 1 |
| `--max-length 4294967296` | narrowed to **0** ("off"), exit 0 | exit 1 |
| `--genome-size 0x10` | read as **0**, enrichment silently off | exit 1 |
| `--min-mapq 300`, `-1`; `--genome-size -1` | accepted | exit 1 |

**Fix:** `include/toolkit/cli_args.hpp`, `parse_bounded_int` — the whole field a
base-10 integer inside the option's own range (length 0–2^31-1, MAPQ 0–255,
genome size ≥ 0), else refused by name. **Tests:** `gtk_cli_args` (unit),
`test_cli.sh` (the binary).

### G5 — the count column failed open (Medium)

```
chr1 100 300 BC 0 / BC -3 / BC abc / BC 99999999999999
sizes -> "fragments":1000003, "malformed":0, exit 0
```

`0`, `-3` and `abc` were each one copy; `99999999999999` was clamped to
1,000,000. Same class as `TORTURE_2026-08-17.md` T2 (a corrupt BEDPE MAPQ
passed as Q255), in the sibling column — and the fuzz suite's
`count column is bounded` case **pinned the fail-open behaviour as correct**.
**Fix:** a present column must be an integer in [1, 1,000,000] or the record is
malformed — since G10, only when the column is read at all (`--with-counts`). That test now asserts the opposite, with the history written above
it.

### G6 — 64-bit accumulators (Medium)

`FripCounter` kept `bases_total`/`bases_in_peaks` as `int64_t` (signed overflow
is UB); `SizeHistogram` kept the length sum as `uint64_t` (wraps silently, and
`mean()` then reports a plausible length). The count column makes both
reachable from a small file: one line can be a million 2 Gbp fragments
(`--max-length 0`).

| run (Release, after the fix) | lines | fragments | result |
|---|---|---|---|
| `frip --keep-dups --max-length 0` | 4,700 | 4.7 × 10^9 | **exit 1 after 142 s**: "FRiP total fragment bases exceeds the 64-bit range" |
| `sizes --max-length 0` | 9,300 | 9.3 × 10^9 | **exit 1 after 41 s**: "summed fragment length exceeds the 64-bit range" |
| `sizes --max-length 0` (control) | 4,700 | 4.7 × 10^9 | exit 0, exact mean — 9.4 × 10^18 fits in 64 bits |

The thresholds are 4,612 lines (signed) and 9,224 (unsigned). The pre-fix
binary was not run to completion on these inputs, so this record does not claim
what it printed; the defect is the unchecked arithmetic, and the line counts
above are where it is reached.

**Fix:** `checked_add` in `types.hpp`, throwing `std::overflow_error`. FRiP
checks once per batch, outside the OpenMP region (a batch of 4,096 max-length
fragments is under 2^44 and cannot overflow inside it; an exception may not
leave a parallel region); `SizeHistogram::merge` checks before modifying
anything. **Tests:** `checked_add` at both boundaries and in both directions
(`gtk_frip`); the histogram driven to the exact overflowing merge — the 34th
doubling of one 2^31−1 bp fragment — with the histogram unchanged after the
refusal (`gtk_fragment_stream`). The FRiP *wiring* is covered by the measured
run above, not by ctest: a unit test would need 4.3 × 10^9 fragments.

### G7 — BEDPE (Low)

`chr1 150 100 chr1 250 300` became part of a 150 bp fragment with 0 malformed;
min/max built a plausible span from an inverted mate. A MAPQ of `-7` was counted
in `dropped_mapq`, reporting a corrupt record as a quality decision; `256`
passed any threshold. **Fix:** both are malformed. **Tests:**
`gtk_fragment_stream`, `gtk_fuzz_parsers` (a `256` case).

### G8 — nothing parsed, exit 0 (Low)

A fragment file of `chr1 300 100` and `chr1 abc 200` reported 0 fragments,
2 malformed, exit 0 — indistinguishable from an empty library. A tolerated
handful of bad lines in a large file is deliberate (`read_peaks_bed`'s comment
gives the asymmetry); *every* record being malformed is the wrong file or the
wrong delimiter. **Fix:** all three subcommands exit 2 after the report when
`records_read > 0` and every record was malformed. An empty file is not that
case and still exits 0. **Test:** `test_cli.sh`.

### G9 — zero-length peaks (Low)

A peak file of only `start == end` intervals exited 2 with "1 stream contigs
(e.g. 'chr1') are absent from the peak file" — about a file naming chr1.
`PeakSet::build` drops zero-length intervals before interning their contig.
**Fix:** `frip` refuses an empty merged cover by name, with the 0-based,
half-open reminder. **Tests:** `gtk_peak_set` (the library reports an empty set
and counts its inputs), `test_cli.sh` (the CLI refuses).

### G10 — a BED score column multiplied every statistic (High)

`stream_fragment_bed` read column 5 of every fragment BED as the 10x
read-pair count. In a BED5/BED6 that column is a score — `bedtools bamtobed`
writes MAPQ there:

```
chr1 100 300 read1 60 +     ->  sizes: "fragments":60, exit 0
```

Every size, duplicate and FRiP count was multiplied by the score, with nothing
in the report to show it. Recorded as an open residual after the first pass
(G5 had made a `.` score malformed but left numeric scores multiplying), then
fixed on instruction. **Fix:** column 5 is read only under `--with-counts`
(`FragmentFilter::use_count_column`); by default each interval is one fragment
and column 5 is not parsed. With the flag, G5's strict count rule applies.
**Breaking for anyone running `genomic_toolkit` directly on a 10x
`fragments.tsv`**: without `--with-counts` the statistics count unique
fragments, not read pairs. No in-repo caller is affected — `pto-cloud`, the
Nextflow modules, the wheel API and `tests/validation` pass no fragment BED
with a count column (checked by search). **Tests:** `gtk_fragment_stream`
(BED5/BED6 scores 60, 0, `.` and 1000 are one fragment each by default, through
both `stream_fragment_bed` and `stream_fragments`; the same file with the flag
holds column 5 to the count rule), `gtk_fuzz_parsers` (the count cases with and
without the flag; the fragment-BED fuzz loop draws the flag at random),
`test_cli.sh` (`sizes` and `frip` on a BED6, and the count cases with and without
the flag).

## What did not break

- **Peak BED:** CRLF, inverted, negative, non-numeric, partially numeric,
  two-column, past int32 (named), empty (named), overlapping and unsorted
  records all behaved correctly. `read_peaks_bed` is fatal on malformed input by
  design.
- **Fragment BED / BEDPE:** CRLF, trailing junk (`300abc`), negative, hex,
  scientific notation, leading sign or space, the int32 ceiling, and zero-length
  fragments (dropped by the length filter, not counted as malformed).
- **BAM:** a mid-block cut and a flipped byte (CRC32 mismatch) exit 1;
  header-only, empty and non-BAM files exit 1 with a named error.
- **Sanitizers:** 10/10 suites clean after the fixes, including the new
  `gtk_cli_args` and `gtk_cli`, with the fuzz suite at ~26 s under ASan.
- **FRiP correctness on real data:** unchanged and exact against samtools +
  bedtools (below).

## Behaviour changes a user will see

| input | before | now |
|---|---|---|
| truncated BGZF/BAM | partial statistics, exit 0 | exit 1, nothing reported |
| `frip` on unsorted input, marking duplicates | exit 0 | exit 2 after the report |
| every record malformed (`sizes`/`markdup`/`frip`) | exit 0 | exit 2 after the report |
| malformed numeric option | exit 134, or silently reinterpreted | exit 1, named |
| peak file of only zero-length intervals | exit 2, naming message | exit 1, named |
| fragment BED with a column 5 (a BED5/BED6 score, or a 10x count) | emitted column-5 times | one fragment per interval; pass `--with-counts` for 10x read-pair counts |
| `--with-counts`, column 5 not an integer in 1–1,000,000 | 1 copy, or clamped | record malformed |

`pto-cloud`'s parameter bounds (`common/models.py`: lengths 0–1,000,000, MAPQ
0–255, genome size ≤ 10^11) sit inside the new ranges, so no job the API
accepts is refused by the binary; the Nextflow modules pass none of these
options; the wheel API passes values through unbounded by design.

## Residuals

- **`--with-counts` is a claim the user makes about the file.** Nothing checks
  that column 5 of a flagged file is really a count; a BED6 passed with the flag
  is multiplied by its scores exactly as every fragment BED was before G10.
- **A non-seekable BGZF source is not checked for its EOF marker**
  (`bgzf_check_EOF` returns 2). The binary dispatches on file extension, so this
  is not a path the CLI reaches today; a future stdin mode would reopen G2.
- **CRAM truncation was not exercised** — no CRAM fixture with a local
  reference. CRAM carries its own EOF container and is excluded from G2's check.
- **FRiP's overflow wiring is verified by one measured run, not by ctest** (see
  G6).
- **Unit tests for new APIs (`parse_bounded_int`, `checked_add`) cannot compile
  against the pre-fix tree**, so their pre-fix evidence is the probe results and
  the 37 `test_cli.sh` failures, not a red unit run.
- **`std::bad_alloc` resilience is still untested** (`TORTURE_2026-08-17.md`:
  macOS enforces no `RLIMIT_AS`).
- **Sanitizers ran on macOS arm64 / Apple clang only**, not Linux/GCC.
- **`cuttag_profiler`'s BED/GTF reader and its matrix writer** remain un-audited
  (out of scope above).

## Verification

| Configuration | Result |
|---|---|
| Release ctest `gtk_*` | 10/10 |
| ASan+UBSan ctest (`-fno-sanitize-recover=all`) | 10/10, no sanitizer report |
| ASan+UBSan ctest, the 8 pre-existing suites before any change | 8/8 |
| `scripts/test_cli.sh`, fixed / pre-fix binary | 56/56 / **37 failures** |
| `validate_against_bedtools.sh`, ENCFF121ZQX chr22 | counts agree exactly (22,473 / 66,467) |
| `frip` genome-wide, ENCFF121ZQX × ENCFF187JYM vs samtools+bedtools | 1,005,639 / 2,886,752, exact |
| `tests/integration/run_pipeline_workflow.sh` | exit 0, 20 assertions |
| `nextflow/tests`, rebuilt engines on `PATH` | 76 passed |
| `pto-core/packaging/tests` | 68 passed, 1 skipped |
| Crafted overflow inputs | `frip` exit 1 at 142 s; `sizes` exit 1 at 41 s |
| `shellcheck -S style` on the new scripts | clean |
