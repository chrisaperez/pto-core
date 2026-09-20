# Audit — `cuttag_profiler`, 2026-09-11

Audit target 4: `cuttag_profiler profile` — the regions × bins matrix builder,
its TSV serializers (`--out-matrix`, `--out-profile`), the BED/GTF region
readers, the BAM binning, CLI option handling, and the same options as the HTTP
`/api/profile` handler receives them. Built as shipped (Release, htslib 1.22),
under AddressSanitizer + UndefinedBehaviorSanitizer, and with
`-DPTO_CLOUD_BUILD=ON` (the worker image's configuration), on Apple M4 Pro /
Apple clang. Started after `main` (`6a434f5`) was confirmed green on both CI
workflows.

Every finding was reproduced against the pre-fix binary before it was patched.
The new `scripts/test_cli.sh` fails **19 of 24** checks on the pre-fix binary
and passes 24/24 after.

**Not re-covered.** The HTTP server's own security surface (all nine findings of
`SECURITY_HTTP_2026-08-15.md` are closed) and the output-overwrites-input and
write-failure cases closed in `TORTURE_2026-09-10.md`.

## Method

1. **Sanitizer baseline first.** The 10 existing suites under
   `-fsanitize=address,undefined -fno-sanitize-recover=all -O1 -g`: **10/10
   clean**. Every finding below came from hostile inputs driven at the binary,
   none from the existing suites.
2. **Black-box probes** on in-tree `demo.bam`: option values past `int`, negative
   and out of range; regions on contigs absent from the BAM, alone and mixed with
   real ones; a window at the size limits; RPKM output precision; a GTF whose
   attributes contain the name key as a substring.
3. **Confirmed under UBSan** where the defect was arithmetic.
4. Fixed on the option *types* (`ProfileOptions`, `FilterOptions`) so the CLI and
   HTTP paths share one set of guards — the lesson of `REVIEW_2026-08-15`
   finding 2 — and pinned with `test_signal_calc` cases plus
   `scripts/test_cli.sh` for what lives in `src/main.cpp`.

## Findings

| # | Severity | Defect | Status |
|---|---|---|---|
| C1 | **High** | Skipped regions (contig absent from the BAM) are averaged into the meta-profile as missing values, shrinking it in proportion, exit 0 | Fixed |
| C2 | **High** | Every region skipped produces an all-zero profile with exit 0 | Fixed |
| C3 | Medium | Signed overflow (UB) in bin geometry for a window the 10M-bin ceiling allows; Release writes wrong bin offsets, exit 0 | Fixed |
| C4 | Medium | `--min-mapq`/`--threads` narrowed to `int` before any check; MAPQ, fragment length, read extension and thread count unvalidated — CLI and HTTP | Fixed |
| C5 | Low | Skipped regions' matrix rows have an empty name | Fixed |
| C6 | Low | GTF name attribute matched as a substring: `havana_gene_name` answers for `gene_name` | Fixed |

### C1 — skipped regions dilute the meta-profile (High)

`demo.bam`, 50 promoter regions, the second 25 moved to a contig not in the BAM:

```
good 25 only:           column_mean  5.000363  5.200377 ... 61.154434
the same 25 + 25 absent: column_mean  2.500181  2.600189 ... 30.577217   exit 0
```

`compute_matrix` filled a skipped row with `missing_value` (0) and divided the
column sums by *all* rows. The only signal was "25 skipped" at the end of the
stderr summary. deepTools, the benchmark's parity reference, drops such
regions. **Fix:** skipped rows are excluded from the column mean (divided by
kept rows), from BPM's summed signal, and from scaling; they stay in the matrix,
named (C5). Off-contig bins inside a *kept* window are unchanged — they hold the
missing value and count as such, the parity `tests/validation` pins with
deepTools' `--missingDataAsZero`. **Tests:** `test_signal_calc` (raw and BPM,
the mixed set's mean equals the kept set's); `test_cli.sh` (the mixed and good
profiles are byte-identical).

### C2 — every region skipped (High)

All 50 regions on an absent contig: an 80-bin profile of `0.000000`, exit 0 — a
chromosome-naming mismatch (`chr1` vs `1`) indistinguishable from a dead
library, the trap `genomic_toolkit frip` already refuses. **Fix:** the column
mean is the missing value (no region contributed, so it is undefined), and
`profile` writes its outputs and then exits **2**, naming an absent contig.
**Tests:** `test_signal_calc` (mean equals a distinguishable missing value);
`test_cli.sh` (exit 2, message, outputs still on disk).

### C3 — bin geometry overflow (Medium)

`validate()` bounded the bin *count* at 10M but not the window, and every bin
boundary is `width * b / nbins`. Upstream 4e18, downstream 1, bin size 1e12 is
4,000,001 bins — accepted:

```
UBSan: signal_calc.cpp:196:40: runtime error: signed integer overflow:
       4000000000000000001 * 3 cannot be represented in type 'std::int64_t'
Release: exit 0; bin 2 offset -3999999805843057753, correct -3999997500000625000
```

The same products appear in `BamReader::query_bins` and `add_bases`. The window
is reachable from the CLI and from `/api/profile` JSON (behind the session
token); `pto-cloud` caps windows at 1,000,000 bp and cannot reach it. **Fix:**
`ProfileOptions::kMaxWindow = INT64_MAX / kMaxBins` (~922 Gbp), enforced in
`validate()`, so every `width × bin` product fits. **Tests:** `test_signal_calc`
(the reproducing window refused by `validate()` and by `compute_matrix`; the cap
accepted exactly and refused one base past it); `test_cli.sh`.

### C4 — option values (Medium)

| input | before |
|---|---|
| `--min-mapq 4294967296` | narrowed to **0**, exit 0 |
| `--min-mapq 300` | accepted: every alignment filtered out, exit 0 |
| `--min-mapq -1` | accepted |
| `--threads 4294967297` | narrowed to **1** |
| `--threads -3` | accepted, read as "all cores" |
| `--max-fragment -5` | accepted: the fragment gate silently off |
| `--extend-reads -100` | accepted |

The HTTP handler read `threads` and `minMapq` from JSON as `int` — the same
narrowing — and never validated the filters. Partial numbers (`--bin-size 50x`)
were already refused. **Fix:** `narrow_option()` range-checks a 64-bit value
before narrowing, used by the CLI (`--min-mapq`, `--threads`, `serve --port`) and
the HTTP handler; new `FilterOptions::validate()` (MAPQ 0–255, fragment length
and extension 0–2^40), called by both entry points and again by
`compute_matrix`; `ProfileOptions::validate()` refuses a thread count outside
0–4096. **Tests:** `test_signal_calc`, `test_cli.sh`.

### C5 — unnamed skipped rows (Low)

`row_names[i]` was assigned only on the success path, so the 25 skipped rows
above were written to `--out-matrix` with an empty first field. **Fix:** every
row is named before the skip test (column 4, or `chrom:start-end`). **Tests:**
both suites.

### C6 — GTF name attribute (Low)

`gene_id "G1"; havana_gene_name "WRONG"; gene_name "RIGHT";` named the region
`WRONG`: `attrs.find("gene_name")` matched inside the longer key. **Fix:** the key
must start the field or follow a separator and be followed by the value
delimiter. **Tests:** `test_signal_calc` (GTF and GFF3 forms; a record with only
the longer key gets no name); `test_cli.sh`.

## What did not break

- **Partial numbers** in CLI options were already refused (`to_int` checks
  `from_chars` consumed the whole field).
- **Coordinates:** negative, inverted and GTF-0 coordinates refused; the centre
  anchor and window edges use checked or saturating arithmetic (earlier
  regressions, still holding).
- **Output precision:** fixed 6-decimal output lost nothing on the demo data under
  RPKM — every nonzero cell was ≥ 125, and the 422 printed zeros are true zeros.
- **Determinism and sanitizers:** 11/11 suites clean under ASan+UBSan after the
  fixes, including the threaded skip tracking (each row's flag is written only by
  the worker that owns the row).
- **Cloud build:** with the HTTP server compiled out, 13/13 including the new CLI
  suite and the three no-server assertions.

## Behaviour changes a user will see

| input | before | now |
|---|---|---|
| some regions on contigs absent from the BAM | profile shrunk in proportion | profile over the regions that exist; skipped rows kept and named |
| every region skipped | all-zero profile, exit 0 | missing-value profile written, exit **2** |
| window over ~922 Gbp | UB; wrong offsets in Release | exit 1, named limit |
| `--min-mapq` outside 0–255, `--threads` outside 0–4096, negative `--max-fragment`/`--extend-reads` | narrowed or accepted | exit 1, named range |
| the same values in `/api/profile` JSON | narrowed or accepted | 400 with the message |
| GTF with a longer key containing `gene_name` | wrong name | the `gene_name` value |

`pto-cloud`'s `CuttagParams` bounds (windows ≤ 1,000,000 bp, bin 1–100,000, MAPQ
0–255, non-negative fragment and extension) and the Nextflow module's
`--threads ${task.cpus}` all sit inside the new ranges.

## Residuals

- **Parity with deepTools was not re-measured.** `computeMatrix`/`bamCoverage`
  are not installed on the audit machine. C1 changes the mean only when a region
  is skipped; the recorded `bench-gold` case had none, but "unchanged" is
  inferred, not observed.
- **`/api/profile` does not refuse an all-skipped request**: it returns 200 with
  `regionsSkipped == regions` and a missing-value profile. The CLI's exit 2 has no
  HTTP counterpart; whether that should be a 4xx is a dashboard decision.
- **Fixed 6-decimal TSV output can print a nonzero value below 5e-7 as
  `0.000000`.** Not reproduced on realistic data (above); changing the format
  would perturb every consumer, so it is recorded, not changed.
- **The GTF key match can still land inside a quoted value** such as
  `note "see gene_name here"`, where the key is preceded by a space and followed
  by one. No such annotation was found; the fix covers the observed case.
- **A skipped row's cells are still `missing_value` in `--out-matrix`**, so a
  consumer that averages the TSV itself, rather than reading `--out-profile`,
  reproduces C1 unless it drops those rows.
- **Sanitizers ran on macOS arm64 / Apple clang only.**

## Verification

| Configuration | Result |
|---|---|
| ASan+UBSan ctest, the 10 pre-existing suites before any change | 10/10 |
| C3 reproducer under UBSan, before | `signed integer overflow` at `signal_calc.cpp:196` |
| `scripts/test_cli.sh`, pre-fix / fixed binary | **19 failures** of 24 / 24/24 |
| Release ctest (11 suites incl. `profiler_cli`, `server_linkage`, `http_security`) | 11/11 |
| ASan+UBSan ctest after the fixes | 11/11, no sanitizer report |
| `-DPTO_CLOUD_BUILD=ON` build + ctest | 13/13 |
| `tests/integration/run_pipeline_workflow.sh` (profile stage) | exit 0, 20 assertions, profile unchanged (centre 20.3× flanks) |
| `nextflow/tests`, rebuilt engines on `PATH` | 76 passed |
| `pto-core/packaging/tests` | 68 passed, 1 skipped |
| deepTools parity benchmark | **not run** (deepTools not installed) |
