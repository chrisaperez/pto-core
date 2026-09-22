# Audit — `fastq_stream`, 2026-09-11

Adversarial pass over `modules/fastq_stream`, the engine behind the cloud's
`fastq_qc` and `fastq_trim` job types. It follows
`docs/TORTURE_2026-09-10.md`, which already closed input truncation (F1),
BGZF EOF (F2), gzip truncation (F3), non-printable SEQ/QUAL (F7), the missing
final newline (F11) and CRLF handling. This pass was aimed past that record, at
four vectors: malformed parsing, quality encoding, block-boundary truncation,
and accumulator overflow.

Built as shipped (`FQ_ENABLE_HTS=ON`, `FQ_NATIVE=OFF`) on Apple M4 Pro, plus a
sanitizer build with `-fsanitize=address,undefined -fno-sanitize-recover=all`.
The baseline (commit `17ade94`) passed all 3 suites under that build before any
change. Every finding was reproduced against the frozen baseline binary first.

## Findings

| # | Severity | Defect | Status |
|---|---|---|---|
| A1 | **High** | Phred+64 input silently processed as Phred+33: wrong QC, no quality trimming, exit 0 | Fixed |
| A2 | Medium | Control and non-ASCII bytes in the header and `+` lines written to the output | Fixed |
| A3 | Medium | Whether a record over the 64 KiB limit is accepted depends on its byte position and the input's framing | Fixed |
| A4 | Low | Reads over the 1,024-cycle table silently missing from every per-cycle array | Fixed |
| A5 | Low | JSON carries Q20/Q30/GC/mean-quality only as 6-digit rates, never their numerators | Fixed |
| A6 | Low | A BAM QUAL value above 93 is reported as a "malformed FASTQ record" | Fixed |
| A7 | Low | CR-only line endings reported as "reads longer than this build supports", or as a truncated record | Fixed |
| A8 | Info | The overflow comment in `phred_calculator.hpp` argued from a typical chunk, not the worst case | Fixed |

### A1 — Phred+64 read as Phred+33 (High)

Nothing in the module knew any offset but 33. The same 2,000 reads, encoded
both ways:

| | Phred+33 | Same reads, Phred+64 (baseline) |
|---|---|---|
| mean quality | Q19.85 | **Q50.85** |
| Q30 rate | 0.265 | **1.000** |
| quality-trimmed reads | 2,000 | **0** |
| reads kept | 52 | **2,000** |
| exit | 0 | **0** |

**Fix.**
- `--phred-offset 33|64`. At 64, each quality byte is validated against
  `'@'..'~'` and rebased to +33 in place, in the worker, before anything reads
  it. Every kernel, filter and report then sees +33, and the output FASTQ is
  +33, as with fastp's `--phred64`. A byte below `'@'` (Solexa, or a +33 file
  declared wrong) is refused, not wrapped.
- When the offset is not given, `looks_like_phred64()` reads the run's own
  quality histogram at no per-base cost. Phred+64 read at +33 has three
  properties together: nothing below Q31 (no byte under `'@'`), nothing above
  Q72 (no byte over `'i'`), and at least 1% of bases at Q61 or higher
  (Phred+64 Q30+). No Phred+33 short-read platform writes Q61, and HiFi's
  ceiling is Q93. The histogram was widened internally to one bin per
  printable value (Q0–Q93) so the ceiling is exact; the JSON keeps its 64-bin
  schema.
  If it fires, the report is still written and the run exits **2**, the
  repository's code for "the input's shape made the answer wrong".
  `--phred-offset 33` overrides it.
- `--phred-offset` with BAM/CRAM input is refused, because those store
  qualities as numbers.

**Revised after review.** The first detection rule was "nothing below Q31,
something above Q41". On 3,000-read Phred+33 fixtures it flagged AVITI-like
(Q31–Q50, and a Q44 plateau) and HiFi-like (Q31–Q93) runs with exit 2. A
second draft put the bulk line at Q55 and still flagged a Q42–Q60 capped
long-read profile, since 31% of its bases were at Q55+. A running per-base
maximum, tried for the Q72 ceiling, cost 4% on a single-worker run and was
replaced by the widened histogram, which costs nothing. Final rule, same
fixtures:

| Profile (3,000 reads) | first rule | final rule |
|---|---|---|
| AVITI-like Phred+33, Q31–Q50 | exit 2 | **exit 0** |
| AVITI-like Phred+33, Q44 plateau | exit 2 | **exit 0** |
| HiFi-like Phred+33, Q31–Q93 | exit 2 | **exit 0** |
| capped long-read Phred+33, Q42–Q60 | exit 2 | **exit 0** |
| NovaSeq binned Phred+33 | exit 0 | exit 0 |
| Phred+64, Illumina 1.3 (Q0–Q40) | exit 2 | exit 2 |
| Phred+64, Illumina 1.5 with `B` (Q2) | exit 2 | exit 2 |
| Phred+64, Q30–Q41 | exit 2 | exit 2 |
| Phred+64, Q2–Q20 only | exit 2 | exit 0 (stated miss) |

After the fix, undeclared Phred+64 exits 2, and `--phred-offset 64` reproduces
the Phred+33 row exactly: Q19.85, 2,000 trimmed, 52 kept, byte-identical
output.

### A2 — header and `+` lines unvalidated (Medium)

F7 validated SEQ and QUAL and stopped there. That is the sibling-line
recurrence `pto-core/CLAUDE.md` warns about. Measured on the baseline, each of
these exited 0 with the bytes copied into the output:
- `@r1\0evil`
- `@r1\e[2J\a`
- `+\x01\x02`

A NUL truncates a read name in every C string API an aligner uses, so two
distinct reads can come out under one name.

**Fix.** `valid_header_line()` allows printable ASCII plus TAB, since Illumina
and 10x headers carry spaces and tabs. It is checked per record before
anything else. A trailing CR is still stripped as a CRLF ending first. Bytes
at or above 0x80 are refused too (see Residuals).

### A3 — the record limit depended on position (Medium)

`buffer.hpp` documents a 64 KiB record limit (`kPrefix`). It was enforced only
on the partial record carried from one chunk to the next. A record over the
limit that fell wholly inside a 256 KiB chunk was never checked.

One 80,010-byte record (a 40 kbp read), identical everywhere except its
offset, on the baseline:

| offset | raw | gzip | BGZF |
|---|---|---|---|
| 0 | exit 0 | exit 0 | exit 0 |
| 100,000 | exit 0 | exit 0 | exit 0 |
| 192,144 (straddles) | exit 1 | exit 1 | exit 1 |

So a long-read file passed or failed depending on the bytes before each read.
Raw/gzip and BGZF also cut chunks at different offsets (262,144 vs 261,120), so
framing alone could flip the outcome.

**Fix.** Every record is checked in the worker against the same limit,
including its final newline. A record that passes can therefore never trip the
carry check, and the two checks agree. `test_record_limit_does_not_depend_on_position`
covers 4 offsets × 3 framings × over/under the limit.

**Contract change:** a record between 64 and 256 KiB that used to pass at a
friendly offset is now always refused. That was already the documented limit.

### A4 — per-cycle table truncation was silent (Low)

`kMaxCycles` is 1,024. With 50 reads of 2 kb, the JSON said `"cycles": 1024`,
and the per-cycle base counts summed to 51,200 of 100,000 bases, with no
indication of the gap.

**Fix.** `QcStats::reads_over_cycle_limit` (one compare per read), the JSON
fields `per_cycle.cycle_limit` and `per_cycle.reads_longer_than_cycle_limit`,
and a text-report line when the count is non-zero.

### A5 — rates without numerators (Low)

The summary printed `mean_quality`, `q20_rate`, `q30_rate` and `gc_rate` with
`%.6g`. Their counts cannot be recovered from six significant digits, so
shard results cannot be merged exactly. That is the reason
`genomic_toolkit`'s frip report carries `bases_in_peaks` beside `base_frip`.

A downstream merge across shards needs those counters to recompute a rate
correctly rather than averaging six-digit rates, which is not the same number.

**Fix.** Additive fields `quality_sum`, `q20_bases`, `q30_bases` and `gc_bases`.
No field was renamed. `input.phred_offset` was added as well.

### A6 — BAM QUAL > 93 (Low)

A BAM with QUAL 94 at one base failed closed, but as
`malformed FASTQ record '@q1': sequence or quality line contains a byte outside printable ASCII`.
The input was a BAM the user never converted by hand.

**Fix.** `HtsInput::decode_qual` names the read, the value and the base:
`read 'q1': QUAL value 94 at base 8 is outside the SAM range 0..93`. 93 still
round-trips as `~` (`test_qual_past_93_is_refused_by_name`).

### A7 — CR-only line endings (Low)

A file with classic Mac line endings has no `\n`. Over a chunk it was refused
as "a FASTQ record exceeds the 65536 byte carry limit; input has reads longer
than this build supports". Under a chunk it was refused as "input ends
mid-record".

**Fix.** Both paths now say it looks like CR-only line endings and give the
`tr` command.

### A8 — the overflow bound's comment (Info)

The header claimed a 256 KiB chunk "holds ~650 reads", which is a 150 bp
chunk. The worst case is the smallest record: 8 bytes for one base, in a chunk
of at most 320 KiB with its prefix. That is ~43,000 records, and a per-cycle
quality sum below 5.4e6, still far inside `uint32`. The conclusion was right;
the argument is corrected.

## What did not break

- **Mutation fuzz.** 900 mutants of raw, gzip and BGZF input (bit flips, byte
  replacement, deletion, insertion, duplication, truncation) through the ASan +
  UBSan no-recover binary at 1, 3 and 8 threads, each under a 30 s timeout:
  - No sanitizer report, crash or hang.
  - No compressed mutant that decoded to different data.
  - Every corrupted gzip (287) and BGZF (311) mutant exited 1.
  - Raw mutants: 283 exited 1 and 19 exited 0, as mutated text can still be
    valid FASTQ.
- **Truncation at block boundaries** (G2's lesson, i.e. strict EOF, applied as
  a sweep):

  | Framing | Cut at | Outcome |
  |---|---|---|
  | BGZF, record-aligned blocks | all 80 member boundaries | 80/80 error |
  | BGZF, 4 KiB blocks | all 245 member boundaries | 245/245 error |
  | BGZF, 4 KiB blocks | 150 random offsets | 150/150 error |
  | gzip, single member | 150 random offsets | 150/150 error |
  | gzip, 8 members | 150 random offsets | 150/150 error |
  | gzip, 8 members | the 7 member boundaries | **7/7 silent, fewer reads** (see Residuals) |
  | raw FASTQ | record boundaries | silent (see Residuals) |
  | raw FASTQ | 150 random offsets | 149 error; 1 accepted, where the cut removed only the final newline (F11, correct) |

  `test_bgzf_cut_at_every_member_boundary_fails` pins the BGZF sweep.
- **CRLF.** CRLF and mixed LF/CRLF inputs across several chunks gave
  byte-identical output and identical `summary`, `filtering`,
  `quality_histogram` and `per_cycle` JSON, with adapter and quality trimming
  on. Pinned by `test_crlf_and_mixed_line_endings_match_lf`.
- **Accumulators at depth.** 32,000,000 records (64,000,000 bases, 320 MB)
  through a pipe at 12 threads, built from four record types with hand-computed
  totals. Every count matched exactly: `reads_in`, `bases_in`, `n_bases`, all
  64 histogram bins, per-cycle base counts, per-cycle mean quality, and the GC,
  Q20, Q30 and mean-quality rates. Pinned at 1.6M records by
  `test_accumulators_are_exact_over_many_chunks`. Every total is `uint64`, and
  the only 32-bit counters are per chunk (A8's bound).
- **Sequence/quality length mismatch, invalid SEQ/QUAL bytes, truncation
  mid-record.** Still refused, as F7/F11 left them.

## Residuals

- **Plain gzip and raw FASTQ have no end-of-stream marker.** A multi-member
  gzip cut exactly between members, or a raw file cut exactly between records,
  is a valid shorter input. No parser can tell it from a complete one. BGZF's
  EOF block is what makes F2's check possible, and nothing equivalent exists
  for these two. Only an external checksum or manifest would catch it.
- **Phred+64 detection is a heuristic.**
  - False positive: Phred+33 data confined to Q31–Q72 with at least 1% of
    bases at Q61+ and nothing past Q72. No short-read platform or HiFi run
    fits; a quality-filtered subset might. `--phred-offset 33` is the override.
  - False negative: a Phred+64 file with no base at Q30 or above (+64). It
    never reaches Q61 when read at +33 and is not detected. The unit test pins
    this miss rather than hiding it.
  - It fires after the run: the output FASTQ has already been written with
    +33 trimming when exit 2 arrives.
- **Records over 64 KiB (reads over ~32 kbp) are unsupported, now
  consistently.** Raising the limit means raising `kPrefix` (buffer memory per
  lane). That is a capability decision, not an audit fix.
- **Header bytes at or above 0x80 are refused.** A UTF-8 sample name in a read
  header now fails the run. FASTQ headers are ASCII in every sequencer and SRA
  output seen, but this is stricter than some tools.
- **Accepted, deliberately:**
  - A `+` line naming a different read (`+r2` under `@r1`), as seqtk and fastp
    also accept.
  - Empty records under an explicit `-l 0`.
- **Not done, though planned:**
  - A throughput probe for BGZF with very small blocks. The reader batches
    4 members per chunk regardless of their size.
  - TSan was not rerun. The changes add no shared state (per-record,
    worker-local logic only), but that is an argument, not a run.
  - No new code is vectorised, so the x86 kernels were not re-executed.

## Verification

| Check | Result |
|---|---|
| Baseline, ASan + UBSan no-recover ctest | 3/3 |
| After, Release ctest (trimmer, CLI, hts_input) | 3/3; CLI 39/39; 0 first-party warnings |
| After, ASan + UBSan no-recover ctest | 3/3; `test_trimmer` 39,657/39,657 |
| New `test_cli.sh` against the baseline binary | 7 failures (A1 ×4, A2, A7, A3), the pre-fix evidence. The new unit tests use APIs the baseline lacks, so they cannot compile against it |
| Output vs baseline, 425 MB BGZF trim | byte-identical FASTQ; all pre-existing JSON fields identical |
| Throughput vs baseline, 7 interleaved runs (final build) | 1 worker raw trim and qc-only: 0.500 → 0.500 s; BGZF -t 12: 0.300 → 0.300 s |
| Phred+64 detector, 9 platform profiles | table under A1 |
| `shellcheck -S style tests/test_cli.sh` | clean |
