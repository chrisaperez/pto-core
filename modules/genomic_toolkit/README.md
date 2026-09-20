# genomic_toolkit

Fragment-level QC for ATAC-seq, CUT&Tag and CUT&RUN: size distribution,
duplicate rate and FRiP in a single streaming pass, with no intermediate files
and no sort.

```bash
genomic_toolkit frip sample.bam --peaks peaks.bed --min-length 100
```

The shell pipeline this replaces reads the data three times, writes two
intermediate BEDs, and spends most of its wall clock in a name sort and a JVM
heap:

```bash
samtools sort -n sample.bam -o ns.bam                 # spills to /tmp
bedtools bamtobed -bedpe -i ns.bam | awk '...' \
  | sort -k1,1 -k2,2n > frags.bed                     # 8 GB of text at scale
bedtools intersect -u -a frags.bed -b peaks.bed | wc -l
java -jar picard.jar MarkDuplicates I=sample.bam ...  # another full pass
```

That pipeline is the structural critique in miniature. Nothing in the
computation requires a sort or a temporary file. Both exist because the stages
are separate processes that can only speak to each other through the filesystem,
and the memory profile of the whole thing is inherited from that arrangement
rather than designed.

## What it does

| Subcommand | Replaces | Output |
|:--|:--|:--|
| `sizes` | `bedtools bamtobed` plus `awk` plus an R histogram | mean, median, percentiles, sub-nucleosomal and mononucleosome fractions, optional TSV histogram |
| `markdup` | Picard `MarkDuplicates`, for counting only | duplicate rate, distinct fragments, Lander-Waterman library-size estimate |
| `frip` | `bedtools intersect -u \| wc -l` | fragment FRiP, base FRiP, enrichment over uniform, peak-merge diagnostics |

Input is a coordinate-sorted BAM or CRAM, a BEDPE, or a fragment BED or TSV with
three or more columns. Each fragment-BED line is one fragment. `.gz` works
wherever the build has htslib.

## Measured

Against the repository's `cuttag_profiler/data/demo.bam` (800k reads, 400k
pairs, 2,000 promoter peaks) on an M-series laptop:

| | Wall clock | Intermediate bytes written | Fragments in peaks |
|:--|--:|--:|--:|
| `samtools sort -n` + `bamtobed` + `sort` + `intersect` + `wc` | 1.42 s | 15.9 MB | 148,680 |
| `genomic_toolkit frip --keep-dups` | 0.05 s | 0 | 148,680 |

The counts are exactly equal, which is the only reason the 28x is worth quoting.
A speedup over a subtly different computation is worthless.
`scripts/validate_against_bedtools.sh` reproduces both sides, and it is the
check to run before believing any change to `peak_set.hpp`.

The repository-level harness measures the same subcommand against a fuller
reference pipeline at a larger tier: 38.80x at 23.5 MB peak RSS against
1,681.2 MB, with fragment and in-peak counts exact. See the
[root README](../../README.md).

## How it avoids the sort

Three observations, one per bottleneck.

**Fragments never become text.** Everything downstream of the reader speaks a
16-byte coordinate-only `Fragment` struct from `include/toolkit/types.hpp`,
handed to consumers in 4,096-record batches. That is 64 KiB per batch, and
therefore L2-resident, invalidated as soon as the callback returns. All three
consumers run off the same batch, so `frip` marks duplicates and accumulates the
size histogram in the pass it was going to make anyway.

**Duplicate detection needs no global state on sorted input.** Two fragments are
duplicates when their 5' ends agree: same contig, start, end and orientation. On
a coordinate-sorted file every fragment sharing a start arrives contiguously, so
the only state required is a set keyed on (end, orientation) that lives for
exactly one start position. This means memory tracks the depth of the deepest
pile-up, which is 3 entries on the demo data and a few hundred on real coverage,
and never the size of the file. There is nothing to spill. That is the entire
reason this does not need Picard's heap.

If the input is not coordinate-sorted, `--unsorted` switches to a 64-bit
fingerprint table over all distinct fragments, open-addressed with linear
probing. Memory becomes O(distinct fragments). It exists so that unsorted input
is a slower answer rather than no answer.

Out-of-order input is detected in the default mode and exits 2. That includes
the interleaved-contig case, which has ascending coordinates throughout and
passes a naive monotonicity check while silently undercounting. The check holds
for `markdup` and, since 2026-09-11, for `frip` whenever it is marking
duplicates. Before then `frip` ran the same marker, did not look at its verdict,
and counted the duplicates it had missed, at exit 0. `frip --keep-dups` needs no
order and is unaffected.

**FRiP is a reduction, not an intersection.** Peaks are merged into a disjoint
ascending cover at load time and queried by binary search, in
`include/toolkit/peak_set.hpp`. Merging is not an optimisation. It is the
definition of the statistic: a fragment spanning three overlapping peaks counts
once, where the naive `bedtools intersect -wa | wc -l` over an unmerged peak file
counts it three times and can report a FRiP above 1.

This module has its own interval structure rather than reusing
`cuttag_profiler`'s, and not only because the repository has no shared layer. It
is a different structure for a different question. The profiler's implicit
interval tree enumerates the set of annotations overlapping a window. FRiP needs
only whether a fragment touches anything and by how many bases, which a merged
cover answers in one binary search without materializing hits.

## Trade-offs and things it is not

Stated plainly, because each of these changes a number someone may publish.

* **`markdup` is not a drop-in `MarkDuplicates`.** Identity uses the aligned
  fragment span where Picard uses the unclipped 5' position, so a duplicate pair
  whose copies were soft-clipped differently is one duplicate to Picard and two
  distinct fragments here. Rates agree closely for ATAC, CUT&Tag and WGS, where
  clipping is rare. They will not agree for primer-clipped amplicon panels.
* **No optical-duplicate distinction and no UMI awareness.** Read names are not
  parsed for tile coordinates, so `estimate_library_size()` treats every
  duplicate as a PCR duplicate and is therefore a lower bound on complexity.
* **It reports, and it does not rewrite.** No marked BAM is emitted. The
  representative kept at a position is the first seen rather than the
  highest-quality one, which is immaterial for counting and not for variant
  calling, so the tool declines to produce an output that would invite the
  latter.
* **int32 coordinates.** Halving the hot arrays is the point of the layout.
  Contigs at or beyond 2^31 bp, which includes some conifer and amphibian
  assemblies, are rejected by name rather than truncated.
* **The size histogram tracks up to 2,000 bp.** Longer fragments go to an
  overflow bucket and are reported separately. Percentiles that fall in it
  saturate at 2,000 rather than being invented, and the mean uses true lengths.
* **Single-end input is opt-in** via `--single-end`, and its fragments are
  flagged as inferred spans, because a read length is not an insert size.
* **FRiP with peaks called on the same sample is circular.** This tool cannot
  detect that. Neither can bedtools.

## Input expectations and failure modes

Input this tool cannot trust is refused rather than summarised. Each of these is
a case where the previous behaviour was to produce a number.

| Input | What happens |
|:--|:--|
| A BGZF file (BAM or `.gz`) with no end-of-file marker | Exit 1 before reporting anything. Cut at a block boundary it decodes into a clean prefix, which is how a truncated BAM once reported half a library at exit 0. |
| A fragment file in which every record is malformed | Exit 2 after its report, rather than reading as an empty library. |
| A count column present but not an integer in 1 to 1,000,000, under `--with-counts` | That record is malformed. It is never read as one copy and never clamped. |
| A corrupt MAPQ field | Malformed record. It is not a maximum-confidence one, which is what a discarded parse result once made it. |
| `--min-length` or `--max-length` outside 0 to 2^31-1 | Refused by name, exit 1. |
| `--min-mapq` outside 0 to 255 | Refused by name, exit 1. |
| `--genome-size` below 0 | Refused by name, exit 1. |
| `--histogram` naming an input, its index, or the peak file | Refused. |
| CRAM without a reachable reference | CRAM needs `REF_PATH` pinned locally, or it hangs rather than failing. |

Details in
[`AUDIT_2026-09-11_genomic_toolkit.md`](../../docs/AUDIT_2026-09-11_genomic_toolkit.md).
This module's record is 12 confirmed defects, 4 of them High or Critical.

**`--with-counts` is the flag most likely to change your numbers.** Column 5 of
a fragment BED is a read-pair count only when you pass it. By default every
interval is one fragment and column 5 is not read at all, because in a BED5 or
BED6 that column is a *score*: until 2026-09-11 a record scored `60` counted as
60 fragments, at exit 0, multiplying every statistic in the report. A 10x
`fragments.tsv` needs `--with-counts` to count read pairs. Without it, the
statistics count unique fragments instead. Both are defensible answers. Only one
of them is the one you asked for.

## Chromosome naming

A peak file saying `1` against a BAM saying `chr1` produces a FRiP of exactly
zero with no error anywhere. It is the most common wrong number in this corner
of bioinformatics, and it is wrong in the worst way: zero is a value a reviewer
will accept.

`frip` exits 2 and names the offending contigs when nothing overlapped anything
*and* the naming does not line up. A few unmatched contigs, such as decoys,
scaffolds and chrM, are normal and are reported without being an error.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Standalone or from the top-level aggregate; both are supported and both are
tested. There are two optional dependencies and neither is a hard failure.

**htslib** (`-DGTK_WITH_HTSLIB=AUTO|ON|OFF`, default `AUTO`) enables BAM, CRAM
and `.gz` input. Without it the BEDPE and fragment-BED paths still build and the
entire test suite still passes, because every test drives the text paths on
purpose. This is a deliberate departure from `cuttag_profiler`, where htslib is
the whole premise. Here it is the one dependency that makes the module hard to
build on a locked-down cluster, and a fragment file is what an scATAC user has
anyway. `ON` turns a missing htslib into a configure error.

**OpenMP** parallelises the FRiP reduction. Without it the same loop runs
serially with bit-identical results, because every accumulator is an integer.
Unlike `scrna_matrix`, a missing OpenMP is not a hard failure here: the queries
are latency-bound binary searches that stop scaling around four threads, so
requiring the toolchain would be a fabricated dependency. There is no
`--threads` flag for the same reason; set `OMP_NUM_THREADS`.

## Tests

Dependency-free, in the house style. A bare compiler is enough.

```bash
ctest --test-dir build -R gtk_ --output-on-failure
```

| Suite | Covers |
|:--|:--|
| `gtk_peak_set` | merging (overlapping, nested, book-ended), half-open boundaries, zero-length peaks, coordinate ceiling, and a differential test against a per-base bitmap over randomised interval sets |
| `gtk_dup_marker` | both modes agreeing on the same fragment set, the O(pile-up) memory claim asserted as a number, out-of-order and contig-revisit detection, library-size estimate |
| `gtk_fragment_stream` | BEDPE pairing, filters, malformed-line accounting, CRLF, missing final newline, a line longer than the 1 MiB read buffer, coordinate ceiling, batching boundaries |
| `gtk_frip` | end-to-end counts, the overlapping-peak double-count trap, duplicate exclusion, the naming mismatch, enrichment, unsorted and overlapping and CRLF peak files, `checked_add` at the 64-bit boundaries |
| `gtk_fuzz_parsers` | seeded mutation fuzzer over BEDPE and fragment BED; the MAPQ and count columns failing closed; coordinate boundaries; degenerate files |
| `gtk_file_identity` | `--histogram` may not name an input, its index, or the peak file |
| `gtk_cli_args` | the checked integer parse behind every numeric option: partial numbers, int32 narrowing, int64 overflow, per-option ranges |
| `gtk_cli` | the binary as a black box: exit codes and refusals that live in `src/main.cpp`, plus truncated-BAM cases, which skip without htslib or samtools |
| `gtk_region_partition` | htslib only: `--region` tiles partition a file exactly |
| `gtk_bam_contig_order` | htslib only: FRiP on a multi-contig BAM does not depend on the peak file's contig order |

The peak-set differential test is the important one. Its reference states the
overlap rule as a bitmap over positions, so it compares two independent
formulations rather than one expression against a restatement of itself. A test
that restates the implementation proves only that you can write the same thing
twice.

One overflow threshold from the audit is worth knowing, because it was measured
rather than reasoned about: with a signed base total, a crafted file of 4,612
lines was enough to overflow. Unsigned doubles that to 9,224. Those correspond
to 4.7e9 and 9.3e9 fragments respectively.

`src/main.cpp`'s behaviour is tested through `scripts/test_cli.sh`, registered
as the `gtk_cli` suite, because no unit suite links it.

## License

MIT, under the root [`../../LICENSE`](../../LICENSE).
