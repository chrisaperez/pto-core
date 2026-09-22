# peaks (`pto-peaks`)

A streaming epigenomic peak caller: single-pass pileup, a multi-scale rolling
background, SIMD Poisson tail probabilities, Benjamini-Hochberg FDR control, and
narrowPeak out. No sort, no intermediate bedGraph.

```bash
pto-peaks --chrom-sizes hg38.chrom.sizes fragments.bed > peaks.narrowPeak
```

```bash
pto-peaks -o peaks.narrowPeak atac.bam
```

A BAM carries its own contig table, so `--chrom-sizes` is required only for BED
input.

**Read [Concordance with MACS3](#concordance-with-macs3) before treating this as
a drop-in replacement.** It is a real, working peak caller with a
MACS3-equivalent summit detector. It is not a numerical match for MACS3's peak
boundaries, and that gap is measured rather than estimated.

## What it does

| | |
|:--|:--|
| Input | A coordinate-sorted fragment BED (`chrom`, `start`, `end`), or a coordinate-sorted BAM or BGZF file. Format is detected from the first byte, and either works on stdin. |
| Output | narrowPeak (BED6+4), the format MACS2 writes and every genome browser reads. Column 10 is the summit offset from the peak start, which is MACS2's own convention. |
| Background | Three centred sliding windows (1 kb, 5 kb, 10 kb) over the fragment pileup, evaluated per base. There is no control or input track. |
| Significance | Poisson tail probability against the local background, then a Benjamini-Hochberg step-up over all candidates. |

A BAM yields one fragment per proper pair, spanning TLEN from the leftmost mate.
Never both mates, which would double every pileup. A gzip- or BGZF-compressed
fragment BED such as a 10x `fragments.tsv.gz` is not read directly; pipe it
through `zcat` into stdin. Records arriving out of order are refused rather than
mis-placed, and the run exits 2.

## CLI reference

```
usage: pto-peaks --chrom-sizes FILE [options] [fragments.bed]

  --chrom-sizes FILE   two columns, name and length. Required for
                       BED input; a BAM carries its own
  --min-mapq N         BAM only: skip alignments below this (default 0)
  -o, --out FILE       narrowPeak output (default: stdout). May not
                       be the input or the --chrom-sizes file
  --cutoff X           -log10 p to open a peak (default 2.0)
  --qvalue X           Benjamini-Hochberg FDR level (default 0.05)
  --min-length N       discard calls shorter than N bp (default 50)
  --max-gap N          sub-threshold bases tolerated inside a peak
                       (default 30)
  --extend-peaks N     pad every reported peak by N bp on each side
                       (default 0). 0 keeps pto-peaks' sharp,
                       high-resolution bounds; a positive value is
                       a compatibility mode for the wider spans
                       MACS2/MACS3 produce. The summit position and
                       the -log10 p / -log10 q columns are
                       unchanged; the 3' edge is clamped to the
                       contig length
  --lambda-bg X        genome-wide background, fragments per base.
                       Used as a floor under the local estimate;
                       0 (the default) means local scales only
  -h, --help           this text
```

This text is transcribed from the binary's own parser in `src/main.cpp`.
`--help` is always the authority if the two ever disagree, because this file is
plain prose and can go stale.

### Exit codes

```
0  peaks called (possibly none)
1  usage or I/O failure
2  the input's shape made the answer wrong: unsorted records, a contig
   missing from the sizes file, a malformed field. Deterministic, so
   retrying reaches the same answer; a valid report is printed first.
```

Exit 2 is never retried: any caller that retries failures automatically -- a
workflow manager, a job scheduler -- should treat it as terminal, because a
retry of a deterministic wrong-shape failure spends the job's timeout arriving
at the same answer.

### Failure modes on malformed input

Every one of these is a refusal rather than a default, and each exists because
the alternative was found to produce a plausible number:

| Input | What happens |
|:--|:--|
| Records out of coordinate order | Refused at the offending record, exit 2. |
| A contig in the data absent from `--chrom-sizes` | Refused, exit 2. |
| A contig name that is empty or duplicated, in a BED or a BAM header | Refused. The BAM header is held to the same rules as the chrom.sizes loader, so a genome is not valid or invalid depending on its file format. |
| A contig length below 1 | Refused. |
| A corrupt BGZF block | Refused on the per-block CRC32. A corrupt block that inflates into plausible alignments is exactly the silent wrong answer this module exists to refuse. |
| A BGZF stream missing its end-of-file marker | Refused as a truncation, rather than treated as a clean EOF. |
| An empty BGZF block mid-stream | Transparent, and this is deliberate. Empty blocks are legal anywhere, the EOF marker is one, and concatenating BGZF streams puts them mid-file. Stopping at the first one dropped every record after it and reported success. |
| A degenerate background lambda (zero, negative, NaN, infinite) | Scores p = 1. See below. |
| An output path naming an input | Refused, including the chrom-sizes and index files. |

The degenerate-lambda rule is the one worth understanding rather than just
obeying. The naive Poisson evaluation tends toward p = 0 as lambda falls, so it
reports every centromere, assembly gap and unmeasured window as the most
significant peak in the genome. A caller that ranks assembly gaps above real
biology is not slightly miscalibrated. It is inverted.

### `--extend-peaks`: matching MACS-style spans on request

`pto-peaks` reports sharp, high-resolution intervals by default. A positive
value pads every reported interval symmetrically and changes nothing else: the
summit position and the `-log10 p` and `-log10 q` columns are byte-identical to
a sharp run, padding is applied after the `--min-length` gate, and the 3' edge
is clamped to the contig length.

On the ENCODE fixture in [`CONCORDANCE.md`](CONCORDANCE.md), `--extend-peaks 90`
takes the median call width from 96 bp to about 276 bp, landing on MACS3's own
275 bp median. It does not change the caller's judgment about where a peak is,
only how wide the reported box around it is.

## Concordance with MACS3

**The 0.98 parity gates this project set for itself are not met, and I am not
going to gloss over that here.** The full record, methodology and reproduction
steps are in [`CONCORDANCE.md`](CONCORDANCE.md). Measured against real MACS3
3.0.4 on real ENCODE K562 ATAC-seq (`ENCSR217QAB`, 2,889,130 fragments):

| Metric | Value | Reading |
|:--|:--|:--|
| Median summit distance, matched peaks | **1 bp** | MACS3-equivalent |
| Summit within 50 bp, matched peaks | 97.6% | MACS3-equivalent |
| Pearson r on peak strength (-log10 p) | **0.990** | MACS3-equivalent |
| Peak count | 31,835 against MACS3's 31,607 | agree to 0.7% |
| Summit within 50 bp, all matched pairs | 0.9616 | short of the 0.98 gate |
| Spearman r on -log10 p | 0.8568 | short of the 0.98 gate |
| Mean IoU | **0.458** | short, and the reason is below |
| Median peak width | 96 bp against MACS3's 275 bp | **roughly 3x narrower** |

`pto-peaks` finds the same peaks, in the same places, ranked in the same order,
with a summit detector that agrees with MACS3 to within a base pair. What it
does not do is match MACS3's interval width. It draws systematically tighter
boundaries, and that single difference is what drags every interval-overlap
metric below the bar: a 96 bp call sitting inside a 275 bp call is about 0.35
IoU at best, however well the two agree about where the peak is.

Since 2026-09-11 `make validate-peaks` exits 0 against regression floors of 0.95
on summit agreement and 0.85 on strength correlation, chosen just under the
measured values as a drift tripwire. That exit code is not a parity claim. The
0.98 gates are still computed, printed, and reported as not met.

**Practical guidance.** Treat `pto-peaks` output as summit-accurate and
rank-accurate, and as boundary-tight rather than a numerical drop-in for MACS3
narrowPeak calls. If your downstream analysis depends on MACS-style interval
widths, for instance comparing peak widths across tools or using feature sets
tuned to MACS's typical span, use `--extend-peaks` or keep MACS3 for that step.
If it depends on summit location or relative peak strength, the two callers
agree closely.

## Known limits

Each of these is a capability boundary I am stating rather than a defect.

* **No `--broad` mode**, and no two-level FDR stitching. Histone domains such as
  H3K27me3 and H3K36me3 have nothing to swap to here.
* **Paired-end only.** A single-end BAM has TLEN 0 on every record, so it yields
  no fragments and therefore **no peaks at exit 0**. The count shows on stderr
  as "not the leftmost mate". Convert to a fragment BED with an intended
  fragment length first.
* **No control or input channel.** Enrichment is scored against the local
  multi-scale background only, which is closer to MACS2 `--nolambda` than to its
  treatment-versus-control default. A pipeline currently passing `-c input.bam`
  is not doing the same comparison here.
* **Peak splitting is minor and one-directional.** `pto-peaks` occasionally
  reports one MACS3 peak as two. It never merges two MACS3 peaks into one.
* **AVX-512 has never run on hardware that natively supports it.** The kernel
  compiles and is numerically correct against a scalar reference. No throughput
  claim is made for it.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Standalone or from the top-level aggregate; both are supported and tested. There
is no hard dependency, and the fragment-BED path plus the whole test suite work
with no optional dependency present.

| Option | Default | Effect |
|:--|:--|:--|
| `PEAKS_WITH_BGZF` | `AUTO` | BAM and BGZF input: `libdeflate` if found, else `zlib`, else compiled out with the fragment-BED path and tests unaffected. `ON` turns a missing decompressor into a configure error; `OFF` skips both discovery calls. Takes no htslib dependency at all, because BAM records are decoded directly. |
| `PEAKS_ENABLE_AVX512` | `OFF` | Compile the AVX-512 Poisson kernel. Never enable this for a binary that will run on a machine other than the one that built it: an unsupported instruction is a SIGILL partway through someone's batch. |
| `PEAKS_BUILD_TESTS` | `ON` | Build the six test binaries below. |

## Tests

Dependency-free, in the house style, with no GoogleTest and no Catch2.

```bash
ctest --test-dir build -R peaks_ --output-on-failure
```

| Suite | Covers |
|:--|:--|
| `peaks_poisson_model` | Poisson tails against mpmath golden values, the degenerate-background fail-closed rule, the finite-lower-sum fast-path gate, convergence at the iteration ceiling |
| `peaks_sliding_window` | the centred multi-scale background rings, contig-end truncation, the coverage-desert skip, and back pressure |
| `peaks_caller` | the state machine, the float-comparison summit-plateau fix, Benjamini-Hochberg's cutoff-relative floor, batch-size invariance |
| `peaks_bam_streamer` | native BGZF and BAM decoding by byte shift, per-block CRC32, an empty BGZF block mid-stream staying transparent |
| `peaks_file_identity` | an output path may not name an input, including chrom-sizes and index files |
| `peaks_cli` | the binary as a black box: argument parsing, the BED reader's fail-closed checks, the narrowPeak writer, and the exit codes the job contract depends on |

Golden Poisson values come from mpmath at 80 digits via
`scripts/gen_poisson_golden.py`, which is not a build or test dependency. The
values are checked in and the suite needs only a compiler.

Why the state machine cannot be vectorized past the Poisson evaluation is
documented in the header comments of
[`include/peaks/caller.hpp`](include/peaks/caller.hpp); the centred
multi-scale background windows are documented the same way in
[`include/peaks/sliding_window.hpp`](include/peaks/sliding_window.hpp) — this
module has no separate architecture document.

## License

MIT, under the root [`../../LICENSE`](../../LICENSE). This module takes no
htslib dependency and vendors nothing, so there is no third-party grant to
reproduce here. See
[`../../THIRD_PARTY_LICENSES.md`](../../THIRD_PARTY_LICENSES.md) for what the
other modules vendor.
