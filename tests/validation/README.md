# Public-data end-to-end validation

The C++ suites under `modules/*/tests/` prove each module does what its author
intended on fixtures its author wrote. This directory answers the different
question: **on real public sequencing runs, does pto-core produce the same
answers as the established tools it is meant to replace?**

Those are not the same question, and the second one is the one a reviewer
actually cares about. A synthetic FASTQ has no adapter dimers, no quality
cliff at cycle 140, and no duplicate structure; a synthetic BAM has no
soft-clipping and no unplaced contigs. Every interesting disagreement between
two tools lives in exactly that material.

```
run_validation.py list                       # what exists, what can run here
run_validation.py fetch --allow-download     # stage the inputs
run_validation.py run                        # compare against the gold standards
run_validation.py run --case gtk/sizes-vs-samtools
```

## What is compared, against what

| case | Peak Throughput Omics | gold standard | assertion |
|---|---|---|---|
| `cuttag/profile-vs-deeptools` | `cuttag_profiler profile` | `bamCoverage` + `computeMatrix` | TSS meta-profile Pearson r ≥ 0.99, Spearman ρ ≥ 0.99 |
| `cuttag/enriched-regions-vs-bedtools` | `cuttag_profiler profile` | `bedtools intersect` | ≥ 95% of the top signal decile falls in the published peaks |
| `gtk/frip-vs-bedtools` | `genomic_toolkit frip` | `samtools` + `bedtools` | **exact** in-peak fragment count |
| `gtk/sizes-vs-samtools` | `genomic_toolkit sizes` | `samtools view` TLEN | KS distance ≤ 0.02, median within 1 bp |
| `gtk/markdup-vs-picard` | `genomic_toolkit markdup` | `picard MarkDuplicates` | duplicate rate within 5% (see below) |
| `fastq/counts-vs-seqkit` | `fastq_stream --qc-only` | `seqkit stats` | **exact** read and base counts |
| `fastq/trimming-vs-fastp` | `fastq_stream` | `fastp` | surviving reads/bases within 2%, Q30 within 0.05 |
| `fastq/alignment-rate-vs-bowtie2` | `fastq_stream` | `bowtie2` | trimming must not lower the alignment rate |
| `scrna/hnsw-recall-vs-exact` | HNSW k-NN | its own exact path | mean recall ≥ 0.90 at k=15 |
| `scrna/cosine-vs-longhand` | SIMD cosine kernel | scalar Python | max abs error ≤ 1e-5 |
| `scrna/graph-agreement-vs-scanpy` | exact k-NN | `scanpy.pp.neighbors` | informational only |

Thresholds all live in one place — `pto_validation/config.py` — so the
suite's entire tolerance budget is readable on one screen. Each carries the
reason it sits where it does. Loosening one is a real decision; record why.

### Where exact equality is *not* expected, and why

Two of these deliberately do not assert equality, and the reasons are
documented rather than absorbed into a comfortable tolerance:

* **`gtk/markdup-vs-picard`.** `genomic_toolkit` keys a duplicate on the
  aligned span; Picard uses the unclipped 5′ end, and additionally does
  optical-duplicate and UMI handling that `markdup` does not attempt. The two
  therefore disagree on soft-clipped reads *by construction* — this is
  recorded in the top-level `README.md` under "Trade-offs". The 5% bound
  exists to catch that gap **changing**, so both rates are printed on every
  run whether it passes or not.
* **`cuttag/profile-vs-deeptools`.** deeptools assigns a fragment to the bin
  containing its midpoint; `cuttag_profiler`'s `reads` mode adds it to every
  bin its span touches. On 50 bp bins that is a mild smoothing difference
  worth about 0.005 of Pearson r, concentrated at the window flanks. Both
  sides run with `raw` normalization so the comparison tests the binning
  rather than two different scale factors.

Everything else is an integer coming off the same input file, and is asserted
exactly.

## Discrepancies produce diffs, not just failures

A metric that misses its threshold writes an artefact next to the result:

* `write_curve_diff` — per-bin comparison of two signal curves, **sorted by
  largest divergence first**, with both absolute and relative deltas.
* `write_interval_diff` — intervals present on exactly one side, capped at
  5,000 per side (a `chr1`-vs-`1` naming mismatch would otherwise produce a
  three-million-line file that diagnoses nothing, and the cap is stated in the
  file so a truncated report is never mistaken for a complete one).
* `write_histogram_diff` — two distributions as normalised fractions, so a
  difference in upstream filtering does not masquerade as a difference in
  shape.

Artefact paths are listed in `validation.md` and `validation.json`.

## Datasets

`datasets.json` pins four public runs. Every accession is real and was
resolved from the ENCODE and ENA APIs; the checksums are the ones those
portals publish.

| key | what | size |
|---|---|---|
| `pbmc3k` | 10x PBMC 3k, already in the repo | — |
| `cuttag-k562-igg` | CUT&Tag K562 IgG control, SRR8754612 | 48 MB |
| `encode-h3k4me3-gm12878` | ENCODE H3K4me3 ChIP-seq, ENCSR000AKA | 2.8 GB |
| `cuttag-k562-h3k4me3-deep` | CUT&Tag K562 H3K4me3, SRR8754610 | 1.5 GB |

Start with `cuttag-k562-igg`. It is 48 MB and exercises every FASTQ path; the
IgG control carries no real enrichment, which makes it a poor signal test and
an excellent QC test. `*-deep` datasets are excluded unless you pass `--all`.

H3K4me3 is the signal-comparison dataset on purpose: it is a sharp,
promoter-proximal mark, so a TSS profile has a strong, narrowly located peak.
A broad mark produces a flat curve on which two tools can differ substantially
while still correlating at 0.99 — the comparison would prove nothing.

Downloads are cached in `~/.cache/pto-validation` (override with
`PTO_VALIDATION_DATA`), outside the worktree, because these files are
gigabytes and putting them under the repo makes every `git status` slow.

### Network access is opt-in

Nothing here touches the network without `--allow-download`. These tools are
built for hosts that hold the sequencing data, and a validation harness that
silently reaches for the internet is the wrong default on exactly that
machine. Without the flag, a missing input is a clean skip carrying the
command needed to stage it.

Downloads are atomic and checksum-verified: bytes land in a `.part` file that
is renamed into place only after the digest matches. A half-written BAM that
looks complete is the failure mode that costs a day.

ENCODE files are pinned **by accession**, not by URL — ENCODE hrefs are
versioned paths that change when a file is re-released, while the accession is
the stable identity. The portal API supplies the href and md5 at fetch time.

To add a dataset with a plain URL, fetch it once, confirm it is the file you
meant, then `run_validation.py fetch --record-checksums` to pin the digest,
and commit the result. A checksum mismatch on a later run is reported and the
pin is **left alone** — overwriting it is precisely the mistake the pin exists
to prevent.

## Skips are not failures

The exit code is 0 unless something actually failed or errored. A host without
deeptools, or an air-gapped host with nothing staged, is an expected state;
failing there would train everyone to ignore the exit code. Instead the
summary reports the skip count prominently, so a run that quietly validated
nothing is visible rather than green.

Consequently a case that produces **no metrics at all** is reported as an
error, not a pass. Silently succeeding because every comparison was skipped is
the most dangerous outcome a validation harness can have, and it is what you
get by default if you only check for failures.

## Dependencies

Standard library only, matching the C++ suites' dependency-free rule — these
tools run on clusters that have a compiler and an interpreter and nothing
else, and the harness that proves them correct has to run in the same place.
Pearson, Spearman, Kolmogorov–Smirnov and the interval Jaccard are all written
out longhand in `metrics.py`.

The `scrna/*` cases are the exception: they need `numpy`, and `anndata` to
read the `.h5ad` fixture. They skip cleanly when those are absent rather than
making them a suite-wide requirement. Point `PYTHONPATH` at the built
extension:

```bash
PYTHONPATH=build/modules/scrna_matrix python tests/validation/run_validation.py run --datasets pbmc3k
```

Reference tools (`samtools`, `bedtools`, `bowtie2`, `picard`, `deeptools`,
`seqkit`, `fastp`) are discovered on `PATH`. Each missing one skips the cases
that need it and nothing else.

## Adding a case

```python
@case("gtk/my-check", module="genomic_toolkit", reference=["samtools"])
def my_check(ctx: Context) -> list[Metric]:
    """One-line summary, shown by `list`."""
    binary = ctx.binary("genomic_toolkit")
    bam = ctx.input("bam")            # stages it, or raises Unavailable
    ...
    return [Metric("my.metric", observed, threshold, Direction.MAX)]
```

Cases return metrics; they do not assert, print, or exit. Raise `Unavailable`
to skip. Anything else that escapes becomes an error carrying its traceback,
and the remaining cases still run — a five-hour suite that aborts on the first
missing index file is a suite nobody runs twice.

Register the case name in the relevant dataset's `cases` list in
`datasets.json`.
