# pto-core

Cache-aware C++20 engines for FASTQ QC, fragment QC and CUT&Tag/CUT&RUN signal
profiling, with a Python API over them.

These tools are built to run **on the machine that holds the sequencing data**:
no service dependency, no data leaving the host, no intermediate files written
to shared scratch. This package is the zero-friction way to get them onto a
login node or a laptop and produce a real number.

```bash
pip install pto-core
pto-core doctor  # which engines this wheel resolves, where from, and their versions
```

```python
import pto

qc = pto.run_qc("sample_R1.fastq.gz")
print(qc["summary"]["reads_in"], qc["summary"]["q30_rate"])
print(f"{qc.wall_seconds:.2f}s, peak RSS {qc.peak_rss_mb} MB")

sizes = pto.fragment_sizes("atac.bam")
enrich = pto.frip("atac.bam", "peaks.narrowPeak", genome_size=2_913_022_398)
print(sizes["median"], enrich["frip"])
```

Every call returns a `Result`: the tool's own JSON report under `.metrics`
(indexable directly, as above), whatever files it wrote under `.artifacts`, and
what the run cost under `.wall_seconds` / `.peak_rss_mb`.

## What is in the wheel

| Engine | What it does |
|---|---|
| `fastq_stream` | Streaming FASTQ QC and adapter/quality trimming; pipes clean reads into an aligner without materialising an intermediate FASTQ |
| `genomic_toolkit` | Streaming fragment QC — size distribution, duplicate rate and FRiP in one pass, no sort and no temp files. **`markdup` reports a rate; it does not write a BAM — see below.** |
| `cuttag_profiler` | CUT&Tag / CUT&RUN reference-point signal matrices from an indexed BAM |
| `pto-peaks` | Streaming peak caller — fragment BED or coordinate-sorted BAM to narrowPeak in one pass, no intermediate bedGraph. **Narrow and paired-end only; see below.** |

### `pto-peaks` is narrow-only, paired-end-only, and not yet MACS-validated

Worth knowing before it replaces anything:

- **No `--broad`.** Histone domains (H3K27me3, H3K36me3) have no mode here.
- **Paired-end only.** One fragment per proper pair, spanning TLEN from the
  leftmost mate. A single-end BAM has TLEN 0 on every record, so it produces no
  fragments and therefore **no peaks, at exit 0** — the count shows as `not the
  leftmost mate` on `.stderr`. Convert to a fragment BED first.
- **No control/input.** Enrichment is scored against the local multi-scale
  background (1kb/5kb/10kb) only, which is nearer MACS2 `--nolambda` than its
  treatment-versus-control default.
- **Measured against MACS3, and the concordance gate is not yet met.** The
  comparison harness (`make validate-peaks`) has been run against real ENCODE
  K562 ATAC-seq data (2026-09-10): summit placement and significance ranking
  are MACS3-equivalent (median summit offset 1 bp, r=0.99 on signal strength),
  but reported intervals run roughly 3x narrower than MACS3's defaults, for a
  mean IoU of 0.46 against this project's own 0.98 gate. `--extend-peaks N` is
  an opt-in width-matching mode for pipelines built around MACS-style spans.
  Full numbers: [`modules/peaks/CONCORDANCE.md`](../modules/peaks/CONCORDANCE.md).

### `genomic_toolkit markdup` reports a rate; it does not write a deduplicated BAM

`pto.mark_duplicates(...)` replaces Picard `MarkDuplicates`' **metrics file**,
not the step that produces a `*.markdup.bam`. Keep whatever you already use to
deduplicate reads; this gives you the same kind of numbers, faster, in the
same streaming pass as the rest of `genomic_toolkit`. Two differences are
deliberate, not bugs:

- **Duplicate identity uses the aligned fragment span**, where Picard uses the
  unclipped 5' end. A pair whose copies were soft-clipped differently is one
  duplicate to Picard and two here. Rates agree closely for ATAC/CUT&Tag/WGS,
  where clipping is rare; they will not for primer-clipped amplicon panels.
- **No optical-duplicate distinction and no UMI awareness.** Read names are not
  parsed for tile coordinates, so the library-size estimate treats every
  duplicate as a PCR duplicate and is therefore a lower bound on true library
  complexity.

The rate is comparable to Picard's on the same library, not identical to it.
Treat a swap as a change of method and say so in your methods section. Full
detail: [`modules/genomic_toolkit/README.md`](../modules/genomic_toolkit/README.md#trade-offs-and-things-it-is-not).

`scrna_matrix` — the cache-aligned sparse-matrix engine and SIMD k-NN — is a
compiled extension module rather than a CLI and is published separately:

```bash
pip install scrna-matrix
```

It is not bundled here because its configure step is a **hard failure** without
OpenMP, by design. Vendoring it would mean a macOS host without `libomp` could
not install a wheel it only wanted a FRiP number from.

## No runtime dependencies

`pip install pto-core` pulls in nothing. The API returns the engines' own JSON as
plain dicts, so there is no validation library to resolve against whatever numpy
or pydantic a core facility's shared interpreter is already pinned to.

## Python API

| Function | Engine |
|---|---|
| `pto.run_qc(reads, ...)` | `fastq_stream --qc-only` |
| `pto.trim(reads, out, ...)` | `fastq_stream` |
| `pto.fragment_sizes(alignment, ...)` | `genomic_toolkit sizes` |
| `pto.mark_duplicates(alignment, ...)` | `genomic_toolkit markdup` |
| `pto.frip(alignment, peaks, ...)` | `genomic_toolkit frip` |
| `pto.profile(alignment, regions, out_dir, ...)` | `cuttag_profiler profile` |
| `pto.call_peaks(fragments, output_path, ...)` | `pto-peaks` |

**Every keyword argument defaults to `None`, and a `None` is not passed to the
binary at all** — the engine applies its own documented default. That is
deliberate: the ranges and defaults live in the CLI parsers, and a copy of them
in these signatures would be one more place to drift out of step. `--help` on
any engine is the authority, and `pto-core exec genomic_toolkit --help` reaches
it without hunting for the binary.

### Exit codes are part of the API

`genomic_toolkit` exits 2 when the *input's own shape* makes the answer wrong,
and prints a valid report first. That surfaces as `pto.InvalidInput`, with the
report attached:

```python
try:
    result = pto.frip("atac.bam", "peaks.bed")
except pto.InvalidInput as exc:
    # No fragment overlapped a peak AND the contigs are absent from the peak
    # file -- almost always `chr1` vs `1` naming. Deterministic: retrying
    # reaches the same answer. The report is still worth showing.
    print(exc.result.metrics)
```

A genuine failure — unreadable input, a rejected option — raises `pto.ToolError`
instead. Retry that; never retry `InvalidInput`.

### Threads

Thread count is always passed explicitly, from `pto.available_cpus()`, which
takes the **minimum** of CPU affinity, the cgroup CPU quota and `os.cpu_count()`.
The engines otherwise default to `std::thread::hardware_concurrency()`, which
reports the machine — so a SLURM task granted `--cpus-per-task=4` on a 128-core
node would start 128 workers and spend its life in scheduler contention. Pass
`threads=` to override.

## Command line

```bash
pto-core doctor                  # which engines were found, where, and their versions
pto-core man                     # every command, one line each
pto-core man frip                # ...or one command's required/optional parameters
pto-core qc reads.fastq.gz       # a real answer in one command
pto-core trim reads.fastq.gz clean.fastq.gz --min-mean-q 20
pto-core sizes atac.bam --json
pto-core markdup atac.bam --json
pto-core frip atac.bam peaks.bed --genome-size 2913022398
pto-core profile atac.bam tss.bed out/ --reference-point TSS
pto-core peaks frags.bed peaks.narrowPeak --chrom-sizes hg38.chrom.sizes
pto-core peaks atac.bam peaks.narrowPeak      # a BAM carries its own contigs
pto-core exec genomic_toolkit --help    # hand off to an engine directly
```

`pto-core man` (alias `pto-core glossary`) reads its output off this package's
own `argparse` parser rather than a hand-maintained copy, so a flag it reports
is a flag that actually exists.

The script is `pto-core`, not `pto`: `pto` belongs to the `pto-cli` distribution,
the client for the managed service.

### Putting the engines on `$PATH`

For Nextflow, Snakemake, or any shell pipeline that wants the bare tool names
with no Python in the call path:

```bash
export PATH="$(pto-core bin-dir):$PATH"
```

The binaries are vendored *inside* the package rather than installed onto the
environment's `bin/` on purpose — `fastq_stream` and `genomic_toolkit` are
generic enough names that shadowing a site's existing, differently-built copy
would be a real hazard. This makes it opt-in and explicit.

## Platforms

| Platform | Wheel | htslib |
|---|---|---|
| Linux x86_64 / aarch64 | `manylinux_2_28` | built from source, **statically** linked |
| macOS arm64 (14+) | native | Homebrew, vendored by `delocate` |
| macOS x86_64 (13.3+) | native | Homebrew, vendored by `delocate` |

**There is no `universal2` macOS wheel, and that is a constraint rather than an
omission.** Homebrew's htslib is per-architecture, so an arm64 runner cannot link
an x86_64 `libhts.dylib`; the cross build either fails at link time or produces a
wheel `dyld` rejects on the user's machine. Both architectures are published as
native wheels instead, built one runner per arch — the same resolution
`wheels.yml` already uses for `libomp` and `scrna-matrix`.

Linux wheels link htslib statically and carry no `libhts`-family shared-object
dependency, which is the configuration `ci.yml`'s `static-link` job proves and
asserts. `libcurl` is disabled in that build: it is what drags nghttp2 and brotli
into a static link, and it is also htslib's remote CRAM reference fetcher, which
behind a proxy or on an air-gapped node hangs rather than failing.

## Building from source

The sdist does not ship prebuilt binaries — it compiles pto-core with CMake:

```bash
pip install pto-core --no-binary pto-core
```

That needs a C++20 compiler, CMake ≥ 3.20, and **htslib** (`pkg-config`
discoverable). htslib is required because `cuttag_profiler` cannot be built
without it. To build the other two engines on a host with no htslib:

```bash
pip install ./pto-core -C cmake.define.PTO_BUILD_CUTTAG_PROFILER=OFF \
                       -C cmake.define.FQ_ENABLE_HTS=OFF
```

### Pointing at a tree you built yourself

A site that compiles pto-core against its own htslib does not need to reinstall
anything:

```bash
cmake -S pto-core -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
export PTO_BIN_DIR=$PWD/build/modules/genomic_toolkit
```

`PTO_BIN_DIR` is searched before the vendored directory and before `$PATH`.

## Reproducibility

Put both numbers in a methods section — the wrapper and the engines are versioned
separately:

```bash
pto-core versions
```

```python
import pto
pto.__version__        # this distribution
pto.tool_versions()    # {'fastq_stream': 'fastq_stream 0.1.0', ...}
```

## The dashboard is not here

`cuttag_profiler serve` — the interactive HTML dashboard — is deliberately not
reachable from this package. Its authentication is a session token bound to one
operator's browser, which is sound on a workstation and the wrong shape on a
shared cluster node. Use `pto.profile(...)` and render `profile.tsv`.

## Licence

MIT.
