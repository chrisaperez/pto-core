# cuttag_profiler

Air-gapped CUT&Tag and CUT&RUN signal profiling, with a dashboard that lives
inside the binary.

One statically linked C++20 executable reads an indexed BAM, computes a
reference-point signal matrix over your loci, and, if you want to look at it,
serves an interactive dashboard from an embedded HTTP server bound to
`127.0.0.1`. No Python environment, no coverage-track intermediate, no bytes
leaving the machine.

```bash
cuttag_profiler profile --bam sample.bam --regions promoters.bed \
    --upstream 2000 --downstream 2000 --bin-size 50 \
    --out-matrix matrix.tsv --out-profile profile.tsv
```

```bash
cuttag_profiler serve --data-root /data/experiments
```

## Why this exists

Quantifying normalised coverage over a large locus set is the read-out step for
almost every CUT&Tag experiment. The standard route through deeptools needs two
stages, `bamCoverage` then `computeMatrix`, and materializes a genome-wide
coverage track even when you care about 0.1% of the genome.

That is the structural problem rather than a speed complaint. The cost of the
deeptools route scales with the size of the genome, because it builds a
genome-wide object first. The cost of a direct indexed query scales with the
loci you asked about. In other words, the intermediate is not a slow step on the
way to the answer. It is a different and much larger computation that happens to
contain the answer.

The interactive tools that fix the usability half are mostly hosted services,
and plenty of labs cannot send aligned human sequence data to someone else's
infrastructure. This tool takes the other path: all of the interactivity, none
of the network.

## Measured performance

100,000 promoter windows at plus and minus 2 kb, 50 bp bins, 8M alignments in an
87 MB BAM, 12 threads, mean of 3 replicates on an Apple M-series workstation:

| Stage | Time | Peak RSS |
|:--|--:|--:|
| **`cuttag_profiler`**, BAM to matrix | **1.42 s** | **84.7 MB** |
| `deeptools bamCoverage`, BAM to bigWig | 7.11 s | 320.7 MB |
| `deeptools computeMatrix`, bigWig to matrix | 43.06 s | 432.4 MB |
| deeptools pipeline total | 50.17 s | 432.4 MB |

That is 35.4x faster end to end and 30.4x against `computeMatrix` alone, at 5.1x
lower peak memory. `computeMatrix` consumes a bigWig rather than a BAM, so the
comparable unit of work is the two-stage pipeline an analyst actually runs. Both
stages are timed, and the stage-only figure is given for anyone who already has
a cached coverage track.

The repository-level harness measures a different case on a different tier and
reports 38.18x at 89.3 MB against 447.1 MB, with a meta-profile Pearson r of
0.999977. Those are separate runs in separate contexts, not two versions of one
number. Both are in the [root README](../../README.md).

Expect the advantage to narrow as your locus set approaches whole-genome
coverage, or if you reuse one cached bigWig across many locus sets. The
advantage is structural, so it disappears exactly where the structure stops
applying.

Reproduce it yourself. Every number in this file comes out of these two
commands, and nothing is quoted from memory:

```bash
python3 scripts/make_demo_data.py --out-dir data --prefix bench \
    --regions 100000 --reads 4000000 --chrom-length 250000000
scripts/benchmark.sh --bam data/bench.bam --regions data/bench_promoters.bed
```

`scripts/benchmark.sh` writes `results/summary.txt` for timing, memory and
speedup, and `results/accuracy.txt` for concordance against deeptools. It skips
the baseline arm gracefully when deeptools is not installed.

Concordance with deeptools on the same matrices: meta-profile Pearson r of
0.99999, maximum residual 0.4% of peak after fitting a single global scale
factor, per-locus median r of 0.988, with 99.8% of loci above r = 0.95. The two
tools normalise against different denominators, so the test is whether the
signal *shape* agrees up to one global scale factor.

**One honest gap.** deeptools parity has not been re-measured since the
skipped-region fix of 2026-09-11, because deeptools was not installed on the
audit machine at the time. The fix changes the reported mean only when a region
is skipped, and the benchmarked case had none, so "unchanged" is inferred from
the logic of the fix rather than independently observed.

## Build

Requires CMake 3.20 or newer, a C++20 compiler, and htslib.

```bash
brew install htslib cmake          # macOS
# apt install libhts-dev cmake     # Debian/Ubuntu

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

`cpp-httplib` and `nlohmann/json` are vendored under `third_party/` on purpose:
a build on an air-gapped machine must not reach a package registry.

### Fully removing the network stack

A stock htslib archive contains `hfile_libcurl.o` and `hfile_s3.o`, so the
linked binary carries an HTTP and S3 client even though this tool never
constructs a remote URL. If "the code path is unreachable" is not a strong
enough claim for your compliance review, rebuild htslib without it. This also
drops the libcurl and OpenSSL link dependencies entirely:

```bash
scripts/build_htslib_minimal.sh --prefix "$PWD/third_party/htslib-install"
cmake -S . -B build -DCMAKE_PREFIX_PATH="$PWD/third_party/htslib-install"
```

The script verifies that no `curl_*` symbols survive in the archive.

## Usage

### `profile`, the headless path

```
--bam FILE               indexed BAM/CRAM (a .bai/.csi must exist)
--regions FILE           BED or GTF of loci (.gz accepted)
--out-matrix FILE        write the regions x bins matrix as TSV
--out-profile FILE       write the meta-profile curve as TSV
--reference-point WHAT   TSS | center | TES        [TSS]
--upstream BP            window before the anchor   [2000]
--downstream BP          window after the anchor    [2000]
--bin-size BP            bin width                  [50]
--normalization WHAT     CPM | RPKM | BPM | raw     [CPM]
--count-mode WHAT        reads | depth              [reads]
--min-mapq N             discard alignments below   [0]
--max-fragment BP        drop longer fragments      [1000]
--extend-reads BP        single-end extension       [0 = off]
--no-extend-fragment     score aligned spans, not fragments
--proper-pairs-only      require the SAM proper-pair flag
--ignore-strand          do not orient by the BED strand column
-p, --threads N          worker threads             [all cores]
-q, --quiet              suppress the progress meter
```

With neither `--out-matrix` nor `--out-profile`, the meta-profile goes to stdout
as TSV so it drops into a shell pipeline.

### `serve`, the dashboard

```
--port N                 listen port                [8080]
--host ADDR              bind address               [127.0.0.1]
--data-root DIR          confine file access to DIR (REQUIRED)
--no-confinement         run with no data root; single-user hosts only
--no-browser             do not launch a browser
--no-auth                do not require the session token; local scripts only
--allow-remote           permit a non-loopback bind (off by default)
```

**`--data-root` is required.** Paths are canonicalised and rejected if they
escape it. It used to be optional, and omitting it disabled confinement
entirely: an unauthenticated read of every file the invoking user could open,
which was the shipped default. Running unconfined now takes `--no-confinement`
in as many words, which is greppable and shows up in `ps`.

**Every `/api/` request needs this run's session token.** The server mints one
at startup from the platform CSPRNG and prints it inside the URL it opens:

```
cuttag_profiler dashboard ready at http://127.0.0.1:8080/?t=cd60e58c...
```

Opening that link is the whole authentication step. The dashboard reads the
token out of the URL, forwards it on every call through a single helper, and
strips it from the address bar. For a script, pass it as `?t=` or as
`Authorization: Bearer <token>`, and prefer the header, since a query parameter
can end up in a log. `--no-auth` turns the requirement off for local scripted
use and cannot be combined with `--allow-remote`.

A loopback bind is not a security boundary on its own. It restricts which
network interfaces reach the socket, not which web pages can make your browser
issue requests to it. `Host` and `Origin` validation is what stops a
DNS-rebinding page; the token is what stops another process on the same machine.
The token is compared with `constant_time_equals` and nothing else, because `==`
short-circuits and leaks the length of a correct prefix. All nine findings from
[`SECURITY_HTTP_2026-08-15.md`](../../docs/SECURITY_HTTP_2026-08-15.md) are
closed, with one stated residual: the TOCTOU between path resolution and open is
fully closed for region files and only partly closed for BAMs, because
`sam_index_load` needs a path rather than a descriptor.

Building with `-DPTO_CLOUD_BUILD=ON` removes the server from the binary
altogether: no listener, no `serve` subcommand, and no `socket`, `bind` or
`listen` symbols in the image. That is the right configuration for anything
running the batch `profile` path on someone else's data, and the reasoning is
that a per-run token bound to one operator's browser is the wrong isolation
primitive for a multi-tenant service. In other words, the safest way to secure a
listening port there is to ensure it cannot be compiled into existence.

## How the signal is computed

**Windowing.** For a locus with anchor *a*, the window is
`[a - upstream, a + downstream)`, split into
`ceil((upstream + downstream) / bin_size)` bins. The anchor is the region start
for `TSS`, the midpoint for `center`, or the end for `TES`. For minus-strand
features the window is mirrored and the row reversed, so bin index always runs
5' to 3'.

Bin geometry is bounded so that every `width x bin_index` product fits in
`int64_t`, which gives a window ceiling of roughly 922 Gbp at the 10,000,000-bin
limit, enforced identically at the CLI and at the `/api/profile` handler. That
bound exists because an unbounded window was found to overflow signed 64-bit
arithmetic and write wrong bin offsets in a Release build, with no diagnostic
whatsoever.

**What counts.** By default a properly paired alignment contributes its full
sequenced fragment from `TLEN`, counted once via the leftmost mate, which is the
right footprint for CUT&Tag where the fragment marks the protected region.
`--no-extend-fragment` scores aligned spans instead, and `--extend-reads` covers
single-end libraries. Because a fragment can start before the window and still
intersect it, the iterator query is widened to the left by the maximum fragment
length; skipping that step silently truncates signal at window edges.

**Normalisation.** With *N* the mapped library size (read from index metadata
rather than by streaming the file), *C* the fragments overlapping a bin, and *w*
the bin width:

| Mode | Definition |
|:--|:--|
| `raw` | *C*, or mean per-base depth under `--count-mode depth` |
| `CPM` | *C* x 10^6 / *N* |
| `RPKM` | *C* x 10^9 / (*N* x *w*) |
| `BPM` | per-bp rate divided by the summed rate over the supplied loci, x 10^6 |

`BPM` is computed over the locus set you supply and not genome-wide, because no
genome-wide pass is performed. If you need a genome-wide BPM denominator, use
CPM or RPKM.

`--count-mode depth` reports mean per-base depth and honours the CIGAR string,
so deletions and reference skips do not inflate coverage. The default, `reads`,
matches what deeptools reports.

## Input expectations and failure modes

**Skipped regions are the behaviour most likely to change a figure.** A region
is skipped when its contig is not in the BAM, or when its window falls entirely
off the contig. Its row stays in the matrix, named by column 4 or by
`chrom:start-end`, but it is left out of everything computed across regions: the
meta-profile mean, BPM's total, and the scaling.

That fix landed on 2026-09-11 and it matters. Before it, a skipped region was
averaged in as a missing value of zero, so a region set half on an absent contig
produced a profile half as tall. The curve was smooth, plausible, and wrong by a
constant factor the user had no way to see.

| Input | What happens |
|:--|:--|
| A region on a contig absent from the BAM | Skipped, named, and excluded from every cross-region computation. |
| Every region skipped | Outputs are written and the run exits 2, naming a missing contig. That is almost always a `chr1` against `1` naming mismatch rather than an empty profile. |
| A window straddling a contig edge | Kept. Out-of-bounds bins are set to the missing value and counted as such in the mean, which is the parity the validation suite pins against deeptools' `--missingDataAsZero`. |
| An option value out of range, at the CLI or over HTTP | Refused. Guards live on `ProfileOptions` and `FilterOptions` and in `narrow_option()`, never in one entry point, because they were once validated after being narrowed to `int`. |
| A window beyond the bin-geometry ceiling | Refused, rather than overflowing into wrong bin offsets. |
| An output path naming an input | Refused. |
| A GTF attribute key that is a prefix of a longer one | Handled. A name was once taken from the longer key. |

This module's record is 9 confirmed defects, 5 of them High or Critical, across
the reviews and the 2026-09-11 audit. Details in
[`AUDIT_2026-09-11_cuttag_profiler.md`](../../docs/AUDIT_2026-09-11_cuttag_profiler.md).

## HTTP API

The dashboard is a client of these, and so can your scripts be.

| Endpoint | Method | Purpose |
|:--|:--|:--|
| `/api/health` | GET | version, core count |
| `/api/chroms?bam=...` | GET | header contigs, mapped-read total |
| `/api/profile` | POST | compute a matrix; returns profile, heatmap, stats |
| `/api/matrix?token=...` | GET | full matrix as TSV, keyed by a run token |

```bash
curl -s -X POST http://127.0.0.1:8080/api/profile \
  -H 'Content-Type: application/json' \
  -H 'Authorization: Bearer <session-token>' \
  -d '{"bam":"/data/sample.bam","regions":"/data/promoters.bed",
       "upstream":2000,"downstream":2000,"binSize":50,
       "normalization":"CPM","referencePoint":"TSS"}'
```

Full matrices are fetched separately as TSV rather than serialised into every
interactive response. 100k by 80 doubles is not worth pushing through JSON on
each click.

The `/api/profile` response does not echo `countMode`. Do not render fields the
server does not send.

## Threading notes

Work is distributed one locus per task, block partitioned so each worker writes
a contiguous stripe of the matrix and the write path needs no synchronization.
htslib's index and header are immutable once loaded and shared read-only. Its
`htsFile` handles are not, because decompression state is mutable per handle, so
each worker leases a private one.

One counter-intuitive result is baked into the defaults. htslib offers a thread
pool for BGZF inflation, and the obvious move is one pool sized to the machine.
For per-region parallelism that is actively harmful, because it funnels every
block through the pool's queue and adds a hand-off per block. Measured on the
benchmark above:

| Configuration | Wall | User | Sys | Peak RSS |
|:--|--:|--:|--:|--:|
| Shared pool, 12 threads | 5.95 s | 17.9 s | 34.6 s | 133 MB |
| Shared pool, 1 thread | 6.26 s | 7.6 s | 2.9 s | 89 MB |
| **No pool, inline inflation** | **1.38 s** | **5.3 s** | **0.9 s** | **85 MB** |

So the default is no shared pool, and each handle inflates on the thread that
owns it: `choose_bgzf_thread_count()` returns 0 unconditionally
(`bam_reader.hpp:65-67`), so no pool is built for any caller, single-threaded
included -- inline inflation with no queue hand-off wins there too. Rather
than merely adding threads wherever a library offers them, the profiler places
parallelism at the level where the work actually divides.

One more threading hazard is worth naming because it is invisible from the API.
`sam_hdr_name2tid()` takes a non-const `sam_hdr_t*` because on its first call it
builds the header's `hrecs` index in place, unsynchronised. Calling it from a
worker pool had several threads each build a full index and store it over one
another: a data race plus roughly 160 KB leaked across 22 allocations. The fix
is a name-to-tid table built once in the constructor, with the htslib call kept
only for `@SQ AN:` aliases and serialised behind a mutex. Anything taking a
non-const `sam_hdr_t*` is suspect.

## Dashboard design rules

These erode one change at a time, so they are written down.

Colour carries exactly one meaning, which is status. Panels, headings, rules and
the sidebar are neutrals; the accent appears in the signal curve, the focus ring
and the primary button and nowhere else. That is what makes an amber "regions
skipped" pill readable at a glance.

Numbers use `font-variant-numeric: tabular-nums` throughout so metric columns
align. The heatmap ramp is single-hue between two CSS custom properties, which
the canvas reads so it re-derives on a theme flip; it replaced a multi-hue ramp
that invented category edges in continuous data and whose darkest stop was
hard-coded to the dark theme's background, producing a near-black slab in light
mode. Everything is inlined, with no CDN, no webfont and no build step, because
the Content-Security-Policy and the air gap both forbid a second request.

## Layout

```
CMakeLists.txt            static htslib resolution, asset embedding, tests
cmake/embed_assets.cmake  turns web/ into a C++ byte-array translation unit
include/profiler/         public headers (no htslib leakage)
src/                      bam_reader, bed_reader, signal_calc, http_server, main
web/index.html            the dashboard: one file, no build step, no dependencies
tests/                    synthetic-BAM suite; no fixtures, no reference genome
scripts/                  demo data, benchmark harness, minimal htslib build
docs/APPLICATION_NOTE.tex manuscript; every figure is a macro fed from results/
```

Tests synthesise their own BAMs through htslib's write API at run time, so the
suite needs no checked-in fixture files and no reference genome.

Path handling opens once and asks `fstat` and `S_ISREG` of the descriptor rather
than of the name, and region files are parsed straight from that descriptor. Do
not reintroduce an exists-then-open shape; the `O_NONBLOCK` in that open is what
stops a FIFO hanging the handler thread, and it is not a tuning choice.

## License

MIT, under the root [`../../LICENSE`](../../LICENSE). Vendored `cpp-httplib` and
`nlohmann/json` keep their own MIT grants, reproduced in
[`../../THIRD_PARTY_LICENSES.md`](../../THIRD_PARTY_LICENSES.md).
