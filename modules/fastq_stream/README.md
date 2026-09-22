# fastq_stream

Low-latency, in-memory FASTQ streaming, QC and trimming.

It reads gzip or BGZF FASTQ, computes single-pass QC, trims adapters and
low-quality tails, and streams clean reads straight into an aligner through a
named pipe, without ever materializing an intermediate compressed FASTQ on disk.

```bash
fastq_stream -i sample.fq.gz --stdout | bwa-mem2 mem ref.fa - > sample.sam
```

Pure C++20, depending only on `libdeflate` and `zlib`. No network, no service,
no runtime beyond libc++ or libstdc++.

The point of eliminating that intermediate is not disk space. A conventional
decompress, trim and recompress pipeline writes a file and then reads it back,
which is two full passes over the data and a compression round-trip that exists
only because the stages are separate processes.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
```

Requires CMake 3.20 or newer and a C++20 compiler. `libdeflate` is found from
the system if installed, otherwise fetched and built statically. On macOS with
Homebrew, add `-DCMAKE_PREFIX_PATH=/opt/homebrew`.

| Option | Default | Meaning |
|:--|:--|:--|
| `FQ_NATIVE` | ON | `-march=native` or `-mcpu=native` |
| `FQ_BUILD_TESTS` | ON | build tests and microbenchmarks |
| `FQ_STATIC_LINK` | OFF | static libstdc++ and libgcc, Linux only |
| `FQ_SANITIZE` | OFF | ASan and UBSan |
| `FQ_ENABLE_AVX512` | OFF | compile the AVX-512BW kernels |

Turn `FQ_NATIVE` **off** for any redistributable binary. A `-march=native` build
faults on older hosts, and the fault is a SIGILL partway through a run rather
than a diagnostic at startup. The SIMD kernels are runtime-dispatched through
CPUID regardless, so turning it off costs nothing in kernel selection.

```bash
./build/test_trimmer     # 2,415 assertions
./build/benchmark        # component and end-to-end microbenchmarks
```

## Usage

```bash
fastq_stream -i in.fq.gz --qc-only -j qc.json
```

```bash
fastq_stream -i in.fq.gz -o clean.fq --window 4 --window-mean 20 -l 36
```

```bash
fastq_stream -i in.fq.gz -o /tmp/clean.fq --mkfifo &
bwa-mem2 mem ref.fa /tmp/clean.fq > out.sam
```

That last form is the intended production mode. `fastq_stream --help` lists
every option and is the authority if it and this file disagree.

The options worth knowing: `--window` and `--window-mean` implement
Trimmomatic's `SLIDINGWINDOW` semantics; `--trim-front` and `--trim-tail` are
fixed-length cuts; `-l`/`--min-len`, `--max-n-rate` and `--min-mean-q` are the
read-level filters; `-a`/`--adapter` and `--adapter-overlap` control adapter
removal; `-t`/`--threads` is a total budget rather than a worker count.

## Input framing matters more than anything else

This is the most consequential finding in the module, and it is not about the
QC arithmetic at all.

`libdeflate` exposes no streaming interface. It inflates a complete gzip member
into a caller-sized buffer, which is the source of its speed rather than an API
oversight. A single-member plain gzip stream, which is what `bcl2fastq`, DRAGEN
and `gzip` itself produce, has exactly one member spanning the whole file, so no
interior DEFLATE block boundary can be located without inflating everything
before it. Parallel decompression of such a stream is structurally impossible,
for any tool, no matter how many cores it is given.

BGZF, the SAM/BAM block-gzip framing, is a sequence of independent members that
each declare their own compressed size and expand to at most 64 KiB. Members can
be located by header arithmetic and inflated concurrently.

Both paths are implemented, and the tool reports which one it used.

| Framing | Compressed throughput | Inflater threads |
|:--|--:|--:|
| BGZF | 914 MB/s | 3 |
| Plain gzip | 372 MB/s | 0 |

Measured on 200k reads of 150 bp, Apple M-series, 12 cores, byte-identical
output from both paths.

That 2.5x costs nothing to obtain. It follows from running `bgzip` rather than
`gzip` once, at the moment the FASTQ is generated, and it then benefits every
downstream consumer of that file permanently. If you control how FASTQs are
written, this is the single highest-value change available to you. The
bottleneck was never the arithmetic. It was the container.

## Measured performance

Apple M-series, 12 cores, NEON, 150 bp reads. Rows are labelled by actual worker
count, since three threads always go to the reader, assembler and writer.

| Config | Throughput | Peak pool |
|:--|--:|--:|
| `--threads 4`, 1 worker | 645 MB/s | 14 MB |
| `--threads 6`, 3 workers | 1,863 MB/s | 25 MB |
| `--threads 8`, 5 workers | 3,063 MB/s | 34 MB |
| `--threads 12`, 9 workers | 4,304 MB/s | 83 MB |

**Peak RSS on a 200k-read run is about 10 MB, and it does not grow with input
size.** The ceiling is structural: resident memory is bounded by
(lanes x slots + in-flight) x 256 KiB, which is a property of the topology and
not of the file. A 2 GB FASTQ and a 2 TB FASTQ have the same footprint.

Component figures, best of 5: SPSC hand-off 14 ns at 70 M items per second;
`sum_u8` 55 GB/s; `count_ge` 49 GB/s; `add_widen` 43 GB/s; sliding-window trim
89 GB/s; adapter seed-and-extend 0.95 GB/s.

Adapter matching is roughly 50x slower than every other kernel and is the
pipeline's compute bottleneck. `--no-adapter` is substantially faster when you
do not need it.

These are single-host numbers on a laptop with hybrid cores. Consecutive runs of
an unchanged kernel varied by up to 2x, which is why every figure is a
best-of-N. Treat them as ranking evidence rather than absolutes.

### The comparison against fastp and Trimmomatic

The module-level head-to-head has not been run, because neither tool was
available on the development host. Those cells are marked TODO rather than
filled with estimates.

The repository-level harness does have a measured `fastq/qc` case against fastp:
7.48x faster at 4.7 MB peak RSS against fastp's 1,163 MB, with read and base
counts exact against the reference. That is in the [root README](../../README.md).

To produce the module-level numbers yourself:

```bash
./scripts/download_giab_data.sh --subset 20000000
python3 scripts/run_benchmarks.py \
    --input data/giab/HG002.novaseq.pcr-free.30x.R1.subset.fastq.gz \
    --bgzf-input data/giab/HG002.novaseq.pcr-free.30x.R1.subset.bgzf.gz \
    --threads 16 --replicates 3
```

The harness runs each measurement in an isolated child process, because
`ru_maxrss` over `RUSAGE_CHILDREN` is a monotonic high-water mark and otherwise
leaks between tools. It reports medians over replicates, and it measures
`fastq_stream` on **both** framings so that the comparison against fastp's
plain-gzip input is fair rather than flattering.

## Adapter trimming: pick your overlap

On 200k random reads with 40,000 planted TruSeq adapters:

| `--adapter-overlap` | Trimmed | Spurious | Rate |
|--:|--:|--:|--:|
| 4 (default) | 43,528 | 3,528 | 2.2% |
| 8 | 40,248 | 248 | 0.15% |

The default of 4 matches the cutadapt and fastp convention, and it spuriously
trims about 2% of reads on random sequence. **`--adapter-overlap 8` is what I
recommend for production.** It removes roughly 93% of false calls and only
misses adapter remnants under 8 bp.

## Input expectations and failure modes

Every entry here is a refusal rather than a default, and most of them exist
because the alternative was found producing a plausible number.

| Input | What happens |
|:--|:--|
| Qualities in Phred+64 with no declaration | Exits 2 after printing the report. Detection reads the quality histogram rather than the per-base path. |
| `--phred-offset 64` | Validated and rebased to Phred+33 in the worker, before anything reads a quality byte. |
| A record over 64 KiB | Refused. Every record is held to the limit, not just a carried partial one: checking only the partial made acceptance depend on byte position and framing. |
| A header or `+` line with non-printable bytes | Refused. These are validated exactly like SEQ and QUAL, because an earlier version stopped at SEQ/QUAL and let NUL and ESC in read names reach the output. |
| A truncated gzip or BGZF stream | Refused. A BGZF stream must end in its end-of-file marker. |
| `-i X -o X` | Refused. This once truncated the input. |
| A FIFO as input | Works. The BAM probe no longer consumes the first 32 KiB. |
| Thread-spawn failure | Reported, rather than `std::terminate`. |

**The Phred+64 case is the one to understand.** Read as Phred+33, a Phred+64
file moves mean Q from 19.85 to 50.85 and the Q30 rate from 0.265 to 1.000. The
run completes, exits 0 and prints a QC report showing a flawless library. That
is the failure mode this whole project is arranged against: the output is not an
error, it is a number, and the number is wrong.

One residual is open and stated rather than papered over. Plain gzip and raw
FASTQ carry no end-of-stream marker analogous to BGZF's, so a file truncated
exactly between records or between gzip members is, to any parser, a valid
shorter input rather than a detectable truncation. Only an external checksum or
manifest would catch it.

## Limitations

* **Single-end only.** There is no paired-end support.
* **Plain gzip is serial.** The parallel-inflate claim is BGZF-only, and that is
  a property of the format rather than of this implementation.
* **Maximum record 64 KiB**, roughly 32 kbp. Illumina yes, ONT and PacBio no.
  Raising it means raising per-lane buffer memory, which is a capacity decision
  and not a bug fix.
* **Output is uncompressed by design.** The target is a FIFO into an aligner.
* **No BAM input.** FASTQ only.
* `--threads N` is a total budget, so `N` at 4 or below yields one worker.
* **AVX-512 is compiled but has never executed on hardware that natively
  supports it**, and the AVX2 kernels have run on x86-64 only under Rosetta 2,
  which executes AVX2 without reporting it in CPUID. Selection of that path was
  forced rather than observed, and no throughput claim is made for either ISA.

## Correctness

2,415 assertions, no external test framework. SIMD paths are asserted
bit-identical to a scalar double-precision reference at every read length from 0
to 300, a range chosen to cover the sub-vector tails where such kernels
typically diverge. Raw, gzip and BGZF inputs must produce byte-identical output
invariant under thread count. Malformed input must raise rather than silently
truncate.

The module is clean under ThreadSanitizer and under ASan plus UBSan with
`-fno-sanitize-recover=all`:

```bash
cmake -S . -B build-tsan -DCMAKE_CXX_FLAGS=-fsanitize=thread \
      -DCMAKE_EXE_LINKER_FLAGS=-fsanitize=thread
cmake --build build-tsan -j && ./build-tsan/test_trimmer
```

This module carries the largest share of the repository's defect record: **27
confirmed defects, 9 of them High or Critical**, across the 2026-09-10 torture
pass and the 2026-09-11 audit. Every one has a regression test added in the
commit that fixed it. The 2026-09-11 pass additionally drove a 900-mutant corpus
(bit flips, byte replacement, deletion, insertion, duplication and truncation
across raw, gzip and BGZF framing) through the sanitizer build at 1, 3 and 8
threads with no crash, hang or sanitizer report. Details in
[`AUDIT_2026-09-11_fastq_stream.md`](../../docs/AUDIT_2026-09-11_fastq_stream.md)
and [`TORTURE_2026-09-10.md`](../../docs/TORTURE_2026-09-10.md).

One rule from that record is worth repeating here because it looks like a
tuning knob and is not: the pipeline's `Backoff` must reach an actual sleep.
Yielding forever cost about 11 cores of system time, and shorter sleeps were
measured and are worse.

## Layout

```
include/fastq_stream/   simd, lockfree_queue, buffer, record,
                        phred_calculator, adapter_trimmer, reader,
                        pipeline, report
src/                    reader.cpp, pipeline.cpp, main.cpp
tests/                  test_trimmer.cpp, benchmark.cpp
scripts/                download_giab_data.sh, run_benchmarks.py
docs/                   ARCHITECTURE.md, PAPER_DRAFT.tex
```

The lock-free ring, its cache-line separation and the order-preserving lane
dispatch are explained in
[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md).

MIT licensed.
