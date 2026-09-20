# fastq_stream — Architecture

This document describes what the system actually does, including the places
where the original design brief could not be implemented as written. Those are
called out explicitly in §3 and §9 rather than hidden, because they change what
the tool can honestly claim.

---

## 1. Problem statement

An NGS core facility's preprocessing step is usually a shell pipeline:

```
zcat in.fq.gz | trimmomatic ... | gzip > out.fq.gz && bwa-mem2 mem ref.fa out.fq.gz
```

Three costs dominate, and none of them is the actual QC arithmetic:

1. **Serial decompression.** `zlib` inflate runs at roughly 200–400 MB/s on one
   core and cannot be parallelised for a single-member gzip stream.
2. **Intermediate materialisation.** The trimmed FASTQ is compressed, written to
   disk, and read back by the aligner — paying deflate (far more expensive than
   inflate), then disk write, then disk read, then inflate again.
3. **Pipe and interpreter overhead.** Each stage copies through kernel pipe
   buffers, often wrapped in Python or Bash glue.

`fastq_stream` removes (2) entirely, reduces (1) where the input framing permits
it, and replaces (3) with in-process lock-free hand-offs.

---

## 2. Pipeline topology

```
                    ┌──────────── BGZF input only ────────────┐
                    │                                          │
 reader ──[SPSC]×D──┴─▶ inflaters ──[SPSC]×D──▶ assembler      │
   │                                                │          │
   └────────────── plain gzip / raw ────[SPSC]──────┘          │
                                                    │
                                              [SPSC]×W
                                                    ▼
                            writer ◀──[SPSC]×W──◀ workers
```

| Stage | Threads | Responsibility |
|---|---|---|
| reader | 1 | Frame the input. BGZF: slice whole members. Plain gzip: streaming zlib inflate. Raw: pass through. |
| inflaters | D | `libdeflate_gzip_decompress` on batches of 4 BGZF members. BGZF only. |
| assembler | 1 | Splice records spanning chunk boundaries; cut only on 4-line multiples; dispatch to workers. |
| workers | W | Parse records, accumulate QC, trim, filter, serialise survivors. |
| writer | 1 | Ordered write to file, FIFO, or stdout. |

Default split with `--threads N`: `spare = max(1, N-3)`, `D = spare/3` (BGZF
only), `W = spare - D`. The reader, assembler and writer always exist, so
`--threads 4` yields one worker. The text report always prints the real
breakdown — see §9.

### 2.1 Why every queue is genuinely SPSC

The brief specified lock-free SPSC ring buffers. An SPSC ring is only correct
with exactly one producer thread and one consumer thread; a shared queue read
by a pool of workers is MPMC and needs different (slower) algorithms.

Rather than weaken the claim, the topology is built so no ring is ever touched
by more than two threads: the reader owns one lane *per* inflater, each inflater
owns exactly one output lane, the assembler owns one lane per worker, and each
worker owns one output lane. `SpscRing` is therefore literally SPSC everywhere.

### 2.2 Ordering without a reorder buffer

Chunk *k* is dispatched to lane `k mod D`. Each lane preserves order internally.
A collector reading lanes round-robin in the same order therefore observes the
original stream order — no sequence numbers to sort, no reorder buffer, no
priority queue.

This imposes one requirement: **a worker must emit exactly one output buffer per
input chunk**, even an empty one. Otherwise the collector cannot know how many
buffers belong to chunk *k*. That is why output buffers are drawn from a pool
sized `kChunkBytes + kPrefix` — large enough that one input chunk's survivors
always fit in one output buffer (trimming only ever shrinks a record, and the
id/`+` lines are copied verbatim).

Termination uses a per-lane `closed` flag. A lane that is drained *and* closed
means no further chunks exist anywhere, because dispatch was strictly
round-robin from lane 0.

### 2.3 Ring implementation

`include/fastq_stream/lockfree_queue.hpp`:

- Capacity is a power of two; one slot is sacrificed so `head == tail`
  unambiguously means empty.
- `head_` and `tail_` sit on separate cache lines
  (`std::hardware_destructive_interference_size`, falling back to 128 bytes on
  arm64 and 64 elsewhere).
- Each side caches the opposite index in thread-private storage, so a
  steady-state push/pop never reads the other core's line. The cached value is
  refreshed only when the ring appears full (producer) or empty (consumer).
- `try_push`/`try_pop` are release/acquire respectively; there is no CAS and no
  fence stronger than needed.

Measured: ~14 ns median hand-off, ~70 M items/s between two threads.
Clean under ThreadSanitizer.

---

## 3. The libdeflate constraint (important)

**The brief specified "multi-threaded block reader leveraging `libdeflate` to
decompress raw gzipped blocks". This is not possible for arbitrary `.fastq.gz`
input, and the reason is structural, not an implementation shortcut.**

`libdeflate` deliberately exposes no streaming API. `libdeflate_gzip_decompress`
takes a complete member and a caller-provided output buffer whose size must be
known up front. That design is what makes it fast, and it means:

- A **single-member** gzip stream — which is what `bcl2fastq`, DRAGEN and plain
  `gzip` produce — has exactly one member spanning the entire file. Feeding it
  to libdeflate requires the whole 50 GB compressed file resident in RAM, and
  there is no way to find an interior DEFLATE block boundary without inflating
  everything before it. Neither parallel decompression nor bounded memory is
  achievable for this framing, by libdeflate or by anything else.
- **BGZF** (the SAM/BAM/`bgzip` framing) is a sequence of independent gzip
  members, each declaring its compressed size in the header and expanding to at
  most 64 KiB. Members can be located by header arithmetic alone and inflated
  concurrently.

So the tool implements both and reports which path it took:

| Input | Path | Parallel inflate | Memory |
|---|---|---|---|
| BGZF | libdeflate, D threads | yes | O(1) |
| plain gzip (single or multi-member) | streaming zlib inflate, 1 thread | no | O(1) |
| uncompressed | pass-through | n/a | O(1) |

Measured on 200k×150 bp reads (Apple M-series, 12 cores), byte-identical output
from both paths:

| Framing | Compressed throughput | Inflater threads |
|---|---:|---:|
| BGZF | 914 MB/s | 3 |
| plain gzip | 372 MB/s | 0 |

**2.5× purely from framing.** The practical consequence for a facility: running
`bgzip` instead of `gzip` at the point of FASTQ generation, once, unlocks
parallel decompression for every downstream consumer forever. Any published
comparison must state which framing was used — quoting the BGZF number against a
competitor's plain-gzip number would be a rigged benchmark, which is why
`scripts/run_benchmarks.py` measures both.

---

## 4. Memory model

- Payload bytes are **never copied between stages**; only `Buffer` handles move
  through the rings.
- Buffers are 64-byte aligned, 256 KiB payload, drawn from a recycling pool.
- Each allocation carries `kPrefix` = 64 KiB of headroom *in front of* `data()`,
  so the assembler can splice a straddling record onto the front of the next
  chunk in place instead of memmoving the whole payload. This also bounds the
  largest supported FASTQ record at 64 KiB (~32 kbp reads); longer inputs are
  rejected with an explicit error, never silently corrupted.
- The pool free list is guarded by a mutex. This is deliberate and not a
  contradiction of "lock-free": it is touched once per 256 KiB of stream
  (~200k times for a 50 GB dataset) against hundreds of millions of lock-free
  ring operations. The *data path* is lock-free; buffer recycling is not, and
  does not need to be.

Measured peak RSS on a 200k-read dataset: **~10 MB**, against the 512 MB budget.
The ceiling is structural — `(lanes × ring slots + in-flight) × 256 KiB` — and
does not grow with input size.

---

## 5. SIMD kernels

`include/fastq_stream/simd.hpp` dispatches at compile time: AVX-512BW → AVX2 →
NEON → scalar. There is no runtime dispatch; build once per target ISA
(`-DFQ_NATIVE=ON` uses `-march=native` / `-mcpu=native`).

| Kernel | Purpose | x86-64 | arm64 |
|---|---|---|---|
| `sum_u8` | mean Phred | `_mm256_sad_epu8` | `vpadalq_u16(vpaddlq_u8(...))` |
| `count_ge` | Q20/Q30 fractions | `max_epu8`+`cmpeq` (avoids signed-compare trap) | `vcgeq_u8` |
| `add_widen_u8_to_u32` | per-cycle quality | `_mm256_cvtepu8_epi32` | `vmovl_u8`/`vmovl_u16` |
| `mismatches` | adapter extension | `cmpeq`+`sad`, early exit | `vceqq_u8`+`vbicq_u8`, early exit |

**Every SIMD path is asserted bit-identical to a scalar reference at every
length from 0 to 300**, including the sub-vector tails where hand-written
kernels usually break, and `add_widen` is checked for writes past `n`.

Per-cycle quality sums accumulate into a 32-bit scratch folded into 64-bit
totals once per chunk. A 256 KiB chunk holds ~650 reads; overflow at one cycle
would need ~58 M reads, so the inner loop stays in 32-bit lanes with no risk.

Measured (NEON, best of 5, 150 bp reads): `sum_u8` 55 GB/s, `count_ge` 49 GB/s,
`add_widen` 43 GB/s, `trim_sliding_window` 89 GB/s.

---

## 6. Trimming semantics

- **`--trim-front/--trim-tail Q`** — hard 5'/3' trim below Phred Q.
- **`--window N --window-mean Q`** — Trimmomatic `SLIDINGWINDOW` semantics:
  scan 5'→3' with a rolling sum and truncate at the **start** of the first
  window whose mean falls below Q. This sacrifices up to `N-1` good bases, which
  is Trimmomatic's documented behaviour and is asserted in the tests. It is not
  a bug; tools differ here, and a comparison must use matched parameters.
- **Adapters** — seed-and-extend (§7).
- **Filters** — `--min-len`, `--max-n-rate`, `--min-mean-q`.

Order: front → tail → window → adapter → filters.

---

## 7. Adapter matching, and its measured specificity cost

A naive all-offsets scan costs O(read × adapter) and dominates everything else.
Instead each adapter contributes an exact seed (its first `--adapter-seed`
bases); only exact seed hits are extended, using the early-exit SIMD Hamming
kernel. A 4-byte seed is compared as a single unaligned 32-bit load — replacing
a variable-length `memcmp` here was worth **3.4×** (0.25 → 0.95 GB/s).

**The trade-off is real and measurable.** On 200k random 150 bp reads with
40,000 planted TruSeq adapters:

| `--adapter-overlap` | Reads trimmed | Spurious | Spurious rate |
|---:|---:|---:|---:|
| 4 (default) | 43,528 | 3,528 | 2.2% of reads |
| 6 | 40,968 | 968 | 0.6% |
| 8 | 40,248 | 248 | 0.15% |
| 10 | 40,181 | 181 | 0.11% |

At the default overlap of 4, a 4-mer matches random sequence roughly every 256
positions and the mismatch budget at overlap 4 is zero, so short spurious 3'
trims are common. The default is kept at 4 for consistency with the
cutadapt/fastp family, but **`--adapter-overlap 8` is recommended** for
production use: it removes ~93% of spurious calls while only missing adapter
remnants shorter than 8 bp, which have negligible effect on alignment.

A second, inherent limitation: a read whose adapter copy carries a sequencing
error *inside the seed window* is missed entirely. This is the same compromise
fastp makes. Lower `--adapter-seed` to trade throughput for sensitivity.

---

## 8. Correctness

`tests/test_trimmer.cpp` — 2,415 assertions, no external test framework (the
project's premise is that it drops onto a cluster with only a compiler):

- SIMD vs scalar at every length 0–300.
- SPSC ring: full/empty disambiguation, wrap-around, and a 2 M-item threaded
  hand-off asserting exact ordering and checksum.
- `'@'` as a quality character — proving records cannot be found by scanning
  for `@` — plus CRLF handling and truncated-record detection.
- **End-to-end equivalence: raw, plain gzip and real BGZF inputs must produce
  byte-identical output and identical statistics**, and the result must be
  invariant under `--threads 1, 2, 3, 16`. This is what actually guards the
  parallel paths.
- Malformed input (truncated record, seq/qual length mismatch) must raise, not
  silently truncate.

Clean under **ThreadSanitizer** and **ASan+UBSan**.

---

## 9. Known limitations

1. **Plain gzip cannot be inflated in parallel** (§3). The headline
   parallel-decompression claim applies to BGZF input only.
2. **No paired-end support.** Single-end only. Proper PE requires synchronised
   R1/R2 record streams and paired filtering decisions; the current dispatch
   would need to key on read pairs rather than byte chunks.
3. **`--threads N` is a total budget, not a worker count.** Three threads are
   always spent on reader/assembler/writer, so `N ≤ 4` gives one worker and
   `-t 1` still spawns four threads. The report prints the real breakdown; the
   benchmark labels rows by actual worker count for this reason.
4. **Max record size 64 KiB** (§4) — fine for Illumina, not for ONT/PacBio.
5. **Output is uncompressed.** By design: the target is a FIFO into an aligner.
   Writing `.gz` output would reintroduce the deflate cost the tool exists to
   avoid.
6. **Scaling saturates.** Measured 645 / 1863 / 3063 / 4304 MB/s at 1 / 3 / 5 / 9
   workers — 6.7× from 9× the workers. Beyond core count, oversubscription
   reverses the gains.
7. **No BAM input** despite the project title; FASTQ only. BAM would need a
   record-level parser on top of the existing BGZF layer.
8. **Adapter specificity at default settings** (§7).
