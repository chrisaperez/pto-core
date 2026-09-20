# pto-core — design for four feature extensions

Date: 2026-08-16
Status: **implemented**, except where §7 says otherwise. §§1–6 are the design as
written before the work; §7 records what was actually built, where the design
turned out to be wrong, and what is deliberately still open. Where the two
disagree, §7 and the code are authoritative.

Ground rules this design is written against (from `CLAUDE.md` / `README.md`):

* Modules share no code and no runtime. Nothing here introduces a cross-module
  header, a shared `include/`, or a build dependency on the top-level
  `CMakeLists.txt`. Each module stays standalone-buildable.
* Strict local execution. Two features below have a live path that would reach
  the network by default (CRAM reference resolution, §1.4; PyPI-time
  toolchain fetches, §4). Both are closed explicitly, not by convention.
* Tests stay framework-free (`assert` + a `main`), matching
  `tests/test_matrix_ops.cpp` and `tests/test_http_security.cpp`.
* New allocation-driving inputs get a validated ceiling, in the same place for
  every entry path. That is what `ProfileOptions::validate()` exists for, and
  §3.5 and §2.4 each add one for the same reason.

---

## 1. BAM/CRAM streaming ingestion for `fastq_stream`

### 1.1 Where the seam goes

`fastq_stream`'s pipeline (`src/pipeline.cpp`) is five ranks: reader →
inflaters → assembler → workers → writer. Everything from the assembler
rightwards is defined over *bytes of FASTQ text* — `last_record_boundary()`
counts newlines, `RecordIterator` splits on `\n`, `RecordProcessor` writes
`@id/seq/+/qual` back out. Threading a `bam1_t`-shaped record type through
those ranks would mean a second `RecordProcessor`, a second `RecordIterator`,
and a variant `Buffer` payload — a large blast radius for no gain.

**Decision: BAM/CRAM decode replaces the reader rank only.** `HtsInput`
produces `Buffer`s already containing whole FASTQ records; the assembler,
workers and writer are untouched, and the SPSC lane topology is unchanged.

```
                 (existing)  reader ──▶ inflaters ──▶ assembler ──▶ workers ──▶ writer
   (BAM/CRAM)  hts_reader ──────────────────────────▶ assembler ──▶ workers ──▶ writer
```

Because `HtsInput` only ever ends a chunk on a record boundary, the
assembler's `carry` is permanently empty in this mode — it costs one `memchr`
pass at memory bandwidth and buys us a single code path. Do not add a
"pre-framed, skip the assembler" flag for this; the measurement to justify it
does not exist yet.

On the "zero-copy" requirement: honestly stated, BAM ingestion cannot be
zero-copy at the record level. `bam1_t` stores SEQ 4-bit packed and QUAL as raw
phred values, neither of which is the byte layout the SIMD trim kernels read.
The nibble expansion and the `+33` are *required transforms*, not copies we
chose. What the design preserves is the property that actually matters: the
expansion writes **directly into a pooled `Buffer`**, and from that point to
the writer the payload moves by pointer only — same as the gzip path.

### 1.2 New files

```
modules/fastq_stream/include/fastq_stream/hts_input.hpp   (new)
modules/fastq_stream/src/hts_input.cpp                    (new)
modules/fastq_stream/tests/test_hts_input.cpp             (new)
modules/fastq_stream/CMakeLists.txt                       (+ FQ_ENABLE_HTS)
modules/fastq_stream/include/fastq_stream/pipeline.hpp    (+ Config fields)
modules/fastq_stream/src/pipeline.cpp                     (reader rank switch)
modules/fastq_stream/src/main.cpp                         (+ CLI flags)
```

htslib is a heavy dependency for a module whose entire current dependency set
is zlib + libdeflate + pthreads, and clusters that run `fastq_stream` on raw
Illumina output do not need it. **`FQ_ENABLE_HTS` defaults to `OFF`**; with it
off, `hts_input.cpp` is not compiled and `--input *.bam` fails with a build-time
diagnostic rather than a missing-symbol link error. (`cuttag_profiler` requires
htslib unconditionally; that is correct there and wrong here.)

### 1.3 Header

```cpp
// SPDX-License-Identifier: MIT
// fastq_stream — BAM/CRAM ingestion.
//
// Decodes alignment records straight into pooled FASTQ-text Buffers at the
// reader rank, so nothing downstream of the assembler knows the input was not
// a FASTQ file. See docs/DESIGN_EXTENSIONS_2026-08-16.md §1.
#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include "fastq_stream/buffer.hpp"

namespace fq {

// Longest read this path will emit. Bounds the flush watermark below and, with
// kPrefix, the largest record the assembler could ever be asked to carry.
inline constexpr std::size_t kHtsMaxReadLen = 32 * 1024;

enum class HtsLayout {
  kSingle,        // one record per read; mates, if any, emitted independently
  kInterleaved,   // /1 then /2 back to back; requires name-grouped input
};

struct HtsOptions {
  HtsLayout layout = HtsLayout::kSingle;
  // Local FASTA for CRAM reference resolution. Empty + CRAM input is a hard
  // error, never a silent fetch from the EBI reference registry (§1.4).
  std::string reference_fasta;
  // Threads handed to htslib's own BGZF decode pool. The pipeline reassigns
  // its inflater rank's budget here (§1.6).
  int decode_threads = 0;
  // Substituted when QUAL is absent (0xff). Default 'I' == Q40.
  char missing_qual = 'I';
  // Refuse to buffer more than this many unpaired records while waiting for a
  // mate in kInterleaved mode. See §1.5.
  std::size_t max_pending_mates = 1u << 20;
};

struct HtsStats {
  uint64_t records_read = 0;
  uint64_t emitted = 0;
  uint64_t skipped_secondary = 0;    // FLAG 0x100
  uint64_t skipped_supplementary = 0;// FLAG 0x800
  uint64_t skipped_no_seq = 0;       // SEQ == '*'
  uint64_t synthesized_qual = 0;     // QUAL == 0xff
  uint64_t reverse_complemented = 0; // FLAG 0x10
  uint64_t singletons = 0;           // kInterleaved: mate never arrived
};

// Reader-rank source. Mirrors GzipInput's shape so pipeline.cpp can hold
// either behind one `next(Buffer&)` call.
class HtsInput {
 public:
  HtsInput(const std::string& path, BufferPool& pool, const HtsOptions& opts);
  ~HtsInput();
  HtsInput(const HtsInput&) = delete;
  HtsInput& operator=(const HtsInput&) = delete;

  // Fills `out` with whole FASTQ records only. False at end of stream.
  bool next(Buffer& out);

  const HtsStats& stats() const noexcept;
  uint64_t compressed_bytes() const noexcept;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

// True for .bam/.cram/.sam by htslib format probe, not by file extension.
bool looks_like_alignment(const std::string& path);

}  // namespace fq
```

### 1.4 Decode core — the three correctness traps

```cpp
// src/hts_input.cpp (core of the fill loop)

// SEQ is stored in alignment orientation. A read that aligned to the minus
// strand must be reverse-complemented and its QUAL reversed to recover the
// sequence the instrument actually produced. Getting this wrong is silent:
// the FASTQ is well-formed, the QC report looks plausible, and every
// downstream realignment is wrong for ~half the reads.
static constexpr char kNt16[] = "=ACMGRSVTWYHKDBN";
static constexpr char kNt16Comp[] = "=TGKCYSBAWRDMHVN";

void decode_seq(const bam1_t* b, char* out) {
  const uint8_t* packed = bam_get_seq(b);
  const int32_t n = b->core.l_qseq;
  if ((b->core.flag & BAM_FREVERSE) != 0) {
    for (int32_t i = 0; i < n; ++i) {
      out[n - 1 - i] = kNt16Comp[bam_seqi(packed, i)];
    }
  } else {
    for (int32_t i = 0; i < n; ++i) out[i] = kNt16[bam_seqi(packed, i)];
  }
}

void decode_qual(const bam1_t* b, char* out, char missing, HtsStats& st) {
  const uint8_t* q = bam_get_qual(b);
  const int32_t n = b->core.l_qseq;
  if (q[0] == 0xff) {                      // QUAL absent
    std::memset(out, missing, static_cast<size_t>(n));
    ++st.synthesized_qual;
    return;
  }
  if ((b->core.flag & BAM_FREVERSE) != 0) {
    for (int32_t i = 0; i < n; ++i) out[n - 1 - i] = static_cast<char>(q[i] + 33);
  } else {
    for (int32_t i = 0; i < n; ++i) out[i] = static_cast<char>(q[i] + 33);
  }
}
```

**Trap 1 — secondary and supplementary alignments.** One read can appear in
many records. Emitting all of them duplicates reads (secondary, `0x100`) and
emits hard-clipped fragments as if they were whole reads (supplementary,
`0x800`). Both are dropped unconditionally and counted; this is not a
user-tunable filter, because keeping them cannot produce a correct FASTQ.

**Trap 2 — CRAM reaches the network.** htslib resolves CRAM reference sequences
via `REF_PATH`/`REF_CACHE`, whose compiled-in default includes
`https://www.ebi.ac.uk/...`. A CRAM decode on a stock htslib will therefore
attempt an outbound request, sending the file's reference MD5s off-host. That
is a direct violation of this repo's central constraint, and it happens with no
diagnostic. Closed on open:

```cpp
if (fmt == cram) {
  if (opts.reference_fasta.empty()) {
    throw std::runtime_error(
        "CRAM input requires --reference <genome.fa>: without it htslib would "
        "resolve reference sequences over the network, which this tool does "
        "not do. Point it at the FASTA the CRAM was compressed against.");
  }
  // Belt and braces: pin the resolver to the local file *and* blank the
  // remote registry for this process, so a stale REF_PATH in the environment
  // cannot reintroduce the fetch.
  ::setenv("REF_PATH", opts.reference_fasta.c_str(), /*overwrite=*/1);
  if (hts_set_fai_filename(fp, (opts.reference_fasta + ".fai").c_str()) < 0) {
    throw std::runtime_error("cannot load FASTA index " +
                             opts.reference_fasta + ".fai (run: samtools faidx)");
  }
}
```

**Trap 3 — buffer watermark.** A record is written as a unit; a chunk must never
be cut mid-record. Flush when the remaining capacity drops below the worst-case
record size:

```cpp
// id + 2 newlines-worth of framing, seq, '+', qual, 4 newlines.
constexpr std::size_t kWorstRecord = 1024 + 2 * kHtsMaxReadLen + 8;
static_assert(kWorstRecord <= kChunkBytes, "chunk cannot hold one record");

if (buf.capacity() - used < kWorstRecord) { buf.set_size(used); return true; }
```

Reads longer than `kHtsMaxReadLen` are rejected by name with an explicit error,
matching how the assembler already handles over-long records rather than
silently truncating.

### 1.5 Pairing

Coordinate-sorted BAM places mates arbitrarily far apart — on human WGS, up to
the whole file. A name-keyed pending map is therefore unbounded, which is
exactly the memory-exhaustion shape this repo has closed twice elsewhere.

Design:

* Read `@HD SO:` from the header. `SO:coordinate` + `kInterleaved` is a **hard
  error at open time**, with the remedy in the message
  (`samtools collate -u -O in.bam | fastq_stream --input -`). Failing at open
  costs nothing; failing at 40 GB of resident pending records costs a node.
* For `SO:queryname` / `SO:unsorted` / unspecified, pair with a single-slot
  hold: keep the previous record, and if the next one shares its QNAME and is
  the opposite `0x40`/`0x80`, emit both. Name-grouped input needs nothing more
  than this, and it is O(1) memory.
* If input is name-*grouped* but not strictly adjacent, the bounded map is the
  fallback, capped by `max_pending_mates`; exceeding the cap is an error naming
  the flag, not a silent drop.
* Unmatched reads at EOF are emitted as singletons and counted (`st.singletons`),
  reported in the run summary. Silently discarding them would change the read
  count without saying so.

`kSingle` layout has none of these problems and is the default.

### 1.6 Pipeline and thread-budget changes

```cpp
// pipeline.hpp — Config additions
bool hts_input = false;          // set by main.cpp from looks_like_alignment()
HtsOptions hts;                  // ignored unless hts_input
```

```cpp
// pipeline.cpp — rank sizing
const int spare = std::max(1, hw - 3);
int n_inflate = 0, n_worker = 0, n_hts_decode = 0;
if (cfg.hts_input) {
  // Our libdeflate BGZF fast path cannot be used: BAM framing is record-level,
  // not line-level, so members must be inflated by htslib itself in order.
  // The inflater budget moves into htslib's own decode pool.
  n_hts_decode = std::max(1, spare / 3);
  n_worker = std::max(1, spare - n_hts_decode);
} else {
  n_inflate = (fmt == GzipFormat::kBgzf) ? std::max(1, spare / 3) : 0;
  n_worker = std::max(1, spare - n_inflate);
}
```

The reader thread body becomes a two-line switch over a
`std::variant<GzipInput, HtsInput>` (or a tiny `struct Source` with a
`next(Buffer&)` virtual — the call is once per 256 KiB, so the indirect call is
free). `parallel_inflate` is false in HTS mode, so the existing `direct_lane`
path is reused verbatim.

`RunReport` gains `HtsStats hts_stats;` and the summary prints the skip/RC/
singleton counts. These numbers are the only evidence a user has that the
conversion was faithful; they are not optional output.

### 1.7 Tests (`tests/test_hts_input.cpp`, framework-free)

1. Synthesise a 4-record BAM in-process with `sam_hdr_parse` + `bam_set1`
   (no external fixture file, matching how the gzip tests build their inputs).
2. **Reverse-strand round-trip**: a record with `FLAG 0x10` must decode to the
   reverse complement of the stored SEQ, with QUAL reversed. Assert byte-exact
   against a hand-written expectation, not against a second implementation.
3. Secondary + supplementary records are dropped; counters match.
4. `QUAL == 0xff` yields `missing_qual` for every base and increments the
   counter.
5. `SEQ == '*'` is skipped, not emitted as a zero-length record.
6. Chunk-boundary: enough records to force ≥3 flushes; concatenated output
   parses cleanly with `RecordIterator` and every record is `well_formed()`.
7. Coordinate-sorted header + `kInterleaved` throws at construction.
8. End-to-end: `run_pipeline` over the synthetic BAM produces the same QC stats
   as `run_pipeline` over the equivalent FASTQ text. This is the test that
   proves the seam is behaviour-preserving.

### 1.8 CMake

```cmake
option(FQ_ENABLE_HTS "Enable BAM/CRAM ingestion via htslib" OFF)

if(FQ_ENABLE_HTS)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(HTSLIB REQUIRED htslib)
  target_sources(fastq_stream_core PRIVATE src/hts_input.cpp)
  target_compile_definitions(fastq_stream_core PUBLIC FQ_ENABLE_HTS=1)
  target_include_directories(fastq_stream_core SYSTEM PRIVATE ${HTSLIB_INCLUDE_DIRS})
  target_link_libraries(fastq_stream_core PRIVATE ${HTSLIB_LINK_LIBRARIES})
  target_link_directories(fastq_stream_core PRIVATE ${HTSLIB_LIBRARY_DIRS})
endif()
```

Note `SYSTEM PRIVATE`: htslib's headers do not survive `-Wall -Wextra
-Wpedantic`, and CI's warning grep is scoped to `modules/*/(src|include|tests)/`
— keeping htslib out of the public include set keeps that grep honest. Add a CI
job that configures `fastq_stream` standalone with `FQ_ENABLE_HTS=OFF` to
assert the dependency stays optional.

---

## 2. HNSW index serialization for `scrna_matrix`

### 2.1 What has to change first

`hnsw_build_and_search{,_dense}()` construct `hnswlib::HierarchicalNSW` as a
**local** and destroy it on return. There is no index object to save. The
feature is therefore primarily a refactor: lift the index into a first-class
type, then express both existing functions in terms of it so there is exactly
one construction path.

```
modules/scrna_matrix/include/matrix/hnsw_index.hpp   (new, header-only)
modules/scrna_matrix/include/matrix/knn_graph.hpp    (build/search delegate to it)
modules/scrna_matrix/src/python_bindings.cpp         (+ HnswIndex class)
modules/scrna_matrix/tests/test_hnsw_persist.cpp     (new)
modules/scrna_matrix/tests/test_python_bindings.py   (+ save/load round-trip)
```

### 2.2 The file format, and why hnswlib's own is not enough

`HierarchicalNSW::saveIndex()` writes the graph, the vectors and the level
structure. It does **not** write the space. `loadIndex()` requires you to hand
it a `SpaceInterface*` and a `max_elements` — and it trusts them. Load a
50-dim index with a space constructed for 2000 dims and hnswlib computes
offsets from the wrong `data_size_`; the result is out-of-bounds reads inside
the loaded blob, not an error.

So the format is **our header plus hnswlib's payload**, written as two files
that must agree:

```
atlas.hnsw        hnswlib payload, written by saveIndex()
atlas.hnsw.meta   our fixed-layout binary header, below
```

Two files rather than one concatenated file because the vendored hnswlib's
`loadIndex` takes a path, not a stream, in the version pinned here — confirm
against `third_party/hnswlib/hnswalg.h` before landing; if the stream overload
is present, prefer a single self-describing file and drop the sidecar.

```cpp
// hnsw_index.hpp
struct HnswMeta {
  char     magic[8];        // "SCRNAHNW"
  uint32_t format_version;  // 1
  uint32_t dim;
  uint64_t element_count;
  uint32_t M;
  uint32_t ef_construction;
  uint32_t metric;          // 0 = ScrnaCosineSpace, 1 = InnerProductSpace
  uint32_t normalized;      // 1 = vectors were unit-normalised on insert
  uint64_t payload_bytes;   // expected size of the .hnsw file
  uint64_t payload_crc32;   // over the .hnsw file
  uint32_t scalar_layout;   // 0 = IEEE754 binary32 little-endian
  uint32_t reserved[5];
};
static_assert(sizeof(HnswMeta) == 88, "meta layout is part of the file format");
```

`metric` and `normalized` are load-bearing, not decoration: an index built
under `ScrnaCosineSpace` (which normalises on insert and returns `1 - cos`) and
queried under `InnerProductSpace` gives numerically plausible but wrong
neighbours. Mismatch is a refusal.

`scalar_layout` exists so a big-endian or non-IEEE754 host refuses the file
instead of reading garbage floats. x86-64 ↔ aarch64 interchange is fine and is
the case that actually matters (build the atlas on a cluster node, query on a
laptop).

### 2.3 API

```cpp
class HnswIndex {
 public:
  HnswIndex(std::size_t dim, std::size_t capacity, const HnswParams& p = {});

  void add_dense(const float* data, std::size_t n, int num_threads = 0);
  void add_csr(const Block_CSR<float>& mat, int num_threads = 0);

  // Query rows that need not be members of the index — the atlas-mapping case.
  // `exclude_self_label`: pass npos when querying foreign cells.
  [[nodiscard]] KnnGraph query_dense(const float* q, std::size_t n,
                                     std::size_t k, int num_threads = 0) const;
  [[nodiscard]] KnnGraph query_csr(const Block_CSR<float>& q, std::size_t k,
                                   int num_threads = 0) const;

  void save_index(const std::string& path) const;
  [[nodiscard]] static HnswIndex load_index(const std::string& path,
                                            const HnswParams& p = {});

  [[nodiscard]] std::size_t dim() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;
  void set_ef(std::size_t ef);

  // Optional cell barcodes, label-indexed. Serialized alongside the meta so a
  // loaded atlas can answer with names rather than integer labels.
  void set_labels(std::vector<std::string> names);
  [[nodiscard]] const std::vector<std::string>& labels() const noexcept;
 private:
  std::unique_ptr<hnswlib::SpaceInterface<float>> space_;   // must outlive index_
  std::unique_ptr<hnswlib::HierarchicalNSW<float>> index_;  // declared second
  HnswMeta meta_{};
  std::vector<std::string> labels_;
};
```

Declaration order is deliberate and worth a comment in the header:
`HierarchicalNSW` keeps a raw `SpaceInterface*` and calls through it in its
destructor. `space_` must be declared **before** `index_` so it is destroyed
after. This is the same class of lifetime bug the existing code avoids by
keeping both locals in one scope; moving them into a class makes ordering
explicit and easy to break in a later edit.

The two existing free functions keep their signatures and become thin wrappers
(`HnswIndex idx(...); idx.add_*(); return idx.query_*();`), so the graph-building
behaviour, the `max_index_bytes` guard and every existing test are unchanged.

### 2.4 Loading is a trust boundary

`load_index` takes a path and, from bytes in that file, decides how much memory
to allocate and how to index into a buffer. hnswlib's loader does not validate
its header. Given this repo's review trail, that gets treated as untrusted
input:

```cpp
HnswIndex HnswIndex::load_index(const std::string& path, const HnswParams& p) {
  HnswMeta m = read_meta(path + ".meta");                    // exact size read
  if (std::memcmp(m.magic, "SCRNAHNW", 8) != 0) throw ...;
  if (m.format_version != 1) throw ...;
  if (m.scalar_layout != 0) throw ...;                        // endianness
  if (m.dim == 0 || m.dim > kMaxDim) throw ...;
  // Same ceiling that gates a *built* index, applied to a *loaded* one.
  const std::size_t need = hnsw_index_bytes(m.element_count, m.dim, m.M);
  if (need > p.max_index_bytes) throw std::invalid_argument(...);
  if (file_size(path) != m.payload_bytes) throw ...;          // truncation
  if (crc32(path) != m.payload_crc32) throw ...;              // corruption
  ...
}
```

**Residual, to be recorded in `modules/scrna_matrix/docs/AUDIT.md`:** these
checks catch truncation, corruption, version skew and metric mismatch. They do
**not** make `load_index` safe against a deliberately malicious `.hnsw`
payload — the graph body is consumed by vendored hnswlib code that we do not
audit, and a hostile link array can send it out of bounds regardless of our
header. The honest statement for the docs is: *loading an index file is
equivalent to trusting whoever produced it.* Do not let this API become
reachable from the `cuttag_profiler` HTTP surface or any other unauthenticated
entry point.

### 2.5 Python bindings

```cpp
py::class_<HnswIndex>(m, "HnswIndex")
    .def(py::init<std::size_t, std::size_t, const HnswParams&>(),
         py::arg("dim"), py::arg("capacity"), py::arg("params") = HnswParams{})
    .def("add_dense",
         [](HnswIndex& self, py::array_t<float, py::array::c_style | py::array::forcecast> a,
            int threads) {
           auto info = a.request();
           if (info.ndim != 2) throw std::invalid_argument("expected a 2-D array");
           // GIL released: construction is the multi-threaded, minutes-long part.
           py::gil_scoped_release rel;
           self.add_dense(static_cast<const float*>(info.ptr),
                          static_cast<std::size_t>(info.shape[0]), threads);
         },
         py::arg("data"), py::arg("threads") = 0)
    .def("query_dense", ..., py::call_guard<py::gil_scoped_release>())
    .def("save_index", &HnswIndex::save_index, py::arg("path"))
    .def_static("load_index", &HnswIndex::load_index,
                py::arg("path"), py::arg("params") = HnswParams{})
    .def_property_readonly("dim", &HnswIndex::dim)
    .def("__len__", &HnswIndex::size);
```

The whole class sits inside the existing `#if defined(SCRNA_ENABLE_HNSW)` guard,
and `has_hnsw()` stays the runtime probe.

Target workflow this unlocks:

```python
idx = scrna_matrix.HnswIndex(dim=50, capacity=1_300_000)
idx.add_dense(pbmc_pca)            # once, offline: minutes
idx.save_index("pbmc1M.hnsw")

idx = scrna_matrix.HnswIndex.load_index("pbmc1M.hnsw")   # seconds, mmap-free
g = idx.query_dense(new_cells_pca, k=15)                 # milliseconds
```

### 2.6 Tests

* Round-trip: build over a fixed seed, save, load, query — neighbour indices and
  similarities **byte-identical** to querying the in-memory index. Not
  "approximately equal": serialization must be lossless, so any drift is a bug.
* Truncated payload → throws, does not crash. Flip a byte → CRC throws.
* `dim` mismatch between meta and requested space → throws before any allocation.
* Metric mismatch (cosine index, inner-product load) → throws.
* `element_count` inflated to 1e12 in the meta → refused by the `max_index_bytes`
  guard, not attempted.
* Cross-arch: CI stores an index built on x86-64 as an artifact and loads it in
  the aarch64 job (and the reverse). This is the only way the portability claim
  gets tested.

---

## 3. Genomic interval indexing for `cuttag_profiler`

### 3.1 Two problems, two structures

"Find all annotations overlapping this query interval" and "intersect a
coordinate-sorted BAM against an annotation set" are different problems and
deserve different answers:

| Access pattern | Structure | Cost |
|---|---|---|
| Random point/range queries (HTTP dashboard, enhancer–promoter links, gene bodies) | implicit augmented interval array (§3.2) | O(log n + m) per query, one contiguous allocation |
| Streaming intersection against a coordinate-sorted BAM | sweep line with an active set (§3.4) | O(n + m), no index at all |

Building an interval tree to walk a sorted BAM would be strictly worse than the
sweep — the BAM is *already sorted*, which is the sweep's entire precondition.
Both are in scope; do not implement one and use it for both.

### 3.2 `IntervalIndex` — layout

Pointer-based red-black interval trees are the textbook answer and the wrong
one here: a node per interval is a cache miss per level. Use the implicit
augmented-array form (the cgranges / Li layout) — sort by start, augment each
position with the maximum end in its implicit subtree, navigate by index
arithmetic. One allocation per contig, no pointers, and the hot arrays are
struct-of-arrays so a descent touches starts and max_ends without dragging
names and payloads through cache.

```cpp
// include/profiler/interval_index.hpp — header-only, no htslib include.
#pragma once
#include <cstdint>
#include <span>
#include <string>
#include <vector>
#include "profiler/types.hpp"

namespace profiler {

// Per-contig coordinates are stored as int32. Human chr1 is 2.5e8, three
// orders under the limit, and halving the hot arrays is the point of this
// structure. Contigs at or beyond INT32_MAX (some amphibian and conifer
// assemblies) are rejected by name at build time rather than silently
// truncated -- build() throws with the contig and its length.
using Coord = std::int32_t;

class ContigIntervals {
 public:
  // Sorted-by-start SoA. `max_end_[i]` is the maximum end over the implicit
  // subtree rooted at i, which is what turns the prefix scan into a descent.
  alignas(kCacheLineSize) std::vector<Coord>    starts_;
  alignas(kCacheLineSize) std::vector<Coord>    ends_;
  alignas(kCacheLineSize) std::vector<Coord>    max_end_;
  alignas(kCacheLineSize) std::vector<std::uint32_t> payload_;  // index into the owner's Region vector

  // Appends the payload ids of every interval overlapping [qs, qe) to `out`.
  // `out` is caller-owned scratch, reused across queries -- the query itself
  // allocates nothing. Returns the number appended.
  std::size_t overlap(Coord qs, Coord qe, std::vector<std::uint32_t>& out) const;
};

class IntervalIndex {
 public:
  // `regions` is borrowed and must outlive the index; payload ids are indices
  // into it. Immutable after build, therefore safe to share const across the
  // profiler's worker threads -- same contract as BamReader's index/header.
  static IntervalIndex build(const std::vector<Region>& regions);

  // Contig ids are resolved once, against the BAM header's ordering, so the
  // hot path never hashes a chromosome name.
  void bind_contigs(const std::vector<std::string>& header_order);

  [[nodiscard]] std::span<const std::uint32_t>
  overlap(int tid, Coord qs, Coord qe, std::vector<std::uint32_t>& scratch) const;

  [[nodiscard]] std::size_t size() const noexcept;
 private:
  std::vector<ContigIntervals> by_contig_;   // indexed by bound tid
  std::vector<int> tid_of_local_;
  const std::vector<Region>* regions_ = nullptr;
};

}  // namespace profiler
```

The zero-copy requirement is met by the `scratch`/`std::span` pair: a query
returns a view over caller-owned storage holding **payload ids**, never copies
of `Region` (which own `std::string`s and would allocate per hit). Callers
resolve `regions[id]` only for the hits they actually use.

### 3.3 Query

```cpp
// Half-open overlap: [as, ae) meets [bs, be) iff as < be && bs < ae.
// Stated once, here, because the off-by-one in this predicate is the single
// most common defect in interval code and BED is already half-open.
std::size_t ContigIntervals::overlap(Coord qs, Coord qe,
                                     std::vector<std::uint32_t>& out) const {
  if (qs >= qe || starts_.empty()) return 0;
  const std::size_t before = out.size();
  // Descend the implicit tree; prune any subtree whose max_end_ <= qs, since
  // no interval within it can reach the query.
  // ... (index-arithmetic descent, no recursion, no allocation)
  return out.size() - before;
}
```

Two properties to assert in tests rather than assume:

* **Duplicate starts** and **fully nested intervals** (`[100,200)` containing
  `[120,130)`) are the cases a naive "binary search then scan forward" gets
  wrong. The `max_end_` augmentation exists precisely for the nested case.
* An interval with `start == end` (zero-length; legal in BED for insertion
  points) overlaps nothing under the half-open predicate. That is correct, and
  should be an explicit test so nobody "fixes" it later.

### 3.4 Sweep line for BAM intersection

```cpp
// src/interval_sweep.cpp
// Preconditions, checked not assumed: the BAM is coordinate-sorted (@HD SO:)
// and `index` was built from the same contig ordering. Both are verified at
// entry; a coordinate-unsorted BAM makes the active set unbounded, which is
// the memory-exhaustion shape ProfileOptions::kMaxBins already guards against
// elsewhere in this module.
//
// Active set: a min-heap keyed by interval end. Each alignment start pops
// everything that ended before it, then every remaining member overlaps.
// O((n + m) log a) with a = max simultaneous overlap depth, one pass, no index.
```

This is the path that answers "signal over gene bodies of arbitrary length" and
"multi-region intersections", which fixed-width reference-point windowing cannot
express.

### 3.5 Wiring, and a warning about the HTTP surface

New CLI: `--annotation <bed/gtf>` alongside the existing `--regions`, plus
`--mode overlap|reference-point`. `read_regions()` already parses both formats;
`IntervalIndex::build` consumes its output directly, so no new parser.

If an `/api/overlap` endpoint is added to the dashboard, note that
`docs/SECURITY_HTTP_2026-08-15.md` finding **H1c (no authentication) is still
open**. A new endpoint that accepts caller-supplied intervals and returns
per-hit results is an unauthenticated, result-set-sized allocation — the exact
class `ProfileOptions::kMaxBins` was introduced to close. Therefore:

* Cap hits per query and total response bytes, validated in **one** place that
  every entry path calls, following the `ProfileOptions::validate()` precedent.
* Reject `qs >= qe`, negative coordinates, and unknown contig names explicitly.
* Add the endpoint to `tests/test_http_security.cpp` in the same pass, not later.
* Extend `docs/SECURITY_HTTP_2026-08-15.md` with the new surface — it is a
  living record, per `CLAUDE.md`.

Recommendation: land §3.2–3.4 with CLI access only, and treat the HTTP endpoint
as a separate change reviewed on its own.

### 3.6 Tests (`tests/test_interval_index.cpp`)

Nested, identical, adjacent (`[0,10)` vs `[10,20)` — must **not** overlap),
zero-length, duplicate starts, single interval, empty index, query entirely
left/right of all data, contig with no intervals, unknown tid. Plus a
randomised differential test against an O(n) brute-force scan over 10⁴ random
intervals — cheap, and it is what actually catches descent bugs.

---

## 4. Zero-dependency wheels for `scrna_matrix`

### 4.1 Layout

`pyproject.toml` goes in **`modules/scrna_matrix/`**, not the repo root. The
module must stay self-contained and independently buildable; a root-level
`pyproject.toml` referring down into `modules/` would make the wheel build
depend on the aggregate, which is the thing `CLAUDE.md` forbids.

```
modules/scrna_matrix/pyproject.toml            (new)
modules/scrna_matrix/python/scrna_matrix/__init__.py  (new, thin re-export)
modules/scrna_matrix/README.md                 (used as the PyPI long description)
.github/workflows/wheels.yml                   (new)
```

### 4.2 `pyproject.toml`

```toml
[build-system]
requires = ["scikit-build-core>=0.9", "pybind11>=2.12"]
build-backend = "scikit_build_core.build"

[project]
name = "scrna-matrix"
version = "0.1.0"
requires-python = ">=3.9"
description = "Cache-aligned sparse matrix engine and SIMD k-NN for single-cell data"
readme = "README.md"
license = { text = "MIT" }
dependencies = ["numpy>=1.21"]

[tool.scikit-build]
cmake.version = ">=3.20"
cmake.source-dir = "."
wheel.packages = ["python/scrna_matrix"]
build-dir = "build/{wheel_tag}"

[tool.scikit-build.cmake.define]
SCRNA_BUILD_PYTHON = "ON"
SCRNA_BUILD_TESTS  = "OFF"
SCRNA_ENABLE_HNSW  = "ON"
# Non-negotiable for a redistributable wheel: -march=native and
# SCRNA_HNSW_NATIVE_ISA bake the *build host's* ISA into the binary. The wheel
# is installed on machines we have never seen; an AVX-512 instruction reaching
# a Skylake-SP-less CPU is a SIGILL at import time. The runtime dispatch in
# simd_math.hpp (per-function target attributes) is what delivers AVX-512
# *without* this, and is precisely why it was written that way.
SCRNA_NATIVE = "OFF"
SCRNA_HNSW_NATIVE_ISA = "OFF"
```

numpy is a **runtime** dependency only — `pybind11/numpy.h` binds the NumPy C
API at import, so no numpy headers are needed at build time and there is no
oldest-supported-numpy pin to manage.

### 4.3 The OpenMP problem

`scrna_matrix` **fails configuration without OpenMP**, by design, and CI asserts
it stays a hard error. That interacts with wheels in a way that decides the
whole workflow:

* **Linux**: manylinux images ship `libgomp`. `auditwheel repair` vendors
  `libgomp.so.1` into the wheel. Works out of the box.
* **macOS**: Apple Clang ships no OpenMP runtime. `libomp` must be installed on
  the runner (`cmake/PtoOpenMP.cmake` already resolves Homebrew's keg-only
  layout — reuse it, do not write a second resolver). `delocate` then vendors
  `libomp.dylib` into the wheel.
* **macOS cross-compilation is the trap**: Homebrew's `libomp` is per-arch. An
  `arm64` runner cannot link an `x86_64` `libomp.dylib`, so `universal2` and
  `ARCHFLAGS`-style cross builds fail or, worse, produce a wheel that dyld
  rejects on the user's machine. **Build macOS wheels natively per arch** —
  `macos-13` for `x86_64`, `macos-14` for `arm64` — and ship two wheels, not a
  universal2.
* **Linux aarch64**: use a native ARM runner (`ubuntu-24.04-arm`) rather than
  QEMU. QEMU works but turns a 10-minute build into an hour and emulates the
  CPU feature detection that `simd_math.hpp` dispatches on — which is exactly
  the code you least want tested under emulation.

### 4.4 Workflow

```yaml
name: wheels
on:
  push:
    tags: ["v*"]
  workflow_dispatch:

jobs:
  build:
    name: ${{ matrix.os }} ${{ matrix.arch }}
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
        include:
          - { os: ubuntu-latest,    arch: x86_64  }
          - { os: ubuntu-24.04-arm, arch: aarch64 }
          - { os: macos-13,         arch: x86_64  }   # Intel runner
          - { os: macos-14,         arch: arm64   }   # Apple Silicon runner
    steps:
      - uses: actions/checkout@v4
        with: { submodules: recursive }   # third_party/hnswlib

      - name: Install libomp (macOS)
        if: runner.os == 'macOS'
        run: brew install libomp

      - uses: pypa/cibuildwheel@v2.21
        with:
          package-dir: modules/scrna_matrix
        env:
          CIBW_ARCHS: ${{ matrix.arch }}
          CIBW_BUILD: "cp39-* cp310-* cp311-* cp312-* cp313-*"
          CIBW_SKIP: "*-musllinux_i686 pp*"
          # Baseline ISA. Runtime dispatch does the rest; see §4.2.
          CIBW_ENVIRONMENT: >-
            SKBUILD_CMAKE_DEFINE="SCRNA_NATIVE=OFF;SCRNA_HNSW_NATIVE_ISA=OFF"
          CIBW_ENVIRONMENT_MACOS: >-
            MACOSX_DEPLOYMENT_TARGET=${{ matrix.os == 'macos-13' && '13.0' || '14.0' }}
            SKBUILD_CMAKE_DEFINE="SCRNA_NATIVE=OFF;SCRNA_HNSW_NATIVE_ISA=OFF;OpenMP_ROOT=$(brew --prefix libomp)"
          CIBW_TEST_REQUIRES: "pytest numpy scipy"
          CIBW_TEST_COMMAND: >-
            python -c "import scrna_matrix; print(scrna_matrix.simd_isa(), scrna_matrix.has_hnsw())" &&
            pytest {package}/tests/test_python_bindings.py -q

      - uses: actions/upload-artifact@v4
        with:
          name: wheels-${{ matrix.os }}-${{ matrix.arch }}
          path: wheelhouse/*.whl

  publish:
    needs: build
    if: startsWith(github.ref, 'refs/tags/v')
    runs-on: ubuntu-latest
    environment: pypi
    permissions: { id-token: write }   # Trusted Publishing; no API token in the repo
    steps:
      - uses: actions/download-artifact@v4
        with: { pattern: wheels-*, merge-multiple: true, path: dist }
      - uses: pypa/gh-action-pypi-publish@release/v1
```

The `CIBW_TEST_COMMAND` import check is not ceremony — it runs the wheel on a
fresh interpreter and is the thing that catches a missing vendored `libomp` or a
baked-in ISA. Add a dedicated job that runs the x86_64 wheel under
`SCRNA_FORCE_ISA=scalar` to prove the dispatch fallback works on a CPU without
AVX2, since no GitHub runner will give you one.

Offline note: `cibuildwheel` fetches its toolchain, and `fastq_stream`'s CMake
has a `FetchContent` fallback for libdeflate. Neither affects `scrna_matrix`
wheels (no network fetch in its CMake, hnswlib is vendored), but the wheel job
must not be presented as an offline build — it is release tooling, and the
"no network" constraint is about the *runtime* on the sequencing host.

---

## 5. Cross-cutting notes

* §1 and §3 both take a hard dependency on htslib behaviour; only §3's module
  already has it. Keep §1's optional.
* §2 and §3 each add a new place where an external file or a caller-supplied
  bound drives an allocation. Both get a validated ceiling in a single function
  every entry path calls — `ProfileOptions::validate()` is the pattern.
* Per `REVIEW_2026-08-15`'s retrospective: four of five findings were a fixed
  defect recurring in a sibling function. When landing the reverse-complement
  handling (§1.4), the half-open overlap predicate (§3.3) and the meta
  validation (§2.4), grep for the same predicate elsewhere before calling any of
  them done.
* Living documents to extend, not replace: `docs/SECURITY_HTTP_2026-08-15.md`
  (§3.5), `modules/scrna_matrix/docs/AUDIT.md` (§2.4 residual),
  `docs/REVIEW_2026-08-15.md` if any of this reopens a closed finding.

## 6. Suggested landing order

| # | Change | Depends on | Risk |
|---|---|---|---|
| 1 | §2 refactor: lift `HnswIndex` out of the two free functions, no behaviour change | — | low; existing tests are the oracle |
| 2 | §2 save/load + meta validation + Python bindings | 1 | low |
| 3 | §3 `IntervalIndex` + tests, CLI only, no HTTP | — | low; self-contained |
| 4 | §3 sweep-line BAM intersection | 3 | medium; needs sortedness preconditions checked |
| 5 | §4 wheels (independent of everything else; ship early for feedback) | 2 if the index API is to be in the first release | medium; macOS/libomp is the whole risk |
| 6 | §1 `HtsInput` single-end + pipeline seam | — | medium |
| 7 | §1 interleaved pairing + CRAM reference pinning | 6 | medium; §1.4 trap 2 is a correctness/privacy issue, not a feature |
| 8 | §3.5 HTTP `/api/overlap`, reviewed on its own against H1c | 3, 4 | high; unauthenticated surface |

Items 1–3 are safe to run in parallel with 5. Item 8 should not land until H1c
is either closed or explicitly accepted for the new endpoint.

---

## 7. Implementation record (2026-08-16)

Everything in §6 items 1–3, 5, 6 and 7 is landed and tested. What follows is
what changed against the design above and why, plus what is deliberately still
open.

### Landed

| Design | Files | Tests |
|---|---|---|
| §2 HnswIndex + serialization | `modules/scrna_matrix/include/matrix/hnsw_index.hpp`, `knn_graph.hpp` (HNSW free functions moved out), `src/python_bindings.cpp` | `tests/test_hnsw_persist.cpp` (12 cases), `tests/test_python_bindings.py` (+9) |
| §3.2–3.3 IntervalIndex | `modules/cuttag_profiler/include/profiler/interval_index.hpp` | `tests/test_interval_index.cpp` (12 cases incl. a 20k-query differential) |
| §4 wheels | `modules/scrna_matrix/pyproject.toml`, `python/scrna_matrix/__init__.py`, `.github/workflows/wheels.yml` | wheel built and installed locally; `sdist` + `scalar-fallback` jobs |
| §1 BAM/CRAM ingestion | `modules/fastq_stream/include/fastq_stream/hts_input.hpp`, `src/hts_input.cpp`, pipeline/CLI/report seams | `tests/test_hts_input.cpp` (12 cases) |

Full tree: `ctest` 9/9 pass, no warnings under `modules/*/(src|include|tests)/`.

### Corrections to the design above

1. **§2.2 was wrong about metric mismatch being a refusal.** `ScrnaCosineSpace`
   and hnswlib's `InnerProductSpace` both return `1 - <a,b>` over unit-normalised
   vectors — they are the same distance computed by different kernels. Refusing
   a mismatch would have blocked a legitimate and useful choice (load an index
   with the faster kernel available on the query host). `HnswMetric` therefore
   records the distance *family*, and `use_scrna_space` stays a free
   performance choice at load time.

2. **§2.4's validation list was incomplete, and the test caught it.** Every
   check as designed established only that the sidecar was self-consistent and
   *plausible*. A sidecar declaring `dim = 4096` for an 8-dim payload passes all
   of them — and that is exactly the input that makes hnswlib stride out of
   bounds. `validate_payload_header()` was added: it reads the leading 96 bytes
   of the hnswlib payload and requires
   `label_offset_ - offsetData_ == dim * sizeof(float)`, plus agreement on
   `max_elements_`, `cur_element_count` and `M_`. Recorded in `AUDIT.md`.
   This is the REVIEW_2026-08-15 pattern again: the obvious validation was
   written, the sibling path it did not cover was live.

3. **§3.3's half-open predicate was wrong for degenerate intervals**, also
   caught by the differential test rather than by inspection. The familiar
   `as < be && bs < ae` is equivalent to "the intersection is non-empty" only
   when both intervals are non-empty; for a zero-length interval strictly
   inside the query it reports an overlap with a feature covering no bases. The
   implementation now carries the `as < ae` term, the rule is "start >= end is
   empty and overlaps nothing" applied to both sides, and the test's reference
   states it independently as `max(as,bs) < min(ae,be)`.

4. **§4 named a `SCRNA_NATIVE` option that does not exist.** `scrna_matrix`'s
   CMake has no `-march=native` switch at all — only `SCRNA_HNSW_NATIVE_ISA`,
   already OFF by default. The wheel config sets that one and nothing else.

5. **§4 understated the sdist problem.** `modules/scrna_matrix/CMakeLists.txt`
   includes the repo-root `cmake/PtoOpenMP.cmake` with `OPTIONAL`, so it is
   present when building from a checkout but *not* from a published sdist,
   whose root is the module directory. macOS therefore sets `OpenMP_ROOT` from
   `brew --prefix libomp` explicitly, which works either way, and CI has a job
   that installs from the sdist to keep that honest.

6. **§1's `HtsInput` reports `compressed_bytes() == 0`.** htslib exposes no
   byte counter across formats; reporting a number the code cannot stand behind
   would be worse than a zero. The text report now prints the alignment format
   and the conversion counters instead of labelling a BAM "uncompressed".

7. **§1.5's `max_pending_mates` cap was dropped as unnecessary.** With input
   rejected up front unless it is name-grouped, a single-slot hold pairs
   everything a well-formed name-grouped file can present, in O(1) memory by
   construction rather than by a configurable bound nobody would know how to
   set. Records whose mate never arrives are emitted as singletons and counted.

8. **§1.6's rank sizing had to move earlier.** htslib's decode pool is sized at
   open time, so the thread budget is now computed before the source is
   constructed rather than after.

### Verified, not merely designed

- Reverse-strand restoration is asserted byte-exact against a hand-written
  expectation (`test_reverse_strand_is_restored`), and separately by requiring a
  forward and a reverse-strand record of the same read to produce identical
  FASTQ.
- The BAM seam is behaviour-preserving: `test_pipeline_matches_equivalent_fastq`
  runs the full pipeline over a BAM and over the equivalent FASTQ text and
  requires identical QC totals.
- The macOS OpenMP concern is real and confirmed: `otool -L` on a locally built
  wheel shows an absolute `/opt/homebrew/opt/libomp/lib/libomp.dylib`, which is
  what `delocate` must vendor. A wheel shipped without that repair step would
  fail to import on any machine without Homebrew's libomp.
- `FQ_ENABLE_HTS` off and on both build warning-free, and CI asserts both.

### Still open

- **§3.4 sweep-line BAM intersection.** Not started. `IntervalIndex` answers
  random queries; the streaming intersection against a coordinate-sorted BAM is
  a separate, simpler O(n+m) pass that does not need an index at all.
- **§3.5 `/api/overlap`.** Deliberately not landed. `SECURITY_HTTP_2026-08-15`
  finding H1c (no authentication) is still open, and a new endpoint taking
  caller-supplied intervals and returning a result set whose size the caller
  controls is the same memory-exhaustion shape `ProfileOptions::kMaxBins`
  exists to close. It should be reviewed on its own, with a validated hit
  ceiling in a single function every entry path calls, and with
  `tests/test_http_security.cpp` extended in the same change.
- **CLI wiring for `IntervalIndex`.** The structure and its tests are in; the
  `--annotation` / `--mode overlap` flags on `cuttag_profiler` are not.
- **`load_index` fuzzing.** See the residual in `AUDIT.md`: the checks catch
  accident, not a hostile payload.
