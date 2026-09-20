// SPDX-License-Identifier: MIT
//
// A persistable HNSW index.
//
// knn_graph.hpp holds the HNSW *primitives* -- HnswParams, the memory
// estimator, ScrnaCosineSpace and the normalisation helpers. This header holds
// the object those primitives compose into: an index you can build once, save,
// and reload to query new cells against. The two one-shot entry points
// (hnsw_build_and_search / _dense) live here too and are now thin wrappers over
// it, so there is exactly one construction path to reason about.
//
// The motivating workflow is reference mapping. Building a 1.3M-cell atlas
// index takes minutes; querying a few thousand new cells against it takes
// milliseconds. Before this header the two were inseparable, because the index
// was a local inside the search function and died on return.
//
//   HnswIndex idx(/*dim=*/50, /*capacity=*/1'300'000);
//   idx.add_dense(atlas_pca, 1'300'000);
//   idx.save_index("pbmc1M.hnsw");
//   ...
//   auto idx = HnswIndex::load_index("pbmc1M.hnsw");
//   auto g = idx.query_dense(new_cells, n_new, /*k=*/15);
//
// SECURITY: load_index() takes a file whose contents decide how much memory to
// allocate and how to walk a buffer. See the notes on validate_meta() and the
// residual recorded in docs/AUDIT.md -- the header checks here are real, but
// they do not make the vendored graph reader safe against a hostile payload.
#pragma once

#include "matrix/knn_graph.hpp"

#if defined(SCRNA_ENABLE_HNSW)

#include <array>
#include <atomic>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace scrna {

// ---------------------------------------------------------------------------
// On-disk format
// ---------------------------------------------------------------------------
//
// Two files, which must agree:
//
//   <path>        the vendored hnswlib payload, written by saveIndex()
//   <path>.meta   the header below
//
// Two files rather than one because hnswlib v0.9.0's saveIndex/loadIndex take a
// *path*, not a stream (third_party/hnswlib/hnswalg.h:725,755), so we cannot
// prepend our own header inside the same file without reimplementing its
// serialiser. If a future bump adds the stream overload, collapse these into a
// single self-describing file and bump kFormatVersion.
//
// The sidecar is not optional decoration. hnswlib::loadIndex takes the space
// from the caller and trusts it: hand it a space built for the wrong dimension
// and it computes element offsets from the wrong stride, then reads inside the
// loaded blob at those offsets. The result is out-of-bounds reads, not an
// error. `dim` in this header is what makes that unrepresentable.
//
// Every field is written little-endian, field by field, rather than by dumping
// the struct. That costs a few lines and makes the format genuinely defined
// rather than "whatever this compiler laid out today".

inline constexpr std::uint32_t kHnswFormatVersion = 1;
inline constexpr char kHnswMagic[8] = {'S', 'C', 'R', 'N', 'A', 'H', 'N', 'W'};

// Distance family. This records the *semantics*, not which kernel computed
// them: ScrnaCosineSpace and hnswlib's InnerProductSpace both return
// 1 - <a,b> over unit-normalised vectors, so an index built with one is
// correct to query with the other and `use_scrna_space` stays a free
// performance choice at load time. Only a genuinely different family (L2, say)
// would have to be refused, which is why this is an enum with room to grow and
// not a bool.
enum class HnswMetric : std::uint32_t {
  kCosineViaInnerProduct = 0,
};

// Sanity ceiling on the dimension read from a file. 1M floats per vector is
// already three orders of magnitude past any sane embedding; the point is to
// reject a garbage or hostile value before it reaches an allocation.
inline constexpr std::uint32_t kHnswMaxDim = 1u << 20;

struct HnswMeta {
  std::uint32_t format_version = kHnswFormatVersion;
  std::uint32_t dim = 0;
  std::uint64_t element_count = 0;
  std::uint64_t capacity = 0;
  std::uint32_t M = 16;
  std::uint32_t ef_construction = 200;
  std::uint32_t ef_search = 0;
  HnswMetric metric = HnswMetric::kCosineViaInnerProduct;
  // 1 when vectors were unit-normalised on insertion. Always 1 today; recorded
  // so a future raw-vector mode cannot be silently queried as if normalised.
  std::uint32_t normalized = 1;
  std::uint64_t payload_bytes = 0;
  std::uint32_t payload_crc32 = 0;
  // 0 == IEEE-754 binary32, little-endian. A host that does not match refuses
  // the file rather than reading garbage floats. x86-64 <-> aarch64 both match,
  // which is the interchange that actually matters: build the atlas on a
  // cluster node, query it on a laptop.
  std::uint32_t scalar_layout = 0;
  std::uint64_t label_count = 0;
};

namespace detail {

// CRC-32 (IEEE 802.3, reflected). Table built on first use. Present because
// scrna_matrix deliberately has no zlib dependency and a truncated or
// half-written index file is the single most likely way this format goes wrong
// in practice -- an interrupted save on a cluster filesystem.
[[nodiscard]] inline const std::array<std::uint32_t, 256>& crc32_table() {
  static const std::array<std::uint32_t, 256> table = [] {
    std::array<std::uint32_t, 256> t{};
    for (std::uint32_t i = 0; i < 256; ++i) {
      std::uint32_t c = i;
      for (int bit = 0; bit < 8; ++bit) {
        c = (c & 1u) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
      }
      t[i] = c;
    }
    return t;
  }();
  return table;
}

inline std::uint32_t crc32_update(std::uint32_t crc, const unsigned char* p,
                                  std::size_t n) noexcept {
  const auto& t = crc32_table();
  for (std::size_t i = 0; i < n; ++i) {
    crc = t[(crc ^ p[i]) & 0xFFu] ^ (crc >> 8);
  }
  return crc;
}

// CRC and byte count of a file, read in 1 MiB blocks. Throws if unreadable.
inline void crc32_and_size_of_file(const std::string& path, std::uint32_t& crc_out,
                                   std::uint64_t& size_out) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open '" + path + "' for checksumming");
  std::vector<char> buf(1 << 20);
  std::uint32_t crc = 0xFFFFFFFFu;
  std::uint64_t total = 0;
  while (in) {
    in.read(buf.data(), static_cast<std::streamsize>(buf.size()));
    const std::streamsize got = in.gcount();
    if (got <= 0) break;
    crc = crc32_update(crc, reinterpret_cast<const unsigned char*>(buf.data()),
                       static_cast<std::size_t>(got));
    total += static_cast<std::uint64_t>(got);
  }
  if (in.bad()) throw std::runtime_error("read error on '" + path + "'");
  crc_out = crc ^ 0xFFFFFFFFu;
  size_out = total;
}

// Explicit little-endian scalar I/O. Deliberately not memcpy-of-struct.
inline void put_u32(std::ostream& os, std::uint32_t v) {
  unsigned char b[4] = {static_cast<unsigned char>(v & 0xFFu),
                        static_cast<unsigned char>((v >> 8) & 0xFFu),
                        static_cast<unsigned char>((v >> 16) & 0xFFu),
                        static_cast<unsigned char>((v >> 24) & 0xFFu)};
  os.write(reinterpret_cast<const char*>(b), 4);
}

inline void put_u64(std::ostream& os, std::uint64_t v) {
  unsigned char b[8];
  for (int i = 0; i < 8; ++i) b[i] = static_cast<unsigned char>((v >> (8 * i)) & 0xFFu);
  os.write(reinterpret_cast<const char*>(b), 8);
}

inline std::uint32_t get_u32(std::istream& is) {
  unsigned char b[4];
  is.read(reinterpret_cast<char*>(b), 4);
  if (!is) throw std::runtime_error("index metadata is truncated");
  return static_cast<std::uint32_t>(b[0]) | (static_cast<std::uint32_t>(b[1]) << 8) |
         (static_cast<std::uint32_t>(b[2]) << 16) |
         (static_cast<std::uint32_t>(b[3]) << 24);
}

inline std::uint64_t get_u64(std::istream& is) {
  unsigned char b[8];
  is.read(reinterpret_cast<char*>(b), 8);
  if (!is) throw std::runtime_error("index metadata is truncated");
  std::uint64_t v = 0;
  for (int i = 0; i < 8; ++i) v |= static_cast<std::uint64_t>(b[i]) << (8 * i);
  return v;
}

// This build's scalar layout tag. Anything other than IEEE-754 binary32
// little-endian gets a distinct value so it cannot silently read a foreign
// file. Both supported targets (x86-64, aarch64) tag as 0.
[[nodiscard]] inline std::uint32_t host_scalar_layout() noexcept {
  static_assert(sizeof(float) == 4, "float is not 32-bit on this target");
  if constexpr (!std::numeric_limits<float>::is_iec559) return 0xFFFFFFFFu;
  const std::uint32_t probe = 0x01020304u;
  unsigned char bytes[4];
  std::memcpy(bytes, &probe, 4);
  return (bytes[0] == 0x04) ? 0u : 1u;  // 0 = little-endian, 1 = big-endian
}

}  // namespace detail

// ---------------------------------------------------------------------------
// HnswIndex
// ---------------------------------------------------------------------------

class HnswIndex {
 public:
  // Passed as `self_label_base` when the query rows are NOT members of the
  // index, which is the whole point of a saved reference atlas. hnswlib's
  // labeltype cannot take this value for a real element, so the self-exclusion
  // comparison in drain_into() simply never fires.
  static constexpr std::size_t kNoSelfLabel = std::numeric_limits<std::size_t>::max();

  // `capacity` is the maximum number of vectors the index will ever hold;
  // hnswlib allocates for it up front, so it is checked against
  // params.max_index_bytes here rather than after the allocation has already
  // happened.
  HnswIndex(std::size_t dim, std::size_t capacity, const HnswParams& params = {})
      : params_(params) {
    if (dim == 0) throw std::invalid_argument("HnswIndex: dim must be non-zero");
    if (dim > kHnswMaxDim) {
      throw std::invalid_argument("HnswIndex: dim " + std::to_string(dim) +
                                  " exceeds the " + std::to_string(kHnswMaxDim) +
                                  " limit");
    }
    if (capacity == 0) {
      throw std::invalid_argument("HnswIndex: capacity must be non-zero");
    }
    check_capacity(capacity);
    validate_params(params);
    check_budget(capacity, dim, params.M, params.max_index_bytes);

    meta_.dim = static_cast<std::uint32_t>(dim);
    meta_.capacity = capacity;
    meta_.M = static_cast<std::uint32_t>(params.M);
    meta_.ef_construction = static_cast<std::uint32_t>(params.ef_construction);
    meta_.scalar_layout = detail::host_scalar_layout();

    make_space(dim, params.use_scrna_space);
    index_ = std::make_unique<hnswlib::HierarchicalNSW<float>>(
        space_.get(), capacity, params.M, params.ef_construction, params.random_seed);
    set_ef(params.ef_search);
  }

  HnswIndex(const HnswIndex&) = delete;
  HnswIndex& operator=(const HnswIndex&) = delete;
  // Movable: load_index() returns by value, and the unique_ptr members make the
  // defaulted move correct -- index_ holds a raw SpaceInterface* into space_,
  // and moving both pointers together keeps that pointer valid.
  HnswIndex(HnswIndex&&) noexcept = default;
  HnswIndex& operator=(HnswIndex&&) noexcept = default;

  [[nodiscard]] std::size_t dim() const noexcept { return meta_.dim; }
  [[nodiscard]] std::size_t capacity() const noexcept {
    return static_cast<std::size_t>(meta_.capacity);
  }
  [[nodiscard]] std::size_t size() const noexcept {
    if (!index_) return 0;
    // cur_element_count is a std::atomic<size_t> in hnswlib v0.9.0; load it
    // explicitly rather than letting it convert inside a ternary.
    return static_cast<std::size_t>(index_->cur_element_count.load());
  }
  [[nodiscard]] const HnswMeta& meta() const noexcept { return meta_; }
  [[nodiscard]] const std::vector<std::string>& labels() const noexcept {
    return labels_;
  }

  // Optional per-element names (cell barcodes), label-indexed. Serialized with
  // the index so a reloaded atlas can answer with names, not bare integers.
  void set_labels(std::vector<std::string> names) {
    const std::lock_guard<std::mutex> lock(shared_->write_mutex);
    if (!names.empty() && names.size() != size()) {
      throw std::invalid_argument(
          "HnswIndex::set_labels: got " + std::to_string(names.size()) +
          " names for " + std::to_string(size()) + " elements");
    }
    labels_ = std::move(names);
  }

  // 0 selects the same default the one-shot path uses: max(2k, 64) at query
  // time. Stored, not applied, until a k is known -- and since SM9 it is never
  // applied to the index at all: each query passes its own ef to searchKnn, so
  // nothing writes hnswlib's shared ef_ while a search is reading it.
  void set_ef(std::size_t ef) {
    // The header stores ef_search as uint32; a larger value was truncated into
    // it silently, so the saved index queried with a different breadth.
    if (ef > std::numeric_limits<std::uint32_t>::max()) {
      throw std::invalid_argument("HnswIndex: ef_search " + std::to_string(ef) +
                                  " exceeds the uint32 range the index header stores");
    }
    const std::lock_guard<std::mutex> lock(shared_->write_mutex);
    meta_.ef_search = static_cast<std::uint32_t>(ef);
    // Published for the query path to read; queries never write it back.
    shared_->query_ef.store(static_cast<std::uint32_t>(ef), std::memory_order_relaxed);
  }

  // ------------------------------------------------------------------ build

  // Appends `n` rows of a dense row-major (n x dim) matrix. Labels are assigned
  // sequentially from the current size, so repeated calls extend the index.
  void add_dense(const float* data, std::size_t n, int num_threads = 0) {
    if (n == 0) return;
    if (data == nullptr) throw std::invalid_argument("add_dense: null data");
    // Held across base, capacity check and insertion; see write_mutex_.
    const std::lock_guard<std::mutex> lock(shared_->write_mutex);
    const std::size_t base = size();
    require_room(n);
    const std::size_t d = dim();

    // Non-finite input is rejected rather than normalised away: a NaN row
    // yields norm_sq = NaN, fails the `> 0` test, and would be indexed as an
    // all-zero vector that matches nothing -- silently.
    detail::parallel_for_checked(n, resolve_threads(num_threads),
                                 "HNSW index construction", [&](std::size_t i) {
      thread_local std::vector<float> unit;
      unit.resize(d);
      detail::normalize_dense_row(data + i * d, d, unit.data(), /*check_finite=*/true, i);
      index_->addPoint(unit.data(), static_cast<hnswlib::labeltype>(base + i));
    });
    meta_.element_count = size();
  }

  // Appends every row of a Block_CSR, densifying on the fly. The index holds
  // rows x dim floats regardless of input sparsity; see the memory note in
  // knn_graph.hpp before using this on gene-space input.
  void add_csr(const Block_CSR<float>& mat, int num_threads = 0) {
    const std::size_t n = mat.rows();
    if (n == 0) return;
    if (mat.cols() != dim()) {
      throw std::invalid_argument("add_csr: matrix has " + std::to_string(mat.cols()) +
                                  " columns but the index has dim " +
                                  std::to_string(dim()));
    }
    const std::lock_guard<std::mutex> lock(shared_->write_mutex);
    const std::size_t base = size();
    require_room(n);
    const std::size_t d = dim();

    detail::parallel_for_checked(n, resolve_threads(num_threads),
                                 "HNSW index construction", [&](std::size_t i) {
      thread_local std::vector<float> dense;
      if (dense.size() != d) dense.assign(d, 0.0f);
      const auto idx = mat.row_indices(i);
      const auto val = mat.row_values(i);
      detail::densify_normalized(idx, val, dense.data());
      // Scratch is cleared even on throw: an un-cleared buffer would leak this
      // row's values into the next row on the same thread and silently corrupt
      // every later insertion.
      try {
        index_->addPoint(dense.data(), static_cast<hnswlib::labeltype>(base + i));
      } catch (...) {
        detail::clear_dense_row(idx, dense.data());
        throw;
      }
      detail::clear_dense_row(idx, dense.data());
    });
    meta_.element_count = size();
  }

  // ------------------------------------------------------------------ query

  // Queries `n` dense rows against the index.
  //
  // `self_label_base`: when the query rows ARE index members (row i is label
  // self_label_base + i), that label is excluded from its own neighbour list
  // and k+1 candidates are fetched to compensate. Pass kNoSelfLabel -- the
  // default -- when querying foreign cells against a reference atlas, which is
  // the case that must NOT drop a legitimate neighbour just because its label
  // happens to equal the query row's position.
  [[nodiscard]] KnnGraph query_dense(const float* data, std::size_t n, std::size_t k,
                                     int num_threads = 0,
                                     std::size_t self_label_base = kNoSelfLabel) const {
    const std::size_t slots = detail::checked_graph_slots(n, k);  // SM1
    KnnGraph graph;
    graph.k = k;
    graph.neighbor_idx.assign(slots, -1);
    graph.neighbor_sim.assign(slots, 0.0f);
    if (n == 0 || k == 0) return graph;
    if (data == nullptr) throw std::invalid_argument("query_dense: null data");

    if (size() == 0) return graph;  // nothing to find; every slot stays -1

    const std::size_t d = dim();
    const bool exclude_self = self_label_base != kNoSelfLabel;
    const std::size_t ef = query_ef(k);
    const std::size_t fetch = fetch_count(k, exclude_self);

    detail::parallel_for_checked(n, resolve_threads(num_threads), "HNSW query",
                                 [&](std::size_t i) {
      thread_local std::vector<float> unit;
      unit.resize(d);
      detail::normalize_dense_row(data + i * d, d, unit.data(), /*check_finite=*/true, i);
      auto pq = index_->searchKnnWithEf(unit.data(), fetch, ef);
      detail::drain_into(pq, exclude_self ? self_label_base + i : kNoSelfLabel, k,
                         graph.neighbor_idx.data() + i * k,
                         graph.neighbor_sim.data() + i * k);
    });
    return graph;
  }

  [[nodiscard]] KnnGraph query_csr(const Block_CSR<float>& mat, std::size_t k,
                                   int num_threads = 0,
                                   std::size_t self_label_base = kNoSelfLabel) const {
    const std::size_t n = mat.rows();
    const std::size_t slots = detail::checked_graph_slots(n, k);  // SM1
    KnnGraph graph;
    graph.k = k;
    graph.neighbor_idx.assign(slots, -1);
    graph.neighbor_sim.assign(slots, 0.0f);
    if (n == 0 || k == 0) return graph;
    if (mat.cols() != dim()) {
      throw std::invalid_argument("query_csr: matrix has " + std::to_string(mat.cols()) +
                                  " columns but the index has dim " +
                                  std::to_string(dim()));
    }

    if (size() == 0) return graph;  // nothing to find; every slot stays -1

    const std::size_t d = dim();
    const bool exclude_self = self_label_base != kNoSelfLabel;
    const std::size_t ef = query_ef(k);
    const std::size_t fetch = fetch_count(k, exclude_self);

    detail::parallel_for_checked(n, resolve_threads(num_threads), "HNSW query",
                                 [&](std::size_t i) {
      thread_local std::vector<float> dense;
      if (dense.size() != d) dense.assign(d, 0.0f);
      const auto idx = mat.row_indices(i);
      const auto val = mat.row_values(i);
      detail::densify_normalized(idx, val, dense.data());
      try {
        auto pq = index_->searchKnnWithEf(dense.data(), fetch, ef);
        detail::drain_into(pq, exclude_self ? self_label_base + i : kNoSelfLabel, k,
                           graph.neighbor_idx.data() + i * k,
                           graph.neighbor_sim.data() + i * k);
      } catch (...) {
        detail::clear_dense_row(idx, dense.data());
        throw;
      }
      detail::clear_dense_row(idx, dense.data());
    });
    return graph;
  }

  // ------------------------------------------------------------ persistence

  // Writes <path> (hnswlib payload) and <path>.meta (our header).
  //
  // Order matters: the payload is written first, then checksummed, then the
  // meta is written. A crash between the two leaves a payload with no meta,
  // which load_index() refuses -- rather than a meta that describes a file that
  // was never finished.
  void save_index(const std::string& path) const {
    // Serialised with add_*: indexFileSize() below must describe the same
    // graph saveIndex() wrote.
    const std::lock_guard<std::mutex> lock(shared_->write_mutex);

    // Labels set before a later add_*() describe a prefix of the index. The
    // meta used to record them as they were, save returned, and load_index()
    // refused the pair ("carries 10 labels for 20 elements") -- a failure
    // discovered after the build that produced it is gone (SM5). Refuse here,
    // while set_labels() can still fix it.
    if (!labels_.empty() && labels_.size() != size()) {
      throw std::invalid_argument(
          "save_index: the index holds " + std::to_string(size()) + " elements but " +
          std::to_string(labels_.size()) +
          " labels; elements were added after set_labels(). Call set_labels() "
          "with one name per element, or with an empty list, before saving.");
    }

    const std::string meta_path = path + ".meta";
    // A stale sidecar from an earlier save to this path would otherwise sit
    // beside a payload this call fails to finish. Without one, a failure at any
    // point below leaves a payload with no meta, which load_index() refuses.
    std::remove(meta_path.c_str());

    index_->saveIndex(path);

    HnswMeta m = meta_;
    m.element_count = size();
    m.label_count = labels_.size();
    detail::crc32_and_size_of_file(path, m.payload_crc32, m.payload_bytes);

    // hnswlib's saveIndex() writes through an unchecked std::ofstream, so a
    // write that fails partway (disk full, quota, RLIMIT_FSIZE -- which Python
    // turns into EFBIG by ignoring SIGXFSZ) returned normally. The checksum
    // above is of whatever reached the disk, so the meta then SEALED the short
    // payload as intact: save succeeded, and the CRC that was meant to catch
    // truncation matched it (SM4). The byte count hnswlib itself computes for
    // this graph is the independent figure to hold the file to.
    const std::uint64_t expected = index_->indexFileSize();
    if (m.payload_bytes != expected) {
      throw std::runtime_error(
          "save_index: wrote " + std::to_string(m.payload_bytes) + " of " +
          std::to_string(expected) + " bytes to '" + path +
          "' (disk full, quota or file-size limit?). No .meta was written, so "
          "the partial index cannot be loaded by mistake.");
    }

    std::ofstream out(meta_path, std::ios::binary | std::ios::trunc);
    if (!out) throw std::runtime_error("cannot write '" + meta_path + "'");
    out.write(kHnswMagic, sizeof(kHnswMagic));
    detail::put_u32(out, m.format_version);
    detail::put_u32(out, m.dim);
    detail::put_u64(out, m.element_count);
    detail::put_u64(out, m.capacity);
    detail::put_u32(out, m.M);
    detail::put_u32(out, m.ef_construction);
    detail::put_u32(out, m.ef_search);
    detail::put_u32(out, static_cast<std::uint32_t>(m.metric));
    detail::put_u32(out, m.normalized);
    detail::put_u64(out, m.payload_bytes);
    detail::put_u32(out, m.payload_crc32);
    detail::put_u32(out, m.scalar_layout);
    detail::put_u64(out, m.label_count);
    for (int i = 0; i < 4; ++i) detail::put_u32(out, 0);  // reserved
    for (const auto& s : labels_) {
      detail::put_u64(out, s.size());
      out.write(s.data(), static_cast<std::streamsize>(s.size()));
    }
    out.close();
    if (!out) throw std::runtime_error("write error on '" + meta_path + "'");
  }

  // `params` supplies the load-time *choices* -- thread count, ef_search
  // override, the memory ceiling to enforce, and which distance kernel to use.
  // Everything that must match the file (dim, M, metric) comes from the file.
  //
  // `verify_checksum` re-reads the payload to check its CRC. Worth the pass on
  // anything you did not just write; it is the only thing that catches a
  // truncated save on a cluster filesystem. On a multi-GB index it costs a
  // sequential read.
  [[nodiscard]] static HnswIndex load_index(const std::string& path,
                                            const HnswParams& params = {},
                                            bool verify_checksum = true) {
    const HnswMeta m = read_and_validate_meta(path, params, verify_checksum);
    return HnswIndex(path, m, params);
  }

 private:
  // Private ctor used only by load_index, after the meta has been validated.
  HnswIndex(const std::string& path, const HnswMeta& m, const HnswParams& params)
      : params_(params), meta_(m) {
    make_space(m.dim, params.use_scrna_space);
    index_ = std::make_unique<hnswlib::HierarchicalNSW<float>>(
        space_.get(), path, /*nmslib=*/false, static_cast<std::size_t>(m.capacity));
    // A caller-supplied ef_search overrides the stored one; 0 keeps the file's.
    set_ef(params.ef_search ? params.ef_search : m.ef_search);
    labels_ = read_labels(path + ".meta", m);
  }

  void make_space(std::size_t d, bool use_scrna_space) {
    if (use_scrna_space) {
      space_ = std::make_unique<detail::ScrnaCosineSpace>(d);
    } else {
      space_ = std::make_unique<hnswlib::InnerProductSpace>(d);
    }
  }

  // hnswlib's own M bounds. It does not reject anything outside them.
  //   M < 2      mult_ = 1 / log(M) is +inf at M=1, and getRandomLevel()'s
  //              (int)inf is undefined behaviour (UBSan); at M=0, maxM0_ is 0
  //              and mutuallyConnectNewElement() indexes a vector at -4.
  //   M > 10000  capped to 10000 with a stderr warning, while meta_.M and the
  //              memory budget keep the requested value -- so the index built,
  //              saved, and then refused to load: "declares M=20000 but the
  //              payload was built with M=10000".
  // docs/AUDIT_2026-09-11_scrna_matrix.md, SM3.
  static constexpr std::size_t kMinM = 2;
  static constexpr std::size_t kMaxM = 10000;

  static void validate_params(const HnswParams& p) {
    if (p.M < kMinM || p.M > kMaxM) {
      throw std::invalid_argument("HnswIndex: M=" + std::to_string(p.M) +
                                  " is outside [" + std::to_string(kMinM) + ", " +
                                  std::to_string(kMaxM) + "]");
    }
    // Both are stored in the uint32 header; ef_search is checked in set_ef().
    if (p.ef_construction > std::numeric_limits<std::uint32_t>::max()) {
      throw std::invalid_argument("HnswIndex: ef_construction " +
                                  std::to_string(p.ef_construction) +
                                  " exceeds the uint32 range the index header stores");
    }
  }

  // A KnnGraph stores neighbour ids as int32 (drain_into narrows hnswlib's
  // labeltype into it), so an index able to hold more than INT32_MAX elements
  // could answer with a WRAPPED, negative id -- indistinguishable from the -1
  // "no neighbour" sentinel. The memory budget makes that unreachable in
  // practice (~300 GB), which is why it was first left as a residual; it is
  // guarded here instead, because "unreachable" is a property of a default
  // that callers may raise. docs/AUDIT_2026-09-11_scrna_matrix.md, SM10.
  static constexpr std::size_t kMaxCapacity =
      static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max());

  static void check_capacity(std::size_t capacity) {
    if (capacity > kMaxCapacity) {
      throw std::invalid_argument(
          "HnswIndex: capacity " + std::to_string(capacity) + " exceeds " +
          std::to_string(kMaxCapacity) +
          "; neighbour ids are int32 in a KnnGraph, so a larger index could "
          "return a wrapped id. Shard the atlas instead.");
    }
  }

  static void check_budget(std::size_t n, std::size_t d, std::size_t M,
                           std::size_t limit) {
    const std::size_t needed = hnsw_index_bytes(n, d, M);
    if (needed > limit) {
      throw std::invalid_argument(
          "HNSW index would need ~" + std::to_string(needed >> 20) + " MB for " +
          std::to_string(n) + " points x " + std::to_string(d) + " dims, over the " +
          std::to_string(limit >> 20) +
          " MB limit. Reduce dimensionality (PCA/HVG) or raise max_index_bytes.");
    }
  }

  void require_room(std::size_t n) const {
    const std::size_t have = size();
    if (n > capacity() || have > capacity() - n) {
      throw std::invalid_argument(
          "HnswIndex: adding " + std::to_string(n) + " points to an index holding " +
          std::to_string(have) + " would exceed its capacity of " +
          std::to_string(capacity()));
    }
  }

  [[nodiscard]] int resolve_threads(int requested) const {
    int nt = requested > 0 ? requested : params_.num_threads;
#ifdef _OPENMP
    if (nt <= 0) nt = omp_get_max_threads();
#else
    nt = 1;
#endif
    return nt;
  }

  // The search breadth for one query, computed and PASSED rather than stored.
  //
  // This used to be `apply_query_ef(k)`, which called index_->setEf() -- a
  // write to hnswlib's shared ef_ immediately before a search that reads it.
  // Two concurrent queries with different k therefore raced, which TSan
  // reported as a write in query_dense against a read in searchKnn
  // (docs/AUDIT_2026-09-11_scrna_matrix.md, SM9). hnswlib's searchKnn now
  // takes the value as an argument (third_party/hnswlib/VERSION.txt, patch 4),
  // so the query path reads shared state and never writes it.
  //
  // ef must be at least k or the search returns fewer than k candidates. Same
  // default as before: max(2k, 64) unless explicitly set.
  [[nodiscard]] std::size_t query_ef(std::size_t k) const {
    const std::size_t stored = shared_->query_ef.load(std::memory_order_relaxed);
    const std::size_t ef = stored ? stored : std::max<std::size_t>(2 * k, 64);
    return std::max(ef, k < std::numeric_limits<std::size_t>::max() ? k + 1 : k);
  }

  // How many candidates to ask hnswlib for.
  //
  // `k` is the caller's row STRIDE and is deliberately not clamped -- a row
  // with fewer than k available neighbours keeps the -1 sentinel, which is the
  // documented contract (see KnnGraph, test_knn_k_exceeds_rows). What is
  // clamped is the request: asking for more candidates than the index holds
  // makes hnswlib size a heap for a result set that cannot exist. One extra
  // candidate is fetched when the query rows are index members, to pay for
  // dropping each row's own label.
  [[nodiscard]] std::size_t fetch_count(std::size_t k, bool exclude_self) const {
    const std::size_t held = size();
    std::size_t fetch = std::min(k, held);
    if (exclude_self && fetch < held) ++fetch;
    return fetch;
  }

  // Reads the sidecar and rejects anything that would make the payload load
  // unsafe or the results wrong. Ordered cheapest-first, and every check
  // happens BEFORE hnswlib is handed the file.
  static HnswMeta read_and_validate_meta(const std::string& path,
                                         const HnswParams& params,
                                         bool verify_checksum) {
    const std::string meta_path = path + ".meta";
    std::ifstream in(meta_path, std::ios::binary);
    if (!in) {
      throw std::runtime_error(
          "cannot open '" + meta_path +
          "'. An index saved by this library is a pair of files; the .meta "
          "sidecar carries the dimension and metric that make the payload safe "
          "to load.");
    }
    char magic[sizeof(kHnswMagic)];
    in.read(magic, sizeof(magic));
    if (!in || std::memcmp(magic, kHnswMagic, sizeof(magic)) != 0) {
      throw std::runtime_error("'" + meta_path + "' is not a scrna_matrix index");
    }

    HnswMeta m;
    m.format_version = detail::get_u32(in);
    if (m.format_version != kHnswFormatVersion) {
      throw std::runtime_error("index format version " +
                               std::to_string(m.format_version) +
                               " is not supported by this build (expected " +
                               std::to_string(kHnswFormatVersion) + ")");
    }
    m.dim = detail::get_u32(in);
    m.element_count = detail::get_u64(in);
    m.capacity = detail::get_u64(in);
    m.M = detail::get_u32(in);
    m.ef_construction = detail::get_u32(in);
    m.ef_search = detail::get_u32(in);
    const std::uint32_t metric_raw = detail::get_u32(in);
    m.normalized = detail::get_u32(in);
    m.payload_bytes = detail::get_u64(in);
    m.payload_crc32 = detail::get_u32(in);
    m.scalar_layout = detail::get_u32(in);
    m.label_count = detail::get_u64(in);
    for (int i = 0; i < 4; ++i) (void)detail::get_u32(in);  // reserved

    if (metric_raw != static_cast<std::uint32_t>(HnswMetric::kCosineViaInnerProduct)) {
      throw std::runtime_error("index uses distance metric " +
                               std::to_string(metric_raw) +
                               ", which this build does not implement");
    }
    m.metric = HnswMetric::kCosineViaInnerProduct;

    if (m.normalized != 1) {
      throw std::runtime_error(
          "index was not built with unit-normalised vectors; this build only "
          "queries normalised indexes and would return wrong similarities");
    }
    if (m.scalar_layout != detail::host_scalar_layout()) {
      throw std::runtime_error(
          "index was written with a different float layout (file tag " +
          std::to_string(m.scalar_layout) + ", this host " +
          std::to_string(detail::host_scalar_layout()) +
          "); rebuild it on this architecture");
    }
    if (m.dim == 0 || m.dim > kHnswMaxDim) {
      throw std::runtime_error("index declares an implausible dimension of " +
                               std::to_string(m.dim));
    }
    // The payload's own mult_ is loaded rather than recomputed, so an M of 0
    // or 1 here reaches the same undefined level draw on the next add_*().
    if (m.M < kMinM || m.M > kMaxM) {
      throw std::runtime_error("index declares M=" + std::to_string(m.M) +
                               ", outside [" + std::to_string(kMinM) + ", " +
                               std::to_string(kMaxM) + "]");
    }
    // A file may declare a capacity this build cannot represent in a KnnGraph.
    check_capacity(static_cast<std::size_t>(m.capacity));
    if (m.capacity == 0 || m.element_count > m.capacity) {
      throw std::runtime_error("index declares " + std::to_string(m.element_count) +
                               " elements in a capacity of " +
                               std::to_string(m.capacity));
    }
    // The same ceiling that gates a *built* index, applied to a *loaded* one.
    // Without this, a meta claiming 1e12 elements sends the allocation straight
    // to hnswlib.
    check_budget(static_cast<std::size_t>(m.capacity), m.dim, m.M,
                 params.max_index_bytes);

    // Bind the declared dimension to what the payload actually contains. Every
    // check above only established that `dim` is *plausible*; a sidecar saying
    // 4096 for an 8-dim payload passes all of them, and is exactly the input
    // that makes hnswlib stride through the loaded blob out of bounds. This is
    // the check that closes it.
    validate_payload_header(path, m);

    std::uint32_t crc = 0;
    std::uint64_t bytes = 0;
    if (verify_checksum) {
      detail::crc32_and_size_of_file(path, crc, bytes);
      if (bytes != m.payload_bytes) {
        throw std::runtime_error("index payload '" + path + "' is " +
                                 std::to_string(bytes) + " bytes but the metadata "
                                 "says " + std::to_string(m.payload_bytes) +
                                 " (truncated or overwritten)");
      }
      if (crc != m.payload_crc32) {
        throw std::runtime_error("index payload '" + path +
                                 "' fails its checksum (corrupt)");
      }
    }
    return m;
  }

  // Cross-checks our sidecar against the leading fields of the hnswlib payload.
  //
  // hnswlib::saveIndex writes a fixed prefix of scalars before the graph body
  // (third_party/hnswlib/hnswalg.h:725). Its loadIndex reads them back but
  // never checks them against the space the caller supplied -- data_size_ comes
  // from us, size_data_per_element_ comes from the file, and a disagreement is
  // silent. Reading 96 bytes here and insisting the two agree is what turns a
  // mismatched sidecar from an out-of-bounds read into an exception.
  //
  // If the vendored hnswlib is bumped, re-derive these offsets from its
  // saveIndex; that is why each one is spelled out rather than computed.
  static void validate_payload_header(const std::string& path, const HnswMeta& m) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open index payload '" + path + "'");
    unsigned char buf[96];
    in.read(reinterpret_cast<char*>(buf), sizeof(buf));
    if (!in) {
      throw std::runtime_error("index payload '" + path +
                               "' is too short to contain a header");
    }
    const auto u64_at = [&buf](std::size_t off) {
      std::uint64_t v = 0;
      for (int i = 0; i < 8; ++i) v |= static_cast<std::uint64_t>(buf[off + i]) << (8 * i);
      return v;
    };
    const std::uint64_t max_elements = u64_at(8);
    const std::uint64_t cur_count = u64_at(16);
    const std::uint64_t label_offset = u64_at(32);
    const std::uint64_t offset_data = u64_at(40);
    const std::uint64_t M = u64_at(72);

    // label_offset_ == offsetData_ + data_size_, and data_size_ is dim floats.
    // This is the one relation that pins the payload's stride to our dim.
    if (label_offset < offset_data ||
        (label_offset - offset_data) != static_cast<std::uint64_t>(m.dim) * sizeof(float)) {
      throw std::runtime_error(
          "index metadata declares dim " + std::to_string(m.dim) + " (" +
          std::to_string(static_cast<std::uint64_t>(m.dim) * sizeof(float)) +
          " bytes per vector) but the payload stores " +
          std::to_string(label_offset - offset_data) +
          " bytes per vector; the two files do not belong together");
    }
    if (max_elements != m.capacity || cur_count != m.element_count) {
      throw std::runtime_error(
          "index metadata and payload disagree on size (metadata: " +
          std::to_string(m.element_count) + "/" + std::to_string(m.capacity) +
          ", payload: " + std::to_string(cur_count) + "/" +
          std::to_string(max_elements) + ")");
    }
    if (M != m.M) {
      throw std::runtime_error("index metadata declares M=" + std::to_string(m.M) +
                               " but the payload was built with M=" +
                               std::to_string(M));
    }
  }

  static std::vector<std::string> read_labels(const std::string& meta_path,
                                              const HnswMeta& m) {
    std::vector<std::string> out;
    if (m.label_count == 0) return out;
    if (m.label_count != m.element_count) {
      throw std::runtime_error("index carries " + std::to_string(m.label_count) +
                               " labels for " + std::to_string(m.element_count) +
                               " elements");
    }
    std::ifstream in(meta_path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot reopen '" + meta_path + "'");
    // Fixed-size header: magic + 7 u32 + 3 u64 + 2 u32 + 1 u64 + 4 reserved u32.
    constexpr std::streamoff kHeaderBytes =
        sizeof(kHnswMagic) + 4 * 7 + 8 * 3 + 4 * 2 + 8 + 4 * 4;
    in.seekg(kHeaderBytes, std::ios::beg);
    out.reserve(static_cast<std::size_t>(m.label_count));
    for (std::uint64_t i = 0; i < m.label_count; ++i) {
      const std::uint64_t len = detail::get_u64(in);
      // A label length is attacker-controlled in the same sense the rest of the
      // header is; cap it so a garbage value cannot request a huge allocation.
      if (len > (1u << 20)) {
        throw std::runtime_error("index label " + std::to_string(i) +
                                 " declares an implausible length");
      }
      std::string s(static_cast<std::size_t>(len), '\0');
      in.read(s.data(), static_cast<std::streamsize>(len));
      if (!in) throw std::runtime_error("index metadata is truncated in its labels");
      out.push_back(std::move(s));
    }
    return out;
  }

  HnswParams params_;
  HnswMeta meta_{};
  // DECLARATION ORDER IS LOAD-BEARING. HierarchicalNSW keeps a raw
  // SpaceInterface* and dereferences it (clear() -> the dist func params) from
  // its own destructor. space_ must be declared BEFORE index_ so that it is
  // destroyed AFTER it. Reordering these two lines is a use-after-free that
  // will not show up until a destructor runs under a sanitizer.
  std::unique_ptr<hnswlib::SpaceInterface<float>> space_;
  std::unique_ptr<hnswlib::HierarchicalNSW<float>> index_;
  std::vector<std::string> labels_;

  // Mutable shared state, behind one pointer so the class stays movable
  // (load_index returns by value, and neither std::mutex nor std::atomic is).
  struct Shared {
    std::mutex write_mutex;
    // The configured query breadth, or 0 for "max(2k, 64) per query". Written
    // only by set_ef, read by every query -- atomic because a query may run
    // concurrently with a caller changing it.
    std::atomic<std::uint32_t> query_ef{0};
  };
  std::unique_ptr<Shared> shared_ = std::make_unique<Shared>();
  // write_mutex serialises every mutation: add_dense, add_csr, set_labels, save_index.
  // add_*() reads size() as the first new label, and hnswlib's addPoint() on a
  // label that already exists UPDATES that element in place rather than
  // failing. Two concurrent add_dense() calls of 2,000 rows each -- reachable
  // from Python, which releases the GIL for both -- therefore took the same
  // base, and the index held 2,000 elements with no exception (SM2). Queries
  // do not take it: hnswlib's searchKnn is safe against a concurrent addPoint,
  // and serialising them would forfeit parallel query throughput.
};

// ---------------------------------------------------------------------------
// One-shot entry points
// ---------------------------------------------------------------------------
//
// Build an index over the input and immediately query every row of it against
// itself, excluding self-matches. Both are now wrappers over HnswIndex so the
// construction, normalisation and budget-check logic exists once.

// Approximate k-NN over an already-dense row-major matrix (n x dim), e.g. a PCA
// embedding. This is the path that scales to 1M+ cells; see the memory note in
// knn_graph.hpp for why gene-space input usually should not be used at that
// size.
[[nodiscard]] inline KnnGraph hnsw_build_and_search_dense(const float* data, std::size_t n,
                                                          std::size_t dim, std::size_t k,
                                                          const HnswParams& params = {}) {
  const std::size_t slots = detail::checked_graph_slots(n, k);  // SM1 -- overflow check only
  // query_dense() below builds its OWN complete KnnGraph, sentinels included,
  // on every path that reaches it -- so a graph built here would exist only to
  // be thrown away the instant this function returns query_dense()'s. Building
  // and filling n*k slots (two full-size allocations plus O(n*k) writes) just
  // to discard them was exactly that on every call where n, k and dim were all
  // non-zero, i.e. every call that was not this early-exit. Finding 19
  // (2026-09-16 audit).
  if (n == 0 || k == 0 || dim == 0) {
    KnnGraph graph;
    graph.k = k;
    graph.neighbor_idx.assign(slots, -1);
    graph.neighbor_sim.assign(slots, 0.0f);
    return graph;
  }

  HnswIndex index(dim, n, params);
  index.add_dense(data, n);
  return index.query_dense(data, n, k, params.num_threads, /*self_label_base=*/0);
}

// Approximate k-NN directly over a Block_CSR, densifying each row on the fly.
//
// Memory: the index holds n_rows x n_cols floats regardless of input sparsity.
// For wide gene-space matrices this is usually the wrong call -- reduce
// dimensionality first and use hnsw_build_and_search_dense. The max_index_bytes
// guard will refuse rather than let the process thrash into swap.
[[nodiscard]] inline KnnGraph hnsw_build_and_search(const Block_CSR<float>& mat,
                                                    std::size_t k,
                                                    const HnswParams& params = {}) {
  const std::size_t n = mat.rows();
  const std::size_t dim = mat.cols();

  const std::size_t slots = detail::checked_graph_slots(n, k);  // SM1 -- overflow check only
  // See hnsw_build_and_search_dense's comment: query_csr() below builds its
  // own complete KnnGraph on every path that reaches it, so a graph built
  // here would exist only to be discarded on every call that was not this
  // early-exit. Finding 19 (2026-09-16 audit).
  if (n == 0 || k == 0 || dim == 0) {
    KnnGraph graph;
    graph.k = k;
    graph.neighbor_idx.assign(slots, -1);
    graph.neighbor_sim.assign(slots, 0.0f);
    return graph;
  }

  // Preserve the CSR-specific budget message: it is the one that tells a user
  // their *gene-space* input is the problem and that sparsity does not help.
  const std::size_t needed = hnsw_index_bytes(n, dim, params.M);
  if (needed > params.max_index_bytes) {
    throw std::invalid_argument(
        "HNSW index would need ~" + std::to_string(needed >> 20) + " MB for " +
        std::to_string(n) + " cells x " + std::to_string(dim) +
        " genes (hnswlib stores dense vectors, so input sparsity does not help). "
        "Reduce dimensionality (PCA/HVG) and use hnsw_build_and_search_dense, or "
        "raise max_index_bytes past " +
        std::to_string(params.max_index_bytes >> 20) + " MB.");
  }

  HnswIndex index(dim, n, params);
  index.add_csr(mat);
  return index.query_csr(mat, k, params.num_threads, /*self_label_base=*/0);
}

}  // namespace scrna

#endif  // SCRNA_ENABLE_HNSW
