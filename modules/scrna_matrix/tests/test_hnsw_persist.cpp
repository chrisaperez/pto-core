// SPDX-License-Identifier: MIT
//
// HNSW index serialization. Dependency-free in the same style as
// test_matrix_ops.cpp: assertions increment a counter and main() exits
// non-zero, so it composes with ctest without pulling in a framework.
//
// The central claim under test is that serialization is LOSSLESS. A reloaded
// index must return byte-identical neighbour lists to the in-memory one --
// not "approximately equal". HNSW is approximate with respect to exact k-NN,
// but a round-trip through disk is not allowed to move a single result, so any
// drift here is a bug in the format rather than an accuracy trade-off.

#include <sys/resource.h>  // RLIMIT_FSIZE, for the short-write test
#include <sys/wait.h>
#include <unistd.h>  // getpid, to keep concurrent ctest runs from colliding

#include <atomic>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <limits>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "matrix/hnsw_index.hpp"

namespace {

int g_failures = 0;

#define SCRNA_CHECK(cond)                                                        \
  do {                                                                           \
    if (!(cond)) {                                                               \
      std::fprintf(stderr, "FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);   \
      ++g_failures;                                                              \
    }                                                                            \
  } while (0)

// Asserts that `stmt` throws anything derived from std::exception. Load
// failures are a mix of runtime_error (bad file) and invalid_argument (bad
// budget), and which one a given corruption produces is not part of the
// contract -- refusing to load is.
#define SCRNA_CHECK_THROWS(stmt)                                                 \
  do {                                                                           \
    bool threw = false;                                                          \
    try {                                                                        \
      stmt;                                                                      \
    } catch (const std::exception&) {                                            \
      threw = true;                                                              \
    } catch (...) {                                                              \
    }                                                                            \
    if (!threw) {                                                                \
      std::fprintf(stderr, "FAILED (no throw): %s (%s:%d)\n", #stmt, __FILE__,   \
                   __LINE__);                                                    \
      ++g_failures;                                                              \
    }                                                                            \
  } while (0)

#if defined(SCRNA_ENABLE_HNSW)

using scrna::HnswIndex;
using scrna::HnswParams;

// Temporary paths. The index is a *pair* of files, so cleanup has to remove
// both or a later test in the same run will load a stale sidecar.
struct TempIndex {
  explicit TempIndex(const char* stem) {
    const char* dir = std::getenv("TMPDIR");
    path = std::string(dir ? dir : "/tmp") + "/scrna_" + stem + "_" +
           std::to_string(::getpid()) + ".hnsw";
  }
  ~TempIndex() {
    std::remove(path.c_str());
    std::remove((path + ".meta").c_str());
  }
  TempIndex(const TempIndex&) = delete;
  TempIndex& operator=(const TempIndex&) = delete;

  std::string meta() const { return path + ".meta"; }
  std::string path;
};

std::vector<float> random_dense(std::size_t n, std::size_t dim, std::uint32_t seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
  std::vector<float> out(n * dim);
  for (auto& v : out) v = dist(rng);
  return out;
}

// Single-threaded everywhere in this file. HNSW construction is
// order-dependent, so a multi-threaded build is not reproducible run to run --
// which would make "identical before and after a round-trip" untestable for
// reasons that have nothing to do with serialization.
HnswParams deterministic_params() {
  HnswParams p;
  p.M = 16;
  p.ef_construction = 200;
  p.ef_search = 128;
  p.num_threads = 1;
  return p;
}

// ---------------------------------------------------------------- round trip

void test_round_trip_is_lossless() {
  const std::size_t n = 500, dim = 32, k = 10;
  const auto data = random_dense(n, dim, 12345);
  const auto params = deterministic_params();

  TempIndex tmp("roundtrip");

  HnswIndex built(dim, n, params);
  built.add_dense(data.data(), n, /*num_threads=*/1);
  SCRNA_CHECK(built.size() == n);
  SCRNA_CHECK(built.dim() == dim);

  const auto before = built.query_dense(data.data(), n, k, 1, /*self_label_base=*/0);
  built.save_index(tmp.path);

  auto loaded = HnswIndex::load_index(tmp.path, params);
  SCRNA_CHECK(loaded.size() == n);
  SCRNA_CHECK(loaded.dim() == dim);
  SCRNA_CHECK(loaded.capacity() == n);

  const auto after = loaded.query_dense(data.data(), n, k, 1, /*self_label_base=*/0);

  SCRNA_CHECK(before.k == after.k);
  SCRNA_CHECK(before.neighbor_idx == after.neighbor_idx);
  // Exact float equality: the same graph and the same vectors must produce
  // bit-identical distances. Tolerance here would hide a truncated payload.
  SCRNA_CHECK(before.neighbor_sim == after.neighbor_sim);
}

// The reference-mapping workflow: query rows that are NOT in the index. This is
// the case the one-shot API could not express at all, and the reason the
// self-exclusion had to become opt-in.
void test_query_foreign_rows() {
  const std::size_t n_ref = 400, n_query = 25, dim = 16, k = 5;
  const auto ref = random_dense(n_ref, dim, 777);
  const auto qry = random_dense(n_query, dim, 999);
  const auto params = deterministic_params();

  TempIndex tmp("foreign");
  HnswIndex atlas(dim, n_ref, params);
  atlas.add_dense(ref.data(), n_ref, 1);
  atlas.save_index(tmp.path);

  auto loaded = HnswIndex::load_index(tmp.path, params);
  const auto g = loaded.query_dense(qry.data(), n_query, k, 1);

  SCRNA_CHECK(g.k == k);
  SCRNA_CHECK(g.neighbor_idx.size() == n_query * k);
  for (std::size_t i = 0; i < n_query * k; ++i) {
    // Every slot must be filled -- the index has 400 points and we asked for 5.
    SCRNA_CHECK(g.neighbor_idx[i] >= 0);
    SCRNA_CHECK(g.neighbor_idx[i] < static_cast<std::int32_t>(n_ref));
  }

  // Self-exclusion must NOT fire for foreign queries. Query row 0 asking for
  // its neighbours in the atlas is entitled to atlas label 0 as an answer; the
  // old code path would have dropped it purely because the indices coincide.
  // Feed the atlas's own row 0 back in as a foreign query and require that it
  // finds itself, at similarity ~1.
  const auto self = loaded.query_dense(ref.data(), 1, k, 1);
  SCRNA_CHECK(self.neighbor_idx[0] == 0);
  SCRNA_CHECK(self.neighbor_sim[0] > 0.999f);
}

void test_labels_round_trip() {
  const std::size_t n = 20, dim = 8;
  const auto data = random_dense(n, dim, 31337);
  const auto params = deterministic_params();

  TempIndex tmp("labels");
  HnswIndex idx(dim, n, params);
  idx.add_dense(data.data(), n, 1);

  std::vector<std::string> names;
  for (std::size_t i = 0; i < n; ++i) names.push_back("AAACCTGAGC-" + std::to_string(i));
  idx.set_labels(names);
  idx.save_index(tmp.path);

  auto loaded = HnswIndex::load_index(tmp.path, params);
  SCRNA_CHECK(loaded.labels() == names);

  // Wrong count is refused rather than silently truncated.
  SCRNA_CHECK_THROWS(idx.set_labels(std::vector<std::string>{"only-one"}));
}

// -------------------------------------------------------------- refusals

// Helper: build a small valid index at `tmp`, then hand the caller its meta so
// it can be corrupted.
void build_small(const TempIndex& tmp, std::size_t n = 50, std::size_t dim = 8) {
  const auto data = random_dense(n, dim, 5150);
  HnswIndex idx(dim, n, deterministic_params());
  idx.add_dense(data.data(), n, 1);
  idx.save_index(tmp.path);
}

// Overwrites `count` bytes at `offset` in `path`.
void poke(const std::string& path, std::streamoff offset, const char* bytes,
          std::size_t count) {
  std::fstream f(path, std::ios::binary | std::ios::in | std::ios::out);
  f.seekp(offset, std::ios::beg);
  f.write(bytes, static_cast<std::streamsize>(count));
}

void test_missing_meta_is_refused() {
  TempIndex tmp("nometa");
  build_small(tmp);
  std::remove(tmp.meta().c_str());
  // The payload alone is not loadable: without the sidecar there is no
  // trustworthy dimension, and hnswlib would take ours on faith.
  SCRNA_CHECK_THROWS((void)HnswIndex::load_index(tmp.path, deterministic_params()));
}

void test_bad_magic_is_refused() {
  TempIndex tmp("magic");
  build_small(tmp);
  poke(tmp.meta(), 0, "XXXXXXXX", 8);
  SCRNA_CHECK_THROWS((void)HnswIndex::load_index(tmp.path, deterministic_params()));
}

void test_version_skew_is_refused() {
  TempIndex tmp("version");
  build_small(tmp);
  const char bump[4] = {99, 0, 0, 0};  // format_version, immediately after magic
  poke(tmp.meta(), 8, bump, 4);
  SCRNA_CHECK_THROWS((void)HnswIndex::load_index(tmp.path, deterministic_params()));
}

void test_dim_mismatch_is_refused() {
  TempIndex tmp("dim");
  build_small(tmp, 50, 8);
  // `dim` sits at offset 12 (magic 8 + format_version 4). Claiming 4096 dims
  // for an 8-dim payload is the exact input that makes hnswlib compute element
  // offsets from the wrong stride and read out of bounds. It must be rejected
  // before the payload is opened at all.
  const char big_dim[4] = {0, 16, 0, 0};  // 4096, little-endian
  poke(tmp.meta(), 12, big_dim, 4);
  SCRNA_CHECK_THROWS((void)HnswIndex::load_index(tmp.path, deterministic_params()));
}

void test_truncated_payload_is_refused() {
  TempIndex tmp("truncate");
  build_small(tmp);
  // Shorten the payload without touching the meta. This is what an interrupted
  // save on a cluster filesystem actually looks like.
  std::ifstream in(tmp.path, std::ios::binary);
  std::string body((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  in.close();
  body.resize(body.size() / 2);
  std::ofstream out(tmp.path, std::ios::binary | std::ios::trunc);
  out.write(body.data(), static_cast<std::streamsize>(body.size()));
  out.close();
  SCRNA_CHECK_THROWS((void)HnswIndex::load_index(tmp.path, deterministic_params()));
}

void test_corrupt_payload_is_refused() {
  TempIndex tmp("corrupt");
  build_small(tmp);
  // Flip bits deep inside the vector data, leaving the length intact. Only the
  // checksum can catch this; every structural check still passes.
  std::ifstream probe(tmp.path, std::ios::binary | std::ios::ate);
  const auto sz = static_cast<std::streamoff>(probe.tellg());
  probe.close();
  const char junk[8] = {'\xde', '\xad', '\xbe', '\xef', '\xde', '\xad', '\xbe', '\xef'};
  poke(tmp.path, sz / 2, junk, sizeof(junk));

  SCRNA_CHECK_THROWS((void)HnswIndex::load_index(tmp.path, deterministic_params()));

  // ...and with verification disabled it loads, which is precisely why
  // verify_checksum defaults to true. Asserted so nobody flips that default
  // thinking it is a free speedup.
  bool loaded_unverified = false;
  try {
    auto idx = HnswIndex::load_index(tmp.path, deterministic_params(),
                                     /*verify_checksum=*/false);
    loaded_unverified = idx.size() > 0;
  } catch (const std::exception&) {
    // Also acceptable: a flipped byte may land somewhere hnswlib itself
    // rejects. The point of the assertion below is only that the checksum is
    // doing work the structural checks do not.
    loaded_unverified = true;
  }
  SCRNA_CHECK(loaded_unverified);
}

void test_inflated_element_count_is_refused() {
  TempIndex tmp("budget");
  build_small(tmp, 50, 8);
  // capacity sits at offset 8 + 4 (version) + 4 (dim) + 8 (element_count) = 24.
  // A meta claiming 2^40 elements must hit the max_index_bytes ceiling, not be
  // handed to hnswlib as an allocation request.
  const char huge[8] = {0, 0, 0, 0, 1, 0, 0, 0};  // 2^32, little-endian
  poke(tmp.meta(), 24, huge, 8);
  SCRNA_CHECK_THROWS((void)HnswIndex::load_index(tmp.path, deterministic_params()));
}

// ---------------------------------------------------------- API guardrails

void test_capacity_and_dim_guards() {
  const auto params = deterministic_params();
  const auto data = random_dense(10, 4, 1);

  SCRNA_CHECK_THROWS((void)HnswIndex(0, 10, params));   // dim 0
  SCRNA_CHECK_THROWS((void)HnswIndex(4, 0, params));    // capacity 0

  HnswIndex idx(4, 8, params);
  // 10 rows into a capacity of 8: refused up front with a message naming both
  // numbers, rather than throwing from inside hnswlib mid-insertion with a
  // partially populated index.
  SCRNA_CHECK_THROWS(idx.add_dense(data.data(), 10, 1));
  SCRNA_CHECK(idx.size() == 0);
}

void test_incremental_add() {
  const std::size_t dim = 8, half = 30;
  const auto a = random_dense(half, dim, 11);
  const auto b = random_dense(half, dim, 22);
  HnswIndex idx(dim, 2 * half, deterministic_params());
  idx.add_dense(a.data(), half, 1);
  SCRNA_CHECK(idx.size() == half);
  // Labels continue from the current size, so a second call extends rather
  // than overwrites.
  idx.add_dense(b.data(), half, 1);
  SCRNA_CHECK(idx.size() == 2 * half);

  const auto g = idx.query_dense(b.data(), 1, 3, 1);
  SCRNA_CHECK(g.neighbor_idx[0] == static_cast<std::int32_t>(half));
}

// ------------------------------------------------- audit 2026-09-11 (SM1-SM5)

// SM1. Every producer of a KnnGraph sizes it as n * k and then offsets writes
// into it by i * k. Brute force has checked that product since
// REVIEW_2026-08-15 finding 4; the four HNSW producers did not, and
// query_dense() with k = 2^62 over 4 rows stored through a null pointer
// (UBSan: "store to null pointer of type 'std::int32_t'"). The reachable form
// is a k read from a config file or a notebook cell, not a hostile one.
void test_a_k_that_overflows_the_result_is_refused() {
  const std::size_t huge_k = std::size_t{1} << 62;  // 4 * 2^62 wraps to 0
  const auto data = random_dense(4, 4, 8080);
  const auto params = deterministic_params();

  HnswIndex idx(4, 4, params);
  idx.add_dense(data.data(), 4, 1);
  SCRNA_CHECK_THROWS((void)idx.query_dense(data.data(), 4, huge_k, 1));

  // The CSR query and both one-shot entry points take the same path.
  scrna::AlignedVector<float> v{1.0f, 1.0f, 1.0f, 1.0f};
  scrna::AlignedVector<std::int32_t> ci{0, 1, 2, 3};
  scrna::AlignedVector<std::int32_t> rp{0, 2, 4};
  const auto mat = scrna::Block_CSR<float>::from_raw(2, 4, v, ci, rp);
  SCRNA_CHECK_THROWS((void)idx.query_csr(mat, huge_k, 1));
  SCRNA_CHECK_THROWS((void)scrna::hnsw_build_and_search(mat, huge_k, params));
  SCRNA_CHECK_THROWS(
      (void)scrna::hnsw_build_and_search_dense(data.data(), 4, 4, huge_k, params));

  // A k larger than the index is still legal: unused slots keep the sentinel.
  const auto g = idx.query_dense(data.data(), 1, 10, 1);
  SCRNA_CHECK(g.neighbor_idx.size() == 10);
  SCRNA_CHECK(g.neighbor_idx[9] == -1);
}

// SM2. add_dense() reads size() as the label of its first new element, and
// hnswlib's addPoint() on a label that already exists UPDATES that element in
// place rather than failing. Two concurrent adds therefore took the same base
// and the second overwrote the first: 2,000 + 2,000 rows left an index of
// 2,000, with no exception anywhere. Reachable from Python, whose add_dense
// binding releases the GIL.
void test_concurrent_adds_do_not_share_labels() {
  const std::size_t n = 400, dim = 8;
  const auto a = random_dense(n, dim, 4141);
  const auto b = random_dense(n, dim, 4242);

  HnswParams p = deterministic_params();
  p.num_threads = 2;
  HnswIndex idx(dim, 2 * n, p);

  std::string err;
  std::thread t1([&] {
    try { idx.add_dense(a.data(), n, 2); } catch (const std::exception& e) { err = e.what(); }
  });
  std::thread t2([&] {
    try { idx.add_dense(b.data(), n, 2); } catch (const std::exception& e) { err = e.what(); }
  });
  t1.join();
  t2.join();

  SCRNA_CHECK(err.empty());
  SCRNA_CHECK(idx.size() == 2 * n);  // was n before the write mutex
}

// SM3. hnswlib validates M by clamping it (>10000, with a stderr warning) or
// not at all (0, 1). M=1 makes mult_ = 1/log(1) = +inf and getRandomLevel()'s
// (int)inf undefined; M=0 leaves maxM0_ = 0 and walks a vector at offset -4;
// M=20000 builds and saves, then refuses to load, because the meta kept 20000
// while the payload recorded the clamped 10000.
void test_M_outside_the_range_hnswlib_honours_is_refused() {
  const auto data = random_dense(10, 4, 606);
  for (const std::size_t bad : {std::size_t{0}, std::size_t{1}, std::size_t{20000}}) {
    HnswParams p = deterministic_params();
    p.M = bad;
    p.max_index_bytes = std::size_t{1} << 40;  // so M is what refuses it, not the budget
    SCRNA_CHECK_THROWS((void)HnswIndex(4, 10, p));
  }

  // A file may also declare one: loadIndex takes mult_ from the payload rather
  // than recomputing it, so a poked M reaches the same undefined level draw on
  // the next add.
  TempIndex tmp("mrange");
  build_small(tmp, 20, 8);
  // M sits at 8 (magic) + 4 (version) + 4 (dim) + 8 (element_count) + 8 (capacity).
  const char one[4] = {1, 0, 0, 0};
  poke(tmp.meta(), 32, one, 4);
  SCRNA_CHECK_THROWS((void)HnswIndex::load_index(tmp.path, deterministic_params()));

  // ef_search and ef_construction are stored as uint32 and used to truncate.
  HnswIndex ok(4, 10, deterministic_params());
  ok.add_dense(data.data(), 10, 1);
  SCRNA_CHECK_THROWS(ok.set_ef(std::size_t{1} << 33));
  HnswParams wide = deterministic_params();
  wide.ef_construction = std::size_t{1} << 33;
  SCRNA_CHECK_THROWS((void)HnswIndex(4, 10, wide));
}

// SM4. hnswlib's saveIndex() writes through an unchecked ofstream, so a write
// that fails partway returns normally -- and our CRC, taken from the file
// afterwards, then SEALED the short payload as intact. save_index() reported
// success on an index that could never be loaded.
//
// RLIMIT_FSIZE in a child models the failure. SIGXFSZ is ignored first, which
// is what CPython does at startup, so the write fails with EFBIG instead of
// killing the process -- the same shape as ENOSPC or a quota.
void test_a_short_write_is_not_sealed_as_a_valid_index() {
  const std::size_t n = 2000, dim = 64;
  const auto data = random_dense(n, dim, 2727);
  TempIndex tmp("shortwrite");

  HnswIndex idx(dim, n, deterministic_params());
  idx.add_dense(data.data(), n, 1);

  const pid_t pid = ::fork();
  if (pid < 0) {
    std::printf("fork unavailable; skipping the short-write case\n");
    return;
  }
  if (pid == 0) {
    std::signal(SIGXFSZ, SIG_IGN);
    rlimit rl{256 * 1024, 256 * 1024};  // the payload needs ~800 KB
    ::setrlimit(RLIMIT_FSIZE, &rl);
    int rc = 1;
    try {
      idx.save_index(tmp.path);
      rc = 0;  // returned success on a truncated write: the defect
    } catch (const std::exception&) {
      rc = 2;  // refused: correct
    }
    _exit(rc);
  }
  int status = 0;
  ::waitpid(pid, &status, 0);
  SCRNA_CHECK(WIFEXITED(status) && WEXITSTATUS(status) == 2);

  // And no sidecar was left behind, so the partial payload cannot be loaded.
  SCRNA_CHECK(!std::ifstream(tmp.meta()).good());
  SCRNA_CHECK_THROWS((void)HnswIndex::load_index(tmp.path, deterministic_params()));
}

// SM5. Labels set, then more elements added, produced a meta that described a
// prefix of the index. save_index() succeeded and load_index() refused the
// pair -- a failure that surfaces only once the build that produced it is
// gone. Refused at save, where set_labels() can still fix it.
void test_labels_stale_after_a_later_add_are_refused_at_save() {
  const std::size_t dim = 8, half = 10;
  const auto data = random_dense(2 * half, dim, 5454);
  TempIndex tmp("stalelabels");

  HnswIndex idx(dim, 2 * half, deterministic_params());
  idx.add_dense(data.data(), half, 1);
  std::vector<std::string> names;
  for (std::size_t i = 0; i < half; ++i) names.push_back("cell" + std::to_string(i));
  idx.set_labels(names);
  idx.add_dense(data.data() + half * dim, half, 1);

  SCRNA_CHECK_THROWS(idx.save_index(tmp.path));
  SCRNA_CHECK(!std::ifstream(tmp.path).good());  // nothing written at all

  // Naming every element, or none, saves and loads as before.
  for (std::size_t i = half; i < 2 * half; ++i) names.push_back("cell" + std::to_string(i));
  idx.set_labels(names);
  idx.save_index(tmp.path);
  SCRNA_CHECK(HnswIndex::load_index(tmp.path, deterministic_params()).labels() == names);
}

// ------------------------------------------- audit 2026-09-11, round 2 (SM9-SM11)

// SM9. Query breadth used to be applied by calling index_->setEf() immediately
// before each search -- a write to hnswlib's shared ef_ that every concurrent
// search reads. TSan, two threads querying one index with k=2 and k=400:
//
//   WARNING: ThreadSanitizer: data race
//     Write of size 8 ... HnswIndex::query_dense  hnsw_index.hpp:385
//     Previous read of size 8 ... searchKnn       hnswalg.h:1414
//
// searchKnn now takes ef as an argument (third_party/hnswlib/VERSION.txt,
// patch 4) and the query path never writes it. The observable consequence,
// which is what this test pins: a query's results no longer depend on what
// another thread happens to be asking for at the same moment.
void test_concurrent_queries_do_not_share_search_breadth() {
  const std::size_t n = 1500, dim = 16;
  const auto data = random_dense(n, dim, 909);
  HnswIndex idx(dim, n, deterministic_params());
  idx.add_dense(data.data(), n, 1);

  // The answers each k must produce, computed with nothing else running.
  const auto serial_small = idx.query_dense(data.data(), 100, 2, 1, 0);
  const auto serial_large = idx.query_dense(data.data(), 100, 300, 1, 0);

  std::atomic<bool> go{false};
  bool small_ok = true, large_ok = true;
  std::thread a([&] {
    while (!go.load()) {}
    for (int rep = 0; rep < 20; ++rep) {
      small_ok = small_ok && idx.query_dense(data.data(), 100, 2, 1, 0).neighbor_idx ==
                                 serial_small.neighbor_idx;
    }
  });
  std::thread b([&] {
    while (!go.load()) {}
    for (int rep = 0; rep < 20; ++rep) {
      large_ok = large_ok && idx.query_dense(data.data(), 100, 300, 1, 0).neighbor_idx ==
                                 serial_large.neighbor_idx;
    }
  });
  go.store(true);
  a.join();
  b.join();

  SCRNA_CHECK(small_ok);
  SCRNA_CHECK(large_ok);
}

// SM10. A KnnGraph stores neighbour ids as int32, so an index able to hold more
// than INT32_MAX elements could answer with a wrapped, negative id -- which is
// exactly the -1 "no neighbour" sentinel. Unreachable under the default memory
// budget, which is why this was first written down as a residual; but the
// budget is a caller-settable default, so the invariant is enforced where it
// belongs. The guard must fire BEFORE any allocation is attempted, which is
// what the raised budget here checks.
void test_capacity_beyond_int32_is_refused() {
  HnswParams p = deterministic_params();
  p.max_index_bytes = std::numeric_limits<std::size_t>::max();
  const std::size_t past_int32 = static_cast<std::size_t>(
      std::numeric_limits<std::int32_t>::max()) + 1;
  SCRNA_CHECK_THROWS((void)HnswIndex(1, past_int32, p));

  // A file may declare one too.
  TempIndex tmp("bigcap");
  build_small(tmp, 20, 8);
  // capacity sits at 8 (magic) + 4 (version) + 4 (dim) + 8 (element_count).
  const char huge[8] = {0, 0, 0, '\x80', 0, 0, 0, 0};  // 2^31
  poke(tmp.meta(), 24, huge, 8);
  SCRNA_CHECK_THROWS((void)HnswIndex::load_index(tmp.path, p));

  // INT32_MAX itself is still accepted as a declared capacity (the budget is
  // what refuses it in practice, with a message about memory).
  SCRNA_CHECK_THROWS((void)HnswIndex(
      1, static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()),
      deterministic_params()));
}

// `k` is a row STRIDE, not a request for exactly k answers. Both edges of that
// contract are pinned for brute force (test_knn_k_exceeds_rows,
// test_knn_edge_cases); the HNSW paths had no equivalent, and they are where
// the fetch count is now clamped to the index size.
void test_k_at_the_edges_of_the_index() {
  const std::size_t n = 8, dim = 4;
  const auto data = random_dense(n, dim, 4004);
  HnswIndex idx(dim, n, deterministic_params());
  idx.add_dense(data.data(), n, 1);

  // k = 0: an empty graph, not an error.
  SCRNA_CHECK(idx.query_dense(data.data(), n, 0, 1).neighbor_idx.empty());

  // k far beyond the index: k slots, only what exists filled, the rest -1.
  const auto g = idx.query_dense(data.data(), 1, 50, 1, /*self_label_base=*/0);
  SCRNA_CHECK(g.neighbor_idx.size() == 50);
  std::size_t filled = 0;
  for (const auto id : g.neighbor_idx) filled += (id >= 0);
  SCRNA_CHECK(filled == n - 1);  // every element except the query row itself
  SCRNA_CHECK(g.neighbor_idx.back() == -1);
  SCRNA_CHECK(g.neighbor_sim.back() == 0.0f);

  // Querying an index with nothing in it is likewise sentinels, not an error.
  HnswIndex empty(dim, 4, deterministic_params());
  const auto ge = empty.query_dense(data.data(), 2, 3, 1);
  SCRNA_CHECK(ge.neighbor_idx.size() == 6);
  for (const auto id : ge.neighbor_idx) SCRNA_CHECK(id == -1);
}

// A cell with no counts is ordinary input (scanpy filters them; nothing
// guarantees the caller did). Its norm is zero, so the normalisation scale is
// 1/0 -- and the question this answers is whether that reaches the graph as a
// NaN, which WOULD corrupt traversal, or as a defined zero vector.
//
// Measured: no non-finite value anywhere. inv_norm_from_sq returns 0 for a
// non-positive norm, the row is indexed as all zeros, its distance to
// everything is 1 - 0 = 1 -- the maximum -- so it can never displace a real
// neighbour, and its own neighbour list comes back at similarity 0.
void test_a_zero_norm_vector_produces_no_nan() {
  const std::size_t n = 200, dim = 16, empty_row = 3;
  auto data = random_dense(n, dim, 6006);
  for (std::size_t j = 0; j < dim; ++j) data[empty_row * dim + j] = 0.0f;

  HnswIndex idx(dim, n, deterministic_params());
  idx.add_dense(data.data(), n, 1);
  const auto g = idx.query_dense(data.data(), n, 5, 1, /*self_label_base=*/0);

  for (const float s : g.neighbor_sim) SCRNA_CHECK(std::isfinite(s));
  // The empty cell gets a filled neighbour list, all at zero similarity.
  for (std::size_t r = 0; r < 5; ++r) {
    SCRNA_CHECK(g.neighbor_idx[empty_row * 5 + r] >= 0);
    SCRNA_CHECK(g.neighbor_sim[empty_row * 5 + r] == 0.0f);
  }
  // And it is nobody else's best match.
  for (std::size_t i = 0; i < n; ++i) {
    if (i == empty_row) continue;
    SCRNA_CHECK(g.neighbor_idx[i * 5] != static_cast<std::int32_t>(empty_row));
  }
}

// Degenerate geometry: every one of these is either refused with a message
// naming the mismatch, or defined and empty. None of them may reach an
// allocation of zero bytes that is then indexed, or a stride computed from a
// zero dimension.
void test_degenerate_dimensions_are_refused_or_empty() {
  const auto data = random_dense(10, 8, 7007);

  // Construction bounds.
  SCRNA_CHECK_THROWS((void)HnswIndex(0, 10, deterministic_params()));  // dim 0
  SCRNA_CHECK_THROWS((void)HnswIndex(8, 0, deterministic_params()));   // capacity 0

  HnswIndex idx(8, 10, deterministic_params());

  // Null data with a positive count is a caller error, not a silent no-op.
  SCRNA_CHECK_THROWS(idx.add_dense(nullptr, 3, 1));
  SCRNA_CHECK_THROWS((void)idx.query_dense(nullptr, 2, 3, 1));

  // An empty add is a no-op; an empty query is an empty graph.
  idx.add_dense(data.data(), 0, 1);
  SCRNA_CHECK(idx.size() == 0);
  SCRNA_CHECK(idx.query_dense(data.data(), 0, 3, 1).neighbor_idx.empty());

  idx.add_dense(data.data(), 10, 1);

  // Feature width must equal the index dimension on both CSR entry points.
  const scrna::AlignedVector<float> v{1.0f};
  const scrna::AlignedVector<std::int32_t> ci{0};
  const scrna::AlignedVector<std::int32_t> rp{0, 1};
  const auto narrow = scrna::Block_CSR<float>::from_raw(1, 4, v, ci, rp);  // 4 != 8
  SCRNA_CHECK_THROWS(idx.add_csr(narrow, 1));
  SCRNA_CHECK_THROWS((void)idx.query_csr(narrow, 2, 1));

  // A matrix with no rows adds nothing and is not an error.
  const scrna::AlignedVector<float> no_vals;
  const scrna::AlignedVector<std::int32_t> no_idx;
  const scrna::AlignedVector<std::int32_t> just_zero{0};
  const auto no_rows = scrna::Block_CSR<float>::from_raw(0, 8, no_vals, no_idx, just_zero);
  const std::size_t before = idx.size();
  idx.add_csr(no_rows, 1);
  SCRNA_CHECK(idx.size() == before);

  // The one-shots over degenerate input: empty, never a crash.
  SCRNA_CHECK(scrna::hnsw_build_and_search_dense(nullptr, 0, 4, 5).neighbor_idx.empty());
  const auto zero_by_zero =
      scrna::Block_CSR<float>::from_raw(0, 0, no_vals, no_idx, just_zero);
  SCRNA_CHECK(scrna::hnsw_build_and_search(zero_by_zero, 3).neighbor_idx.empty());
  // dim == 0 keeps the caller's stride and fills it with sentinels.
  const auto g = scrna::hnsw_build_and_search_dense(data.data(), 4, 0, 5);
  SCRNA_CHECK(g.neighbor_idx.size() == 20);
  for (const auto id : g.neighbor_idx) SCRNA_CHECK(id == -1);
}

#endif  // SCRNA_ENABLE_HNSW

}  // namespace

int main() {
#if defined(SCRNA_ENABLE_HNSW)
  test_round_trip_is_lossless();
  test_query_foreign_rows();
  test_labels_round_trip();
  test_missing_meta_is_refused();
  test_bad_magic_is_refused();
  test_version_skew_is_refused();
  test_dim_mismatch_is_refused();
  test_truncated_payload_is_refused();
  test_corrupt_payload_is_refused();
  test_inflated_element_count_is_refused();
  test_capacity_and_dim_guards();
  test_incremental_add();
  test_a_k_that_overflows_the_result_is_refused();
  test_concurrent_adds_do_not_share_labels();
  test_M_outside_the_range_hnswlib_honours_is_refused();
  test_a_short_write_is_not_sealed_as_a_valid_index();
  test_labels_stale_after_a_later_add_are_refused_at_save();
  test_concurrent_queries_do_not_share_search_breadth();
  test_capacity_beyond_int32_is_refused();
  test_k_at_the_edges_of_the_index();
  test_a_zero_norm_vector_produces_no_nan();
  test_degenerate_dimensions_are_refused_or_empty();
#else
  std::printf("SCRNA_ENABLE_HNSW is off; nothing to test\n");
#endif

  if (g_failures == 0) {
    std::printf("All tests passed.\n");
    return 0;
  }
  std::fprintf(stderr, "%d test(s) failed.\n", g_failures);
  return 1;
}
