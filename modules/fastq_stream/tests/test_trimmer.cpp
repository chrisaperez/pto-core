// SPDX-License-Identifier: MIT
// fastq_stream — unit tests.
//
// Deliberately dependency-free: no GoogleTest, no Catch2. The project's whole
// premise is that it drops into a cluster with nothing but a compiler, and the
// test suite honours the same constraint.
#include <libdeflate.h>
#include <zlib.h>

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <functional>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <limits>
#include <vector>

#include "fastq_stream/adapter_trimmer.hpp"
#include "fastq_stream/lockfree_queue.hpp"
#include "fastq_stream/phred_calculator.hpp"
#include "fastq_stream/pipeline.hpp"
#include "fastq_stream/record.hpp"
#include "fastq_stream/report.hpp"
#include "fastq_stream/simd.hpp"

namespace {

int g_failures = 0;
int g_checks = 0;

#define CHECK(cond)                                                            \
  do {                                                                         \
    ++g_checks;                                                                \
    if (!(cond)) {                                                             \
      std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << "  " << #cond      \
                << "\n";                                                       \
      ++g_failures;                                                            \
    }                                                                          \
  } while (0)

#define CHECK_EQ(a, b)                                                         \
  do {                                                                         \
    ++g_checks;                                                                \
    const auto va_ = (a);                                                      \
    const auto vb_ = (b);                                                      \
    if (!(va_ == vb_)) {                                                       \
      std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << "  " << #a         \
                << " (" << va_ << ") != " << #b << " (" << vb_ << ")\n";        \
      ++g_failures;                                                            \
    }                                                                          \
  } while (0)

void section(const char* name) { std::cerr << "-- " << name << "\n"; }

// --------------------------------------------------------------- SIMD kernels
// Scalar references. The SIMD paths must agree exactly, at every length,
// including the sub-vector tails where hand-written kernels usually break.

uint64_t ref_sum(const uint8_t* p, std::size_t n) {
  uint64_t s = 0;
  for (std::size_t i = 0; i < n; ++i) s += p[i];
  return s;
}
std::size_t ref_count_ge(const uint8_t* p, std::size_t n, uint8_t t) {
  std::size_t c = 0;
  for (std::size_t i = 0; i < n; ++i) c += (p[i] >= t);
  return c;
}
std::size_t ref_mismatch(const uint8_t* a, const uint8_t* b, std::size_t n) {
  std::size_t d = 0;
  for (std::size_t i = 0; i < n; ++i) d += (a[i] != b[i]);
  return d;
}
// Independent of scalar::find_seed4 (which every ISA including scalar itself
// is checked against below) -- this is the "restate the loop adapter_trimmer
// used before find_seed4 existed" oracle, byte-by-byte with memcmp rather
// than a word load, so a bug shared between the production scalar reference
// and a vector kernel would still show up as a disagreement here.
std::size_t ref_find_seed4(const uint8_t* p, std::size_t from, std::size_t to_inclusive,
                          uint32_t seed) {
  for (std::size_t s = from; s <= to_inclusive; ++s) {
    if (std::memcmp(p + s, &seed, 4) == 0) return s;
  }
  return fq::simd::kNotFound;
}

// Every ISA this CPU can execute, through the explicit-ISA entry points, plus
// whatever the dispatcher selected. The second matters because FQ_FORCE_ISA can
// select a path CPUID does not report -- Rosetta 2 executes AVX2 and advertises
// none of it -- and CI's per-ISA sweep reaches kernels only that way.
std::vector<fq::simd::Isa> isas_under_test() {
  using fq::simd::Isa;
  std::vector<Isa> out;
  for (const Isa isa : {Isa::kScalar, Isa::kNeon, Isa::kAvx2, Isa::kAvx512}) {
    if (fq::simd::isa_supported(isa)) out.push_back(isa);
  }
  if (std::find(out.begin(), out.end(), fq::simd::isa_level()) == out.end()) {
    out.push_back(fq::simd::isa_level());
  }
  return out;
}

void test_simd() {
  section("simd kernels vs scalar reference, every executable ISA");
  std::mt19937 rng(12345);
  std::uniform_int_distribution<int> byte(0, 255);

  // Offsets 0..7 into the buffers so no kernel only ever sees an aligned start.
  std::vector<uint8_t> a_buf(320), b_buf(320);
  for (auto& v : a_buf) v = static_cast<uint8_t>(byte(rng));
  for (auto& v : b_buf) v = static_cast<uint8_t>(byte(rng));
  // Force a run of equal bytes so mismatch counting is exercised both ways.
  std::memcpy(b_buf.data() + 40, a_buf.data() + 40, 90);

  for (const fq::simd::Isa isa : isas_under_test()) {
    std::cerr << "   isa " << fq::simd::isa_name(isa) << "\n";
    for (std::size_t off = 0; off < 8; ++off) {
      const uint8_t* a = a_buf.data() + off;
      const uint8_t* b = b_buf.data() + off;
      for (std::size_t n = 0; n + off <= 308; ++n) {
        CHECK_EQ(fq::simd::sum_u8_isa(isa, a, n), ref_sum(a, n));
        for (uint8_t thr : {uint8_t(0), uint8_t(53), uint8_t(63), uint8_t(255)}) {
          CHECK_EQ(fq::simd::count_ge_isa(isa, a, n, thr), ref_count_ge(a, n, thr));
        }
        // cap = n means "never saturate", so the result must be exact.
        const std::size_t exact = ref_mismatch(a, b, n);
        CHECK_EQ(fq::simd::mismatches_isa(isa, a, b, n, n), exact);
        // Below the exact count, only "> cap" is promised.
        if (exact > 0) CHECK(fq::simd::mismatches_isa(isa, a, b, n, exact - 1) > exact - 1);

        std::vector<uint32_t> dst(n + 20, 7);
        fq::simd::add_widen_u8_to_u32_isa(isa, a, n, dst.data());
        bool widen_ok = true;
        for (std::size_t i = 0; i < n; ++i) {
          if (dst[i] != 7u + a[i]) widen_ok = false;
        }
        for (std::size_t i = n; i < dst.size(); ++i) {
          if (dst[i] != 7u) widen_ok = false;  // must not write past n
        }
        CHECK(widen_ok);
      }
    }
    // The early-exit path must saturate rather than lie about the count.
    CHECK(fq::simd::mismatches_isa(isa, a_buf.data(), b_buf.data(), 300, 2) > 2);
  }

  // The dispatching entry points are the ones production calls.
  CHECK_EQ(fq::simd::sum_u8(a_buf.data(), 300), ref_sum(a_buf.data(), 300));
  CHECK_EQ(fq::simd::count_ge(a_buf.data(), 300, 53), ref_count_ge(a_buf.data(), 300, 53));
  CHECK(fq::simd::isa_supported(fq::simd::Isa::kScalar));
  CHECK(fq::simd::isa_compiled(fq::simd::isa_level()));
}

// find_seed4's block widths (16/32/64 for NEON/AVX2/AVX-512) mean the
// interesting failures live at block boundaries, not in the middle of one --
// so unlike test_simd's length sweep above, this plants a match at every
// single offset across several full blocks of the widest kernel, rather than
// sampling lengths.
void test_find_seed4() {
  section("find_seed4 vs scalar reference, every executable ISA, every offset");
  using fq::simd::kNotFound;

  const uint32_t seed = 0x44434241;  // bytes 'A','B','C','D' in memory order
  uint8_t seed_bytes[4];
  std::memcpy(seed_bytes, &seed, 4);

  // Background bytes outside 'A'..'D' so they can never accidentally form the
  // seed word themselves -- a "no match" case that happened to match would
  // hide a real bug rather than exercise the no-match path.
  std::mt19937 rng(20260916);
  std::uniform_int_distribution<int> bg(int('E'), int('Z'));

  constexpr std::size_t kRangeLen = 200;  // > 3 full 64-wide (AVX-512) blocks
  std::vector<uint8_t> base(kRangeLen + 3);
  for (auto& v : base) v = static_cast<uint8_t>(bg(rng));
  const std::size_t to = kRangeLen - 4;  // last offset with 4 bytes in range

  for (const fq::simd::Isa isa : isas_under_test()) {
    std::cerr << "   isa " << fq::simd::isa_name(isa) << "\n";

    // 1. No match anywhere, at every `from`.
    for (std::size_t from = 0; from <= to; from += 5) {
      CHECK_EQ(fq::simd::find_seed4_isa(isa, base.data(), from, to, seed), kNotFound);
    }

    // 2. A single planted match, at every offset -- this is the one that
    //    walks every block-boundary residue for every compiled ISA.
    for (std::size_t pos = 0; pos <= to; ++pos) {
      std::vector<uint8_t> planted = base;
      std::memcpy(planted.data() + pos, seed_bytes, 4);
      const std::size_t got = fq::simd::find_seed4_isa(isa, planted.data(), 0, to, seed);
      const std::size_t want = ref_find_seed4(planted.data(), 0, to, seed);
      CHECK_EQ(got, want);
      CHECK_EQ(got, pos);
    }

    // 3. Two matches: the 5'-most must win, including pairs that straddle a
    //    block boundary from opposite sides (a kernel that reduced blocks
    //    out of order, or picked the wrong lane within one, would return the
    //    later position instead).
    for (std::size_t first : {std::size_t(0), std::size_t(1), std::size_t(15), std::size_t(16),
                              std::size_t(17), std::size_t(31), std::size_t(32),
                              std::size_t(63), std::size_t(64), std::size_t(100)}) {
      // delta >= 4: a smaller gap would make the second 4-byte plant
      // overwrite part of the first one, corrupting both.
      for (std::size_t delta : {std::size_t(4), std::size_t(20), std::size_t(63),
                                std::size_t(100)}) {
        const std::size_t second = first + delta;
        if (second > to) continue;
        std::vector<uint8_t> planted = base;
        std::memcpy(planted.data() + first, seed_bytes, 4);
        std::memcpy(planted.data() + second, seed_bytes, 4);
        CHECK_EQ(fq::simd::find_seed4_isa(isa, planted.data(), 0, to, seed), first);
      }
    }

    // 4. `from` excludes an earlier match; a later one is still found.
    {
      std::vector<uint8_t> planted = base;
      std::memcpy(planted.data() + 10, seed_bytes, 4);
      std::memcpy(planted.data() + 90, seed_bytes, 4);
      CHECK_EQ(fq::simd::find_seed4_isa(isa, planted.data(), 50, to, seed), std::size_t(90));
      CHECK_EQ(fq::simd::find_seed4_isa(isa, planted.data(), 50, 80, seed), kNotFound);
    }

    // 5. Degenerate and single-offset ranges.
    CHECK_EQ(fq::simd::find_seed4_isa(isa, base.data(), 5, 4, seed), kNotFound);
    {
      std::vector<uint8_t> one = base;
      std::memcpy(one.data() + 7, seed_bytes, 4);
      CHECK_EQ(fq::simd::find_seed4_isa(isa, one.data(), 7, 7, seed), std::size_t(7));
      CHECK_EQ(fq::simd::find_seed4_isa(isa, one.data(), 8, 8, seed), kNotFound);
    }
  }

  // The dispatching entry point production calls.
  {
    std::vector<uint8_t> one = base;
    std::memcpy(one.data() + 50, seed_bytes, 4);
    CHECK_EQ(fq::simd::find_seed4(one.data(), 0, to, seed), std::size_t(50));
  }
}

// ----------------------------------------------------------------- SPSC ring

void test_ring_basic() {
  section("spsc ring, single-threaded semantics");
  fq::SpscRing<int, 8> ring;
  CHECK(ring.empty());
  CHECK(!ring.try_pop().has_value());

  for (int i = 0; i < 7; ++i) CHECK(ring.try_push(int(i)));
  CHECK(!ring.try_push(99));  // one slot is reserved to disambiguate full/empty
  CHECK_EQ(ring.size_approx(), std::size_t(7));

  for (int i = 0; i < 7; ++i) {
    auto v = ring.try_pop();
    CHECK(v.has_value());
    CHECK_EQ(*v, i);
  }
  CHECK(ring.empty());

  // Wrap around several times.
  for (int round = 0; round < 100; ++round) {
    CHECK(ring.try_push(int(round)));
    auto v = ring.try_pop();
    CHECK(v.has_value() && *v == round);
  }
}

void test_ring_threaded() {
  section("spsc ring, producer/consumer handoff");
  constexpr int kN = 2'000'000;
  fq::SpscRing<int, 1024> ring;
  std::atomic<bool> done{false};
  uint64_t sum = 0;

  std::thread consumer([&] {
    int expected = 0;
    while (expected < kN) {
      if (auto v = ring.try_pop()) {
        if (*v != expected) {  // order must be exact
          ++g_failures;
          std::cerr << "FAIL ring order: got " << *v << " expected " << expected << "\n";
          return;
        }
        sum += static_cast<uint64_t>(*v);
        ++expected;
      }
    }
    done.store(true);
  });

  for (int i = 0; i < kN; ++i) {
    while (!ring.try_push(int(i))) { /* spin */ }
  }
  consumer.join();

  ++g_checks;
  const uint64_t expect = uint64_t(kN) * (kN - 1) / 2;
  if (sum != expect) {
    ++g_failures;
    std::cerr << "FAIL ring sum " << sum << " != " << expect << "\n";
  }
  CHECK(done.load());
}

// ------------------------------------------------------------------ trimming

void test_quality_trimming() {
  section("quality trimming");
  // Phred 40 = 'I', Phred 2 = '#'.
  const std::string q_all_good(50, 'I');
  const auto* g = reinterpret_cast<const uint8_t*>(q_all_good.data());
  CHECK_EQ(fq::trim_sliding_window(g, 50, 4, 20.0), std::size_t(50));
  CHECK_EQ(fq::trim_back(g, 50, 20), std::size_t(50));
  CHECK_EQ(fq::trim_front(g, 50, 20), std::size_t(0));

  // Good for 30 bases (Q40 'I') then a Q2 '#' tail. Trimmomatic SLIDINGWINDOW
  // semantics cut at the *start* of the first failing window, which sacrifices
  // up to window-1 good bases. Here the window at offset 28 covers
  // (73,73,35,35) = 216 >= 4*33 + 80 and survives; offset 29 covers
  // (73,35,35,35) = 178 and fails, so the read is cut to 29 — not 30.
  std::string q = std::string(30, 'I') + std::string(20, '#');
  const auto* p = reinterpret_cast<const uint8_t*>(q.data());
  const std::size_t kept = fq::trim_sliding_window(p, q.size(), 4, 20.0);
  CHECK_EQ(kept, std::size_t(29));

  // A read that fails from the very first window is fully trimmed.
  const std::string bad(50, '#');
  CHECK_EQ(fq::trim_sliding_window(reinterpret_cast<const uint8_t*>(bad.data()), 50, 4,
                                   20.0),
           std::size_t(0));

  // Window wider than the read: leave it alone.
  CHECK_EQ(fq::trim_sliding_window(p, 3, 4, 20.0), std::size_t(3));
  // Window disabled.
  CHECK_EQ(fq::trim_sliding_window(p, q.size(), 0, 20.0), q.size());

  // Leading/trailing hard trims.
  std::string mixed = "##" + std::string(20, 'I') + "###";
  const auto* m = reinterpret_cast<const uint8_t*>(mixed.data());
  CHECK_EQ(fq::trim_front(m, mixed.size(), 20), std::size_t(2));
  CHECK_EQ(fq::trim_back(m, mixed.size(), 20), std::size_t(22));

  CHECK(fq::mean_phred(g, 50) > 39.9);
  CHECK_EQ(fq::count_n(reinterpret_cast<const uint8_t*>("ACGTNNAC"), 8), std::size_t(2));
}

// ------------------------------------------------------------------ adapters

void test_adapter_trimming() {
  section("adapter trimming");
  fq::AdapterTrimmer t;
  t.add(fq::kTruSeqAdapter);  // AGATCGGAAGAGC

  const std::string insert = "TTGCACAGTTGCACAGTTGCACAGTTGCACAG";  // 32 bp, no seed

  // Full adapter at a known offset.
  {
    const std::string read = insert + std::string(fq::kTruSeqAdapter) + "TTTTTT";
    CHECK_EQ(t.trim(reinterpret_cast<const uint8_t*>(read.data()), read.size()),
             insert.size());
  }

  // Adapter runs off the 3' end (partial overlap) — still trimmed.
  {
    const std::string read = insert + "AGATCGG";
    CHECK_EQ(t.trim(reinterpret_cast<const uint8_t*>(read.data()), read.size()),
             insert.size());
  }

  // Overlap shorter than --adapter-overlap is left alone.
  {
    t.set_min_overlap(6);
    const std::string read = insert + "AGAT";
    CHECK_EQ(t.trim(reinterpret_cast<const uint8_t*>(read.data()), read.size()),
             read.size());
    t.set_min_overlap(4);
  }

  // One mismatch outside the seed is tolerated at the default 20% rate.
  {
    std::string ad(fq::kTruSeqAdapter);
    ad[8] = 'T';  // AGATCGGA[A->T]GAGC
    const std::string read = insert + ad + "TTTT";
    CHECK_EQ(t.trim(reinterpret_cast<const uint8_t*>(read.data()), read.size()),
             insert.size());
  }

  // No adapter: a random-ish read must survive untouched.
  {
    const std::string read =
        "TTGCACAGTTGCACAGTTGCACAGTTGCACAGTTGCACAGTTGCACAGTTGCACAG";
    CHECK_EQ(t.trim(reinterpret_cast<const uint8_t*>(read.data()), read.size()),
             read.size());
  }

  // Empty trimmer is a no-op.
  {
    fq::AdapterTrimmer none;
    CHECK(none.empty());
    CHECK_EQ(none.trim(reinterpret_cast<const uint8_t*>(insert.data()), insert.size()),
             insert.size());
  }
}

// Restates AdapterTrimmer::trim's pre-find_seed4 algorithm exactly -- a
// byte-at-a-time seed test rather than a bulk vector search -- independent of
// AdapterTrimmer's own implementation, as the oracle for the fuzz below. Any
// divergence between this and the shipped trim() is a bug in the find_seed4
// integration (the loop restructuring in adapter_trimmer.hpp), since
// find_seed4 itself is already checked directly, exhaustively, against
// scalar::find_seed4 in test_find_seed4 above.
std::size_t ref_trim(const std::vector<std::string>& adapters, std::size_t min_overlap,
                     double max_mismatch_rate, std::size_t seed_len, const uint8_t* seq,
                     std::size_t len) {
  if (adapters.empty() || len == 0 || len < min_overlap) return len;
  std::size_t best = len;
  for (const std::string& ad : adapters) {
    const auto* a = reinterpret_cast<const uint8_t*>(ad.data());
    const std::size_t alen = ad.size();
    const std::size_t seed = std::min({seed_len, min_overlap, alen});
    const std::size_t limit = len - min_overlap;
    for (std::size_t s = 0; s <= limit && s < best; ++s) {
      if (std::memcmp(seq + s, a, seed) != 0) continue;
      const std::size_t overlap = std::min(alen, len - s);
      if (overlap < min_overlap) continue;
      const std::size_t cap =
          static_cast<std::size_t>(max_mismatch_rate * static_cast<double>(overlap));
      if (ref_mismatch(seq + s, a, overlap) <= cap) {
        best = s;
        break;
      }
    }
  }
  return best;
}

void test_adapter_trim_matches_reference_after_find_seed4() {
  section("adapter trim vs an independent reference, 20000 random reads");
  std::mt19937 rng(20260916);
  std::uniform_int_distribution<int> base_pick(0, 3);
  const char bases[4] = {'A', 'C', 'G', 'T'};

  const std::vector<std::vector<std::string>> adapter_sets = {
      {std::string(fq::kTruSeqAdapter)},
      {std::string(fq::kTruSeqAdapter), std::string(fq::kNexteraAdapter)},
      {"AAAA"},
      {"ACGT"},
      {"GATTACA"},
  };

  for (int trial = 0; trial < 20000; ++trial) {
    const std::size_t len = 1 + (rng() % 60);
    std::vector<uint8_t> seq(len);
    for (auto& b : seq) b = static_cast<uint8_t>(bases[base_pick(rng)]);

    const auto& adapters = adapter_sets[rng() % adapter_sets.size()];
    // Pure random DNA almost never hits a 4-19 byte seed by chance, which
    // would make the fuzz mostly exercise the empty-result path -- splice in
    // a (possibly mutated) adapter copy about half the time so the seed
    // actually fires and the extend-and-check logic, not just "no match", is
    // under test at every trial.
    if (len >= 4 && (rng() % 2 == 0)) {
      const std::string& ad = adapters[rng() % adapters.size()];
      const std::size_t pos = rng() % len;
      const std::size_t n = std::min(ad.size(), len - pos);
      for (std::size_t i = 0; i < n; ++i) {
        uint8_t c = static_cast<uint8_t>(ad[i]);
        if (rng() % 6 == 0) c = static_cast<uint8_t>(bases[base_pick(rng)]);  // planted mismatch
        seq[pos + i] = c;
      }
    }

    fq::AdapterTrimmer at;
    for (const auto& ad : adapters) at.add(ad);
    at.set_min_overlap(1 + (rng() % 8));
    at.set_max_mismatch_rate((rng() % 6) * 0.1);
    at.set_seed_len(1 + (rng() % 6));

    const std::size_t got = at.trim(seq.data(), len);
    const std::size_t want = ref_trim(adapters, at.min_overlap(), at.max_mismatch_rate(),
                                      at.seed_len(), seq.data(), len);
    CHECK_EQ(got, want);
  }
}

// ------------------------------------------------------------------- records

void test_record_parsing() {
  section("record framing");
  const std::string one = "@r1\nACGT\n+\nIIII\n";
  const std::string two = one + "@r2\nAC\n+\n@@\n";  // '@' as a quality char

  CHECK_EQ(fq::last_record_boundary(one.data(), one.size()), one.size());
  CHECK_EQ(fq::last_record_boundary(two.data(), two.size()), two.size());
  // A truncated second record must not be reported as complete.
  const std::string partial = one + "@r2\nAC\n+\n";
  CHECK_EQ(fq::last_record_boundary(partial.data(), partial.size()), one.size());
  CHECK_EQ(fq::last_record_boundary("@r1\nACGT\n", 9), std::size_t(0));

  fq::RecordIterator it(reinterpret_cast<const uint8_t*>(two.data()), two.size());
  fq::RecordView r;
  CHECK(it.next(r));
  CHECK(r.well_formed());
  CHECK_EQ(r.seq_len, std::size_t(4));
  CHECK_EQ(std::string(reinterpret_cast<const char*>(r.qual), r.qual_len), "IIII");
  CHECK(it.next(r));
  CHECK(r.well_formed());
  // The quality line here is "@@" — proof that '@' cannot delimit records.
  CHECK_EQ(std::string(reinterpret_cast<const char*>(r.qual), r.qual_len), "@@");
  CHECK(!it.next(r));

  // CRLF input.
  const std::string crlf = "@r1\r\nACGT\r\n+\r\nIIII\r\n";
  fq::RecordIterator it2(reinterpret_cast<const uint8_t*>(crlf.data()), crlf.size());
  CHECK(it2.next(r));
  CHECK(r.well_formed());
  CHECK_EQ(r.seq_len, std::size_t(4));
}

// ------------------------------------------------------------------ QC stats

void test_qc_stats() {
  section("QC accumulation and merge");
  fq::QcStats a;
  fq::CycleScratch sc;
  const std::string seq = "ACGTNACGTG";
  const std::string qual = "IIIIIIIIII";  // Phred 40
  fq::accumulate_read(a, sc, reinterpret_cast<const uint8_t*>(seq.data()),
                      reinterpret_cast<const uint8_t*>(qual.data()), seq.size());
  a.fold(sc);

  CHECK_EQ(a.reads_in, uint64_t(1));
  CHECK_EQ(a.bases_in, uint64_t(10));
  CHECK_EQ(a.qsum_in, uint64_t(400));
  CHECK_EQ(a.q30, uint64_t(10));
  CHECK_EQ(a.n_bases, uint64_t(1));
  CHECK_EQ(a.gc, uint64_t(5));  // ACGTNACGTG -> C,G at 1,2; C,G at 6,7; G at 9
  CHECK_EQ(a.max_cycle, std::size_t(10));
  CHECK_EQ(a.pos_qsum[0], uint64_t('I'));
  CHECK_EQ(a.pos_base[4][fq::kN], uint64_t(1));

  fq::QcStats b = a;
  fq::QcStats sum;
  sum.merge(a);
  sum.merge(b);
  CHECK_EQ(sum.reads_in, uint64_t(2));
  CHECK_EQ(sum.qsum_in, uint64_t(800));
  CHECK_EQ(sum.pos_base[4][fq::kN], uint64_t(2));
  CHECK(sum.mean_quality() > 39.9 && sum.mean_quality() < 40.1);

  // The 32-bit scratch must fold cleanly over many reads.
  fq::QcStats big;
  fq::CycleScratch bsc;
  for (int i = 0; i < 5000; ++i) {
    fq::accumulate_read(big, bsc, reinterpret_cast<const uint8_t*>(seq.data()),
                        reinterpret_cast<const uint8_t*>(qual.data()), seq.size());
    if (i % 700 == 0) big.fold(bsc);
  }
  big.fold(bsc);
  CHECK_EQ(big.pos_count[0], uint64_t(5000));
  CHECK_EQ(big.pos_qsum[0], uint64_t('I') * 5000);
}

// qhist is a per-BASE histogram (94 bins, no per-cycle structure), unlike
// pos_qsum/pos_count/pos_base above -- privatising it four ways round-robin
// on the base index (i & 3) means every read length must still land the
// right total in each bin regardless of how its length divides by 4: a
// length that leaves a "partial" group of bases is not a tail case here the
// way it would be for a SIMD-width split, since every single base index maps
// to exactly one of the four sub-histograms and back, but it is exactly the
// property most likely to go wrong in an off-by-one across that mapping.
void test_qhist_is_exact_across_every_read_length_and_fold_point() {
  section("qhist survives a 4-way split at every read length");
  fq::QcStats st;
  fq::CycleScratch sc;
  const uint8_t q10 = fq::kPhredOffset + 10;
  const uint8_t q20 = fq::kPhredOffset + 20;
  uint64_t expect_q10 = 0, expect_q20 = 0;
  // Lengths 1..11 cover every residue class of (length mod 4), and folding at
  // irregular points (not a multiple of 4, not aligned with the reads) checks
  // that a fold mid-way through a length's own accumulation loses nothing --
  // it cannot, since fold happens between accumulate_read() calls, never
  // inside one, but the test is cheap enough to pin that explicitly.
  for (std::size_t len = 1; len <= 11; ++len) {
    std::vector<uint8_t> seq(len, 'A');
    std::vector<uint8_t> qual(len);
    for (std::size_t i = 0; i < len; ++i) {
      const bool even = (i % 2 == 0);
      qual[i] = even ? q10 : q20;
      (even ? expect_q10 : expect_q20) += 1;
    }
    fq::accumulate_read(st, sc, seq.data(), qual.data(), len);
    if (len % 3 == 0) st.fold(sc);
  }
  st.fold(sc);

  CHECK_EQ(st.qhist[10], expect_q10);
  CHECK_EQ(st.qhist[20], expect_q20);
  CHECK_EQ(st.bases_in, expect_q10 + expect_q20);
  for (std::size_t q = 0; q < fq::kQualBins; ++q) {
    if (q != 10 && q != 20) CHECK_EQ(st.qhist[q], uint64_t(0));
  }
}

// ------------------------------------------------------------- fixture files

std::string make_fastq(std::size_t n_reads, bool with_adapters) {
  std::mt19937 rng(999);
  std::uniform_int_distribution<int> pick(0, 3);
  const char bases[] = "ACGT";
  std::string out;
  out.reserve(n_reads * 340);
  for (std::size_t i = 0; i < n_reads; ++i) {
    std::string seq;
    seq.reserve(150);
    for (int j = 0; j < 150; ++j) seq.push_back(bases[pick(rng)]);
    if (with_adapters && (i % 5 == 0)) {
      // Simulate a short insert: adapter read-through from base 100.
      seq.replace(100, std::string(fq::kTruSeqAdapter).size(),
                  std::string(fq::kTruSeqAdapter));
    }
    std::string qual(150, 'I');
    if (i % 7 == 0) qual.replace(130, 20, std::string(20, '#'));  // bad 3' tail
    out += "@read_" + std::to_string(i) + " 1:N:0:ATCG\n";
    out += seq + "\n+\n" + qual + "\n";
  }
  return out;
}

void write_plain_gzip(const std::filesystem::path& p, const std::string& data) {
  gzFile f = gzopen(p.string().c_str(), "wb6");
  if (f == nullptr) throw std::runtime_error("gzopen failed");
  std::size_t off = 0;
  while (off < data.size()) {
    const unsigned n = static_cast<unsigned>(std::min<std::size_t>(1 << 20, data.size() - off));
    if (gzwrite(f, data.data() + off, n) != static_cast<int>(n)) {
      gzclose(f);
      throw std::runtime_error("gzwrite failed");
    }
    off += n;
  }
  gzclose(f);
}

// Emit a real BGZF file so the block-parallel path is genuinely exercised,
// rather than assumed to work because the plain path does.
//
// `eof_marker` appends the 28-byte empty member bgzip always writes. The reader
// now requires it -- a stream without one is indistinguishable from a file cut
// between members -- so only the truncation test turns it off.
void write_bgzf(const std::filesystem::path& p, const std::string& data,
                bool eof_marker = true, std::size_t block = fq::kBgzfMaxBlock) {
  libdeflate_compressor* c = libdeflate_alloc_compressor(6);
  if (c == nullptr) throw std::runtime_error("libdeflate_alloc_compressor failed");
  std::ofstream out(p, std::ios::binary);
  std::vector<uint8_t> tmp(fq::kBgzfMaxBlock + 1024);

  std::size_t off = 0;
  while (off < data.size()) {
    const std::size_t n = std::min<std::size_t>(block, data.size() - off);
    const std::size_t clen =
        libdeflate_deflate_compress(c, data.data() + off, n, tmp.data(), tmp.size());
    if (clen == 0) throw std::runtime_error("deflate failed");

    const std::size_t total = 12 + 6 + clen + 8;
    const uint16_t bsize = static_cast<uint16_t>(total - 1);
    uint8_t hdr[18] = {31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 'B', 'C', 2, 0, 0, 0};
    hdr[16] = static_cast<uint8_t>(bsize & 0xff);
    hdr[17] = static_cast<uint8_t>(bsize >> 8);
    out.write(reinterpret_cast<const char*>(hdr), 18);
    out.write(reinterpret_cast<const char*>(tmp.data()), static_cast<std::streamsize>(clen));

    const uint32_t crc = libdeflate_crc32(0, data.data() + off, n);
    const uint32_t isize = static_cast<uint32_t>(n);
    uint8_t trailer[8];
    for (int i = 0; i < 4; ++i) trailer[i] = static_cast<uint8_t>(crc >> (8 * i));
    for (int i = 0; i < 4; ++i) trailer[4 + i] = static_cast<uint8_t>(isize >> (8 * i));
    out.write(reinterpret_cast<const char*>(trailer), 8);
    off += n;
  }
  if (eof_marker) {
    static const uint8_t kEof[28] = {0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43,
                                     0x02, 0,    0x1b, 0,    3, 0, 0, 0, 0, 0,    0,    0, 0,    0};
    out.write(reinterpret_cast<const char*>(kEof), sizeof(kEof));
  }
  libdeflate_free_compressor(c);
}

std::string slurp(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

std::size_t count_records(const std::string& s) {
  std::size_t nl = 0;
  for (char c : s) nl += (c == '\n');
  return nl / 4;
}

// ------------------------------------------------------------- end-to-end

void test_pipeline_end_to_end() {
  section("end-to-end pipeline (raw, gzip, BGZF)");
  const auto dir = std::filesystem::temp_directory_path() / "fastq_stream_test";
  std::filesystem::create_directories(dir);

  const std::string fastq = make_fastq(4000, /*with_adapters=*/true);
  const auto raw_path = dir / "in.fq";
  const auto gz_path = dir / "in.fq.gz";
  const auto bgzf_path = dir / "in.bgzf.fq.gz";

  { std::ofstream f(raw_path, std::ios::binary); f << fastq; }
  write_plain_gzip(gz_path, fastq);
  write_bgzf(bgzf_path, fastq);

  auto run = [&](const std::filesystem::path& in, const char* tag) {
    fq::Config cfg;
    cfg.input = in.string();
    cfg.output = (dir / (std::string("out_") + tag + ".fq")).string();
    cfg.threads = 6;
    cfg.window = 4;
    cfg.window_mean = 20.0;
    cfg.min_len = 15;
    cfg.trim_adapters = true;
    return std::make_pair(fq::run_pipeline(cfg), cfg.output);
  };

  const auto [rep_raw, out_raw] = run(raw_path, "raw");
  const auto [rep_gz, out_gz] = run(gz_path, "gz");
  const auto [rep_bgzf, out_bgzf] = run(bgzf_path, "bgzf");

  CHECK(rep_raw.format == fq::GzipFormat::kRaw);
  CHECK(rep_gz.format == fq::GzipFormat::kPlain);
  CHECK(rep_bgzf.format == fq::GzipFormat::kBgzf);
  CHECK(rep_bgzf.inflater_threads >= 1);  // the parallel path really ran

  // Every input framing must yield byte-identical output and identical stats.
  const std::string a = slurp(out_raw);
  const std::string b = slurp(out_gz);
  const std::string c = slurp(out_bgzf);
  CHECK(a == b);
  CHECK(b == c);
  CHECK(!a.empty());

  CHECK_EQ(rep_raw.stats.reads_in, uint64_t(4000));
  CHECK_EQ(rep_gz.stats.reads_in, uint64_t(4000));
  CHECK_EQ(rep_bgzf.stats.reads_in, uint64_t(4000));
  CHECK_EQ(rep_raw.stats.reads_out, rep_gz.stats.reads_out);
  CHECK_EQ(rep_raw.stats.bases_out, rep_bgzf.stats.bases_out);
  CHECK_EQ(rep_raw.stats.qsum_in, rep_bgzf.stats.qsum_in);
  CHECK_EQ(count_records(a), rep_raw.stats.reads_out);

  // 1 read in 5 carries adapter, 1 in 7 has a bad tail — both must show up.
  CHECK(rep_raw.stats.adapter_trimmed >= 700);
  CHECK(rep_raw.stats.quality_trimmed >= 500);
  CHECK(rep_raw.stats.bases_out < rep_raw.stats.bases_in);

  // Memory ceiling: the whole point of the design.
  CHECK(rep_bgzf.pool_bytes < 512ull * 1024 * 1024);

  // Order preservation: read_0 first, read_3999 last (if it survived).
  CHECK(a.rfind("@read_0 ", 0) == 0);

  // Thread count must not change the answer.
  for (int t : {1, 2, 3, 16}) {
    fq::Config cfg;
    cfg.input = bgzf_path.string();
    cfg.output = (dir / ("out_t" + std::to_string(t) + ".fq")).string();
    cfg.threads = t;
    cfg.window = 4;
    cfg.window_mean = 20.0;
    cfg.min_len = 15;
    const auto rep = fq::run_pipeline(cfg);
    CHECK_EQ(rep.stats.reads_out, rep_raw.stats.reads_out);
    CHECK(slurp(cfg.output) == a);
  }

  // qc-only writes nothing but still counts everything.
  {
    fq::Config cfg;
    cfg.input = gz_path.string();
    cfg.qc_only = true;
    cfg.threads = 4;
    const auto rep = fq::run_pipeline(cfg);
    CHECK_EQ(rep.stats.reads_in, uint64_t(4000));
    CHECK_EQ(rep.bytes_written, uint64_t(0));
  }

  std::filesystem::remove_all(dir);
}

void test_malformed_input() {
  section("malformed input is rejected, not silently truncated");
  const auto dir = std::filesystem::temp_directory_path() / "fastq_stream_bad";
  std::filesystem::create_directories(dir);
  const auto p = dir / "bad.fq";
  { std::ofstream f(p); f << "@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\n"; }  // truncated

  fq::Config cfg;
  cfg.input = p.string();
  cfg.output = (dir / "out.fq").string();
  cfg.threads = 4;
  bool threw = false;
  try {
    fq::run_pipeline(cfg);
  } catch (const std::exception&) {
    threw = true;
  }
  CHECK(threw);

  // Sequence/quality length mismatch.
  const auto p2 = dir / "bad2.fq";
  { std::ofstream f(p2); f << "@r1\nACGTACGT\n+\nIIII\n"; }
  cfg.input = p2.string();
  threw = false;
  try {
    fq::run_pipeline(cfg);
  } catch (const std::exception&) {
    threw = true;
  }
  CHECK(threw);

  std::filesystem::remove_all(dir);
}

// ------------------------------------------------------- security regressions
// One test per closed finding in docs/REVIEW_2026-08-15.md. These are the cases
// the original suites did not reach, so they are written to fail loudly against
// the pre-fix code rather than to re-cover ground already covered above.

// REVIEW finding 1 (High): the offset loop `s <= limit` is inclusive, so
// clamping `limit` to 0 for len < min_overlap still ran one iteration and the
// 4-byte word-seed load read past the end of the read.
void test_adapter_short_read_overread() {
  section("regression: adapter seed scan on sub-min_overlap reads");
  fq::AdapterTrimmer t;
  t.add_defaults();  // seed 4, min_overlap 4
  CHECK_EQ(t.min_overlap(), std::size_t{4});

  // Heap, exact-size allocations: ASan then has a redzone immediately after the
  // read. A stack array (or the pipeline's oversized chunk buffer) hides the
  // overread, which is precisely why this went unnoticed.
  for (std::size_t len = 1; len <= 5; ++len) {
    auto* buf = static_cast<uint8_t*>(std::malloc(len));
    std::memset(buf, 'A', len);
    CHECK_EQ(t.trim(buf, len), len);  // no match possible, no overread
    std::free(buf);
  }

  // Same shape with a read whose bytes are an adapter prefix: the seed *would*
  // match, so this exercises the extend path's length arithmetic too.
  for (std::size_t len = 1; len <= 3; ++len) {
    auto* buf = static_cast<uint8_t*>(std::malloc(len));
    std::memcpy(buf, fq::kTruSeqAdapter.data(), len);
    CHECK_EQ(t.trim(buf, len), len);
    std::free(buf);
  }

  // A larger min_overlap must widen the guard, not just the seed comparison.
  fq::AdapterTrimmer wide;
  wide.add_defaults();
  wide.set_min_overlap(12);
  for (std::size_t len = 1; len <= 11; ++len) {
    auto* buf = static_cast<uint8_t*>(std::malloc(len));
    std::memcpy(buf, fq::kTruSeqAdapter.data(), len);
    CHECK_EQ(wide.trim(buf, len), len);
    std::free(buf);
  }

  // The guard must not suppress a legitimate match at exactly min_overlap.
  {
    const std::size_t len = 4;
    auto* buf = static_cast<uint8_t*>(std::malloc(len));
    std::memcpy(buf, fq::kTruSeqAdapter.data(), len);
    CHECK_EQ(t.trim(buf, len), std::size_t{0});  // adapter starts at offset 0
    std::free(buf);
  }
}

// REVIEW finding 3 (High): mean_phred subtracted len*33 from an unsigned sum,
// so any byte below '!' wrapped the result to ~1.8e19. That value is not merely
// wrong -- it is larger than any --min-mean-q, so the drop filter inverted and
// kept exactly the corrupt reads it exists to reject.
void test_mean_phred_underflow() {
  section("regression: mean_phred underflow below the Phred offset");

  const uint8_t bad[4] = {1, 2, 3, 4};  // every byte below '!' (33)
  CHECK_EQ(fq::mean_phred(bad, 4), 0.0);

  const uint8_t mixed[4] = {33 + 40, 1, 33 + 40, 2};
  CHECK_EQ(fq::mean_phred(mixed, 4), 20.0);  // (40 + 0 + 40 + 0) / 4

  // Zero bytes are the common real-world case: a truncated or NUL-padded
  // quality line. Before the fix this returned ~1.8e19.
  const uint8_t zeros[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  CHECK_EQ(fq::mean_phred(zeros, 8), 0.0);

  // Valid input must be unchanged by the clamp.
  std::vector<uint8_t> ok(64, static_cast<uint8_t>(33 + 30));
  CHECK_EQ(fq::mean_phred(ok.data(), ok.size()), 30.0);
  CHECK_EQ(fq::mean_phred(ok.data(), 0), 0.0);

  // The inverted-filter consequence, stated as an ordering property: a corrupt
  // read must now sort below any threshold a caller would set.
  CHECK(fq::mean_phred(bad, 4) < 20.0);
}


// Regression: trim_sliding_window's threshold cast was undefined behaviour.
//
// The threshold was computed as
//     static_cast<uint32_t>(min_mean * static_cast<double>(window));
// which is UB whenever the truncated value is not representable as uint32_t.
// All three of these reach it straight from the CLI, and UBSan reported
// "outside the range of representable values of type 'unsigned int'" on each:
//
//     --window-mean -5      ->  -20
//     --window-mean 1e30    ->  4e+30
//     --window-mean nan     ->  nan
//
// main.cpp now range-checks the option, but the kernel is a public header that
// tests and library consumers call directly, so it clamps as well. NaN is
// handled first and explicitly, since every comparison against NaN is false and
// a plain `<= 0` test would fall straight through to the cast.
void test_sliding_window_threshold_domain() {
  // A uniformly high-quality read: 40 bases at Q40 ('I').
  std::vector<uint8_t> good(40, static_cast<uint8_t>(33 + 40));
  // A uniformly low-quality read: 40 bases at Q2.
  std::vector<uint8_t> bad(40, static_cast<uint8_t>(33 + 2));

  // Negative and zero thresholds mean "no bound"; nothing is trimmed.
  for (double m : {-1e30, -20.0, -0.5, 0.0}) {
    CHECK_EQ(fq::trim_sliding_window(good.data(), good.size(), 4, m), good.size());
    CHECK_EQ(fq::trim_sliding_window(bad.data(), bad.size(), 4, m), bad.size());
  }

  // NaN is not a bound either: keep the read rather than invoke the cast.
  const double nan_v = std::numeric_limits<double>::quiet_NaN();
  CHECK_EQ(fq::trim_sliding_window(good.data(), good.size(), 4, nan_v), good.size());
  CHECK_EQ(fq::trim_sliding_window(bad.data(), bad.size(), 4, nan_v), bad.size());

  // Absurdly high thresholds saturate instead of wrapping. No window can meet
  // them, so the cut is at position 0 for every read.
  //
  // The saturation value matters: `limit` and the window's ASCII offset are
  // added together, and clamping `limit` to UINT32_MAX in a uint32_t sum would
  // wrap that addition to `offset - 1` and invert the comparison on exactly
  // this input. The sum is done in uint64_t for that reason.
  for (double m : {1e9, 1e30, 1e300, std::numeric_limits<double>::infinity()}) {
    CHECK_EQ(fq::trim_sliding_window(good.data(), good.size(), 4, m), 0u);
    CHECK_EQ(fq::trim_sliding_window(bad.data(), bad.size(), 4, m), 0u);
  }

  // And the ordinary path still behaves: Q40 survives a Q20 bound, Q2 does not.
  CHECK_EQ(fq::trim_sliding_window(good.data(), good.size(), 4, 20.0), good.size());
  CHECK_EQ(fq::trim_sliding_window(bad.data(), bad.size(), 4, 20.0), 0u);

  // Window wider than the read is a no-op regardless of the threshold.
  CHECK_EQ(fq::trim_sliding_window(good.data(), 3, 4, 1e30), 3u);
  CHECK_EQ(fq::trim_sliding_window(good.data(), good.size(), 0, 1e30), good.size());
}

// ------------------------------------------ torture regressions, 2026-09-10
// One test per defect in docs/TORTURE_2026-09-10.md. Each was run against the
// pre-fix code and fails there.

void write_bytes(const std::filesystem::path& p, const std::string& s) {
  std::ofstream f(p, std::ios::binary);
  f << s;
}

// Runs the pipeline, returning the report or an exception message.
struct RunResult {
  bool threw = false;
  std::string message;
  fq::RunReport rep;
};
RunResult try_run(const fq::Config& cfg) {
  RunResult r;
  try {
    r.rep = fq::run_pipeline(cfg);
  } catch (const std::exception& e) {
    r.threw = true;
    r.message = e.what();
  }
  return r;
}

fq::Config plain_cfg(const std::filesystem::path& in, const std::filesystem::path& out) {
  fq::Config cfg;
  cfg.input = in.string();
  cfg.output = out.string();
  cfg.threads = 4;
  cfg.window = 0;
  cfg.min_len = 0;
  cfg.trim_adapters = false;
  return cfg;
}

// Truncated compressed input decoded "successfully". A gzip missing its CRC
// trailer, a gzip cut to its header, and a BGZF file cut between members all
// exited 0 -- the last with half the reads.
void test_truncated_compressed_input() {
  section("regression: truncated gzip/BGZF input is an error, not fewer reads");
  const auto dir = std::filesystem::temp_directory_path() / "fastq_stream_trunc";
  std::filesystem::create_directories(dir);
  const std::string fastq = make_fastq(3000, false);
  const auto out = dir / "out.fq";

  const auto gz = dir / "full.fq.gz";
  write_plain_gzip(gz, fastq);
  const std::string gz_bytes = slurp(gz);
  {
    const auto r = try_run(plain_cfg(gz, out));
    CHECK(!r.threw);
    CHECK_EQ(r.rep.stats.reads_in, uint64_t(3000));
  }
  for (const std::size_t keep : {gz_bytes.size() - 8, gz_bytes.size() - 1,
                                 gz_bytes.size() / 2, std::size_t(10), std::size_t(2)}) {
    const auto p = dir / "cut.fq.gz";
    write_bytes(p, gz_bytes.substr(0, keep));
    const auto r = try_run(plain_cfg(p, out));
    CHECK(r.threw);
  }
  // Two whole gzip members back to back are one valid stream.
  {
    const auto p = dir / "two.fq.gz";
    write_bytes(p, gz_bytes + gz_bytes);
    const auto r = try_run(plain_cfg(p, out));
    CHECK(!r.threw);
    CHECK_EQ(r.rep.stats.reads_in, uint64_t(6000));
  }

  const auto bgzf_full = dir / "full.bgz";
  const auto bgzf_nomark = dir / "nomark.bgz";
  write_bgzf(bgzf_full, fastq);
  write_bgzf(bgzf_nomark, fastq, /*eof_marker=*/false);
  const std::string full = slurp(bgzf_full);
  {
    const auto r = try_run(plain_cfg(bgzf_full, out));
    CHECK(!r.threw);
    CHECK_EQ(r.rep.stats.reads_in, uint64_t(3000));
  }
  {
    // Every member intact, the marker missing: the truncated-between-members
    // shape. Pre-fix this passed.
    const auto r = try_run(plain_cfg(bgzf_nomark, out));
    CHECK(r.threw);
    CHECK(r.message.find("end-of-file marker") != std::string::npos);
  }
  for (const std::size_t junk : {1, 17}) {
    const auto p = dir / "junk.bgz";
    write_bytes(p, full + std::string(junk, '\0'));
    CHECK(try_run(plain_cfg(p, out)).threw);
  }
  {
    // Concatenated BGZF files put an empty member mid-stream. Legal, and it
    // must not be mistaken for the end.
    const auto p = dir / "cat.bgz";
    write_bytes(p, full + full);
    const auto r = try_run(plain_cfg(p, out));
    CHECK(!r.threw);
    CHECK_EQ(r.rep.stats.reads_in, uint64_t(6000));
  }
  {
    // A header claiming BSIZE=0 describes a member smaller than itself.
    //
    // On a later member that is a malformed BGZF stream and must throw. On the
    // FIRST member it is not: detection then sees no valid BGZF header and
    // takes the plain-gzip path, and a BGZF file is a valid concatenated gzip
    // stream whatever its extra field says -- zlib inflates it with every CRC
    // checked. The property there is "no silently wrong answer": the run
    // either fails or reports every read.
    const std::size_t first =
        (static_cast<unsigned char>(full[16]) | (static_cast<unsigned char>(full[17]) << 8)) + 1u;
    std::string later = full;
    later[first + 16] = 0;
    later[first + 17] = 0;
    const auto p = dir / "bsize0.bgz";
    write_bytes(p, later);
    CHECK(try_run(plain_cfg(p, out)).threw);

    std::string head = full;
    head[16] = 0;
    head[17] = 0;
    write_bytes(p, head);
    const auto r = try_run(plain_cfg(p, out));
    CHECK(r.threw || r.rep.stats.reads_in == 3000);
  }
  std::filesystem::remove_all(dir);
}

// `-i x.fq -o x.fq` truncated the input before reading it.
void test_output_aliasing_input() {
  section("regression: output that is the input is refused before it is truncated");
  const auto dir = std::filesystem::temp_directory_path() / "fastq_stream_alias";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);
  const std::string fastq = make_fastq(3000, false);  // > the reader's 1 MiB stage
  const auto in = dir / "in.fq";
  write_bytes(in, fastq);
  std::filesystem::create_hard_link(in, dir / "hard.fq");
  std::filesystem::create_symlink(in, dir / "soft.fq");

  for (const char* name : {"in.fq", "hard.fq", "soft.fq"}) {
    const auto r = try_run(plain_cfg(in, dir / name));
    CHECK(r.threw);
    CHECK(r.message.find("same file") != std::string::npos);
    CHECK(slurp(in) == fastq);
  }
  // --qc-only writes nothing, so the same path is harmless there.
  fq::Config qc = plain_cfg(in, in);
  qc.qc_only = true;
  CHECK(!try_run(qc).threw);
  CHECK(slurp(in) == fastq);
  std::filesystem::remove_all(dir);
}

// A last record with no final newline failed the whole run.
void test_final_record_without_newline() {
  section("regression: last record without a trailing newline is accepted");
  const auto dir = std::filesystem::temp_directory_path() / "fastq_stream_eol";
  std::filesystem::create_directories(dir);
  const auto in = dir / "in.fq";
  const auto out = dir / "out.fq";
  const std::string r1 = "@r1\nACGTACGT\n+\nIIIIIIII\n";

  write_bytes(in, r1 + "@r2\nGGGG\n+\nIIII");
  auto r = try_run(plain_cfg(in, out));
  CHECK(!r.threw);
  CHECK_EQ(r.rep.stats.reads_in, uint64_t(2));
  CHECK(slurp(out) == r1 + "@r2\nGGGG\n+\nIIII\n");

  // The same across a chunk boundary: the record that lacks a newline is also
  // the one carried over from the previous chunk.
  const std::string big = make_fastq(2000, false);
  write_bytes(in, big.substr(0, big.size() - 1));
  r = try_run(plain_cfg(in, out));
  CHECK(!r.threw);
  CHECK_EQ(r.rep.stats.reads_in, uint64_t(2000));
  CHECK(slurp(out) == big);

  write_bytes(in, "@r1\r\nACGT\r\n+\r\nIIII\r");  // CRLF file, final LF missing
  r = try_run(plain_cfg(in, out));
  CHECK(!r.threw);
  CHECK_EQ(r.rep.stats.reads_in, uint64_t(1));

  write_bytes(in, r1 + "\n\n\r\n");  // trailing blank lines
  r = try_run(plain_cfg(in, out));
  CHECK(!r.threw);
  CHECK_EQ(r.rep.stats.reads_in, uint64_t(1));

  // Still rejected: a record genuinely cut short, with or without a newline.
  for (const char* tail : {"@r2\nGGGG\n+\n", "@r2\nGGGG\n+", "@r2\nGGGG\n+\nII", "@r2"}) {
    write_bytes(in, r1 + tail);
    CHECK(try_run(plain_cfg(in, out)).threw);
  }
  std::filesystem::remove_all(dir);
}

// NUL, control and high-bit bytes in SEQ/QUAL were clamped to Q0 / counted as
// N and the read was trimmed, filtered or emitted as if it were data.
void test_non_printable_bytes_rejected() {
  section("regression: non-printable SEQ/QUAL bytes fail the run");
  const auto dir = std::filesystem::temp_directory_path() / "fastq_stream_bytes";
  std::filesystem::create_directories(dir);
  const auto in = dir / "in.fq";
  const auto out = dir / "out.fq";

  const std::string seq = "ACGTACGTAC";
  const std::string qual = "IIIIIIIIII";
  const auto rec = [](const std::string& s, const std::string& q) {
    return "@r\n" + s + "\n+\n" + q + "\n";
  };
  for (const char bad : {'\0', '\x01', '\x1f', ' ', '\x7f', '\x80', '\xff'}) {
    std::string s = seq, q = qual;
    q[5] = bad;
    write_bytes(in, rec(seq, q));
    CHECK(try_run(plain_cfg(in, out)).threw);
    s[5] = bad;
    write_bytes(in, rec(s, qual));
    CHECK(try_run(plain_cfg(in, out)).threw);
  }
  // The printable extremes are data.
  write_bytes(in, rec("NNACGT.acg", "!!!!!~~~~~"));
  CHECK(!try_run(plain_cfg(in, out)).threw);

  fq::QcStats st;
  fq::CycleScratch sc;
  const uint8_t good_s[4] = {'A', 'N', 'G', 'x'};
  const uint8_t good_q[4] = {'!', 'I', '~', '5'};
  const uint8_t bad_q[4] = {'I', 0, 'I', 'I'};
  const auto a = fq::accumulate_read(st, sc, good_s, good_q, 4);
  CHECK(a.valid);
  CHECK_EQ(a.n_bases, std::size_t(2));  // 'N' and 'x'
  CHECK(!fq::accumulate_read(st, sc, good_s, bad_q, 4).valid);

  // The kernel behind `valid`, against a scalar reference: every length across
  // the 16/32-byte vector boundaries, every byte value at every position.
  int kernel_mismatches = 0;
  for (std::size_t n = 0; n <= 70; ++n) {
    std::vector<uint8_t> buf(n + 1, 'I');
    kernel_mismatches += !fq::simd::all_printable(buf.data(), n);
    for (std::size_t pos = 0; pos < n; ++pos) {
      for (int v = 0; v < 256; ++v) {
        buf[pos] = static_cast<uint8_t>(v);
        const bool want = v >= 33 && v <= 126;
        kernel_mismatches += fq::simd::all_printable(buf.data(), n) != want;
      }
      buf[pos] = 'I';
    }
    buf[n] = 0;  // a bad byte just past the range must not be read
    kernel_mismatches += !fq::simd::all_printable(buf.data(), n);
  }
  CHECK_EQ(kernel_mismatches, 0);
  std::filesystem::remove_all(dir);
}

// accumulate_read counts Q20/Q30 on RAW bytes while qsum and the histogram use
// the clamped Phred value; the report treats them as one distribution. That
// holds only because `qual >= 33 + t` and the clamped `q >= t` agree for every
// byte value, which is checked here rather than assumed -- including bytes the
// pipeline now rejects, since the header is callable without that check.
void test_q20_q30_exact_over_all_bytes() {
  section("q20/q30 from the clamped loop equal count_ge on raw bytes");
  std::vector<uint8_t> qual(256), seq(256, 'A');
  for (int i = 0; i < 256; ++i) qual[static_cast<std::size_t>(i)] = static_cast<uint8_t>(i);
  fq::QcStats st;
  fq::CycleScratch sc;
  fq::accumulate_read(st, sc, seq.data(), qual.data(), qual.size());
  CHECK_EQ(st.q20, uint64_t(ref_count_ge(qual.data(), qual.size(), 33 + 20)));
  CHECK_EQ(st.q30, uint64_t(ref_count_ge(qual.data(), qual.size(), 33 + 30)));
}

// With the default --window 4, --trim-front/--trim-tail trims were not counted.
void test_quality_trim_accounting() {
  section("regression: quality_trimmed counts front/tail trims with the window on");
  const auto dir = std::filesystem::temp_directory_path() / "fastq_stream_qtrim";
  std::filesystem::create_directories(dir);
  const auto in = dir / "in.fq";
  const std::string good(30, 'I');
  write_bytes(in, "@front\n" + std::string(40, 'A') + "\n+\n" + std::string(10, '#') + good +
                      "\n@tail\n" + std::string(40, 'A') + "\n+\n" + good +
                      std::string(10, '#') + "\n@clean\n" + std::string(40, 'A') + "\n+\n" +
                      std::string(40, 'I') + "\n");
  fq::Config cfg = plain_cfg(in, dir / "out.fq");
  cfg.front_q = 20;
  cfg.tail_q = 20;
  cfg.window = 4;
  cfg.window_mean = 20.0;
  const auto r = try_run(cfg);
  CHECK(!r.threw);
  CHECK_EQ(r.rep.stats.quality_trimmed, uint64_t(2));
  std::filesystem::remove_all(dir);
}

// The JSON report embeds the input path; only \n and \t were escaped.
void test_json_escapes_control_characters() {
  section("regression: JSON report escapes every control character");
  fq::Config cfg;
  cfg.input = std::string("a\rb\x01") + "c\x1f.fq";
  fq::RunReport rep;
  std::ostringstream os;
  fq::write_json_report(os, cfg, rep);
  const std::string json = os.str();
  bool raw_control = false;
  for (const char c : json) raw_control |= (static_cast<unsigned char>(c) < 0x20 && c != '\n');
  CHECK(!raw_control);
  CHECK(json.find("a\\u000db\\u0001c\\u001f.fq") != std::string::npos);
}

// --------------------------------------- audit regressions, 2026-09-11
// One test per finding in docs/AUDIT_2026-09-11_fastq_stream.md, plus pins for
// the attacks that found nothing, so the clean results stay clean.

std::string fq_record(const std::string& name, const std::string& seq, const std::string& qual) {
  return "@" + name + "\n" + seq + "\n+\n" + qual + "\n";
}

std::string random_bases(std::mt19937& rng, std::size_t n) {
  static const char kBases[] = "ACGT";
  std::string s(n, 'A');
  for (auto& c : s) c = kBases[rng() % 4];
  return s;
}

// A1: Phred+64 was read as Phred+33 -- mean Q +31, Q30 ~100%, no trimming, exit 0.
void test_phred64_is_converted_and_detected() {
  section("audit A1: Phred+64 input is converted, and detected when not declared");
  const auto dir = std::filesystem::temp_directory_path() / "fastq_stream_p64";
  std::filesystem::create_directories(dir);
  std::mt19937 rng(64);
  std::string p33, p64;
  for (int i = 0; i < 600; ++i) {
    const std::string seq = random_bases(rng, 80);
    std::string q33, q64;
    for (int j = 0; j < 80; ++j) {
      const int q = (i % 5 == 0 && j >= 55) ? 2 : 2 + static_cast<int>(rng() % 39);
      q33 += static_cast<char>(33 + q);
      q64 += static_cast<char>(64 + q);
    }
    p33 += fq_record("r" + std::to_string(i), seq, q33);
    p64 += fq_record("r" + std::to_string(i), seq, q64);
  }
  write_bytes(dir / "p33.fq", p33);
  write_bytes(dir / "p64.fq", p64);

  auto trimming = [&](const char* in, const char* out) {
    fq::Config cfg = plain_cfg(dir / in, dir / out);
    cfg.window = 4;
    cfg.window_mean = 20.0;
    cfg.min_len = 15;
    return cfg;
  };
  const auto r33 = try_run(trimming("p33.fq", "out33.fq"));
  fq::Config c64 = trimming("p64.fq", "out64.fq");
  c64.phred_offset = 64;
  const auto r64 = try_run(c64);
  CHECK(!r33.threw);
  CHECK(!r64.threw);
  // Converted in place, so everything downstream sees the +33 answer and the
  // output FASTQ is byte-identical to the +33 run's.
  CHECK(slurp(dir / "out33.fq") == slurp(dir / "out64.fq"));
  CHECK_EQ(r64.rep.stats.qsum_in, r33.rep.stats.qsum_in);
  CHECK_EQ(r64.rep.stats.q30, r33.rep.stats.q30);
  CHECK_EQ(r64.rep.stats.quality_trimmed, r33.rep.stats.quality_trimmed);
  CHECK_EQ(r64.rep.stats.reads_out, r33.rep.stats.reads_out);
  CHECK(r33.rep.stats.quality_trimmed > 100);

  // Undeclared, it runs at +33 -- and the detector names it.
  const auto wrong = try_run(trimming("p64.fq", "wrong.fq"));
  CHECK(!wrong.threw);
  CHECK(fq::looks_like_phred64(wrong.rep.stats));
  CHECK(!fq::looks_like_phred64(r33.rep.stats));
  CHECK(wrong.rep.stats.quality_trimmed < r33.rep.stats.quality_trimmed);

  CHECK(!fq::looks_like_phred64(fq::QcStats{}));

  // A byte below '@' is not Phred+64 (Solexa's negative scores, or a +33 file
  // declared wrong), and is refused rather than wrapped. '~' is still in range.
  write_bytes(dir / "bad64.fq", fq_record("r", "ACGT", "II5I"));
  fq::Config bad = plain_cfg(dir / "bad64.fq", dir / "bad.out");
  bad.phred_offset = 64;
  const auto rb = try_run(bad);
  CHECK(rb.threw);
  CHECK(rb.message.find("not Phred+64") != std::string::npos);
  write_bytes(dir / "top64.fq", fq_record("r", "ACGT", "@@~~"));
  bad.input = (dir / "top64.fq").string();
  CHECK(!try_run(bad).threw);
  CHECK(slurp(dir / "bad.out") == fq_record("r", "ACGT", "!!__"));
  std::filesystem::remove_all(dir);
}

// The Phred+64 detector against real platform profiles. Its first version --
// nothing below Q31, something above Q41 -- flagged AVITI and HiFi runs with no
// base under Q31 (both exited 2). Every Phred+33 profile here must pass; every
// Phred+64 one must be caught, except the low-quality case the header states
// as a known miss.
void test_phred64_detector_spares_high_accuracy_phred33() {
  section("audit A1: the Phred+64 detector spares AVITI/HiFi Phred+33");
  const auto dir = std::filesystem::temp_directory_path() / "fastq_stream_p64profiles";
  std::filesystem::create_directories(dir);
  std::mt19937 rng(42);
  auto uniform = [&](int lo, int hi) { return lo + static_cast<int>(rng() % (hi - lo + 1)); };
  auto flagged = [&](const char* name, const std::function<int()>& phred33) {
    std::string data;
    for (int i = 0; i < 400; ++i) {
      std::string q;
      for (int j = 0; j < 150; ++j) q += static_cast<char>(33 + phred33());
      data += fq_record(std::string(name) + std::to_string(i), random_bases(rng, 150), q);
    }
    write_bytes(dir / "p.fq", data);
    fq::Config cfg = plain_cfg(dir / "p.fq", dir / "out.fq");
    cfg.qc_only = true;
    const auto r = try_run(cfg);
    CHECK(!r.threw);
    return fq::looks_like_phred64(r.rep.stats);
  };
  // Phred+33 with nothing below Q31.
  CHECK(!flagged("aviti", [&] { return uniform(31, 50); }));
  CHECK(!flagged("aviti44", [&] { return rng() % 10 ? 44 : uniform(38, 50); }));
  CHECK(!flagged("hifi", [&] { return rng() % 5 < 3 ? 93 : uniform(31, 92); }));
  CHECK(!flagged("hificap", [&] { return uniform(42, 60); }));
  CHECK(!flagged("novaseq", [&] { return std::array<int, 4>{2, 12, 23, 37}[rng() % 4]; }));
  // Phred+64, expressed as the +33 value its byte reads as (q + 31).
  CHECK(flagged("ill13", [&] { return 31 + uniform(0, 40); }));
  CHECK(flagged("ill15", [&] { return 31 + (rng() % 10 ? uniform(15, 41) : 2); }));
  CHECK(flagged("ill_hq", [&] { return 31 + uniform(30, 41); }));
  // The stated miss: Phred+64 with no base at Q30+ never reaches Q61 at 33.
  CHECK(!flagged("ill_poor", [&] { return 31 + uniform(2, 29); }));
  std::filesystem::remove_all(dir);
}

// A2: control and non-ASCII bytes in the header and '+' lines reached the output.
void test_header_and_plus_lines_are_validated() {
  section("audit A2: header and '+' lines are validated");
  const auto dir = std::filesystem::temp_directory_path() / "fastq_stream_hdr";
  std::filesystem::create_directories(dir);
  const auto in = dir / "in.fq";
  const auto out = dir / "out.fq";
  for (const char bad : {'\0', '\x01', '\x07', '\x1b', '\r', '\x7f', '\x80', '\xff'}) {
    write_bytes(in, "@r" + std::string(1, bad) + "x\nACGT\n+\nIIII\n");
    const auto r = try_run(plain_cfg(in, out));
    CHECK(r.threw);
    CHECK(r.message.find("header or '+' line") != std::string::npos);
    // A CR as the LAST byte of a line is a CRLF ending and is stripped; the
    // header case above puts it mid-line, where it is not.
    if (bad == '\r') continue;
    write_bytes(in, "@r\nACGT\n+" + std::string(1, bad) + "\nIIII\n");
    CHECK(try_run(plain_cfg(in, out)).threw);
  }
  // What real headers carry: spaces, tabs, colons, and a '+' that repeats the
  // name. CRLF endings are stripped before the check.
  const std::string good = "@A00123:8:HXXX:1:1101:1000:2000 1:N:0:ACGT\tBC:Z:x\nACGT\n+A00123 x\nIIII\n";
  write_bytes(in, good);
  auto r = try_run(plain_cfg(in, out));
  CHECK(!r.threw);
  CHECK(slurp(out) == good);
  write_bytes(in, "@r1\r\nACGT\r\n+\r\nIIII\r\n");
  r = try_run(plain_cfg(in, out));
  CHECK(!r.threw);
  CHECK(slurp(out) == "@r1\nACGT\n+\nIIII\n");
  std::filesystem::remove_all(dir);
}

// A3: whether a record over the carry limit was accepted depended on where it
// fell relative to a chunk boundary, and so on the input's framing.
void test_record_limit_does_not_depend_on_position() {
  section("audit A3: the record limit is the same at every position and framing");
  const auto dir = std::filesystem::temp_directory_path() / "fastq_stream_reclimit";
  std::filesystem::create_directories(dir);
  std::mt19937 rng(3);

  // Short records totalling exactly `nbytes`: the last one is sized to land.
  auto filler = [&](std::size_t nbytes) {
    std::string out;
    int i = 0;
    const std::size_t unit = fq_record("s00000", std::string(100, 'A'), std::string(100, 'I')).size();
    while (nbytes - out.size() >= unit + 13) {
      char name[32];
      std::snprintf(name, sizeof(name), "s%05d", i++);
      out += fq_record(name, random_bases(rng, 100), std::string(100, 'I'));
    }
    const std::size_t rem = nbytes - out.size();           // 13 .. unit+12
    const std::size_t name_len = (rem % 2 == 0) ? 6 : 5;   // "@" name "\n" b "\n+\n" b "\n"
    const std::size_t body = (rem - name_len - 6) / 2;
    out += fq_record(std::string(name_len, 'p'), random_bases(rng, body), std::string(body, 'I'));
    CHECK_EQ(out.size(), nbytes);
    return out;
  };
  std::string tail;
  for (int i = 0; i < 1500; ++i) tail += fq_record("t" + std::to_string(i), random_bases(rng, 100), std::string(100, 'I'));

  const std::string over = fq_record("long", random_bases(rng, 40000), std::string(40000, 'I'));   // 80,010 B
  const std::string under = fq_record("long", random_bases(rng, 32000), std::string(32000, 'I'));  // 64,010 B
  CHECK(over.size() > fq::kPrefix);
  CHECK(under.size() <= fq::kPrefix);

  const std::size_t positions[] = {0, 100000, fq::kChunkBytes - 70000, fq::kChunkBytes - 30000};
  for (const std::size_t pos : positions) {
    for (const bool is_over : {true, false}) {
      const std::string data = (pos ? filler(pos) : std::string()) + (is_over ? over : under) + tail;
      write_bytes(dir / "r.fq", data);
      write_plain_gzip(dir / "r.fq.gz", data);
      write_bgzf(dir / "r.fq.bgz", data);
      for (const char* name : {"r.fq", "r.fq.gz", "r.fq.bgz"}) {
        fq::Config cfg = plain_cfg(dir / name, dir / "out.fq");
        cfg.qc_only = true;
        const auto r = try_run(cfg);
        const std::string where = std::string(name) + " at " + std::to_string(pos);
        if (is_over) {
          if (!r.threw) std::cerr << "   over-limit record accepted: " << where << "\n";
          CHECK(r.threw);
        } else {
          if (r.threw) std::cerr << "   under-limit record refused: " << where << ": " << r.message << "\n";
          CHECK(!r.threw);
        }
      }
    }
  }
  std::filesystem::remove_all(dir);
}

// A4: reads past the per-cycle table were silent.
void test_reads_over_the_cycle_limit_are_counted() {
  section("audit A4: reads longer than the cycle table are counted");
  const auto dir = std::filesystem::temp_directory_path() / "fastq_stream_cycles";
  std::filesystem::create_directories(dir);
  std::mt19937 rng(4);
  std::string data;
  for (const std::size_t len : {100u, 1024u, 1025u, 2000u}) {
    data += fq_record("L" + std::to_string(len), random_bases(rng, len), std::string(len, 'I'));
  }
  write_bytes(dir / "in.fq", data);
  fq::Config cfg = plain_cfg(dir / "in.fq", dir / "out.fq");
  cfg.qc_only = true;
  const auto r = try_run(cfg);
  CHECK(!r.threw);
  CHECK_EQ(r.rep.stats.reads_over_cycle_limit, uint64_t(2));
  CHECK_EQ(r.rep.stats.max_cycle, fq::kMaxCycles);
  std::ostringstream json;
  fq::write_json_report(json, cfg, r.rep);
  CHECK(json.str().find("\"reads_longer_than_cycle_limit\": 2") != std::string::npos);
  CHECK(json.str().find("\"cycle_limit\": 1024") != std::string::npos);
  std::filesystem::remove_all(dir);
}

// A5 (clean in the audit, pinned): CRLF and mixed endings give the LF answer,
// across chunk boundaries and with every trim on. CR-only is named as such.
void test_crlf_and_mixed_line_endings_match_lf() {
  section("audit: CRLF and mixed line endings match LF; CR-only is named");
  const auto dir = std::filesystem::temp_directory_path() / "fastq_stream_eol2";
  std::filesystem::create_directories(dir);
  const std::string lf = make_fastq(3000, true);   // ~1 MB: several chunks
  std::string crlf, mixed, cr_only;
  std::size_t line = 0;
  for (const char c : lf) {
    if (c == '\n') {
      crlf += "\r\n";
      mixed += (line++ % 3 == 0) ? "\r\n" : "\n";
      cr_only += '\r';
    } else {
      crlf += c;
      mixed += c;
      cr_only += c;
    }
  }
  auto run = [&](const std::string& data, const char* tag) {
    write_bytes(dir / (std::string(tag) + ".fq"), data);
    fq::Config cfg = plain_cfg(dir / (std::string(tag) + ".fq"), dir / (std::string(tag) + ".out"));
    cfg.window = 4;
    cfg.window_mean = 20.0;
    cfg.min_len = 15;
    cfg.trim_adapters = true;
    return try_run(cfg);
  };
  const auto a = run(lf, "lf");
  const auto b = run(crlf, "crlf");
  const auto c = run(mixed, "mixed");
  CHECK(!a.threw && !b.threw && !c.threw);
  CHECK(slurp(dir / "lf.out") == slurp(dir / "crlf.out"));
  CHECK(slurp(dir / "lf.out") == slurp(dir / "mixed.out"));
  for (const auto* r : {&b, &c}) {
    CHECK_EQ(r->rep.stats.qsum_in, a.rep.stats.qsum_in);
    CHECK_EQ(r->rep.stats.adapter_trimmed, a.rep.stats.adapter_trimmed);
    CHECK_EQ(r->rep.stats.quality_trimmed, a.rep.stats.quality_trimmed);
    CHECK_EQ(r->rep.stats.bases_out, a.rep.stats.bases_out);
  }
  // CR-only, over a chunk (the carry-limit path) and under one (end of stream).
  for (const std::size_t size : {cr_only.size(), std::size_t(4000)}) {
    const auto r = run(cr_only.substr(0, size), "cr");
    CHECK(r.threw);
    CHECK(r.message.find("CR-only") != std::string::npos);
  }
  std::filesystem::remove_all(dir);
}

// A6 (clean in the audit, pinned): every total exact over many chunks. The
// audit ran 32M records through a pipe; this is the same arithmetic at a size
// a unit suite can afford, still spanning dozens of chunks.
void test_accumulators_are_exact_over_many_chunks() {
  section("audit: QC accumulators are exact over many chunks");
  const auto dir = std::filesystem::temp_directory_path() / "fastq_stream_deep";
  std::filesystem::create_directories(dir);
  // GC/"5?" = Q20,Q30   AT/"!~" = Q0,Q93   NN/"II" = Q40,Q40   GA/"++" = Q10,Q10
  const std::string unit = "@\nGC\n+\n5?\n@\nAT\n+\n!~\n@\nNN\n+\nII\n@\nGA\n+\n++\n";
  constexpr std::uint64_t kUnits = 400000;   // 1.6M records, 16 MB
  {
    std::ofstream f(dir / "deep.fq", std::ios::binary);
    std::string block;
    for (int i = 0; i < 1000; ++i) block += unit;
    for (std::uint64_t i = 0; i < kUnits / 1000; ++i) f << block;
  }
  fq::Config cfg = plain_cfg(dir / "deep.fq", dir / "out.fq");
  cfg.qc_only = true;
  cfg.max_n_rate = -1.0;
  cfg.threads = 8;
  const auto r = try_run(cfg);
  CHECK(!r.threw);
  const fq::QcStats& s = r.rep.stats;
  CHECK_EQ(s.reads_in, 4 * kUnits);
  CHECK_EQ(s.bases_in, 8 * kUnits);
  CHECK_EQ(s.gc, 3 * kUnits);
  CHECK_EQ(s.n_bases, 2 * kUnits);
  CHECK_EQ(s.q20, 5 * kUnits);
  CHECK_EQ(s.q30, 4 * kUnits);
  CHECK_EQ(s.qsum_in, 243 * kUnits);
  CHECK_EQ(s.qhist[0], kUnits);
  CHECK_EQ(s.qhist[10], 2 * kUnits);
  CHECK_EQ(s.qhist[20], kUnits);
  CHECK_EQ(s.qhist[30], kUnits);
  CHECK_EQ(s.qhist[40], 2 * kUnits);
  CHECK_EQ(s.qhist[93], kUnits);   // '~' has its own bin now; the report folds it into 63
  std::ostringstream json;
  fq::write_json_report(json, cfg, r.rep);
  const std::string js = json.str();
  const std::size_t h = js.find("\"quality_histogram\": [");
  const std::string hist = js.substr(h, js.find(']', h) - h);
  CHECK_EQ(std::count(hist.begin(), hist.end(), ','), 63);   // still 64 bins
  CHECK(hist.rfind("," + std::to_string(kUnits)) != std::string::npos);
  CHECK_EQ(s.pos_base[0][fq::kG], 2 * kUnits);
  CHECK_EQ(s.pos_base[1][fq::kC], kUnits);
  CHECK_EQ(s.pos_qsum[0], uint64_t('5' + '!' + 'I' + '+') * kUnits);
  CHECK_EQ(s.pos_qsum[1], uint64_t('?' + '~' + 'I' + '+') * kUnits);
  std::filesystem::remove_all(dir);
}

// G2's lesson, swept (clean in the audit, pinned): a BGZF file cut at ANY
// member boundary fails, not only at the one the F2 test happened to pick.
void test_bgzf_cut_at_every_member_boundary_fails() {
  section("audit: BGZF cut at every member boundary fails");
  const auto dir = std::filesystem::temp_directory_path() / "fastq_stream_bgzfsweep";
  std::filesystem::create_directories(dir);
  const std::string data = make_fastq(300, false);
  write_bgzf(dir / "full.bgz", data, /*eof_marker=*/true, /*block=*/3000);
  const std::string full = slurp(dir / "full.bgz");
  {
    const auto r = try_run(plain_cfg(dir / "full.bgz", dir / "out.fq"));
    CHECK(!r.threw);
    CHECK_EQ(r.rep.stats.reads_in, uint64_t(300));
  }
  std::size_t off = 0;
  int boundaries = 0;
  int accepted = 0;
  while (off + 18 <= full.size()) {
    const std::size_t bsize = (static_cast<unsigned char>(full[off + 16]) |
                               (static_cast<unsigned char>(full[off + 17]) << 8)) + 1u;
    off += bsize;
    if (off >= full.size()) break;   // the complete file
    write_bytes(dir / "cut.bgz", full.substr(0, off));
    fq::Config cfg = plain_cfg(dir / "cut.bgz", dir / "out.fq");
    cfg.qc_only = true;
    accepted += try_run(cfg).threw ? 0 : 1;
    ++boundaries;
  }
  CHECK(boundaries > 30);
  CHECK_EQ(accepted, 0);
  std::filesystem::remove_all(dir);
}

}  // namespace

int main() {
  std::cerr << "fastq_stream tests [" << fq::simd::backend() << "]\n";
  try {
    test_simd();
    test_find_seed4();
    test_ring_basic();
    test_ring_threaded();
    test_quality_trimming();
    test_adapter_trimming();
    test_adapter_trim_matches_reference_after_find_seed4();
    test_record_parsing();
    test_qc_stats();
    test_qhist_is_exact_across_every_read_length_and_fold_point();
    test_pipeline_end_to_end();
    test_malformed_input();
    test_adapter_short_read_overread();
    test_mean_phred_underflow();
    test_sliding_window_threshold_domain();
    test_truncated_compressed_input();
    test_output_aliasing_input();
    test_final_record_without_newline();
    test_non_printable_bytes_rejected();
    test_q20_q30_exact_over_all_bytes();
    test_quality_trim_accounting();
    test_json_escapes_control_characters();
    test_phred64_is_converted_and_detected();
    test_phred64_detector_spares_high_accuracy_phred33();
    test_header_and_plus_lines_are_validated();
    test_record_limit_does_not_depend_on_position();
    test_reads_over_the_cycle_limit_are_counted();
    test_crlf_and_mixed_line_endings_match_lf();
    test_accumulators_are_exact_over_many_chunks();
    test_bgzf_cut_at_every_member_boundary_fails();
  } catch (const std::exception& e) {
    std::cerr << "EXCEPTION: " << e.what() << "\n";
    return 1;
  }
  std::cerr << (g_failures == 0 ? "PASS " : "FAIL ") << (g_checks - g_failures) << "/"
            << g_checks << " checks\n";
  return g_failures == 0 ? 0 : 1;
}
