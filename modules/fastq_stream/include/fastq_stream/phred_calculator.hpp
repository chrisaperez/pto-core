// SPDX-License-Identifier: MIT
// fastq_stream — single-pass QC accumulation and quality trimming.
//
// All statistics are computed in one pass over each read while it is still hot
// in L1 from the decompressor. Per-cycle sums accumulate into a 32-bit scratch
// that is folded into the 64-bit totals once per chunk, which keeps the inner
// loop in 32-bit SIMD lanes without any overflow risk. The bound is the
// SMALLEST record, not a typical one: a chunk is at most 320 KiB with its
// carried prefix, a one-base record is 8 bytes, so a chunk holds at most
// ~43,000 non-empty records and a per-cycle quality sum stays below
// 43,000 x 126 = 5.4e6. (This comment used to say "~650 reads", which is a
// 150 bp chunk, not the worst case.) docs/AUDIT_2026-09-11_fastq_stream.md
// checked every total exactly over 32M 2-base records.
#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>

#include "fastq_stream/simd.hpp"

namespace fq {

inline constexpr std::size_t kMaxCycles = 1024;
// The histogram keeps one bin per printable Phred value, Q0-Q93, so the highest
// quality present is exact -- looks_like_phred64 needs to tell Phred+64's Q72
// ceiling from HiFi's Q93, and a clamped top bin cannot. The JSON report still
// prints kReportQualBins, folding Q63+ into the last one; that schema is what
// scripts and pto-cloud's merge read.
inline constexpr std::size_t kQualBins = 94;
inline constexpr std::size_t kReportQualBins = 64;
inline constexpr uint8_t kPhredOffset = 33;  // Sanger / Illumina 1.8+

enum BaseCode : uint8_t { kA = 0, kC = 1, kG = 2, kT = 3, kN = 4 };

inline constexpr std::array<uint8_t, 256> make_base_table() {
  std::array<uint8_t, 256> t{};
  for (auto& v : t) v = kN;
  t['A'] = t['a'] = kA;
  t['C'] = t['c'] = kC;
  t['G'] = t['g'] = kG;
  t['T'] = t['t'] = kT;
  return t;
}
inline constexpr std::array<uint8_t, 256> kBaseTable = make_base_table();

// Per-chunk 32-bit accumulators, folded into QcStats at chunk end.
struct CycleScratch {
  std::array<uint32_t, kMaxCycles> qsum{};
  std::array<uint32_t, kMaxCycles> count{};
  std::array<std::array<uint32_t, 5>, kMaxCycles> base{};
  // Four independent quality histograms, indexed round-robin on the base's
  // position within its read (i & 3) rather than one shared array. Illumina
  // quality is commonly binned to ~4 distinct values, so a run of bases at
  // the same quality is the common case, not a corner case -- incrementing
  // one shared bin makes every such increment wait on the store from the one
  // before it (store-to-load forwarding), while four independent arrays give
  // the CPU four separate dependency chains to overlap. Summed back into one
  // 94-bin histogram in fold() below; the total per bin is the same integer
  // either way, exactly, since this only changes which of four counters an
  // increment lands in, never how many increments happen or what they sum to.
  std::array<std::array<uint32_t, kQualBins>, 4> qh{};
  std::size_t max_cycle = 0;

  void clear() noexcept {
    if (max_cycle == 0) return;
    const std::size_t n = std::min(max_cycle, kMaxCycles);
    std::memset(qsum.data(), 0, n * sizeof(uint32_t));
    std::memset(count.data(), 0, n * sizeof(uint32_t));
    std::memset(base.data(), 0, n * sizeof(base[0]));
    // qh has no per-cycle structure (it is a histogram over quality VALUES,
    // not positions), so unlike the three arrays above -- which only clear
    // the cycles a read could actually have touched -- every one of its
    // 4 x kQualBins entries must be cleared whenever anything was
    // accumulated this chunk, not just the first `n` of them.
    std::memset(qh.data(), 0, sizeof(qh));
    max_cycle = 0;
  }
};

struct QcStats {
  uint64_t reads_in = 0;
  uint64_t bases_in = 0;
  uint64_t reads_out = 0;
  uint64_t bases_out = 0;
  uint64_t qsum_in = 0;  // Phred sum, offset already removed
  uint64_t q20 = 0;
  uint64_t q30 = 0;
  uint64_t gc = 0;
  uint64_t n_bases = 0;
  uint64_t adapter_trimmed = 0;
  uint64_t quality_trimmed = 0;
  uint64_t dropped_short = 0;
  uint64_t dropped_n = 0;
  // Reads longer than kMaxCycles. Their bases past the table count in every
  // total but in no per-cycle array, and before this counter the JSON report
  // gave no sign of it: a 2 kb run reported "cycles": 1024 over half its bases.
  uint64_t reads_over_cycle_limit = 0;

  std::array<uint64_t, kQualBins> qhist{};
  std::array<uint64_t, kMaxCycles> pos_qsum{};
  std::array<uint64_t, kMaxCycles> pos_count{};
  std::array<std::array<uint64_t, 5>, kMaxCycles> pos_base{};
  std::size_t max_cycle = 0;

  void fold(CycleScratch& s) noexcept {
    const std::size_t n = std::min(s.max_cycle, kMaxCycles);
    for (std::size_t i = 0; i < n; ++i) {
      pos_qsum[i] += s.qsum[i];
      pos_count[i] += s.count[i];
      for (int b = 0; b < 5; ++b) pos_base[i][b] += s.base[i][b];
    }
    // The four privatised sub-histograms sum back into one: this is plain
    // integer addition, so the total in each bin is exactly what a single
    // shared histogram would have accumulated, regardless of which of the
    // four counters any given base landed in.
    for (int p = 0; p < 4; ++p) {
      for (std::size_t q = 0; q < kQualBins; ++q) qhist[q] += s.qh[p][q];
    }
    max_cycle = std::max(max_cycle, n);
    s.clear();
  }

  void merge(const QcStats& o) noexcept {
    reads_in += o.reads_in;
    bases_in += o.bases_in;
    reads_out += o.reads_out;
    bases_out += o.bases_out;
    qsum_in += o.qsum_in;
    q20 += o.q20;
    q30 += o.q30;
    gc += o.gc;
    n_bases += o.n_bases;
    adapter_trimmed += o.adapter_trimmed;
    quality_trimmed += o.quality_trimmed;
    dropped_short += o.dropped_short;
    dropped_n += o.dropped_n;
    reads_over_cycle_limit += o.reads_over_cycle_limit;
    for (std::size_t i = 0; i < kQualBins; ++i) qhist[i] += o.qhist[i];
    const std::size_t n = std::min(o.max_cycle, kMaxCycles);
    for (std::size_t i = 0; i < n; ++i) {
      pos_qsum[i] += o.pos_qsum[i];
      pos_count[i] += o.pos_count[i];
      for (int b = 0; b < 5; ++b) pos_base[i][b] += o.pos_base[i][b];
    }
    max_cycle = std::max(max_cycle, o.max_cycle);
  }

  double mean_quality() const noexcept {
    return bases_in ? static_cast<double>(qsum_in) / static_cast<double>(bases_in) : 0.0;
  }
  double q30_rate() const noexcept {
    return bases_in ? static_cast<double>(q30) / static_cast<double>(bases_in) : 0.0;
  }
  double gc_rate() const noexcept {
    return bases_in ? static_cast<double>(gc) / static_cast<double>(bases_in) : 0.0;
  }
};

// Byte class for SEQ: bit 0 = G or C, bit 1 = not ACGT (counted as N). One
// table load answers both totals the QC pass keeps.
inline constexpr std::array<uint8_t, 256> make_class_table() {
  std::array<uint8_t, 256> t{};
  for (std::size_t i = 0; i < 256; ++i) {
    const uint8_t code = kBaseTable[i];
    t[i] = static_cast<uint8_t>((code == kC || code == kG ? 1u : 0u) | (code == kN ? 2u : 0u));
  }
  return t;
}
inline constexpr std::array<uint8_t, 256> kClassTable = make_class_table();

// What accumulate_read learned about a read beyond the statistics it folded in.
struct ReadScan {
  std::size_t n_bases = 0;  // non-ACGT bases over the whole, untrimmed read
  bool valid = true;        // every SEQ and QUAL byte is printable ASCII, '!'..'~'
};

// Accumulate the untrimmed read into stats + per-cycle scratch.
//
// Also reports whether every SEQ/QUAL byte is printable. The statistics below
// still clamp a sub-offset quality byte to Q0 so that a library caller which
// ignores `valid` gets a bounded answer, but the pipeline rejects such reads:
// a NUL is corruption, not a low-quality base.
inline ReadScan accumulate_read(QcStats& st, CycleScratch& sc, const uint8_t* seq,
                                const uint8_t* qual, std::size_t len) noexcept {
  const std::size_t n = std::min(len, kMaxCycles);
  st.reads_over_cycle_limit += (len > kMaxCycles);
  simd::add_widen_u8_to_u32(qual, n, sc.qsum.data());
  for (std::size_t i = 0; i < n; ++i) {
    sc.count[i] += 1;
    const uint8_t code = kBaseTable[seq[i]];
    sc.base[i][code] += 1;
  }
  sc.max_cycle = std::max(sc.max_cycle, n);

  // qsum_in accumulates here, per base, from the same clamped value the
  // histogram uses.
  //
  // It was previously computed as
  //     qsum_in += simd::sum_u8(qual, len) - len * kPhredOffset;
  // which underflows whenever a quality byte falls below the +33 offset --
  // legal to encounter in a corrupt or mis-encoded file, and not rejected by
  // the parser. Four bytes of value 1..4 give sum_u8 = 10 against len*33 = 132,
  // and the uint64_t subtraction wraps to ~1.8e19; the reported mean quality
  // was Q4611686018427387904. The clamp three lines below already guarded the
  // histogram against exactly this input, so only the aggregate was exposed.
  //
  // Folding the sum into the loop that already walks every base costs one add
  // per base and removes the separate SIMD pass, so this is not a slower path.
  //
  //
  // The totals live in locals and are stored once. `seq`/`qual` are uint8_t
  // pointers and uint8_t may alias anything, so a `st.gc += ...` in the loop
  // obliges a load and store per base. That is hygiene rather than a measured
  // win: end to end it was within noise on Apple M4.
  //
  // What stays OUT of this loop IS measured. Q20/Q30 and the printable-byte
  // check were first folded in as per-base compares, on the theory that one
  // pass beats three. With one worker on 425 MB that ran 0.57 s against 0.49 s
  // for the vector passes, which cover 16 or 32 bytes per compare. Anything
  // that can be a whole-read SIMD reduction stays one.
  const bool valid = simd::all_printable(qual, len) && simd::all_printable(seq, len);
  st.q20 += simd::count_ge(qual, len, kPhredOffset + 20);
  st.q30 += simd::count_ge(qual, len, kPhredOffset + 30);

  uint64_t qsum = 0, gc = 0, ns = 0;
  for (std::size_t i = 0; i < len; ++i) {
    const uint8_t qb = qual[i];
    const unsigned q = qb >= kPhredOffset ? static_cast<unsigned>(qb - kPhredOffset) : 0u;
    qsum += q;
    // i & 3 rather than a shared bin: see CycleScratch::qh. Folded into
    // st.qhist once per chunk in QcStats::fold, not touched directly here.
    sc.qh[i & 3u][std::min<unsigned>(q, kQualBins - 1)] += 1;
    const uint8_t cls = kClassTable[seq[i]];
    gc += cls & 1u;
    ns += cls >> 1;
  }
  st.reads_in += 1;
  st.bases_in += len;
  st.qsum_in += qsum;
  st.gc += gc;
  st.n_bases += ns;
  return ReadScan{static_cast<std::size_t>(ns), valid};
}

// 5' trim: advance past leading bases below min_q.
inline std::size_t trim_front(const uint8_t* qual, std::size_t len,
                              uint8_t min_q) noexcept {
  const uint8_t thr = static_cast<uint8_t>(kPhredOffset + min_q);
  std::size_t i = 0;
  while (i < len && qual[i] < thr) ++i;
  return i;
}

// 3' trim: drop trailing bases below min_q.
inline std::size_t trim_back(const uint8_t* qual, std::size_t len,
                             uint8_t min_q) noexcept {
  const uint8_t thr = static_cast<uint8_t>(kPhredOffset + min_q);
  while (len > 0 && qual[len - 1] < thr) --len;
  return len;
}

// Trimmomatic-style SLIDINGWINDOW: scan 5'->3', truncate at the start of the
// first window whose mean Phred falls below min_mean. Windows are evaluated
// with a rolling sum, so the whole scan is O(len) with no division per base.
inline std::size_t trim_sliding_window(const uint8_t* qual, std::size_t len,
                                       std::size_t window,
                                       double min_mean) noexcept {
  if (window == 0 || len < window) return len;

  // The threshold is clamped into uint32_t range BEFORE the narrowing cast.
  //
  // `static_cast<uint32_t>(d)` is undefined behaviour whenever the truncated
  // value of `d` is not representable -- which covers every negative value,
  // every value past 2^32, and NaN. All three are reachable straight from the
  // CLI (`--window-mean -5`, `--window-mean 1e30`, `--window-mean nan`): the
  // option was parsed with strtod and handed here unvalidated, and UBSan
  // reports "outside the range of representable values of type 'unsigned int'"
  // on each. main.cpp now rejects those arguments, but the guard belongs here
  // too -- this is a public header, the kernel is callable directly by the
  // tests and by any library consumer, and a validated CLI is not a property
  // the function itself can rely on.
  //
  // NaN is handled first and explicitly: every comparison against NaN is
  // false, so a `min_mean <= 0` test alone would fall through to the cast.
  // A NaN threshold means "no meaningful bound", which is the same as
  // disabling the window, so the whole read is kept.
  if (!(min_mean == min_mean)) return len;  // NaN
  const double want = min_mean * static_cast<double>(window);
  if (want <= 0.0) return len;              // no window can fall below it

  // The rolling sum and the threshold are compared in 64-bit.
  //
  // `sum` is a raw ASCII sum, at most window*255, and `offset + limit` was
  // previously computed in uint32_t. Saturating `limit` at UINT32_MAX would
  // then have wrapped that addition around to `offset - 1` and inverted the
  // test on exactly the extreme input the saturation exists to handle. In
  // uint64_t neither term can overflow: window is bounded by the read length,
  // and the ceiling below keeps `limit` inside 32 bits.
  constexpr double kLimitCeiling = 4294967296.0;  // 2^32
  const std::uint64_t limit = want >= kLimitCeiling
                                  ? 0xFFFFFFFFull
                                  : static_cast<std::uint64_t>(want);
  const std::uint64_t offset = static_cast<std::uint64_t>(window) * kPhredOffset;
  const std::uint64_t threshold = offset + limit;

  std::uint64_t sum = 0;
  for (std::size_t i = 0; i < window; ++i) sum += qual[i];
  if (sum < threshold) return 0;
  for (std::size_t s = 1; s + window <= len; ++s) {
    sum += qual[s + window - 1];
    sum -= qual[s - 1];
    if (sum < threshold) return s;
  }
  return len;
}

// A header or '+' line: printable ASCII, plus TAB (Illumina and 10x put
// comments after a space or a tab).
//
// These lines were never checked. F7 (docs/TORTURE_2026-09-10.md) closed NUL,
// control and high-bit bytes in SEQ and QUAL and stopped there, so a read name
// carrying NUL, ESC or BEL went into the output FASTQ at exit 0 -- where a NUL
// truncates the name in every C string API an aligner uses, so two different
// reads can come out under one name. The sibling-function recurrence
// pto-core/CLAUDE.md warns about, in the sibling line.
inline bool valid_header_line(const uint8_t* p, std::size_t n) noexcept {
  bool bad = false;
  for (std::size_t i = 0; i < n; ++i) {
    const uint8_t c = p[i];
    bad |= (c != '\t') && static_cast<uint8_t>(c - 0x20u) > 0x5eu;
  }
  return !bad;
}

// Whether quality bytes read at offset 33 can only have been Phred+64.
//
// Phred+64 spans '@'..'i' -- Q0-Q41 on Illumina 1.3-1.8, with 'B' as 1.5's
// Q2 marker. Read at offset 33 that is Q31-Q72, and three things follow:
//
//   1. nothing below Q31 (there is no byte under '@');
//   2. nothing above Q72 (there is no byte over 'i');
//   3. a real share of bases at Q61 or higher -- Phred+64 Q30+, which every
//      Illumina 1.3-1.7 run of usable quality has in bulk. This test
//      requires at least 1%.
//
// The first version checked only (1) and "something above Q41", and that
// flagged exactly the modern Phred+33 data it must not: AVITI's Q44-Q50,
// PacBio HiFi's Q93, whenever a file happened to have no base under Q31.
// Measured on 3,000-read fixtures, both exited 2. (2) rules out HiFi, whose
// ceiling is Q93, and (3) rules out every short-read platform, none of which
// writes Q61 -- AVITI tops out near Q50 -- and a quality-capped long-read
// subset up to Q60. The line was first drawn at Q55 and a Q42-Q60 HiFi subset
// still tripped it (31% of its bases at Q55+); Q61 clears it.
//
// What can still misfire is Phred+33 data confined to Q31-Q72 with 1% or more
// at Q61+ and nothing past Q72. And a Phred+64 file with no base at Q30 or
// above (+64) is NOT detected. Both
// are stated in docs/AUDIT_2026-09-11_fastq_stream.md; --phred-offset settles
// either explicitly.
//
// Decided from the histogram alone, which now has a bin for every printable
// Phred value, so it costs nothing per base. (A running maximum in the QC loop
// was tried first and cost 4% on a single-worker run.)
inline constexpr unsigned kPhred64MinAt33 = 31;   // '@' read at +33
inline constexpr unsigned kPhred64MaxAt33 = 72;   // 'i' read at +33
inline constexpr unsigned kPhred64BulkAt33 = 61;  // Phred+64 Q30

inline bool looks_like_phred64(const QcStats& s) noexcept {
  if (s.bases_in == 0) return false;
  uint64_t below = 0;
  uint64_t bulk = 0;
  uint64_t above = 0;
  for (std::size_t q = 0; q < kPhred64MinAt33; ++q) below += s.qhist[q];
  for (std::size_t q = kPhred64BulkAt33; q <= kPhred64MaxAt33; ++q) bulk += s.qhist[q];
  for (std::size_t q = kPhred64MaxAt33 + 1; q < kQualBins; ++q) above += s.qhist[q];
  if (above > 0) return false;
  // bulk * 100 cannot overflow: bulk <= bases_in, far below 2^57.
  return below == 0 && bulk * 100 >= s.bases_in;
}

inline std::size_t count_n(const uint8_t* seq, std::size_t len) noexcept {
  std::size_t n = 0;
  for (std::size_t i = 0; i < len; ++i) n += (kBaseTable[seq[i]] == kN);
  return n;
}

inline double mean_phred(const uint8_t* qual, std::size_t len) noexcept {
  if (len == 0) return 0.0;
  // Clamp per byte, exactly as accumulate_read does. The previous form,
  //     simd::sum_u8(qual, len) - len * kPhredOffset
  // is the same unsigned underflow fixed in accumulate_read (AUDIT.md S3): any
  // quality byte below the +33 offset wraps the uint64_t subtraction to ~1.8e19.
  //
  // Here that is worse than a bad statistic. This function gates the
  // --min-mean-q drop in pipeline.cpp, and a wrapped mean is enormous, so the
  // `< min_mean_q` test fails and the read is KEPT. The filter inverted on
  // exactly the corrupt input it exists to reject.
  uint64_t s = 0;
  for (std::size_t i = 0; i < len; ++i) {
    s += qual[i] >= kPhredOffset ? static_cast<uint64_t>(qual[i] - kPhredOffset) : 0u;
  }
  return static_cast<double>(s) / static_cast<double>(len);
}

}  // namespace fq
