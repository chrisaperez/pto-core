// SPDX-License-Identifier: MIT
// fastq_stream — seed-and-extend adapter matching.
//
// A naive all-offsets scan costs O(read_len * adapter_len) per read, roughly an
// order of magnitude more work than everything else in the pipeline combined.
// Instead each adapter contributes an exact seed (its first `seed_len` bases);
// only offsets where the seed hits exactly are extended, and the extension uses
// the early-exit SIMD Hamming kernel.
//
// Sensitivity trade-off, stated plainly because it matters for the paper: a
// read whose adapter copy carries a sequencing error inside the seed window is
// missed. This is the same compromise fastp makes. Lowering --adapter-seed
// raises sensitivity and cost; the benchmark script reports both.
#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>
#include <vector>

#include "fastq_stream/simd.hpp"

namespace fq {

// Illumina TruSeq universal adapter prefix — the standard default.
inline constexpr std::string_view kTruSeqAdapter = "AGATCGGAAGAGC";
// Nextera / Illumina DNA Prep transposase mosaic end.
inline constexpr std::string_view kNexteraAdapter = "CTGTCTCTTATACACATCT";

class AdapterTrimmer {
 public:
  AdapterTrimmer() = default;

  void add(std::string_view adapter) {
    if (adapter.empty()) return;
    std::string upper(adapter);
    for (char& c : upper) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    adapters_.push_back(std::move(upper));
  }

  void add_defaults() {
    add(kTruSeqAdapter);
    add(kNexteraAdapter);
  }

  bool empty() const noexcept { return adapters_.empty(); }
  const std::vector<std::string>& adapters() const noexcept { return adapters_; }

  void set_min_overlap(std::size_t n) noexcept { min_overlap_ = std::max<std::size_t>(1, n); }
  void set_max_mismatch_rate(double r) noexcept {
    max_mismatch_rate_ = std::clamp(r, 0.0, 0.5);
  }
  void set_seed_len(std::size_t n) noexcept { seed_len_ = std::max<std::size_t>(1, n); }

  std::size_t min_overlap() const noexcept { return min_overlap_; }
  double max_mismatch_rate() const noexcept { return max_mismatch_rate_; }
  std::size_t seed_len() const noexcept { return seed_len_; }

  // Returns the trimmed read length: the 5'-most offset at which an adapter
  // match begins, or `len` when no adapter is found.
  std::size_t trim(const uint8_t* seq, std::size_t len) const noexcept {
    if (adapters_.empty() || len == 0) return len;
    // A read shorter than the minimum overlap cannot host a match, and the
    // offset loop below is inclusive of `limit` -- so clamping `limit` to 0 is
    // NOT enough to suppress it. Return before any seed load: at len < 4 the
    // 32-bit seed read would run past the end of the read.
    if (len < min_overlap_) return len;
    std::size_t best = len;

    for (const std::string& ad : adapters_) {
      const auto* a = reinterpret_cast<const uint8_t*>(ad.data());
      const std::size_t alen = ad.size();
      const std::size_t seed = std::min({seed_len_, min_overlap_, alen});

      // Only offsets that leave at least min_overlap bases can match, and we
      // never need to look past an offset we have already beaten. `len >=
      // min_overlap_` is guaranteed above, so this cannot underflow and the
      // inclusive bound is genuinely in range.
      const std::size_t limit = len - min_overlap_;
      // Both seed reads below take `seed` bytes at offset s <= limit. Since
      // seed <= min_overlap_, s + seed <= len holds for every iteration.
      assert(seed <= min_overlap_ && "seed must not exceed min_overlap");

      // Fast path: a 4-byte seed is compared as one unaligned 32-bit load.
      // The variable-length memcmp this replaces was not inlined and cost
      // ~200x more than every other kernel in the pipeline. Reading 4 bytes at
      // offset s is in bounds because s <= len - min_overlap and
      // min_overlap >= seed == 4.
      const bool word_seed = (seed == 4) && (min_overlap_ >= 4);
      uint32_t seed_word = 0;
      if (word_seed) std::memcpy(&seed_word, a, 4);

      // s is advanced two ways below: simd::find_seed4 jumps it straight to
      // the next 4-byte seed match (skipping non-matches in bulk instead of
      // testing one offset at a time -- on random DNA a 4-byte seed matches
      // about 1 in 256 offsets, so this is almost always the whole rest of
      // the read in one call), or, for the seed < 4 fallback, `++s` below
      // steps it by one exactly as the loop's own increment used to.
      // `best - 1` is safe: the outer condition (s < best) that keeps this
      // loop running guarantees best >= 1 whenever this line executes.
      for (std::size_t s = 0; s <= limit && s < best;) {
        if (word_seed) {
          s = simd::find_seed4(seq, s, std::min(limit, best - 1), seed_word);
          if (s == simd::kNotFound) break;  // no seed match left in range
        } else if (std::memcmp(seq + s, a, seed) != 0) {
          ++s;
          continue;
        }
        const std::size_t overlap = std::min(alen, len - s);
        if (overlap < min_overlap_) {
          ++s;
          continue;
        }
        const std::size_t cap =
            static_cast<std::size_t>(max_mismatch_rate_ * static_cast<double>(overlap));
        if (simd::mismatches(seq + s, a, overlap, cap) <= cap) {
          best = s;
          break;  // 5'-most hit for this adapter; nothing earlier can exist
        }
        ++s;
      }
    }
    return best;
  }

 private:
  std::vector<std::string> adapters_;
  std::size_t min_overlap_ = 4;
  std::size_t seed_len_ = 4;
  double max_mismatch_rate_ = 0.2;
};

}  // namespace fq
