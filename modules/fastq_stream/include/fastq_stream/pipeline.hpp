// SPDX-License-Identifier: MIT
// fastq_stream — the streaming pipeline.
//
// Topology (BGZF input; the plain-gzip path drops the inflater rank):
//
//   reader ──▶ [SPSC]×D ──▶ inflaters ──▶ [SPSC]×D ──▶ assembler
//                                                          │
//                                                     [SPSC]×W
//                                                          ▼
//                                       writer ◀── [SPSC]×W ◀── workers
//
// Every arrow is a genuine single-producer / single-consumer ring: the reader
// owns one lane per inflater, each inflater owns exactly one output lane, the
// assembler owns one lane per worker, and each worker owns one output lane.
// No ring is ever touched by more than two threads, so the SPSC claim in the
// paper is literal and not a shared-queue approximation.
//
// Ordering falls out of the topology for free. Chunk k is dispatched to lane
// (k mod D), each lane preserves order internally, so a collector that reads
// lanes round-robin in the same order observes the original stream order with
// no reorder buffer and no sequence-number sorting.
#pragma once

#include <atomic>
#include <chrono>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "fastq_stream/adapter_trimmer.hpp"
#include "fastq_stream/buffer.hpp"
#include "fastq_stream/lockfree_queue.hpp"
#include "fastq_stream/phred_calculator.hpp"
#include "fastq_stream/reader.hpp"
#include "fastq_stream/record.hpp"

#if defined(FQ_ENABLE_HTS)
#include "fastq_stream/hts_input.hpp"
#endif

namespace fq {

inline constexpr std::size_t kRingSlots = 16;  // power of two, per lane

inline void cpu_relax() noexcept {
#if defined(__x86_64__) || defined(__i386__)
  __builtin_ia32_pause();
#elif defined(__aarch64__)
  __asm__ __volatile__("yield" ::: "memory");
#endif
}

// Bounded spin, then a few yields, then short sleeps. Keeps hand-off latency
// low when the pipeline is saturated without burning a core when it is starved.
//
// The sleep stage is what makes the second half of that sentence true. The
// previous form spun 128 times and then called yield() forever, and
// sched_yield() returns immediately when nothing else is runnable -- so every
// stage waiting on a slower one was a busy loop in the kernel. Measured on a
// 450 MB plain gzip at -t 12, where the serial inflate is the only real work:
// 1.74 s wall, 6.1 s user and 13.5 s SYSTEM time, i.e. ~11 cores held to do
// one core's work. On a shared node or inside any CPU-quota'd container that
// is the inflater's own CPU being taken from it.
//
// 100 us is short against the unit of work a lane carries (a 256 KiB chunk
// takes on the order of a millisecond to produce or consume), and a stage
// only reaches it after the spin and yield budget is exhausted, i.e. when its
// neighbour is genuinely behind.
class Backoff {
 public:
  void pause() noexcept {
    if (spins_ < kSpinLimit) {
      ++spins_;
      cpu_relax();
    } else if (spins_ < kSpinLimit + kYieldLimit) {
      ++spins_;
      std::this_thread::yield();
    } else {
      std::this_thread::sleep_for(std::chrono::microseconds(kSleepMicros));
    }
  }
  void reset() noexcept { spins_ = 0; }

 private:
  static constexpr unsigned kSpinLimit = 128;
  static constexpr unsigned kYieldLimit = 16;
  static constexpr unsigned kSleepMicros = 100;
  unsigned spins_ = 0;
};

// A ring plus the flag that tells the consumer no more items will ever arrive.
struct Lane {
  SpscRing<Buffer, kRingSlots> ring;
  std::atomic<bool> closed{false};
};

inline void lane_push(Lane& lane, Buffer&& buf, std::atomic<bool>& abort) {
  Backoff bo;
  while (!lane.ring.try_push(std::move(buf))) {
    if (abort.load(std::memory_order_relaxed)) return;
    bo.pause();
  }
}

// Returns false once the lane is drained *and* closed.
inline bool lane_pop(Lane& lane, Buffer& out, std::atomic<bool>& abort) {
  Backoff bo;
  for (;;) {
    if (auto item = lane.ring.try_pop()) {
      out = std::move(*item);
      return true;
    }
    if (lane.closed.load(std::memory_order_acquire)) {
      // Re-check: a producer may have pushed just before closing.
      if (auto item = lane.ring.try_pop()) {
        out = std::move(*item);
        return true;
      }
      return false;
    }
    if (abort.load(std::memory_order_relaxed)) return false;
    bo.pause();
  }
}

struct Config {
  std::string input = "-";
  std::string output = "-";
  int threads = 0;  // 0 => hardware_concurrency

  // Quality trimming
  uint8_t front_q = 0;      // 0 disables
  uint8_t tail_q = 0;       // 0 disables
  std::size_t window = 4;   // 0 disables sliding window
  double window_mean = 20.0;

  // Filters
  std::size_t min_len = 15;
  double max_n_rate = 0.1;   // <0 disables
  double min_mean_q = 0.0;   // 0 disables

  // Adapters
  bool trim_adapters = true;
  std::vector<std::string> adapters;  // empty => built-in defaults
  std::size_t adapter_seed = 4;
  std::size_t adapter_min_overlap = 4;
  double adapter_mismatch = 0.2;

  // BAM/CRAM ingestion. Set by main.cpp from looks_like_alignment(); the
  // hts_* fields are ignored unless it is true. Requires FQ_ENABLE_HTS.
  // `hts_input` and `reference_fasta` stay declared in every build so main.cpp
  // can detect an alignment file and emit "this binary was built without
  // FQ_ENABLE_HTS" rather than not recognising the flag at all.
  bool hts_input = false;
  std::string reference_fasta;   // required for CRAM, see hts_input.cpp
#if defined(FQ_ENABLE_HTS)
  HtsLayout hts_layout = HtsLayout::kSingle;
  char hts_missing_qual = 'I';
#endif

  // Quality encoding of FASTQ text input: 33 or 64. Phred+64 is validated
  // ('@'..'~') and converted to Phred+33 in place before anything reads it, so
  // every kernel, filter and report sees +33 and the output FASTQ is +33.
  int phred_offset = 33;

  bool qc_only = false;   // compute QC, write nothing
  bool make_fifo = false; // mkfifo(output) if it does not exist
  bool verbose = false;
};

struct RunReport {
  QcStats stats;
#if defined(FQ_ENABLE_HTS)
  // Non-zero only on the BAM/CRAM path. These counters are the only evidence a
  // user has that the conversion was faithful -- how many records were dropped
  // as secondary or supplementary, how many were reverse-complemented back to
  // instrument orientation -- so main.cpp prints them rather than treating them
  // as debug output.
  HtsStats hts_stats;
  bool from_alignment = false;
  const char* alignment_format = "";
#endif
  GzipFormat format = GzipFormat::kRaw;
  uint64_t compressed_bytes = 0;
  uint64_t bytes_written = 0;
  int inflater_threads = 0;
  int worker_threads = 0;
  std::size_t pool_blocks = 0;
  std::size_t pool_bytes = 0;
  double wall_seconds = 0.0;
};

// Runs the pipeline to completion. Throws std::runtime_error on fatal errors.
RunReport run_pipeline(const Config& cfg);

}  // namespace fq
