// SPDX-License-Identifier: MIT
// fastq_stream — component microbenchmarks.
//
// These measure the pieces the manuscript makes claims about, in isolation:
// SPSC hand-off latency, SIMD kernel throughput, and end-to-end pipeline
// throughput on a synthetic dataset. Whole-tool comparisons against fastp and
// Trimmomatic live in scripts/run_benchmarks.py; this binary exists so a
// reviewer can reproduce the component numbers without downloading 50 GB.
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include "fastq_stream/adapter_trimmer.hpp"
#include "fastq_stream/lockfree_queue.hpp"
#include "fastq_stream/phred_calculator.hpp"
#include "fastq_stream/pipeline.hpp"
#include "fastq_stream/simd.hpp"

namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point t0) {
  return std::chrono::duration<double>(Clock::now() - t0).count();
}

void row(const char* name, const char* unit, double value, const char* note = "") {
  std::printf("  %-34s %12.2f %-12s %s\n", name, value, unit, note);
}

// Keeps the optimiser from deleting the kernel under test.
template <typename T>
void sink(T&& v) {
  // Written *and* read through the volatile, so the store cannot be elided and
  // the variable is not "set but never used" (-Wunused-but-set-variable, which
  // GCC raises on the write-only form).
  static volatile std::uint64_t consume;
  consume = static_cast<std::uint64_t>(v);
  (void)consume;
}

// Single-shot microbenchmarks on a laptop are not reproducible: on a hybrid
// CPU the same kernel measured 84 GB/s and 36 GB/s in consecutive runs purely
// from P-core/E-core placement and clock drift. Report the best of several
// trials, which is the least noisy estimator of the kernel's actual cost, and
// print the spread so a reader can see how much to trust it.
constexpr int kTrials = 5;

// `fn` runs the workload once and returns elapsed seconds.
template <typename F>
void bench_best(const char* name, double gb, F&& fn) {
  double best = 1e300, worst = 0.0;
  for (int i = 0; i < kTrials; ++i) {
    const double s = fn();
    best = std::min(best, s);
    worst = std::max(worst, s);
  }
  char note[64];
  std::snprintf(note, sizeof(note), "best of %d, %.0f%% spread", kTrials,
                100.0 * (worst - best) / best);
  row(name, "GB/s", gb / best, note);
}

// --------------------------------------------------------------- SIMD kernels

void bench_kernels() {
  std::printf("\nSIMD kernels [%s]\n", fq::simd::backend());

  constexpr std::size_t kLen = 150;              // one Illumina read
  constexpr std::size_t kReads = 2'000'000;
  std::mt19937 rng(7);
  std::uniform_int_distribution<int> q(35, 74);
  std::vector<uint8_t> qual(kLen * 64);
  for (auto& v : qual) v = static_cast<uint8_t>(q(rng));

  const double gb = static_cast<double>(kLen) * kReads / 1e9;

  bench_best("sum_u8 (mean Phred)", gb, [&] {
    const auto t0 = Clock::now();
    uint64_t acc = 0;
    for (std::size_t i = 0; i < kReads; ++i) {
      acc += fq::simd::sum_u8(qual.data() + (i % 64) * kLen, kLen);
    }
    sink(acc);
    return seconds_since(t0);
  });

  bench_best("count_ge (Q30 fraction)", gb, [&] {
    const auto t0 = Clock::now();
    uint64_t acc = 0;
    for (std::size_t i = 0; i < kReads; ++i) {
      acc += fq::simd::count_ge(qual.data() + (i % 64) * kLen, kLen, 33 + 30);
    }
    sink(acc);
    return seconds_since(t0);
  });

  bench_best("add_widen (per-cycle quality)", gb, [&] {
    std::vector<uint32_t> pos(kLen, 0);
    const auto t0 = Clock::now();
    for (std::size_t i = 0; i < kReads; ++i) {
      fq::simd::add_widen_u8_to_u32(qual.data() + (i % 64) * kLen, kLen, pos.data());
    }
    sink(pos[0]);
    return seconds_since(t0);
  });

  bench_best("trim_sliding_window", gb, [&] {
    const auto t0 = Clock::now();
    uint64_t acc = 0;
    for (std::size_t i = 0; i < kReads; ++i) {
      acc += fq::trim_sliding_window(qual.data() + (i % 64) * kLen, kLen, 4, 20.0);
    }
    sink(acc);
    return seconds_since(t0);
  });
}

void bench_adapter() {
  std::printf("\nAdapter matching\n");
  constexpr std::size_t kLen = 150;
  constexpr std::size_t kReads = 1'000'000;

  std::mt19937 rng(11);
  std::uniform_int_distribution<int> pick(0, 3);
  const char bases[] = "ACGT";
  std::vector<std::string> reads(64);
  for (std::size_t i = 0; i < reads.size(); ++i) {
    std::string s;
    s.reserve(kLen);
    for (std::size_t j = 0; j < kLen; ++j) s.push_back(bases[pick(rng)]);
    if (i % 5 == 0) s.replace(110, std::string(fq::kTruSeqAdapter).size(),
                              std::string(fq::kTruSeqAdapter));
    reads[i] = std::move(s);
  }

  fq::AdapterTrimmer t;
  t.add_defaults();
  const double gb = static_cast<double>(kLen) * kReads / 1e9;
  bench_best("seed+extend, 2 adapters", gb, [&] {
    const auto t0 = Clock::now();
    uint64_t acc = 0;
    for (std::size_t i = 0; i < kReads; ++i) {
      const std::string& r = reads[i % reads.size()];
      acc += t.trim(reinterpret_cast<const uint8_t*>(r.data()), r.size());
    }
    sink(acc);
    return seconds_since(t0);
  });
}

// ------------------------------------------------------------------ SPSC ring

void bench_ring() {
  std::printf("\nLock-free SPSC ring\n");
  constexpr int kN = 20'000'000;
  fq::SpscRing<int, 1024> ring;

  const auto t0 = Clock::now();
  std::thread consumer([&] {
    int got = 0;
    fq::Backoff bo;
    while (got < kN) {
      if (ring.try_pop()) {
        ++got;
        bo.reset();
      } else {
        bo.pause();
      }
    }
  });
  {
    fq::Backoff bo;
    for (int i = 0; i < kN; ++i) {
      while (!ring.try_push(int(i))) bo.pause();
      bo.reset();
    }
  }
  consumer.join();
  const double secs = seconds_since(t0);
  row("throughput", "M items/s", kN / secs / 1e6);
  row("mean hand-off latency", "ns", secs / kN * 1e9);
}

// ------------------------------------------------------------- full pipeline

// Quality must follow a realistic NovaSeq profile, not a uniform draw. With
// uniformly random Q2..Q41 the mean is ~Q21, so a Q20 sliding window discards
// ~96% of reads and the pipeline measures a degenerate no-output workload.
// Here Q starts near 37 and decays quadratically toward ~30 by cycle 150, with
// 5% of reads given a genuinely bad 3' tail — which keeps the trim/filter and
// output-write paths on the critical path where they belong.
std::string synth_fastq(std::size_t n_reads) {
  std::mt19937 rng(3);
  std::uniform_int_distribution<int> pick(0, 3);
  std::normal_distribution<double> noise(0.0, 2.0);
  const char bases[] = "ACGT";
  constexpr int kLen = 150;

  std::string out;
  out.reserve(n_reads * 340);
  for (std::size_t i = 0; i < n_reads; ++i) {
    std::string seq(kLen, 'A'), qual(kLen, 'I');
    const bool bad_tail = (i % 20 == 0);
    for (int j = 0; j < kLen; ++j) {
      seq[static_cast<std::size_t>(j)] = bases[pick(rng)];
      const double frac = static_cast<double>(j) / kLen;
      double q = 37.0 - 7.0 * frac * frac + noise(rng);
      if (bad_tail && j > 120) q -= 25.0;
      const int qi = std::clamp(static_cast<int>(q), 2, 41);
      qual[static_cast<std::size_t>(j)] = static_cast<char>(33 + qi);
    }
    if (i % 5 == 0) {
      seq.replace(110, std::string(fq::kTruSeqAdapter).size(),
                  std::string(fq::kTruSeqAdapter));
    }
    out += "@SIM:1:FC:1:1:" + std::to_string(i) + " 1:N:0:ATCGATCG\n";
    out += seq;
    out += "\n+\n";
    out += qual;
    out += "\n";
  }
  return out;
}

void bench_pipeline(std::size_t n_reads) {
  std::printf("\nEnd-to-end pipeline (%zu reads, 150 bp, uncompressed input)\n", n_reads);
  const auto dir = std::filesystem::temp_directory_path() / "fastq_stream_bench";
  std::filesystem::create_directories(dir);
  const auto in = dir / "bench.fq";
  const auto out = dir / "bench.out.fq";

  const std::string data = synth_fastq(n_reads);
  { std::ofstream f(in, std::ios::binary); f << data; }

  // Rows are labelled with the QC/trim worker count actually used, not the
  // --threads value asked for. The pipeline always also runs a reader, an
  // assembler and a writer, and it floors the worker count at 1, so a request
  // for 1, 2 or 4 threads yields the same single worker. Labelling those rows
  // "1/2/4 threads" would report three different numbers for one configuration.
  const int hw = static_cast<int>(std::thread::hardware_concurrency());
  for (int t : {4, 6, 8, hw}) {
    if (t <= 0) continue;
    fq::Config cfg;
    cfg.input = in.string();
    cfg.output = out.string();
    cfg.threads = t;
    cfg.window = 4;
    cfg.window_mean = 20.0;

    // Best of 3, for the same reason the kernels are: a single run of this
    // configuration varied by ~2x across invocations on a hybrid-core host.
    fq::RunReport rep = fq::run_pipeline(cfg);
    for (int trial = 1; trial < 3; ++trial) {
      const fq::RunReport r = fq::run_pipeline(cfg);
      if (r.wall_seconds < rep.wall_seconds) rep = r;
    }
    const double mb = static_cast<double>(data.size()) / (1024.0 * 1024.0);
    char label[64];
    std::snprintf(label, sizeof(label), "--threads %d (%d worker%s)", t,
                  rep.worker_threads, rep.worker_threads == 1 ? "" : "s");
    char note[96];
    std::snprintf(note, sizeof(note), "%.1f MB pool, %llu/%llu reads kept",
                  static_cast<double>(rep.pool_bytes) / (1024.0 * 1024.0),
                  static_cast<unsigned long long>(rep.stats.reads_out),
                  static_cast<unsigned long long>(rep.stats.reads_in));
    row(label, "MB/s", mb / rep.wall_seconds, note);
  }
  std::filesystem::remove_all(dir);
}

}  // namespace

int main(int argc, char** argv) {
  std::size_t reads = 2'000'000;
  if (argc > 1) reads = static_cast<std::size_t>(std::strtoull(argv[1], nullptr, 10));

  std::printf("fastq_stream benchmarks — %u hardware threads\n",
              std::thread::hardware_concurrency());
  bench_kernels();
  bench_adapter();
  bench_ring();
  bench_pipeline(reads);
  std::printf("\n");
  return 0;
}
