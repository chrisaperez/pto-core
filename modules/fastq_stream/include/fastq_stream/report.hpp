// SPDX-License-Identifier: MIT
// fastq_stream — QC report emission (human summary + machine-readable JSON).
//
// The JSON schema is the one scripts/run_benchmarks.py parses and the one the
// manuscript's figures are generated from; keep field names stable.
#pragma once

#include <algorithm>
#include <array>
#include <cstdio>
#include <ostream>
#include <string>

#include "fastq_stream/pipeline.hpp"

namespace fq {

inline constexpr const char* kVersion = "0.1.0";

inline void write_text_report(std::ostream& os, const Config& cfg,
                              const RunReport& rep) {
  const QcStats& s = rep.stats;
  const double mb = static_cast<double>(rep.compressed_bytes) / (1024.0 * 1024.0);
  char line[512];

  os << "fastq_stream " << kVersion << "  [" << simd::backend() << "]\n";
#if defined(FQ_ENABLE_HTS)
  if (rep.from_alignment) {
    os << "  input        : " << cfg.input << "  (" << rep.alignment_format << ")\n";
  } else
#endif
  {
    os << "  input        : " << cfg.input << "  (" << format_name(rep.format) << ")\n";
  }
  os << "  threads      : " << rep.inflater_threads << " inflater + "
     << rep.worker_threads << " qc/trim (+reader, assembler, writer)\n";

  std::snprintf(line, sizeof(line),
                "  reads        : %llu in -> %llu out (%.2f%% kept)\n",
                static_cast<unsigned long long>(s.reads_in),
                static_cast<unsigned long long>(s.reads_out),
                s.reads_in ? 100.0 * static_cast<double>(s.reads_out) /
                                 static_cast<double>(s.reads_in)
                           : 0.0);
  os << line;
  std::snprintf(line, sizeof(line),
                "  bases        : %llu in -> %llu out (%.2f%% kept)\n",
                static_cast<unsigned long long>(s.bases_in),
                static_cast<unsigned long long>(s.bases_out),
                s.bases_in ? 100.0 * static_cast<double>(s.bases_out) /
                                 static_cast<double>(s.bases_in)
                           : 0.0);
  os << line;
  std::snprintf(line, sizeof(line),
                "  quality      : mean Q%.2f  Q20 %.2f%%  Q30 %.2f%%  GC %.2f%%\n",
                s.mean_quality(),
                s.bases_in ? 100.0 * static_cast<double>(s.q20) /
                                 static_cast<double>(s.bases_in)
                           : 0.0,
                100.0 * s.q30_rate(), 100.0 * s.gc_rate());
  os << line;
  std::snprintf(line, sizeof(line),
                "  trimmed      : %llu adapter, %llu quality\n",
                static_cast<unsigned long long>(s.adapter_trimmed),
                static_cast<unsigned long long>(s.quality_trimmed));
  os << line;
  std::snprintf(line, sizeof(line),
                "  dropped      : %llu too short / low quality, %llu excess N\n",
                static_cast<unsigned long long>(s.dropped_short),
                static_cast<unsigned long long>(s.dropped_n));
  os << line;
  if (s.reads_over_cycle_limit > 0) {
    std::snprintf(line, sizeof(line),
                  "  per-cycle    : first %zu cycles only; %llu reads are longer\n",
                  kMaxCycles, static_cast<unsigned long long>(s.reads_over_cycle_limit));
    os << line;
  }
#if defined(FQ_ENABLE_HTS)
  // Conversion diagnostics. These are the only evidence a user has that the
  // BAM -> FASTQ step was faithful -- how many records were dropped because
  // they could never be reads, and how many were reverse-complemented back to
  // instrument orientation -- so they are reported, not hidden behind --verbose.
  if (rep.from_alignment) {
    const HtsStats& h = rep.hts_stats;
    std::snprintf(line, sizeof(line),
                  "  alignment    : %llu records -> %llu reads "
                  "(%llu revcomp, %llu secondary, %llu supplementary, "
                  "%llu no SEQ, %llu synthesised QUAL)\n",
                  static_cast<unsigned long long>(h.records_read),
                  static_cast<unsigned long long>(h.emitted),
                  static_cast<unsigned long long>(h.reverse_complemented),
                  static_cast<unsigned long long>(h.skipped_secondary),
                  static_cast<unsigned long long>(h.skipped_supplementary),
                  static_cast<unsigned long long>(h.skipped_no_seq),
                  static_cast<unsigned long long>(h.synthesized_qual));
    os << line;
    if (h.pairs > 0 || h.singletons > 0) {
      std::snprintf(line, sizeof(line),
                    "  pairing      : %llu pairs, %llu singletons\n",
                    static_cast<unsigned long long>(h.pairs),
                    static_cast<unsigned long long>(h.singletons));
      os << line;
    }
  }
#endif
  std::snprintf(line, sizeof(line),
                "  performance  : %.2fs wall, %.1f MB/s compressed in, "
                "%.1f MB out, %.1f MB peak buffer pool\n",
                rep.wall_seconds, rep.wall_seconds > 0 ? mb / rep.wall_seconds : 0.0,
                static_cast<double>(rep.bytes_written) / (1024.0 * 1024.0),
                static_cast<double>(rep.pool_bytes) / (1024.0 * 1024.0));
  os << line;
}

namespace detail {

inline void json_u64_array(std::ostream& os, const uint64_t* v, std::size_t n) {
  os << '[';
  for (std::size_t i = 0; i < n; ++i) {
    if (i) os << ',';
    os << v[i];
  }
  os << ']';
}

inline std::string json_escape(const std::string& in) {
  std::string out;
  out.reserve(in.size() + 8);
  for (char c : in) {
    switch (c) {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n"; break;
      case '\t': out += "\\t"; break;
      default:
        // RFC 8259 forbids every unescaped control character, not only the
        // two above. A path is arbitrary bytes, so `-i $'a\rb.fq'` produced a
        // report json.load() rejected -- and the cloud worker parses this file.
        if (static_cast<unsigned char>(c) < 0x20) {
          char esc[8];
          std::snprintf(esc, sizeof(esc), "\\u%04x", static_cast<unsigned>(static_cast<unsigned char>(c)));
          out += esc;
        } else {
          out += c;
        }
    }
  }
  return out;
}

}  // namespace detail

inline void write_json_report(std::ostream& os, const Config& cfg,
                              const RunReport& rep) {
  const QcStats& s = rep.stats;
  const std::size_t cycles = s.max_cycle;
  char num[64];
  auto f = [&](double v) {
    std::snprintf(num, sizeof(num), "%.6g", v);
    return std::string(num);
  };

  os << "{\n";
  os << "  \"tool\": \"fastq_stream\",\n";
  os << "  \"version\": \"" << kVersion << "\",\n";
  os << "  \"simd_backend\": \"" << simd::backend() << "\",\n";
  os << "  \"input\": {\n";
  os << "    \"path\": \"" << detail::json_escape(cfg.input) << "\",\n";
#if defined(FQ_ENABLE_HTS)
  os << "    \"format\": \""
     << (rep.from_alignment ? rep.alignment_format : format_name(rep.format))
     << "\",\n";
#else
  os << "    \"format\": \"" << format_name(rep.format) << "\",\n";
#endif
  os << "    \"block_parallel\": " << (rep.format == GzipFormat::kBgzf ? "true" : "false")
     << ",\n";
  os << "    \"compressed_bytes\": " << rep.compressed_bytes << ",\n";
  os << "    \"phred_offset\": " << cfg.phred_offset << "\n";
  os << "  },\n";
  os << "  \"threads\": {\"inflaters\": " << rep.inflater_threads
     << ", \"workers\": " << rep.worker_threads << "},\n";
  os << "  \"summary\": {\n";
  os << "    \"reads_in\": " << s.reads_in << ",\n";
  os << "    \"bases_in\": " << s.bases_in << ",\n";
  os << "    \"reads_out\": " << s.reads_out << ",\n";
  os << "    \"bases_out\": " << s.bases_out << ",\n";
  os << "    \"mean_quality\": " << f(s.mean_quality()) << ",\n";
  os << "    \"q20_rate\": "
     << f(s.bases_in ? static_cast<double>(s.q20) / static_cast<double>(s.bases_in) : 0.0)
     << ",\n";
  os << "    \"q30_rate\": " << f(s.q30_rate()) << ",\n";
  os << "    \"gc_rate\": " << f(s.gc_rate()) << ",\n";
  os << "    \"n_bases\": " << s.n_bases << ",\n";
  // The integer numerators behind the rates above. A rate is printed with six
  // significant digits, so the count cannot be recovered from it and shards
  // cannot be merged exactly from rates -- the reason genomic_toolkit's frip
  // report carries bases_in_peaks beside base_frip. Added, never renamed: the
  // rate fields keep their names for the scripts that read them.
  os << "    \"quality_sum\": " << s.qsum_in << ",\n";
  os << "    \"q20_bases\": " << s.q20 << ",\n";
  os << "    \"q30_bases\": " << s.q30 << ",\n";
  os << "    \"gc_bases\": " << s.gc << "\n";
  os << "  },\n";
  os << "  \"filtering\": {\n";
  os << "    \"adapter_trimmed\": " << s.adapter_trimmed << ",\n";
  os << "    \"quality_trimmed\": " << s.quality_trimmed << ",\n";
  os << "    \"dropped_short_or_low_quality\": " << s.dropped_short << ",\n";
  os << "    \"dropped_excess_n\": " << s.dropped_n << "\n";
  os << "  },\n";

  os << "  \"quality_histogram\": ";
  // The report's 64-bin schema: Q0-Q62 as they are, Q63 and above folded into
  // the last bin, exactly as the histogram was clamped before it grew.
  std::array<uint64_t, kReportQualBins> report_hist{};
  for (std::size_t q = 0; q < kQualBins; ++q) {
    report_hist[std::min(q, kReportQualBins - 1)] += s.qhist[q];
  }
  detail::json_u64_array(os, report_hist.data(), kReportQualBins);
  os << ",\n";

  os << "  \"per_cycle\": {\n";
  os << "    \"cycles\": " << cycles << ",\n";
  // The per-cycle arrays stop at kMaxCycles. Reads past it count in every
  // total and in none of these arrays; without this count a 2 kb run reported
  // "cycles": 1024 over half its bases and nothing said so.
  os << "    \"cycle_limit\": " << kMaxCycles << ",\n";
  os << "    \"reads_longer_than_cycle_limit\": " << s.reads_over_cycle_limit << ",\n";
  os << "    \"mean_quality\": [";
  for (std::size_t i = 0; i < cycles; ++i) {
    if (i) os << ',';
    const double m = s.pos_count[i]
                         ? static_cast<double>(s.pos_qsum[i]) /
                                   static_cast<double>(s.pos_count[i]) -
                               kPhredOffset
                         : 0.0;
    os << f(m);
  }
  os << "],\n";
  static const char* kBaseKeys[5] = {"a", "c", "g", "t", "n"};
  for (int b = 0; b < 5; ++b) {
    os << "    \"" << kBaseKeys[b] << "\": [";
    for (std::size_t i = 0; i < cycles; ++i) {
      if (i) os << ',';
      os << s.pos_base[i][static_cast<std::size_t>(b)];
    }
    os << (b == 4 ? "]\n" : "],\n");
  }
  os << "  },\n";

  os << "  \"performance\": {\n";
  os << "    \"wall_seconds\": " << f(rep.wall_seconds) << ",\n";
  os << "    \"compressed_mb_per_s\": "
     << f(rep.wall_seconds > 0 ? static_cast<double>(rep.compressed_bytes) /
                                     (1024.0 * 1024.0) / rep.wall_seconds
                               : 0.0)
     << ",\n";
  os << "    \"bytes_written\": " << rep.bytes_written << ",\n";
  os << "    \"buffer_pool_blocks\": " << rep.pool_blocks << ",\n";
  os << "    \"buffer_pool_bytes\": " << rep.pool_bytes << "\n";
  os << "  }\n";
  os << "}\n";
}

}  // namespace fq
