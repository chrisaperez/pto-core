// SPDX-License-Identifier: MIT
// fastq_stream — command-line driver.
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include "fastq_stream/pipeline.hpp"
#include "fastq_stream/report.hpp"

namespace {

void usage(std::ostream& os) {
  os << R"(fastq_stream )" << fq::kVersion << R"( — low-latency in-memory FASTQ QC and trimming

USAGE
  fastq_stream -i <in.fq.gz> -o <out.fq|fifo> [options]
  fastq_stream -i in.fq.gz --stdout | bwa-mem2 mem ref.fa - > out.sam

INPUT / OUTPUT
  -i, --in PATH            input FASTQ (gzip/BGZF/plain), or BAM/CRAM/SAM when
                           built with -DFQ_ENABLE_HTS=ON ("-" = stdin) [-]
                           BAM/CRAM is detected only in a regular file, never
                           through a pipe.
  -o, --out PATH           output FASTQ ("-" = stdout) [-]
      --stdout             shorthand for --out -
      --mkfifo             create the output path as a FIFO if it does not exist
      --qc-only            compute QC only; write no reads
      --phred-offset N     quality encoding of FASTQ input, 33 or 64 [33].
                           64 (Illumina 1.3-1.7) is validated and written out
                           as Phred+33. Without this option, input whose
                           qualities fit only Phred+64 fails with exit 2 after
                           the report; pass 33 to confirm a Phred+33 file.

Alignment input (BAM/CRAM/SAM; detected by content, not by file extension)
      --interleaved        emit mates as name/1 then name/2. Needs name-grouped
                           input; a coordinate-sorted file is rejected with the
                           samtools collate command to fix it.
      --reference PATH     reference FASTA. REQUIRED for CRAM: without it htslib
                           resolves references over the network, which this tool
                           does not do.
      --missing-qual C     quality character substituted when QUAL is absent [I]
  -j, --json PATH          write the JSON QC report to PATH
      --report PATH        write the human-readable summary to PATH [stderr]

THREADS
  -t, --threads N          total threads (default: hardware concurrency)

QUALITY TRIMMING
      --trim-front Q       trim 5' bases below Phred Q (0 = off) [0]
      --trim-tail Q        trim 3' bases below Phred Q (0 = off) [0]
      --window N           sliding window width, 0 = off [4]
      --window-mean Q      cut at the first window with mean Phred < Q [20]

FILTERS
  -l, --min-len N          discard reads shorter than N after trimming [15]
      --max-n-rate F       discard reads with an N fraction above F (<0 = off) [0.1]
      --min-mean-q Q       discard reads with mean Phred below Q (0 = off) [0]

ADAPTERS
      --no-adapter         disable adapter trimming
  -a, --adapter SEQ        adapter sequence, ACGTN only (repeatable; replaces the
                           defaults, which are TruSeq AGATCGGAAGAGC and Nextera
                           CTGTCTCTTATACACATCT)
      --adapter-seed N     exact seed length used to anchor a match; must not
                           exceed --adapter-overlap [4]
      --adapter-overlap N  minimum adapter overlap to trim; must not exceed the
                           shortest adapter [4]
      --adapter-mismatch F mismatch rate tolerated in the overlap, 0..0.5 [0.2]

OTHER
  -v, --verbose            report progress on stderr
  -h, --help               this message
      --version            print the version
)";
}

[[noreturn]] void fail(const std::string& msg) {
  std::cerr << "fastq_stream: " << msg << "\n";
  std::exit(2);
}

// Returns argv[++i], or exits if the option is missing its argument.
const char* value(int argc, char** argv, int& i) {
  if (i + 1 >= argc) fail(std::string("option ") + argv[i] + " requires an argument");
  return argv[++i];
}

long parse_long(const char* s, const char* opt) {
  char* end = nullptr;
  const long v = std::strtol(s, &end, 10);
  if (end == s || *end != '\0') fail(std::string("invalid integer for ") + opt + ": " + s);
  return v;
}

double parse_double(const char* s, const char* opt) {
  char* end = nullptr;
  const double v = std::strtod(s, &end);
  if (end == s || *end != '\0') fail(std::string("invalid number for ") + opt + ": " + s);
  return v;
}

// Range-checked variants.
//
// Every numeric option used to be parsed and then narrowed with an unchecked
// static_cast, which failed in two distinct ways:
//
//   * `--trim-tail 300` truncated through uint8_t to 44 and ran silently at a
//     threshold the user never asked for.
//   * `--window-mean -5` / `1e30` / `nan` reached
//     trim_sliding_window's `static_cast<uint32_t>(min_mean * window)`, which
//     is undefined behaviour for a value the destination cannot represent.
//     UBSan reported all three.
//
// The kernel now clamps defensively as well, but a CLI that accepts a value it
// cannot honour and silently substitutes another is its own defect: refusing
// the argument is the only answer that leaves the user's intent intact.
long parse_long_in(const char* s, const char* opt, long lo, long hi) {
  const long v = parse_long(s, opt);
  if (v < lo || v > hi) {
    fail(std::string("value out of range for ") + opt + ": " + s + " (expected " +
         std::to_string(lo) + ".." + std::to_string(hi) + ")");
  }
  return v;
}

double parse_double_in(const char* s, const char* opt, double lo, double hi) {
  const double v = parse_double(s, opt);
  // Written as a positive test so NaN, for which every comparison is false,
  // is rejected rather than slipping through a negated range check.
  if (!(v >= lo && v <= hi)) {
    fail(std::string("value out of range for ") + opt + ": " + s + " (expected " +
         std::to_string(lo) + ".." + std::to_string(hi) + ")");
  }
  return v;
}

// The adapter options interact, and the trimmer resolved every conflict by
// silently substituting a value -- the defect the comment above describes, in
// the three options it did not reach:
//
//   * `-a ''` was dropped by AdapterTrimmer::add, but it still made the
//     adapter list non-empty, so the defaults were not loaded either: adapter
//     trimming turned itself off with the report still saying it ran.
//   * `--adapter-seed 12` was clamped to --adapter-overlap (4): the user asked
//     for a stricter seed and got the loosest one.
//   * `--adapter-overlap 14` exceeds the 13 bp TruSeq default, and an adapter
//     shorter than the overlap can never be trimmed -- no error, no trims.
//   * `--adapter-mismatch` was accepted up to 1.0 and clamped to 0.5 by
//     set_max_mismatch_rate.
void validate_adapters(const fq::Config& cfg) {
  if (!cfg.trim_adapters) return;
  if (cfg.adapter_seed > cfg.adapter_min_overlap) {
    fail("--adapter-seed " + std::to_string(cfg.adapter_seed) + " exceeds --adapter-overlap " +
         std::to_string(cfg.adapter_min_overlap) +
         ": the seed is part of the overlap, so it cannot be longer");
  }
  const std::vector<std::string_view> defaults = {fq::kTruSeqAdapter, fq::kNexteraAdapter};
  std::vector<std::string_view> list(cfg.adapters.begin(), cfg.adapters.end());
  if (list.empty()) list = defaults;
  for (const std::string_view a : list) {
    if (a.empty()) fail("--adapter takes a non-empty sequence");
    for (const char c : a) {
      if (std::strchr("ACGTNacgtn", c) == nullptr || c == '\0') {
        fail("--adapter '" + std::string(a) + "' contains '" + std::string(1, c) +
             "'; adapter sequences are ACGTN only");
      }
    }
    if (a.size() < cfg.adapter_min_overlap) {
      fail("adapter " + std::string(a) + " (" + std::to_string(a.size()) +
           " bp) is shorter than --adapter-overlap " +
           std::to_string(cfg.adapter_min_overlap) + " and could never be trimmed");
    }
  }
}

// Writes a report file and checks that it was written. An ofstream that hits
// EPIPE, ENOSPC or EDQUOT sets badbit and nothing else; the previous code
// never looked, and exited 0 over a truncated JSON report. Exit 1, not 2: this
// is an environment failure, and 2 means the input was the problem.
void write_file(const std::string& path, const char* what,
                const std::function<void(std::ostream&)>& emit) {
  std::ofstream f(path);
  if (f) {
    emit(f);
    f.close();
  }
  if (!f) {
    std::cerr << "fastq_stream: cannot write " << what << " to " << path << ": "
              << std::strerror(errno) << "\n";
    std::exit(1);
  }
}

}  // namespace

int main(int argc, char** argv) {
  fq::Config cfg;
  std::string json_path;
  std::string report_path;
  bool phred_offset_explicit = false;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto is = [&](const char* a, const char* b = nullptr) {
      return arg == a || (b != nullptr && arg == b);
    };

    if (is("-h", "--help")) {
      usage(std::cout);
      return 0;
    } else if (is("--version")) {
      std::cout << fq::kVersion << "\n";
      return 0;
    } else if (is("-i", "--in")) {
      cfg.input = value(argc, argv, i);
    } else if (is("-o", "--out")) {
      cfg.output = value(argc, argv, i);
    } else if (is("--stdout")) {
      cfg.output = "-";
    } else if (is("--mkfifo")) {
      cfg.make_fifo = true;
    } else if (is("--interleaved")) {
#if defined(FQ_ENABLE_HTS)
      cfg.hts_layout = fq::HtsLayout::kInterleaved;
#else
      fail("--interleaved requires a build with -DFQ_ENABLE_HTS=ON");
#endif
    } else if (is("--reference")) {
      cfg.reference_fasta = value(argc, argv, i);
    } else if (is("--missing-qual")) {
#if defined(FQ_ENABLE_HTS)
      {
        const std::string v = value(argc, argv, i);
        // The substituted character becomes a QUAL byte, and the pipeline
        // rejects any QUAL byte outside '!'..'~'. Refuse it here, by name,
        // rather than as a "malformed FASTQ record" the user never wrote.
        if (v.size() != 1 || v[0] < '!' || v[0] > '~') {
          fail("--missing-qual takes a single printable character ('!'..'~')");
        }
        cfg.hts_missing_qual = v[0];
      }
#else
      fail("--missing-qual requires a build with -DFQ_ENABLE_HTS=ON");
#endif
    } else if (is("--qc-only")) {
      cfg.qc_only = true;
    } else if (is("--phred-offset")) {
      const long offset = parse_long(value(argc, argv, i), "--phred-offset");
      if (offset != 33 && offset != 64) fail("--phred-offset must be 33 or 64");
      cfg.phred_offset = static_cast<int>(offset);
      phred_offset_explicit = true;
    } else if (is("-j", "--json")) {
      json_path = value(argc, argv, i);
    } else if (is("--report")) {
      report_path = value(argc, argv, i);
    } else if (is("-t", "--threads")) {
      cfg.threads = static_cast<int>(parse_long_in(value(argc, argv, i), "--threads", 0, 4096));
    } else if (is("--trim-front")) {
      cfg.front_q = static_cast<uint8_t>(
          parse_long_in(value(argc, argv, i), "--trim-front", 0, 93));
    } else if (is("--trim-tail")) {
      cfg.tail_q = static_cast<uint8_t>(
          parse_long_in(value(argc, argv, i), "--trim-tail", 0, 93));
    } else if (is("--window")) {
      cfg.window = static_cast<std::size_t>(
          parse_long_in(value(argc, argv, i), "--window", 0, 1000000));
    } else if (is("--window-mean")) {
      cfg.window_mean = parse_double_in(value(argc, argv, i), "--window-mean", 0.0, 93.0);
    } else if (is("-l", "--min-len")) {
      cfg.min_len = static_cast<std::size_t>(
          parse_long_in(value(argc, argv, i), "--min-len", 0, 1000000000));
    } else if (is("--max-n-rate")) {
      cfg.max_n_rate = parse_double_in(value(argc, argv, i), "--max-n-rate", -1.0, 1.0);
    } else if (is("--min-mean-q")) {
      cfg.min_mean_q = parse_double_in(value(argc, argv, i), "--min-mean-q", 0.0, 93.0);
    } else if (is("--no-adapter")) {
      cfg.trim_adapters = false;
    } else if (is("-a", "--adapter")) {
      cfg.adapters.emplace_back(value(argc, argv, i));
    } else if (is("--adapter-seed")) {
      cfg.adapter_seed = static_cast<std::size_t>(
          parse_long_in(value(argc, argv, i), "--adapter-seed", 1, 1024));
    } else if (is("--adapter-overlap")) {
      cfg.adapter_min_overlap = static_cast<std::size_t>(
          parse_long_in(value(argc, argv, i), "--adapter-overlap", 1, 1024));
    } else if (is("--adapter-mismatch")) {
      cfg.adapter_mismatch = parse_double_in(value(argc, argv, i), "--adapter-mismatch", 0.0, 0.5);
    } else if (is("-v", "--verbose")) {
      cfg.verbose = true;
    } else if (!arg.empty() && arg[0] == '-' && arg != "-") {
      fail("unknown option: " + arg);
    } else {
      fail("unexpected positional argument: " + arg);
    }
  }

  validate_adapters(cfg);

  // Detect alignment input by content rather than extension: a BAM named
  // .fq.gz is still a BAM, and letting it reach the gzip path produces a
  // "malformed FASTQ record" error several ranks downstream that says nothing
  // about the real cause. Without htslib compiled in, run_pipeline reports the
  // build-configuration problem instead.
#if defined(FQ_ENABLE_HTS)
  cfg.hts_input = fq::looks_like_alignment(cfg.input);
#endif
  // BAM and CRAM store quality values as numbers; there is no text offset to
  // choose, and rebasing them by 31 would corrupt every one.
  if (cfg.hts_input && cfg.phred_offset != 33) {
    fail("--phred-offset applies to FASTQ text; BAM/CRAM stores quality values numerically");
  }

  try {
    const fq::RunReport rep = fq::run_pipeline(cfg);

    if (!report_path.empty()) {
      write_file(report_path, "report",
                 [&](std::ostream& os) { fq::write_text_report(os, cfg, rep); });
    } else {
      fq::write_text_report(std::cerr, cfg, rep);
    }

    if (!json_path.empty()) {
      write_file(json_path, "JSON report",
                 [&](std::ostream& os) { fq::write_json_report(os, cfg, rep); });
    }

    // Exit 2, the repository's "the input's shape made the answer wrong": the
    // report above was computed, and every quality figure and trim in it used
    // an offset of 33 that these bytes cannot have been written with. Checked
    // after the fact because it is read off the run's own quality histogram;
    // an explicit --phred-offset is the user saying they know.
    if (!phred_offset_explicit && !cfg.hts_input && fq::looks_like_phred64(rep.stats)) {
      std::cerr << "fastq_stream: read at offset 33, every quality is between Q31 and Q72 "
                   "and at least 1% are Q61 or higher: that is the range Phred+64 "
                   "(Illumina 1.3-1.7) produces, and no Phred+33 short-read platform "
                   "does. The quality figures and trimming above are wrong for such a "
                   "file. Re-run with --phred-offset 64, or --phred-offset 33 if this "
                   "really is Phred+33.\n";
      return 2;
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "fastq_stream: " << e.what() << "\n";
    return 1;
  }
}
