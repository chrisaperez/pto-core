// SPDX-License-Identifier: MIT
//
// pto-peaks: streaming peak caller.
//
//     pto-peaks --chrom-sizes hg38.chrom.sizes fragments.bed > peaks.narrowPeak
//
// Input is a coordinate-sorted fragment BED (chrom, start, end), which is what
// an scATAC pipeline already produces and what `genomic_toolkit` emits. Text
// rather than BAM on purpose: this module has no htslib dependency, so it
// builds on a cluster with a compiler and nothing else, and every test drives
// the same path the CLI does.
//
// Output is narrowPeak (BED6+4), the format MACS2 writes and every genome
// browser reads.
//
// Exit codes follow the repository's convention, which pto-cloud's
// `job_specs.RETRYABLE_EXIT_CODES` depends on:
//
//     0  peaks called (possibly none)
//     1  usage or I/O failure
//     2  the input's shape made the answer wrong -- unsorted records, a contig
//        missing from the sizes file, a malformed field. Deterministic, so
//        retrying reaches the same answer; a valid report is printed first.
#include <algorithm>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <memory>
#include <new>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "peaks/bam_streamer.hpp"
#include "peaks/caller.hpp"
#include "peaks/file_identity.hpp"
#include "peaks/poisson_model.hpp"
#include "peaks/sliding_window.hpp"

namespace {

using pto::peaks::AddResult;
using pto::peaks::BamFilter;
using pto::peaks::BamFragment;
using pto::peaks::BamStatus;
using pto::peaks::BamStreamer;
using pto::peaks::CandidatePeak;
using pto::peaks::choose_seek_target;
using pto::peaks::Coord;
using pto::peaks::kMaxHalfWidth;
using pto::peaks::MultiScaleWindow;
using pto::peaks::PeakCaller;
using pto::peaks::PeakCallOptions;

constexpr int kExitOk = 0;
constexpr int kExitUsage = 1;
constexpr int kExitInvalidInput = 2;

struct Options {
  std::string chrom_sizes;
  std::string input = "-";
  std::string output = "-";
  double lambda_bg = 0.0;
  int min_mapq = 0;
  PeakCallOptions call;
};

void usage(std::FILE* to) {
  std::fprintf(to,
               "pto-peaks " PEAKS_VERSION " -- streaming peak caller\n"
               "\n"
               "usage: pto-peaks --chrom-sizes FILE [options] [fragments.bed]\n"
               "\n"
               "  --chrom-sizes FILE   two columns, name and length. Required for\n"
               "                       BED input; a BAM carries its own\n"
               "  --min-mapq N         BAM only: skip alignments below this (default 0)\n"
               "  -o, --out FILE       narrowPeak output (default: stdout). May not\n"
               "                       be the input or the --chrom-sizes file\n"
               "  --cutoff X           -log10 p to open a peak (default 2.0)\n"
               "  --qvalue X           Benjamini-Hochberg FDR level (default 0.05)\n"
               "  --min-length N       discard calls shorter than N bp (default 50)\n"
               "  --max-gap N          sub-threshold bases tolerated inside a peak\n"
               "                       (default 30)\n"
               "  --extend-peaks N     pad every reported peak by N bp on each side\n"
               "                       (default 0). 0 keeps pto-peaks' sharp,\n"
               "                       high-resolution bounds; a positive value is\n"
               "                       a compatibility mode for the wider spans\n"
               "                       MACS2/MACS3 produce. The summit position and\n"
               "                       the -log10 p / -log10 q columns are\n"
               "                       unchanged; the 3' edge is clamped to the\n"
               "                       contig length\n"
               "  --lambda-bg X        genome-wide background, fragments per base.\n"
               "                       Used as a floor under the local estimate;\n"
               "                       0 (the default) means local scales only\n"
               "  -h, --help           this text\n"
               "\n"
               "Input is a coordinate-sorted fragment BED, or a coordinate-sorted\n"
               "BAM -- the format is detected from the first byte, so either works\n"
               "on stdin. A BAM yields one fragment per proper pair, spanning TLEN\n"
               "from the leftmost mate. A gzip- or BGZF-compressed fragment BED\n"
               "(10x fragments.tsv.gz) is not read directly: zcat it into stdin.\n"
               "\n"
               "Records that arrive out of order are refused rather than\n"
               "mis-placed, and the run exits 2.\n");
}

// std::from_chars writes the parsed prefix BEFORE reporting trailing junk, so
// a discarded `false` still leaves a plausible value behind. Every field here
// therefore checks the error code AND that the whole field was consumed --
// the rule docs/REVIEW_2026-08-15.md records as `stream_bedpe` letting "XX"
// become MAPQ 255.
bool parse_i64(std::string_view text, std::int64_t& out) {
  if (text.empty()) return false;
  const char* first = text.data();
  const char* last = text.data() + text.size();
  const auto res = std::from_chars(first, last, out);
  return res.ec == std::errc{} && res.ptr == last;
}

bool parse_double(std::string_view text, double& out) {
  if (text.empty()) return false;
  // from_chars for double is not available everywhere libstdc++ ships, so this
  // uses strtod with the same "consumed everything" discipline.
  std::string owned(text);
  char* end = nullptr;
  errno = 0;
  const double v = std::strtod(owned.c_str(), &end);
  if (end == owned.c_str() || *end != '\0' || errno != 0) return false;
  out = v;
  return true;
}

// Every numeric option is a finite number. strtod accepts "nan" and "inf", and
// none of this module's range checks rejected them: `--lambda-bg nan` passed
// `< 0.0` (every comparison with NaN is false), then made every local lambda
// NaN, which log_poisson_sf fails closed to p = 1 -- and the run exited 0 with
// zero peaks. `--lambda-bg inf` and `--cutoff inf` did the same by a
// different route. Failing closed is right for a corrupt VALUE met while
// streaming; for an argument the user typed, the answer is to refuse it.
bool parse_finite(std::string_view text, double& out) {
  double v = 0.0;
  if (!parse_double(text, v) || !std::isfinite(v)) return false;
  out = v;
  return true;
}

std::string_view field(std::string_view line, std::size_t& cursor) {
  if (cursor > line.size()) return {};
  const std::size_t begin = cursor;
  std::size_t end = line.find('\t', begin);
  if (end == std::string_view::npos) end = line.size();
  cursor = end + 1;
  return line.substr(begin, end - begin);
}

// Lines, bounded.
//
// std::getline grows its string to whatever precedes the next newline. A
// fragment BED line is ~30 bytes; an input that is not a BED -- a CRAM, a zstd
// file, 300 MB with no newline -- was read into memory WHOLE before being
// refused as malformed (704 MB resident for a 300 MB line). And where that
// allocation failed, under a container memory limit, getline turned the
// std::bad_alloc into a failed stream, the loop took that for end of input,
// and the run exited 0 reporting "0 fragments, 0 peaks": the input unread and
// success claimed. Measured with an allocator that refuses requests over
// 64 MB.
//
// So lines are cut out of one fixed block with memchr, straight from the
// stream buffer, and a line that outgrows the cap is refused where it is found.
// Nothing is allocated per line and no stream state is involved: the first
// version of this used istream::getline into the same fixed buffer and was 4%
// slower end to end than the std::getline it replaced (5.25 s against 5.05 s
// on 6M fragments); this one reads in 4 MiB blocks. An exception from the
// stream buffer propagates to main's handler instead of becoming a failbit the
// loop could mistake for end of input.
constexpr std::size_t kMaxLineBytes = 1u << 20;
constexpr std::size_t kReadBlockBytes = 4u << 20;
static_assert(kReadBlockBytes > kMaxLineBytes, "a block must hold the longest line");

enum class LineStatus { kLine, kEnd, kTooLong, kStreamError };

class LineReader {
 public:
  explicit LineReader(std::istream& in) : buf_(in.rdbuf()), data_(kReadBlockBytes) {}

  LineStatus next(std::string_view& line) {
    if (buf_ == nullptr) return LineStatus::kStreamError;
    for (;;) {
      const char* start = data_.data() + pos_;
      const std::size_t avail = end_ - pos_;
      if (const void* nl = std::memchr(start, '\n', avail)) {
        const auto len = static_cast<std::size_t>(static_cast<const char*>(nl) - start);
        if (len > kMaxLineBytes) return LineStatus::kTooLong;
        line = std::string_view(start, len);
        pos_ += len + 1;
        ++lines_;
        return LineStatus::kLine;
      }
      if (avail > kMaxLineBytes) return LineStatus::kTooLong;
      if (eof_) {
        if (avail == 0) return LineStatus::kEnd;
        line = std::string_view(start, avail);  // a last line with no newline
        pos_ = end_;
        ++lines_;
        return LineStatus::kLine;
      }
      // Keep the partial line, refill behind it. avail <= kMaxLineBytes here,
      // so the block always has room for more.
      if (pos_ > 0) {
        std::memmove(data_.data(), start, avail);
        pos_ = 0;
        end_ = avail;
      }
      const std::streamsize got = buf_->sgetn(data_.data() + end_,
                                              static_cast<std::streamsize>(data_.size() - end_));
      if (got <= 0) {
        eof_ = true;
      } else {
        end_ += static_cast<std::size_t>(got);
      }
    }
  }

  [[nodiscard]] std::size_t lines() const noexcept { return lines_; }

 private:
  std::streambuf* buf_;
  std::vector<char> data_;
  std::size_t pos_ = 0;
  std::size_t end_ = 0;
  std::size_t lines_ = 0;
  bool eof_ = false;
};

// Strips one trailing CR. Done before any test for an empty line: a CRLF blank
// line is "\r", and testing emptiness first made it a malformed record.
void strip_cr(std::string_view& line) {
  if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
}

struct ContigTable {
  std::vector<std::string> names;
  std::vector<Coord> lengths;
  std::unordered_map<std::string, std::int32_t> index;

  void adopt(const std::vector<std::string>& n, const std::vector<Coord>& l) {
    names = n;
    lengths = l;
    index.clear();
    for (std::size_t i = 0; i < names.size(); ++i) {
      index.emplace(names[i], static_cast<std::int32_t>(i));
    }
  }

  bool load(const std::string& path, std::string& error) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
      error = "cannot open chrom sizes file: " + path;
      return false;
    }
    LineReader reader(in);
    std::string_view line;
    for (;;) {
      const LineStatus st = reader.next(line);
      if (st == LineStatus::kEnd) break;
      const std::string where = path + ":" + std::to_string(reader.lines() + 1);
      if (st == LineStatus::kTooLong) {
        error = where + ": line longer than " + std::to_string(kMaxLineBytes) + " bytes";
        return false;
      }
      if (st == LineStatus::kStreamError) {
        error = where + ": read error";
        return false;
      }
      // A chrom.sizes written on Windows, or by a tool that emits CRLF, used to
      // fail every line: "100000\r" is not an integer.
      strip_cr(line);
      if (line.empty() || line[0] == '#') continue;
      std::size_t cursor = 0;
      const std::string_view name = field(line, cursor);
      const std::string_view len_text = field(line, cursor);
      std::int64_t len = 0;
      const std::string lineref = path + ":" + std::to_string(reader.lines());
      if (name.empty()) {
        error = lineref + ": malformed contig record; empty contig name (field 1)";
        return false;
      }
      if (!parse_i64(len_text, len)) {
        error = lineref + ": malformed contig record; length field '" +
                std::string(len_text) + "' is not an integer";
        return false;
      }
      if (len <= 0) {
        error = lineref + ": malformed contig record; length " + std::to_string(len) +
                " must be positive";
        return false;
      }
      if (len > static_cast<std::int64_t>(std::numeric_limits<Coord>::max())) {
        error = lineref + ": malformed contig record; length " + std::to_string(len) +
                " exceeds the maximum representable contig length (" +
                std::to_string(std::numeric_limits<Coord>::max()) + ")";
        return false;
      }
      if (index.count(std::string(name)) != 0) {
        error = lineref + ": duplicate contig " + std::string(name);
        return false;
      }
      index.emplace(std::string(name), static_cast<std::int32_t>(names.size()));
      names.emplace_back(name);
      lengths.push_back(static_cast<Coord>(len));
    }
    if (names.empty()) {
      error = path + ": no contigs";
      return false;
    }
    return true;
  }
};

// Ties the window and the caller together. The interleaving is the part a
// caller must not get wrong: the ring holds a bounded lookahead, so fragments
// and evaluation have to advance together.
class Driver {
 public:
  Driver(const Options& opts, const ContigTable& contigs)
      : opts_(opts), contigs_(contigs), caller_(opts.call) {}

  bool begin_contig(std::int32_t tid) {
    end_contig();
    window_.reset(tid, contigs_.lengths[static_cast<std::size_t>(tid)]);
    caller_.begin_contig(tid);
    active_ = true;
    return true;
  }

  // Before adding a fragment far past whatever came before it: fast-forward
  // the window across the coverage desert in between, when doing so is
  // provably safe. See choose_seek_target's own comment (sliding_window.hpp)
  // for why safety, not merely a big-enough gap, is the deciding question:
  // seek()'s unconditional ring wipe means this only ever fires for a LEADING
  // gap (nothing added to this contig yet) or when bounded-draining with
  // whatever frontier earlier fragments already provided happens to reach all
  // the way to signal_limit() on its own -- an interior gap between two real,
  // narrow peaks does neither, and is left to drain normally rather than risk
  // discarding a peak's own not-yet-observed signal.
  void maybe_seek_before(Coord start) {
    const std::int64_t live = window_.signal_limit();
    bool safe_to_seek;
    if (live < 0) {
      safe_to_seek = true;  // nothing added to this contig yet: nothing to lose
    } else {
      bounded_drain(live);
      safe_to_seek = (static_cast<std::int64_t>(window_.position()) >= live);
    }
    if (!safe_to_seek) return;
    const std::int64_t candidate = static_cast<std::int64_t>(start) - kMaxHalfWidth;
    const std::int64_t t = choose_seek_target(live, candidate, window_.contig_length(),
                                              window_.position(), opts_.call.score_cutoff);
    if (t >= 0) window_.seek(static_cast<Coord>(t));
  }

  AddResult add(Coord start, Coord end) {
    maybe_seek_before(start);
    AddResult r = window_.add_fragment(start, end);
    while (r == AddResult::kNeedsDrain) {
      if (!window_.step()) break;
      visit();
      r = window_.add_fragment(start, end);
    }
    if (r == AddResult::kAccepted) drain();
    return r;
  }

  void end_contig() {
    if (!active_) return;
    // The tail after the last fragment (or the whole contig, if it never
    // carried one at all -- an unplaced/alt contig, or chrY in a female
    // sample) is the same shape of gap as maybe_seek_before's leading case,
    // just with nothing more coming: finish() unlocks evaluable_limit() to
    // the contig's end regardless of the frontier (which alone is stuck at
    // the last fragment's own start, not its end), a bounded drain reaches
    // signal_limit() properly, and everything past it is provably k=0
    // background -- PeakCaller::flush() finalises an open peak from its own
    // already-stored fields regardless of how much of that background is
    // ever stepped through, so there is nothing left to gain by walking it
    // one base at a time. Still gated the same way: at score_cutoff <= 0
    // every base is "above" regardless of pileup (see
    // detail::provably_below_cutoff), so an open peak keeps growing through
    // the background too, and stopping early would cut it short.
    const std::int64_t live = window_.signal_limit();
    window_.finish();
    const std::int64_t stop_at = choose_seek_target(
        live, window_.contig_length() - 1, window_.contig_length(), window_.position(),
        opts_.call.score_cutoff);
    if (stop_at >= 0) {
      bounded_drain(live);
    } else {
      drain();
    }
    caller_.end_contig();
    active_ = false;
  }

  PeakCaller& caller() noexcept { return caller_; }
  [[nodiscard]] const MultiScaleWindow& window() const noexcept { return window_; }

 private:
  void visit() { caller_.observe(window_.position(), window_.pileup(), lambda()); }
  void drain() {
    while (window_.step()) visit();
  }
  void bounded_drain(std::int64_t up_to) {
    while (window_.position() < up_to && window_.step()) visit();
  }
  [[nodiscard]] double lambda() const { return window_.local_lambda(opts_.lambda_bg); }

  const Options& opts_;
  const ContigTable& contigs_;
  MultiScaleWindow window_;
  PeakCaller caller_;
  bool active_ = false;
};

// Streams a coordinate-sorted BAM through the same driver the BED path uses.
// The contig table comes from the BAM header, which is why --chrom-sizes is
// not required for this input.
int stream_bam_records(const Options& opts, const ContigTable& contigs,
                       Driver& driver, std::int64_t& accepted, BamStreamer& bam) {
  (void)opts;
  std::vector<bool> seen(contigs.names.size(), false);
  std::int32_t current = pto::peaks::kNoTid;
  BamFragment f{};
  BamStatus st = BamStatus::kOk;
  while ((st = bam.next(f)) == BamStatus::kOk) {
    const std::size_t tid = static_cast<std::size_t>(f.tid);
    if (f.tid != current) {
      if (seen[tid]) {
        std::fprintf(stderr,
                     "pto-peaks: contig '%s' reappears after another contig; "
                     "the BAM is not coordinate-sorted\n",
                     contigs.names[tid].c_str());
        return kExitInvalidInput;
      }
      seen[tid] = true;
      driver.begin_contig(f.tid);
      current = f.tid;
    }
    const AddResult r = driver.add(f.start, f.end);
    if (r != AddResult::kAccepted) {
      std::fprintf(stderr,
                   "pto-peaks: %s:%d: record refused (%s). A BAM must be "
                   "coordinate-sorted.\n",
                   contigs.names[tid].c_str(), f.start,
                   r == AddResult::kOutOfOrder ? "out of order" : "unusable span");
      return kExitInvalidInput;
    }
    ++accepted;
  }
  if (st != BamStatus::kEndOfFile) {
    std::fprintf(stderr, "pto-peaks: %s\n", BamStreamer::describe(st));
    return kExitInvalidInput;
  }
  std::fprintf(stderr,
               "pto-peaks: BAM: %lld records, %lld fragments (%lld by flag, "
               "%lld by mapq, %lld by span, %lld not the leftmost mate, "
               "%lld unplaced)\n",
               static_cast<long long>(bam.records_seen()),
               static_cast<long long>(bam.fragments_emitted()),
               static_cast<long long>(bam.skipped_flags()),
               static_cast<long long>(bam.skipped_mapq()),
               static_cast<long long>(bam.skipped_span()),
               static_cast<long long>(bam.skipped_unpaired()),
               static_cast<long long>(bam.skipped_unplaced()));
  return kExitOk;
}

int run_bed(const Options& opts, std::istream& in, const ContigTable& contigs,
            Driver& driver, std::int64_t& accepted) {
  std::vector<bool> seen(contigs.names.size(), false);
  std::int32_t current = pto::peaks::kNoTid;
  LineReader reader(in);
  std::string_view line;

  for (;;) {
    const LineStatus st = reader.next(line);
    if (st == LineStatus::kEnd) break;
    if (st == LineStatus::kTooLong) {
      std::fprintf(stderr,
                   "pto-peaks: %s:%zu: line longer than %zu bytes; this is not a "
                   "fragment BED\n",
                   opts.input.c_str(), reader.lines() + 1, kMaxLineBytes);
      return kExitInvalidInput;
    }
    if (st == LineStatus::kStreamError) {
      std::fprintf(stderr, "pto-peaks: %s: read error after line %zu\n",
                   opts.input.c_str(), reader.lines());
      return kExitUsage;
    }
    const std::size_t lineno = reader.lines();
    strip_cr(line);
    if (line.empty() || line[0] == '#' || line.rfind("track", 0) == 0) continue;

    std::size_t cursor = 0;
    const std::string_view chrom = field(line, cursor);
    const std::string_view start_text = field(line, cursor);
    const std::string_view end_text = field(line, cursor);
    std::int64_t start = 0;
    std::int64_t end = 0;
    if (chrom.empty() || !parse_i64(start_text, start) || !parse_i64(end_text, end)) {
      std::fprintf(stderr, "pto-peaks: %s:%zu: malformed BED record\n",
                   opts.input.c_str(), lineno);
      return kExitInvalidInput;
    }

    const auto it = contigs.index.find(std::string(chrom));
    if (it == contigs.index.end()) {
      std::fprintf(stderr,
                   "pto-peaks: %s:%zu: contig '%.*s' is not in %s. Contig naming "
                   "must match (chr1 vs 1 is the usual cause).\n",
                   opts.input.c_str(), lineno, static_cast<int>(chrom.size()),
                   chrom.data(), opts.chrom_sizes.c_str());
      return kExitInvalidInput;
    }
    const std::int32_t tid = it->second;

    if (tid != current) {
      if (seen[static_cast<std::size_t>(tid)]) {
        std::fprintf(stderr,
                     "pto-peaks: %s:%zu: contig '%.*s' reappears after another "
                     "contig; input must be grouped by contig\n",
                     opts.input.c_str(), lineno, static_cast<int>(chrom.size()),
                     chrom.data());
        return kExitInvalidInput;
      }
      seen[static_cast<std::size_t>(tid)] = true;
      driver.begin_contig(tid);
      current = tid;
    }

    // Each condition gets its own message with the actual field values, rather
    // than one "interval outside the contig" line for three different defects
    // (a negative offset, an inverted interval, and a coordinate past the
    // contig's real length are different bugs in different upstream tools).
    if (start < 0) {
      std::fprintf(stderr,
                   "pto-peaks: %s:%zu: negative start coordinate %lld (field 2)\n",
                   opts.input.c_str(), lineno, static_cast<long long>(start));
      return kExitInvalidInput;
    }
    if (end <= start) {
      std::fprintf(stderr,
                   "pto-peaks: %s:%zu: inverted or empty interval [%lld, %lld) "
                   "(end must be > start; fields 2-3)\n",
                   opts.input.c_str(), lineno, static_cast<long long>(start),
                   static_cast<long long>(end));
      return kExitInvalidInput;
    }
    if (end > contigs.lengths[static_cast<std::size_t>(tid)]) {
      std::fprintf(stderr,
                   "pto-peaks: %s:%zu: end coordinate %lld exceeds contig '%.*s' "
                   "length %d (field 3; check --chrom-sizes matches this input's "
                   "assembly)\n",
                   opts.input.c_str(), lineno, static_cast<long long>(end),
                   static_cast<int>(chrom.size()), chrom.data(),
                   contigs.lengths[static_cast<std::size_t>(tid)]);
      return kExitInvalidInput;
    }

    const AddResult r = driver.add(static_cast<Coord>(start), static_cast<Coord>(end));
    if (r != AddResult::kAccepted) {
      const char* why = r == AddResult::kOutOfOrder ? "out of order"
                        : r == AddResult::kBadSpan  ? "span is empty or too long"
                                                    : "outside the contig";
      std::fprintf(stderr, "pto-peaks: %s:%zu: record refused (%s)\n",
                   opts.input.c_str(), lineno, why);
      return kExitInvalidInput;
    }
    ++accepted;
  }
  return kExitOk;
}

// Shared tail: close the last contig, correct for multiple testing, write
// narrowPeak. Both front ends end here so the output is identical whichever
// format came in.
int emit(const Options& opts, const ContigTable& contigs, Driver& driver,
         std::int64_t accepted) {
  driver.end_contig();
  PeakCaller& caller = driver.caller();
  caller.finalize();

  std::ofstream file_out;
  if (opts.output != "-") {
    file_out.open(opts.output);
    if (!file_out) {
      std::fprintf(stderr, "pto-peaks: cannot write %s\n", opts.output.c_str());
      return kExitUsage;
    }
  }
  std::ostream& out = (opts.output == "-") ? std::cout : file_out;

  const std::vector<CandidatePeak>& peaks = caller.peaks();
  const std::vector<float>& q = caller.qscores();
  // Grown, not fixed. snprintf returns the length the line WOULD have needed,
  // which is not the length it wrote, and a contig name is neither this
  // program's to choose nor bounded by anything useful -- a BAM header may
  // carry 64 KiB of it. Writing that return value out of a fixed buffer reads
  // past the end of it and puts whatever was there into the output file.
  std::vector<char> buf(512);
  for (std::size_t i = 0; i < peaks.size(); ++i) {
    const CandidatePeak& p = peaks[i];
    // narrowPeak column 5 is an integer display score, capped at 1000 by the
    // format; the real numbers are in columns 7-9.
    long display = static_cast<long>(q[i] * 10.0f);
    if (display > 1000) display = 1000;
    if (display < 0) display = 0;
    const char* fmt = "%s\t%d\t%d\tpeak_%zu\t%ld\t.\t%.5f\t%.5f\t%.5f\t%d\n";
    const char* name = contigs.names[static_cast<std::size_t>(p.tid)].c_str();
    // --extend-peaks may have pushed the 3' bound past the contig. The caller
    // clamps the 5' side (it knows 0 is the floor) and leaves this side to the
    // writer, the layer that has the contig length -- a narrowPeak `end` past
    // the chromosome is a file downstream tools reject. Column 10 stays valid:
    // the summit is inside the enrichment core, which is inside the contig.
    const Coord contig_len = contigs.lengths[static_cast<std::size_t>(p.tid)];
    const Coord peak_end = p.end < contig_len ? p.end : contig_len;
    const double fold = static_cast<double>(p.total_signal) / std::max(1, peak_end - p.start);
    int n = std::snprintf(buf.data(), buf.size(), fmt, name, p.start, peak_end, i + 1,
                          display, fold, static_cast<double>(p.summit_score),
                          static_cast<double>(q[i]), p.summit_pos - p.start);
    if (n < 0) {
      std::fprintf(stderr, "pto-peaks: cannot format peak_%zu\n", i + 1);
      return kExitUsage;
    }
    if (static_cast<std::size_t>(n) >= buf.size()) {
      buf.resize(static_cast<std::size_t>(n) + 1);
      n = std::snprintf(buf.data(), buf.size(), fmt, name, p.start, peak_end, i + 1,
                        display, fold, static_cast<double>(p.summit_score),
                        static_cast<double>(q[i]), p.summit_pos - p.start);
      if (n < 0 || static_cast<std::size_t>(n) >= buf.size()) {
        std::fprintf(stderr, "pto-peaks: cannot format peak_%zu\n", i + 1);
        return kExitUsage;
      }
    }
    out.write(buf.data(), n);
  }
  // Checked. An ostream that hits ENOSPC, EDQUOT or EFBIG sets badbit and
  // nothing else, and this used to exit 0 over a narrowPeak cut off mid-row --
  // measured with a 512-byte file size limit: 1024 bytes written, 30 peaks
  // called, success. close() is where a buffered tail is actually written.
  out.flush();
  if (file_out.is_open()) file_out.close();
  const bool written = opts.output == "-" ? static_cast<bool>(std::cout)
                                          : static_cast<bool>(file_out);
  if (!written) {
    std::fprintf(stderr, "pto-peaks: writing %s failed: %s\n",
                 opts.output == "-" ? "stdout" : opts.output.c_str(), std::strerror(errno));
    return kExitUsage;
  }

  std::fprintf(stderr,
               "pto-peaks: %lld fragments, %zu peaks (%lld discarded as shorter "
               "than %d bp)\n",
               static_cast<long long>(accepted), peaks.size(),
               static_cast<long long>(caller.discarded_short()), opts.call.min_length);
  if (caller.saturated()) {
    std::fprintf(stderr, "pto-peaks: warning: a counter saturated; signal values "
                         "at the extremes are pinned, not wrapped\n");
  }
  return kExitOk;
}

int run(const Options& opts) {
  std::ifstream file_in;
  if (opts.input != "-") {
    file_in.open(opts.input, std::ios::binary);
    if (!file_in) {
      std::fprintf(stderr, "pto-peaks: cannot open %s\n", opts.input.c_str());
      return kExitUsage;
    }
  }
  std::istream& in = (opts.input == "-") ? std::cin : file_in;

  // gzip's first byte. A fragment BED begins with a contig name or a comment
  // and never with 0x1f, so one peeked byte tells the two apart on a pipe as
  // well as on a file.
  const bool looks_like_bam = in.peek() == 0x1f;

  ContigTable contigs;
  std::int64_t accepted = 0;

  if (looks_like_bam) {
    if (!BamStreamer::available()) {
      // Same hint and same exit code as the BGZF build gives a gzipped
      // fragment BED below. Without BGZF the input cannot be inflated to tell
      // a BAM from bgzipped text, and the commoner of the two -- 10x
      // fragments.tsv.gz -- needs `zcat`, not a rebuild. This is also the
      // build a bare cluster produces, so it is where that input arrives most.
      // One exit code per input regardless of build flags is what
      // `job_specs.RETRYABLE_EXIT_CODES` and Nextflow's retry rule assume.
      std::fprintf(stderr,
                   "pto-peaks: the input is gzip-compressed, and this build has no "
                   "BGZF support. If it is a compressed fragment BED (such as 10x "
                   "fragments.tsv.gz), decompress it first:\n"
                   "    zcat fragments.tsv.gz | pto-peaks --chrom-sizes FILE\n"
                   "If it is a BAM, rebuild with libdeflate or zlib installed.\n");
      return kExitInvalidInput;
    }
    BamStreamer bam;
    // The contig table is filled from the header, so the driver cannot be
    // built until after open(). Both live to the end of the function.
    BamFilter probe;
    probe.min_mapq = opts.min_mapq;
    const BamStatus opened = bam.open(in, probe);
    if (opened != BamStatus::kOk) {
      std::fprintf(stderr, "pto-peaks: %s\n", BamStreamer::describe(opened));
      // The gzip magic says "compressed", not "BAM". The commonest compressed
      // input to a peak caller is not a BAM at all but the fragment BED every
      // 10x scATAC run ships as fragments.tsv.gz -- bgzipped text, which
      // decodes cleanly and then fails the BAM magic. Saying only "bad magic"
      // sends the operator looking for a corrupt BAM they never had.
      if (opened == BamStatus::kBadMagic || opened == BamStatus::kBadBlock) {
        std::fprintf(stderr,
                     "pto-peaks: the input is gzip-compressed. If it is a compressed "
                     "fragment BED (such as 10x fragments.tsv.gz) rather than a BAM, "
                     "decompress it first:\n"
                     "    zcat fragments.tsv.gz | pto-peaks --chrom-sizes FILE\n");
      }
      return kExitInvalidInput;
    }
    contigs.adopt(bam.contig_names(), bam.contig_lengths());
    if (contigs.names.empty()) {
      std::fprintf(stderr, "pto-peaks: the BAM header names no contigs\n");
      return kExitInvalidInput;
    }
    Driver driver(opts, contigs);
    const int rc = stream_bam_records(opts, contigs, driver, accepted, bam);
    if (rc != kExitOk) return rc;
    return emit(opts, contigs, driver, accepted);
  }

  if (opts.chrom_sizes.empty()) {
    std::fprintf(stderr,
                 "pto-peaks: --chrom-sizes is required for BED input (a BAM "
                 "carries its own contig lengths)\n");
    return kExitUsage;
  }
  std::string error;
  if (!contigs.load(opts.chrom_sizes, error)) {
    std::fprintf(stderr, "pto-peaks: %s\n", error.c_str());
    return kExitUsage;
  }
  Driver driver(opts, contigs);
  const int rc = run_bed(opts, in, contigs, driver, accepted);
  if (rc != kExitOk) return rc;
  return emit(opts, contigs, driver, accepted);
}

}  // namespace

int main(int argc, char** argv) {
  Options opts;
  std::vector<std::string> positional;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto next = [&](const char* what) -> const char* {
      if (i + 1 >= argc) {
        std::fprintf(stderr, "pto-peaks: %s needs a value\n", what);
        std::exit(kExitUsage);
      }
      return argv[++i];
    };
    if (arg == "-h" || arg == "--help") {
      usage(stdout);
      return kExitOk;
    } else if (arg == "--version") {
      std::printf("pto-peaks %s\n", PEAKS_VERSION);
      return kExitOk;
    } else if (arg == "--chrom-sizes") {
      opts.chrom_sizes = next("--chrom-sizes");
    } else if (arg == "-o" || arg == "--out") {
      opts.output = next("--out");
    } else if (arg == "--cutoff") {
      if (!parse_finite(next("--cutoff"), opts.call.score_cutoff)) {
        std::fprintf(stderr, "pto-peaks: --cutoff must be a finite number\n");
        return kExitUsage;
      }
    } else if (arg == "--qvalue") {
      if (!parse_finite(next("--qvalue"), opts.call.qvalue_cutoff)) {
        std::fprintf(stderr, "pto-peaks: --qvalue must be a number\n");
        return kExitUsage;
      }
    } else if (arg == "--lambda-bg") {
      if (!parse_finite(next("--lambda-bg"), opts.lambda_bg) || opts.lambda_bg < 0.0) {
        std::fprintf(stderr,
                     "pto-peaks: --lambda-bg must be a finite, non-negative number\n");
        return kExitUsage;
      }
    } else if (arg == "--min-mapq") {
      std::int64_t v = 0;
      if (!parse_i64(next("--min-mapq"), v) || v < 0 || v > 255) {
        std::fprintf(stderr, "pto-peaks: --min-mapq must be in [0, 255]\n");
        return kExitUsage;
      }
      opts.min_mapq = static_cast<int>(v);
    } else if (arg == "--min-length" || arg == "--max-gap" ||
               arg == "--extend-peaks") {
      std::int64_t v = 0;
      if (!parse_i64(next(arg.c_str()), v) || v < 0 ||
          v > std::numeric_limits<Coord>::max()) {
        std::fprintf(stderr, "pto-peaks: %s must be a non-negative integer\n",
                     arg.c_str());
        return kExitUsage;
      }
      if (arg == "--min-length") {
        opts.call.min_length = static_cast<Coord>(v);
      } else if (arg == "--max-gap") {
        opts.call.max_gap = static_cast<Coord>(v);
      } else {
        opts.call.extend_peaks = static_cast<Coord>(v);
      }
    } else if (!arg.empty() && arg[0] == '-' && arg != "-") {
      std::fprintf(stderr, "pto-peaks: unknown option %s\n", arg.c_str());
      usage(stderr);
      return kExitUsage;
    } else {
      positional.push_back(arg);
    }
  }

  if (positional.size() > 1) {
    std::fprintf(stderr, "pto-peaks: at most one input file\n");
    return kExitUsage;
  }
  if (!positional.empty()) opts.input = positional[0];

  try {
    opts.call.validate();
  } catch (const std::exception& e) {
    std::fprintf(stderr, "pto-peaks: %s\n", e.what());
    return kExitUsage;
  }

  // The output is opened only after the last record is read, so this was never
  // a truncate-before-read. It was a replace-after-read: `-o x.bed x.bed`,
  // `-o hg38.chrom.sizes`, and `-o x.bed < x.bed` each exited 0 with the input
  // overwritten by the narrowPeak computed from it. Refused up front, before
  // an hour of streaming is spent producing the file that would do it.
  if (const auto clash =
          pto::peaks::output_overwrites_input(opts.output, {opts.input, opts.chrom_sizes})) {
    std::fprintf(stderr,
                 "pto-peaks: output '%s' is the same file as input '%s'; writing it "
                 "would destroy the input\n",
                 opts.output.c_str(), clash->c_str());
    return kExitUsage;
  }

  std::ios::sync_with_stdio(false);
  // Nothing on the streaming path is expected to throw, but an exception that
  // reaches main is std::terminate -- SIGABRT and a core file where the job
  // contract expects an exit code.
  try {
    return run(opts);
  } catch (const std::bad_alloc&) {
    std::fprintf(stderr, "pto-peaks: out of memory\n");
    return kExitUsage;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "pto-peaks: %s\n", e.what());
    return kExitUsage;
  }
}
