// SPDX-License-Identifier: MIT
#include "fastq_stream/pipeline.hpp"
#include "fastq_stream/file_identity.hpp"

#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <signal.h>
#include <unistd.h>

#include <cerrno>
#include <chrono>
#include <cstring>
#include <mutex>
#include <optional>
#include <system_error>
#include <utility>

namespace fq {
namespace {

// ---------------------------------------------------------------- utilities

// Collects the first fatal error from any stage and trips the abort flag so
// every other thread unwinds instead of deadlocking on a full/empty ring.
class ErrorGate {
 public:
  void fail(std::string msg) {
    std::lock_guard<std::mutex> lk(mu_);
    if (message_.empty()) message_ = std::move(msg);
    abort_.store(true, std::memory_order_relaxed);
  }
  bool failed() const noexcept { return abort_.load(std::memory_order_relaxed); }
  std::atomic<bool>& flag() noexcept { return abort_; }
  const std::string& message() const noexcept { return message_; }

 private:
  std::atomic<bool> abort_{false};
  std::mutex mu_;
  std::string message_;
};

// `-i reads.fq -o reads.fq` used to open the output with O_TRUNC while the
// reader had consumed only its first 1 MiB stage: the run failed with "input
// ends mid-record" and left a zero-byte file where the input had been. The
// same happened through a hard link or a symlink. Checked before anything is
// opened for writing, because afterwards there is nothing left to protect.
// file_identity.hpp holds the semantics (links, "-", regular files only).
void refuse_output_aliasing_input(const Config& cfg) {
  if (const auto clash = output_overwrites_input(cfg.output, {cfg.input, cfg.reference_fasta})) {
    throw std::runtime_error("output '" + cfg.output + "' is the same file as input '" +
                             *clash +
                             "': writing it would destroy the input before it is read");
  }
}

int open_output(const Config& cfg) {
  if (cfg.output == "-") return STDOUT_FILENO;

  struct stat st{};
  const bool exists = ::stat(cfg.output.c_str(), &st) == 0;

  if (!exists && cfg.make_fifo) {
    if (::mkfifo(cfg.output.c_str(), 0666) != 0) {
      throw std::runtime_error("mkfifo('" + cfg.output + "') failed: " +
                               std::strerror(errno));
    }
    // Blocks until the downstream aligner opens the read end.
    const int fd = ::open(cfg.output.c_str(), O_WRONLY);
    if (fd < 0) {
      throw std::runtime_error("open FIFO '" + cfg.output + "' failed: " +
                               std::strerror(errno));
    }
    return fd;
  }

  const bool is_fifo = exists && S_ISFIFO(st.st_mode);
  const int flags = is_fifo ? O_WRONLY : (O_WRONLY | O_CREAT | O_TRUNC);
  const int fd = ::open(cfg.output.c_str(), flags, 0644);
  if (fd < 0) {
    throw std::runtime_error("open '" + cfg.output + "' failed: " +
                             std::strerror(errno));
  }
  return fd;
}

bool write_all(int fd, const std::byte* p, std::size_t n, std::string* err) {
  while (n > 0) {
    const ssize_t got = ::write(fd, p, n);
    if (got < 0) {
      if (errno == EINTR) continue;
      if (err) *err = std::string("write failed: ") + std::strerror(errno);
      return false;
    }
    p += got;
    n -= static_cast<std::size_t>(got);
  }
  return true;
}

// A read id made safe to quote in an error message: bounded, and with any
// non-printable byte replaced, since the id is exactly the field that may be
// binary garbage when this message fires.
std::string printable_id(const RecordView& r) {
  std::string s;
  const std::size_t n = std::min<std::size_t>(r.id_len, 64);
  for (std::size_t i = 0; i < n; ++i) {
    const unsigned char c = r.id[i];
    s += (c >= 0x20 && c < 0x7f) ? static_cast<char>(c) : '?';
  }
  if (r.id_len > n) s += "...";
  return s;
}

// Why the assembler could not carry a partial record into the next chunk. A
// file with classic Mac CR-only line endings has no '\n' at all and arrives
// here as one enormous record; the message used to tell that user their input
// "has reads longer than this build supports".
bool looks_cr_only(const std::vector<std::byte>& carry) {
  bool cr = false;
  bool lf = false;
  for (const std::byte b : carry) {
    cr = cr || b == std::byte{'\r'};
    lf = lf || b == std::byte{'\n'};
  }
  return cr && !lf;
}

constexpr const char* kCrOnlyMessage =
    "the input has CR bytes and no LF line endings: this looks like classic Mac "
    "(CR-only) line endings. Convert first, e.g. tr '\\r' '\\n'";

std::string carry_limit_message(const std::vector<std::byte>& carry) {
  if (looks_cr_only(carry)) return kCrOnlyMessage;
  return "a FASTQ record exceeds the " + std::to_string(kPrefix) +
         "-byte record limit (reads up to ~32 kbp); this build does not support "
         "longer records";
}

// ------------------------------------------------------------------- worker

enum class Outcome { kOk, kOverflow, kInvalidBytes, kInvalidHeader, kInvalidPhred64 };

class RecordProcessor {
 public:
  RecordProcessor(const Config& cfg, const AdapterTrimmer& trimmer)
      : cfg_(cfg), trimmer_(trimmer) {}

  // Applies trimming and filtering to one record, appending the surviving
  // record to `out`/`used`. kOverflow never happens in practice because the
  // output pool is sized to the input chunk plus its prefix.
  Outcome process(const RecordView& r, QcStats& st, CycleScratch& sc, std::byte* out,
                  std::size_t cap, std::size_t& used) {
    const std::size_t len = r.seq_len;
    if (!valid_header_line(r.id, r.id_len) || !valid_header_line(r.plus, r.plus_len)) {
      return Outcome::kInvalidHeader;
    }
    // Phred+64 is validated and rebased to +33 in place, before anything reads
    // a quality byte. The record points into this worker's own chunk buffer,
    // which is writable, so the const_cast removes nothing the storage had.
    // A byte below '@' is not Phred+64 at all (Solexa's negative scores are
    // the usual one) and is refused rather than wrapped.
    if (cfg_.phred_offset == 64) {
      auto* q = const_cast<uint8_t*>(r.qual);
      bool bad = false;
      for (std::size_t i = 0; i < len; ++i) {
        bad |= static_cast<uint8_t>(q[i] - 64u) > 62u;
        q[i] = static_cast<uint8_t>(q[i] - 31u);
      }
      if (bad) return Outcome::kInvalidPhred64;
    }
    const ReadScan scan = accumulate_read(st, sc, r.seq, r.qual, len);
    // A NUL, a control byte or a high-bit byte in SEQ or QUAL is corruption,
    // not a low-quality base. Clamping it to Q0 let it through as a read that
    // merely trimmed badly -- and a NUL in SEQ went on to the aligner as 'N'.
    // Filters fail closed: reject the run and say which read.
    if (!scan.valid) return Outcome::kInvalidBytes;

    std::size_t start = 0;
    std::size_t end = len;

    if (cfg_.front_q > 0) start = trim_front(r.qual, end, cfg_.front_q);
    if (cfg_.tail_q > 0) end = start + trim_back(r.qual + start, end - start, cfg_.tail_q);
    if (cfg_.window > 0) {
      end = start + trim_sliding_window(r.qual + start, end - start, cfg_.window,
                                        cfg_.window_mean);
    }
    // Counted once, over all three quality trims together. This used to be
    // `if (window) { kept < span } else { start/end moved }`, so with the
    // default --window 4 a read trimmed only by --trim-front/--trim-tail was
    // never counted at all.
    if (start != 0 || end != len) st.quality_trimmed += 1;

    if (cfg_.trim_adapters && !trimmer_.empty()) {
      const std::size_t span = end - start;
      const std::size_t kept = trimmer_.trim(r.seq + start, span);
      if (kept < span) {
        st.adapter_trimmed += 1;
        end = start + kept;
      }
    }

    const std::size_t final_len = end - start;

    if (final_len < cfg_.min_len) {
      st.dropped_short += 1;
      return Outcome::kOk;
    }
    if (cfg_.max_n_rate >= 0.0) {
      // An untrimmed read's N count is already known from the QC pass.
      const std::size_t ns =
          final_len == len ? scan.n_bases : count_n(r.seq + start, final_len);
      if (static_cast<double>(ns) > cfg_.max_n_rate * static_cast<double>(final_len)) {
        st.dropped_n += 1;
        return Outcome::kOk;
      }
    }
    if (cfg_.min_mean_q > 0.0 &&
        mean_phred(r.qual + start, final_len) < cfg_.min_mean_q) {
      st.dropped_short += 1;  // counted as a quality drop in the report
      return Outcome::kOk;
    }

    st.reads_out += 1;
    st.bases_out += final_len;

    if (cfg_.qc_only) return Outcome::kOk;

    const std::size_t need = r.id_len + 1 + final_len + 1 + r.plus_len + 1 + final_len + 1;
    if (used + need > cap) return Outcome::kOverflow;

    std::byte* w = out + used;
    auto put = [&w](const void* src, std::size_t n) {
      std::memcpy(w, src, n);
      w += n;
    };
    const char nl = '\n';
    put(r.id, r.id_len);
    put(&nl, 1);
    put(r.seq + start, final_len);
    put(&nl, 1);
    put(r.plus, r.plus_len);
    put(&nl, 1);
    put(r.qual + start, final_len);
    put(&nl, 1);
    used = static_cast<std::size_t>(w - out);
    return Outcome::kOk;
  }

 private:
  const Config& cfg_;
  const AdapterTrimmer& trimmer_;
};

}  // namespace

// ------------------------------------------------------------- run_pipeline

RunReport run_pipeline(const Config& cfg) {
  const auto t0 = std::chrono::steady_clock::now();

  // A closed downstream FIFO must surface as EPIPE, not a fatal signal.
  ::signal(SIGPIPE, SIG_IGN);

  if (!cfg.qc_only) refuse_output_aliasing_input(cfg);

  int hw = cfg.threads > 0 ? cfg.threads
                           : static_cast<int>(std::thread::hardware_concurrency());
  if (hw <= 0) hw = 4;

  BufferPool in_pool(kChunkBytes);
  // Output chunks can be as large as an input chunk plus a spliced-in record.
  BufferPool out_pool(kChunkBytes + kPrefix);

  // Rank sizing: the reader, assembler and writer each need a thread; what is
  // left is split between the decompression rank and the QC workers. Computed
  // BEFORE the source is constructed, because htslib's decode pool has to be
  // sized at open time.
  const int spare = std::max(1, hw - 3);
  int n_inflate = 0;
  int n_worker = 0;

  // Source selection. Both types expose next(Buffer&); rather than give them a
  // common base for a call made once per 256 KiB, hold whichever one is in play
  // and dispatch through a lambda. The reader thread below is the only consumer.
  std::optional<GzipInput> gz;
#if defined(FQ_ENABLE_HTS)
  std::optional<HtsInput> hts;
#endif

  if (cfg.hts_input) {
#if defined(FQ_ENABLE_HTS)
    // This module's block-parallel libdeflate path cannot be used on a BAM:
    // record framing lives inside the BGZF payload, so members must be inflated
    // by htslib itself and in order. The inflater budget therefore moves into
    // htslib's own decode pool rather than being left idle.
    const int n_decode = std::max(1, spare / 3);
    n_worker = std::max(1, spare - n_decode);

    HtsOptions ho;
    ho.layout = cfg.hts_layout;
    ho.reference_fasta = cfg.reference_fasta;
    ho.missing_qual = cfg.hts_missing_qual;
    ho.decode_threads = n_decode;
    hts.emplace(cfg.input, in_pool, ho);
#else
    throw std::runtime_error(
        "'" + cfg.input +
        "' is a BAM/CRAM file, but this binary was built without htslib. "
        "Reconfigure with -DFQ_ENABLE_HTS=ON, or convert it with "
        "`samtools fastq` first.");
#endif
  } else {
    gz.emplace(cfg.input, in_pool);
    // Only BGZF gets inflaters; for plain gzip the reader thread does the
    // (structurally serial) inflate itself.
    n_inflate = (gz->format() == GzipFormat::kBgzf) ? std::max(1, spare / 3) : 0;
    n_worker = std::max(1, spare - n_inflate);
  }

  const GzipFormat fmt = gz ? gz->format() : GzipFormat::kRaw;
  auto pull = [&](Buffer& b) {
#if defined(FQ_ENABLE_HTS)
    if (hts) return hts->next(b);
#endif
    return gz->next(b);
  };

  AdapterTrimmer trimmer;
  if (cfg.trim_adapters) {
    if (cfg.adapters.empty()) {
      trimmer.add_defaults();
    } else {
      for (const auto& a : cfg.adapters) trimmer.add(a);
    }
    trimmer.set_seed_len(cfg.adapter_seed);
    trimmer.set_min_overlap(cfg.adapter_min_overlap);
    trimmer.set_max_mismatch_rate(cfg.adapter_mismatch);
  }

  ErrorGate gate;
  std::atomic<bool>& abort = gate.flag();

  std::vector<std::unique_ptr<Lane>> inflate_in, inflate_out, work_in, work_out;
  for (int i = 0; i < n_inflate; ++i) {
    inflate_in.push_back(std::make_unique<Lane>());
    inflate_out.push_back(std::make_unique<Lane>());
  }
  for (int i = 0; i < n_worker; ++i) {
    work_in.push_back(std::make_unique<Lane>());
    work_out.push_back(std::make_unique<Lane>());
  }

  // Assembler input: the inflater outputs in BGZF mode, or a single lane fed
  // straight from the reader otherwise.
  Lane direct_lane;
  const bool parallel_inflate = n_inflate > 0;

  std::vector<QcStats> worker_stats(static_cast<std::size_t>(n_worker));
  const int fd_out = cfg.qc_only ? -1 : open_output(cfg);
  std::atomic<uint64_t> bytes_written{0};

  // --------------------------------------------------------------- reader
  auto reader = [&] {
    try {
      uint64_t k = 0;
      Buffer buf;
      while (pull(buf)) {
        if (abort.load(std::memory_order_relaxed)) break;
        if (parallel_inflate) {
          lane_push(*inflate_in[k % n_inflate], std::move(buf), abort);
        } else {
          lane_push(direct_lane, std::move(buf), abort);
        }
        ++k;
      }
    } catch (const std::exception& e) {
      gate.fail(e.what());
    }
    for (auto& l : inflate_in) l->closed.store(true, std::memory_order_release);
    direct_lane.closed.store(true, std::memory_order_release);
  };

  // ------------------------------------------------------------ inflaters
  auto inflater = [&](int i) {
    try {
      BgzfInflater bgzf;
      Buffer in;
      while (lane_pop(*inflate_in[i], in, abort)) {
        Buffer out = in_pool.acquire();
        std::string err;
        if (!bgzf.inflate(in, out, &err)) {
          gate.fail(err);
          break;
        }
        in = Buffer{};
        lane_push(*inflate_out[i], std::move(out), abort);
      }
    } catch (const std::exception& e) {
      gate.fail(e.what());
    }
    inflate_out[i]->closed.store(true, std::memory_order_release);
  };

  // ------------------------------------------------------------ assembler
  // Splices partial records across chunk boundaries and dispatches whole-record
  // chunks to the workers. Single-threaded by design: memchr runs at memory
  // bandwidth, an order of magnitude above the inflate rate it feeds on.
  auto assembler = [&] {
    try {
      std::vector<std::byte> carry;
      uint64_t src_k = 0;   // which upstream lane to read next
      uint64_t dst_k = 0;   // which worker lane to write next
      for (;;) {
        Buffer chunk;
        const bool got =
            parallel_inflate
                ? lane_pop(*inflate_out[src_k % n_inflate], chunk, abort)
                : lane_pop(direct_lane, chunk, abort);
        if (!got) break;
        ++src_k;

        if (!carry.empty()) {
          if (carry.size() > kPrefix) {
            gate.fail(carry_limit_message(carry));
            break;
          }
          std::memcpy(chunk.extend_front(carry.size()), carry.data(), carry.size());
          carry.clear();
        }

        const std::size_t cut = last_record_boundary(chunk.chars(), chunk.size());
        const std::size_t tail = chunk.size() - cut;
        if (tail > 0) {
          carry.assign(chunk.data() + cut, chunk.data() + chunk.size());
        }
        if (cut == 0) continue;  // no complete record yet; keep accumulating
        chunk.set_size(cut);

        lane_push(*work_in[dst_k % n_worker], std::move(chunk), abort);
        ++dst_k;
      }

      // End of stream with bytes left over. Two shapes are not corruption:
      //   * the last record has no final newline (an editor, `printf`, or
      //     `head -c` on a file that ended cleanly) -- exactly three newlines
      //     and a non-empty quality line. That used to fail the whole run.
      //   * trailing blank lines after the last record.
      // Anything else is a record cut short and still fails.
      if (!carry.empty() && !gate.failed()) {
        std::size_t newlines = 0;
        bool blank = true;
        for (const std::byte b : carry) {
          const auto c = static_cast<char>(b);
          newlines += (c == '\n');
          blank = blank && (c == '\n' || c == '\r');
        }
        if (blank) {
          // Nothing to process.
        } else if (newlines == 3 && static_cast<char>(carry.back()) != '\n') {
          Buffer last = in_pool.acquire();
          if (carry.size() + 1 > last.capacity()) {
            gate.fail("FASTQ record exceeds the " + std::to_string(last.capacity()) +
                      " byte chunk; input has reads longer than this build supports");
          } else {
            std::memcpy(last.data(), carry.data(), carry.size());
            last.data()[carry.size()] = std::byte{'\n'};
            last.set_size(carry.size() + 1);
            lane_push(*work_in[dst_k % n_worker], std::move(last), abort);
            ++dst_k;
          }
        } else if (looks_cr_only(carry)) {
          // A CR-only file smaller than a chunk never reaches the carry limit;
          // it arrives here, and "input ends mid-record" was no better a guess.
          gate.fail(kCrOnlyMessage);
        } else {
          gate.fail("input ends mid-record: " + std::to_string(carry.size()) +
                    " trailing bytes do not form a complete FASTQ record");
        }
      }
    } catch (const std::exception& e) {
      gate.fail(e.what());
    }
    for (auto& l : work_in) l->closed.store(true, std::memory_order_release);
  };

  // -------------------------------------------------------------- workers
  auto worker = [&](int i) {
    try {
      RecordProcessor proc(cfg, trimmer);
      QcStats& st = worker_stats[static_cast<std::size_t>(i)];
      CycleScratch sc;
      Buffer chunk;
      while (lane_pop(*work_in[i], chunk, abort)) {
        Buffer out = out_pool.acquire();
        std::size_t used = 0;
        RecordIterator it(chunk.bytes(), chunk.size());
        RecordView r;
        bool ok = true;
        while (it.next(r)) {
          if (!r.well_formed()) {
            gate.fail("malformed FASTQ record (expected '@id/seq/+/qual' with "
                      "matching sequence and quality lengths)");
            ok = false;
            break;
          }
          // The record limit, applied to every record rather than only to one
          // that happens to straddle a chunk. The assembler can carry at most
          // kPrefix bytes of a partial record, so an 80 KiB read used to be
          // accepted when it fell inside a chunk and refused when it straddled
          // one: the same file passed or failed with the bytes before it, and
          // raw, gzip and BGZF cut chunks at different offsets. Measured: one
          // 80,010-byte record at offset 0 or 100 kB, exit 0; at 192 kB, exit 1.
          // A record counted here with its final newline never exceeds what the
          // carry holds, so the two checks now agree.
          const auto record_bytes = static_cast<std::size_t>(r.qual + r.qual_len - r.id) + 1;
          if (record_bytes > kPrefix) {
            gate.fail("FASTQ record '" + printable_id(r) + "' is " +
                      std::to_string(record_bytes) + " bytes, over the " +
                      std::to_string(kPrefix) +
                      "-byte record limit (reads up to ~32 kbp); this build does not "
                      "support longer records");
            ok = false;
            break;
          }
          const Outcome o = proc.process(r, st, sc, out.data(), out.capacity(), used);
          if (o == Outcome::kInvalidHeader) {
            gate.fail("malformed FASTQ record '" + printable_id(r) +
                      "': header or '+' line contains a control or non-ASCII byte");
            ok = false;
            break;
          }
          if (o == Outcome::kInvalidPhred64) {
            gate.fail("FASTQ record '" + printable_id(r) +
                      "': a quality byte is outside '@'..'~', so this input is not "
                      "Phred+64 (--phred-offset 64)");
            ok = false;
            break;
          }
          if (o == Outcome::kInvalidBytes) {
            gate.fail("malformed FASTQ record '" + printable_id(r) +
                      "': sequence or quality line contains a byte outside "
                      "printable ASCII ('!'..'~')");
            ok = false;
            break;
          }
          if (o == Outcome::kOverflow) {
            gate.fail("internal error: output chunk overflow");
            ok = false;
            break;
          }
        }
        st.fold(sc);
        chunk = Buffer{};
        if (!ok) break;
        out.set_size(used);
        lane_push(*work_out[i], std::move(out), abort);
      }
    } catch (const std::exception& e) {
      gate.fail(e.what());
    }
    work_out[i]->closed.store(true, std::memory_order_release);
  };

  // --------------------------------------------------------------- writer
  auto writer = [&] {
    try {
      uint64_t k = 0;
      uint64_t written = 0;
      for (;;) {
        Buffer out;
        if (!lane_pop(*work_out[k % n_worker], out, abort)) break;
        ++k;
        if (!cfg.qc_only && out.size() > 0) {
          std::string err;
          if (!write_all(fd_out, out.data(), out.size(), &err)) {
            gate.fail(err);
            break;
          }
          written += out.size();
        }
      }
      bytes_written.store(written, std::memory_order_relaxed);
    } catch (const std::exception& e) {
      gate.fail(e.what());
    }
  };

  // Every stage is started through this one vector so that a failure to start
  // any of them unwinds the ones already running.
  //
  // std::thread's constructor throws std::system_error when the OS refuses a
  // thread -- a container pids limit, RLIMIT_NPROC on Linux (which counts
  // threads), or --threads 4096 against a per-process cap. Previously the
  // exception propagated past already-running, joinable std::thread objects,
  // whose destructors call std::terminate: SIGABRT and a core dump instead of
  // an error. Tripping the gate makes every running stage leave its ring loop,
  // and the joins below then complete normally.
  const std::size_t n_threads = static_cast<std::size_t>(3 + n_inflate + n_worker);
  std::vector<std::thread> threads;
  threads.reserve(n_threads);
  try {
    threads.emplace_back(reader);
    for (int i = 0; i < n_inflate; ++i) threads.emplace_back(inflater, i);
    threads.emplace_back(assembler);
    for (int i = 0; i < n_worker; ++i) threads.emplace_back(worker, i);
    threads.emplace_back(writer);
  } catch (const std::system_error& e) {
    gate.fail("cannot start pipeline thread " + std::to_string(threads.size() + 1) +
              " of " + std::to_string(n_threads) + ": " + e.what() +
              " (lower --threads)");
  }
  for (auto& t : threads) t.join();

  // close() is where NFS and some FUSE filesystems report a deferred write
  // error (EIO, ENOSPC, EDQUOT); ignoring it reports success over a short file.
  if (fd_out >= 0 && fd_out != STDOUT_FILENO && ::close(fd_out) != 0) {
    gate.fail("close '" + cfg.output + "' failed: " + std::strerror(errno));
  }

  if (gate.failed()) throw std::runtime_error(gate.message());

  RunReport rep;
  for (const auto& s : worker_stats) rep.stats.merge(s);
  rep.format = fmt;
#if defined(FQ_ENABLE_HTS)
  if (hts) {
    rep.compressed_bytes = hts->compressed_bytes();
    rep.hts_stats = hts->stats();
    rep.from_alignment = true;
    rep.alignment_format = hts->format_name();
  } else
#endif
  {
    rep.compressed_bytes = gz->compressed_bytes();
  }
  rep.bytes_written = bytes_written.load(std::memory_order_relaxed);
  rep.inflater_threads = n_inflate;
  rep.worker_threads = n_worker;
  rep.pool_blocks = in_pool.blocks_allocated() + out_pool.blocks_allocated();
  rep.pool_bytes = in_pool.bytes_allocated() + out_pool.bytes_allocated();
  rep.wall_seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return rep;
}

}  // namespace fq
