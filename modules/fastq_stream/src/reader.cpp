// SPDX-License-Identifier: MIT
#include "fastq_stream/reader.hpp"

#include <fcntl.h>
#include <unistd.h>

#include <cerrno>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include <libdeflate.h>
#include <zlib.h>

namespace fq {

const char* format_name(GzipFormat f) noexcept {
  switch (f) {
    case GzipFormat::kRaw:
      return "uncompressed";
    case GzipFormat::kPlain:
      return "gzip (single-stream, serial inflate)";
    case GzipFormat::kBgzf:
      return "BGZF (block-parallel inflate)";
  }
  return "unknown";
}

namespace {

constexpr std::size_t kStageBytes = 1 << 20;

// Length of the BGZF member starting at p, or 0 if p is not a BGZF header.
std::size_t bgzf_member_size(const uint8_t* p, std::size_t n) noexcept {
  if (n < 18) return 0;
  if (p[0] != 31 || p[1] != 139 || p[2] != 8 || (p[3] & 4) == 0) return 0;
  const std::size_t xlen =
      static_cast<std::size_t>(p[10]) | (static_cast<std::size_t>(p[11]) << 8);
  const std::size_t extra_end = 12 + xlen;
  if (extra_end > n) return 0;
  std::size_t i = 12;
  while (i + 4 <= extra_end) {
    const uint8_t si1 = p[i];
    const uint8_t si2 = p[i + 1];
    const std::size_t slen =
        static_cast<std::size_t>(p[i + 2]) | (static_cast<std::size_t>(p[i + 3]) << 8);
    if (si1 == 'B' && si2 == 'C' && slen == 2 && i + 6 <= extra_end) {
      const std::size_t bsize =
          static_cast<std::size_t>(p[i + 4]) | (static_cast<std::size_t>(p[i + 5]) << 8);
      // A member cannot be shorter than its own header plus the 8-byte
      // CRC32/ISIZE trailer. BSIZE=0 claimed a 1-byte member, which the reader
      // then framed as a member and handed on.
      if (bsize + 1 < extra_end + 8) return 0;
      return bsize + 1;
    }
    i += 4 + slen;
  }
  return 0;
}

// Buffered byte source over a file descriptor.
class FdSource {
 public:
  explicit FdSource(const std::string& path) : stage_(kStageBytes) {
    if (path == "-") {
      fd_ = STDIN_FILENO;
      owns_ = false;
    } else {
      fd_ = ::open(path.c_str(), O_RDONLY);
      if (fd_ < 0) {
        throw std::runtime_error("cannot open '" + path + "': " + std::strerror(errno));
      }
      owns_ = true;
#if defined(POSIX_FADV_SEQUENTIAL)
      ::posix_fadvise(fd_, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
    }
  }

  ~FdSource() {
    if (owns_ && fd_ >= 0) ::close(fd_);
  }

  FdSource(const FdSource&) = delete;
  FdSource& operator=(const FdSource&) = delete;

  const uint8_t* ptr() const noexcept { return stage_.data() + pos_; }
  std::size_t avail() const noexcept { return fill_ - pos_; }
  void consume(std::size_t n) noexcept { pos_ += n; }
  uint64_t total() const noexcept { return total_; }

  // Guarantee at least n contiguous bytes, growing the stage if needed.
  // Returns the bytes actually available (< n only at end of file).
  std::size_t ensure(std::size_t n) {
    if (avail() >= n) return avail();
    compact();
    if (stage_.size() < n) stage_.resize(n * 2);
    while (fill_ < n && !eof_) read_more();
    return avail();
  }

  // Pull more bytes without a size guarantee. Returns false at EOF.
  bool fill_more() {
    if (avail() > 0) return true;
    compact();
    if (eof_) return false;
    read_more();
    return avail() > 0;
  }

 private:
  void compact() noexcept {
    if (pos_ == 0) return;
    const std::size_t rem = fill_ - pos_;
    if (rem > 0) std::memmove(stage_.data(), stage_.data() + pos_, rem);
    pos_ = 0;
    fill_ = rem;
  }

  void read_more() {
    if (fill_ == stage_.size()) stage_.resize(stage_.size() * 2);
    ssize_t got;
    do {
      got = ::read(fd_, stage_.data() + fill_, stage_.size() - fill_);
    } while (got < 0 && errno == EINTR);
    if (got < 0) throw std::runtime_error(std::string("read failed: ") + std::strerror(errno));
    if (got == 0) {
      eof_ = true;
      return;
    }
    fill_ += static_cast<std::size_t>(got);
    total_ += static_cast<uint64_t>(got);
  }

  int fd_ = -1;
  bool owns_ = false;
  std::vector<uint8_t> stage_;
  std::size_t pos_ = 0;
  std::size_t fill_ = 0;
  uint64_t total_ = 0;
  bool eof_ = false;
};

}  // namespace

// ---------------------------------------------------------------- GzipInput

struct GzipInput::Impl {
  // Parameter deliberately not named `pool`: it would shadow the member of the
  // same name (-Wshadow under GCC), and in a constructor body that reads the
  // member the shadowing is a live footgun rather than a style point.
  Impl(const std::string& path, BufferPool& buffer_pool) : src(path), pool(buffer_pool) {
    detect();
    if (format == GzipFormat::kPlain) {
      std::memset(&zs, 0, sizeof(zs));
      // 15 window bits + 16 => decode a gzip wrapper.
      if (inflateInit2(&zs, 15 + 16) != Z_OK) {
        throw std::runtime_error("inflateInit2 failed");
      }
      zs_ready = true;
    }
  }

  ~Impl() {
    if (zs_ready) inflateEnd(&zs);
  }

  void detect() {
    const std::size_t have = src.ensure(18);
    const uint8_t* p = src.ptr();
    if (have >= 2 && p[0] == 31 && p[1] == 139) {
      format = (bgzf_member_size(p, have) != 0) ? GzipFormat::kBgzf : GzipFormat::kPlain;
    } else {
      format = GzipFormat::kRaw;
    }
  }

  // End of a BGZF stream. Both of these used to pass silently:
  //
  //   * No end-of-file marker. A BGZF file truncated exactly between members
  //     -- a partial copy, an interrupted upload, a disk that filled -- ends
  //     on a complete member, so neither the member framing nor libdeflate's
  //     CRC notices. Whenever the blocks happen to end on record boundaries
  //     (htslib's writer flushes per record once a block is nearly full) the
  //     assembler has no partial record to trip over either. Measured: a
  //     2,000-read file cut in half reported 1,000 reads and exited 0. The
  //     28-byte empty member bgzip and htslib always write is the only
  //     evidence the stream is whole, which is why `samtools quickcheck` and
  //     modules/peaks both require it. The flag is ASSIGNED per member, never
  //     latched: an empty member is legal mid-stream (concatenated BGZF), and
  //     the question is whether the LAST one was empty.
  //   * 1..17 trailing bytes, too short to be a member header, were dropped
  //     by `if (have < 18) break`.
  void finish_bgzf(std::size_t have) const {
    if (have != 0) {
      throw std::runtime_error(std::to_string(have) +
                               " trailing bytes after the last BGZF member: the input "
                               "is truncated or has data appended");
    }
    if (!last_member_empty) {
      throw std::runtime_error(
          "BGZF input has no end-of-file marker: the file is truncated (a complete "
          "BGZF file ends with the 28-byte empty block bgzip writes)");
    }
  }

  bool next_bgzf(Buffer& out) {
    out = pool.acquire();
    std::size_t size = 0;
    std::size_t members = 0;
    while (members < kBgzfBatchMembers) {
      const std::size_t have = src.ensure(18);
      if (have < 18) {  // end of input
        finish_bgzf(have);
        break;
      }
      const std::size_t msize = bgzf_member_size(src.ptr(), have);
      if (msize == 0) {
        throw std::runtime_error("malformed BGZF member header at offset " +
                                 std::to_string(src.total() - have));
      }
      if (src.ensure(msize) < msize) {
        throw std::runtime_error("truncated BGZF member");
      }
      if (size + msize > out.capacity()) break;
      const uint8_t* m = src.ptr();
      // ISIZE, the member's uncompressed length, is its last four bytes.
      last_member_empty = (m[msize - 4] | m[msize - 3] | m[msize - 2] | m[msize - 1]) == 0;
      std::memcpy(out.data() + size, m, msize);
      size += msize;
      src.consume(msize);
      ++members;
    }
    if (size == 0) return false;
    out.set_size(size);
    out.seq = seq++;
    return true;
  }

  bool next_plain(Buffer& out) {
    if (done) return false;
    out = pool.acquire();
    const std::size_t cap = out.capacity();
    zs.next_out = reinterpret_cast<Bytef*>(out.data());
    zs.avail_out = static_cast<uInt>(cap);

    while (zs.avail_out > 0) {
      if (!src.fill_more()) {
        // End of input inside a member. Previously this was `done = true`,
        // so a gzip missing its 8-byte CRC32/ISIZE trailer decoded "fine"
        // with the CRC never checked, and a file cut to its first 2 or 10
        // bytes reported zero reads and exited 0. gzip -t fails all three.
        if (member_open) {
          throw std::runtime_error(
              "gzip input is truncated: the stream ends inside a compressed member "
              "(data or the CRC32/ISIZE trailer is missing)");
        }
        done = true;
        break;
      }
      zs.next_in = const_cast<Bytef*>(src.ptr());
      const uInt in_before = static_cast<uInt>(src.avail());
      zs.avail_in = in_before;

      member_open = true;  // input exists past the last member boundary
      const int rc = inflate(&zs, Z_NO_FLUSH);
      src.consume(in_before - zs.avail_in);

      if (rc == Z_STREAM_END) {
        member_open = false;
        // Concatenated gzip members: reset and keep going.
        if (inflateReset(&zs) != Z_OK) throw std::runtime_error("inflateReset failed");
        if (src.avail() == 0 && !src.fill_more()) {
          done = true;
          break;
        }
        continue;
      }
      if (rc == Z_BUF_ERROR) continue;  // needs more input or more output room
      if (rc != Z_OK) {
        throw std::runtime_error(std::string("gzip inflate error: ") +
                                 (zs.msg ? zs.msg : std::to_string(rc).c_str()));
      }
    }

    const std::size_t produced = cap - zs.avail_out;
    if (produced == 0) return false;
    out.set_size(produced);
    out.seq = seq++;
    return true;
  }

  bool next_raw(Buffer& out) {
    if (!src.fill_more()) return false;
    out = pool.acquire();
    const std::size_t n = std::min(src.avail(), out.capacity());
    std::memcpy(out.data(), src.ptr(), n);
    src.consume(n);
    out.set_size(n);
    out.seq = seq++;
    return true;
  }

  FdSource src;
  BufferPool& pool;
  GzipFormat format = GzipFormat::kRaw;
  z_stream zs{};
  bool zs_ready = false;
  bool done = false;
  bool member_open = true;         // kPlain: detect() saw a gzip header
  bool last_member_empty = false;  // kBgzf: see finish_bgzf()
  uint64_t seq = 0;
};

GzipInput::GzipInput(const std::string& path, BufferPool& pool)
    : impl_(std::make_unique<Impl>(path, pool)) {}

GzipInput::~GzipInput() = default;

GzipFormat GzipInput::format() const noexcept { return impl_->format; }

uint64_t GzipInput::compressed_bytes() const noexcept { return impl_->src.total(); }

bool GzipInput::next(Buffer& out) {
  switch (impl_->format) {
    case GzipFormat::kBgzf:
      return impl_->next_bgzf(out);
    case GzipFormat::kPlain:
      return impl_->next_plain(out);
    case GzipFormat::kRaw:
      return impl_->next_raw(out);
  }
  return false;
}

// ------------------------------------------------------------ BgzfInflater

struct BgzfInflater::Impl {
  Impl() : d(libdeflate_alloc_decompressor()) {
    if (d == nullptr) throw std::runtime_error("libdeflate_alloc_decompressor failed");
  }
  ~Impl() { libdeflate_free_decompressor(d); }
  libdeflate_decompressor* d;
};

BgzfInflater::BgzfInflater() : impl_(std::make_unique<Impl>()) {}
BgzfInflater::~BgzfInflater() = default;

bool BgzfInflater::inflate(const Buffer& in, Buffer& out, std::string* err) {
  const uint8_t* p = in.bytes();
  std::size_t remaining = in.size();
  std::size_t produced = 0;

  while (remaining >= 18) {
    const std::size_t msize = bgzf_member_size(p, remaining);
    if (msize == 0 || msize > remaining) {
      if (err) *err = "malformed BGZF member in batch";
      return false;
    }
    const std::size_t room = out.capacity() - produced;
    std::size_t actual = 0;
    const auto rc = libdeflate_gzip_decompress(impl_->d, p, msize,
                                               out.data() + produced, room, &actual);
    if (rc != LIBDEFLATE_SUCCESS) {
      if (err) *err = "libdeflate_gzip_decompress failed (code " + std::to_string(rc) + ")";
      return false;
    }
    produced += actual;
    p += msize;
    remaining -= msize;
  }
  if (remaining != 0) {
    if (err) *err = "trailing bytes in BGZF batch";
    return false;
  }
  out.set_size(produced);
  out.seq = in.seq;
  return true;
}

}  // namespace fq
