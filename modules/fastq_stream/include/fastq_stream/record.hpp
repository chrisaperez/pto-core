// SPDX-License-Identifier: MIT
// fastq_stream — zero-copy FASTQ record views.
//
// A FASTQ record is exactly four lines. '@' is a legal quality character, so
// record starts cannot be found by scanning for '@' — the only sound way to
// split a stream is to count newlines from a known record boundary. The stream
// begins at one, and the assembler stage preserves that invariant by only ever
// cutting chunks after a multiple of four lines.
#pragma once

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string_view>

namespace fq {

struct RecordView {
  const uint8_t* id = nullptr;
  std::size_t id_len = 0;
  const uint8_t* seq = nullptr;
  std::size_t seq_len = 0;
  const uint8_t* plus = nullptr;
  std::size_t plus_len = 0;
  const uint8_t* qual = nullptr;
  std::size_t qual_len = 0;

  bool well_formed() const noexcept {
    return id_len > 0 && id[0] == '@' && plus_len > 0 && plus[0] == '+' &&
           seq_len == qual_len;
  }
};

// Byte offset just past the last complete 4-line record in [p, p+n).
// Returns 0 when the range holds no complete record.
inline std::size_t last_record_boundary(const char* p, std::size_t n) noexcept {
  std::size_t boundary = 0;
  std::size_t lines = 0;
  const char* cur = p;
  std::size_t remaining = n;
  while (remaining > 0) {
    const void* nl = std::memchr(cur, '\n', remaining);
    if (nl == nullptr) break;
    const std::size_t consumed = static_cast<std::size_t>(
        static_cast<const char*>(nl) - cur + 1);
    cur += consumed;
    remaining -= consumed;
    if (++lines == 4) {
      boundary = static_cast<std::size_t>(cur - p);
      lines = 0;
    }
  }
  return boundary;
}

// Forward iterator over a buffer containing only whole records.
class RecordIterator {
 public:
  RecordIterator(const uint8_t* data, std::size_t size) noexcept
      : cur_(data), end_(data + size) {}

  bool next(RecordView& out) noexcept {
    if (cur_ >= end_) return false;
    if (!take_line(out.id, out.id_len)) return false;
    if (!take_line(out.seq, out.seq_len)) return false;
    if (!take_line(out.plus, out.plus_len)) return false;
    if (!take_line(out.qual, out.qual_len)) return false;
    return true;
  }

 private:
  bool take_line(const uint8_t*& base, std::size_t& len) noexcept {
    if (cur_ >= end_) return false;
    const auto remaining = static_cast<std::size_t>(end_ - cur_);
    const void* nl = std::memchr(cur_, '\n', remaining);
    if (nl == nullptr) return false;
    base = cur_;
    len = static_cast<std::size_t>(static_cast<const uint8_t*>(nl) - cur_);
    // Tolerate CRLF-terminated files without copying.
    if (len > 0 && base[len - 1] == '\r') --len;
    cur_ = static_cast<const uint8_t*>(nl) + 1;
    return true;
  }

  const uint8_t* cur_;
  const uint8_t* end_;
};

}  // namespace fq
