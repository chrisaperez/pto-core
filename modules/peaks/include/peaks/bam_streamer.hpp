// SPDX-License-Identifier: MIT
//
// pto-peaks: native BGZF/BAM ingestion.
//
// Why this exists rather than a call into htslib
// ----------------------------------------------
// `genomic_toolkit` already links htslib and already reads BAM. Reusing it is
// the obvious move and the repository forbids it: modules share no code, and a
// cross-module dependency here would be the shared layer README.md declares a
// non-goal. Linking htslib into THIS module instead would work, and it would
// cost the property the module was built for -- that it compiles on a cluster
// with a compiler and nothing else.
//
// So this reads BAM directly. It is not a general htslib replacement and does
// not try to be: it decodes coordinate-sorted alignment records far enough to
// recover fragment spans, and refuses everything it does not understand. No
// index, no region queries, no CRAM, no writing.
//
// The decompressor is optional and discovered at configure time: libdeflate if
// it is installed (roughly twice zlib's inflate throughput on BGZF blocks),
// zlib otherwise, and if neither is present the whole file compiles to nothing
// and the fragment-BED path still works.
//
// Zero-copy, and what that does and does not mean
// -----------------------------------------------
// Records are parsed in place, out of the decompressed block buffer. Nothing
// is allocated per record: no std::string for the read name, no object for the
// alignment, no vector growth. `test_streaming_allocates_nothing_per_record`
// counts global operator new to prove it.
//
// What it does NOT mean is `reinterpret_cast`. BAM fields are packed
// little-endian at whatever offset the previous record left off, so a record's
// int32 fields are routinely unaligned; casting a pointer to them is undefined
// behaviour, `-fsanitize=alignment` reports it, and on a strict-alignment
// target it faults. The byte-shift readers below are defined for any alignment
// AND for any host byte order, and clang and gcc both fold them into the same
// single unaligned load the cast would have produced. The cast would buy
// nothing and cost correctness.
//
// Fragments, not reads
// --------------------
// A paired-end fragment is emitted once, from the leftmost mate, using TLEN.
// Emitting both mates would double every pileup -- which reads as twice the
// sequencing depth and moves every p-value with it. Mates with TLEN <= 0 are
// the rightmost of their pair and are skipped by construction.
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <istream>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "peaks/sliding_window.hpp"

#if defined(PEAKS_HAVE_BGZF)
#if defined(PEAKS_BGZF_LIBDEFLATE)
#include <libdeflate.h>
#elif defined(PEAKS_BGZF_ZLIB)
#include <zlib.h>
#endif
#endif

namespace pto::peaks {

struct BamFragment {
  std::int32_t tid = kNoTid;
  Coord start = 0;
  Coord end = 0;
};

enum class BamStatus {
  kOk,
  kEndOfFile,
  kBadMagic,     // not a BAM at all
  kTruncated,    // the file stops in the middle of something, or has no EOF block
  kBadBlock,     // BGZF framing is wrong
  kBadChecksum,  // a block decompressed, and to the wrong bytes
  kBadHeader,    // the BAM header's own lengths do not add up
  kBadRecord,    // an alignment record's lengths do not add up
  kUnsupported,  // structurally valid, beyond what this reader handles
};

// SAM flag bits, named so the filter below reads as prose.
namespace samflags {
inline constexpr std::uint16_t kPaired = 0x1;
inline constexpr std::uint16_t kProperPair = 0x2;
inline constexpr std::uint16_t kUnmapped = 0x4;
inline constexpr std::uint16_t kMateUnmapped = 0x8;
inline constexpr std::uint16_t kSecondary = 0x100;
inline constexpr std::uint16_t kQcFail = 0x200;
inline constexpr std::uint16_t kDuplicate = 0x400;
inline constexpr std::uint16_t kSupplementary = 0x800;
}  // namespace samflags

struct BamFilter {
  int min_mapq = 0;
  // The conventional pipeline's exclusions. kDuplicate IS excluded here, unlike
  // genomic_toolkit's default -- that module marks duplicates itself and must
  // not honour an upstream tool's marks silently; this one does not, so an
  // already-marked duplicate is the caller's stated intent.
  std::uint16_t exclude_flags = samflags::kUnmapped | samflags::kMateUnmapped |
                                samflags::kSecondary | samflags::kQcFail |
                                samflags::kDuplicate | samflags::kSupplementary;
  std::uint16_t require_flags = samflags::kProperPair;
  Coord max_span = kMaxFragmentSpan;
};

#if !defined(PEAKS_HAVE_BGZF)

// Built without a decompressor. The type still exists so callers compile; every
// operation reports kUnsupported, and `available()` is how a CLI decides to say
// so before opening anything.
class BamStreamer {
 public:
  static constexpr bool available() noexcept { return false; }
  BamStatus open(std::istream&, const BamFilter&) { return BamStatus::kUnsupported; }
  BamStatus next(BamFragment&) { return BamStatus::kUnsupported; }
  [[nodiscard]] const std::vector<std::string>& contig_names() const noexcept { return names_; }
  [[nodiscard]] const std::vector<Coord>& contig_lengths() const noexcept { return lengths_; }
  [[nodiscard]] std::int64_t records_seen() const noexcept { return 0; }
  [[nodiscard]] std::int64_t fragments_emitted() const noexcept { return 0; }
  [[nodiscard]] std::int64_t skipped_flags() const noexcept { return 0; }
  [[nodiscard]] std::int64_t skipped_mapq() const noexcept { return 0; }
  [[nodiscard]] std::int64_t skipped_span() const noexcept { return 0; }
  [[nodiscard]] std::int64_t skipped_unpaired() const noexcept { return 0; }
  [[nodiscard]] std::int64_t skipped_unplaced() const noexcept { return 0; }
  static const char* describe(BamStatus) noexcept {
    return "this build has no BGZF support";
  }

 private:
  std::vector<std::string> names_;
  std::vector<Coord> lengths_;
};

#else

namespace detail {

// Little-endian loads that are defined for every alignment and every host byte
// order. Both compilers fold these into a single load; see the header comment
// for why the pointer cast is not used.
[[nodiscard]] inline std::uint16_t le16(const unsigned char* p) noexcept {
  return static_cast<std::uint16_t>(static_cast<std::uint16_t>(p[0]) |
                                    (static_cast<std::uint16_t>(p[1]) << 8));
}

[[nodiscard]] inline std::uint32_t le32(const unsigned char* p) noexcept {
  return static_cast<std::uint32_t>(p[0]) |
         (static_cast<std::uint32_t>(p[1]) << 8) |
         (static_cast<std::uint32_t>(p[2]) << 16) |
         (static_cast<std::uint32_t>(p[3]) << 24);
}

[[nodiscard]] inline std::int32_t lei32(const unsigned char* p) noexcept {
  return static_cast<std::int32_t>(le32(p));
}

// A BGZF block inflates to at most 64 KiB by definition of the format.
inline constexpr std::size_t kMaxBlockSize = 65536;
// Room for one whole block plus whatever partial record preceded it.
inline constexpr std::size_t kDataCapacity = 1u << 20;
inline constexpr std::size_t kMaxRecordSize = kDataCapacity - 2 * kMaxBlockSize;
// A BAM header naming a million contigs is already absurd; past this it is
// corruption being used to make the reader allocate.
inline constexpr std::int32_t kMaxContigs = 1 << 20;
inline constexpr std::int32_t kMaxContigNameLength = 1 << 16;

}  // namespace detail

class BamStreamer {
 public:
  static constexpr bool available() noexcept { return true; }

  BamStreamer() : raw_(detail::kMaxBlockSize), data_(detail::kDataCapacity) {
#if defined(PEAKS_BGZF_LIBDEFLATE)
    inflator_.reset(libdeflate_alloc_decompressor());
#elif defined(PEAKS_BGZF_ZLIB)
    auto z = std::make_unique<z_stream>();
    std::memset(z.get(), 0, sizeof(z_stream));
    // Ownership passes to ZStreamDeleter only on success, so a failed init
    // leaves zs_ null and `usable()` below answers for it. A separate
    // `zs_ok_` flag cannot do that job: it is trivially copied by the
    // defaulted move, which leaves a moved-from object claiming a stream it
    // no longer owns.
    if (inflateInit2(z.get(), -15) == Z_OK) {  // raw DEFLATE, no wrapper
      zs_.reset(z.release());
    }
#endif
  }

  BamStreamer(const BamStreamer&) = delete;
  BamStreamer& operator=(const BamStreamer&) = delete;
  // Defaulted, and safe to default only because every owning member cleans
  // itself up: the decompressor through its deleter, the buffers through
  // vector. There is deliberately no destructor -- an `inflateEnd` written
  // there instead would be skipped by move ASSIGNMENT, which does not run one,
  // and zlib's inflate_state (7 KiB) would leak per assignment.
  BamStreamer(BamStreamer&&) noexcept = default;
  BamStreamer& operator=(BamStreamer&&) noexcept = default;

  // Reads the header. The stream must stay alive until the last next() call.
  BamStatus open(std::istream& in, const BamFilter& filter) {
    if (!usable()) return fail(BamStatus::kUnsupported);
    in_ = &in;
    filter_ = filter;
    return read_header();
  }

  [[nodiscard]] const std::vector<std::string>& contig_names() const noexcept {
    return names_;
  }
  [[nodiscard]] const std::vector<Coord>& contig_lengths() const noexcept {
    return lengths_;
  }

  // Next fragment, or kEndOfFile / an error. Nothing is allocated here.
  BamStatus next(BamFragment& out) {
    if (failed_ != BamStatus::kOk) return failed_;
    // open() was never called, or this object was moved from. Both leave
    // members that every path below assumes are live -- in_, and the two
    // buffers, whose data() is null once a move has stolen them.
    if (in_ == nullptr || !usable()) return fail(BamStatus::kUnsupported);
    for (;;) {
      const unsigned char* rec = nullptr;
      std::size_t rec_len = 0;
      const BamStatus st = next_record(rec, rec_len);
      if (st != BamStatus::kOk) return fail(st);
      ++records_seen_;

      const std::int32_t refid = detail::lei32(rec + 0);
      const std::int32_t pos = detail::lei32(rec + 4);
      const std::uint8_t l_read_name = rec[8];
      const std::uint8_t mapq = rec[9];
      const std::uint16_t n_cigar = detail::le16(rec + 12);
      const std::uint16_t flag = detail::le16(rec + 14);
      const std::int32_t l_seq = detail::lei32(rec + 16);
      const std::int32_t tlen = detail::lei32(rec + 28);

      // The record's own lengths must add up inside the record. A corrupt
      // n_cigar or l_seq is how a truncated download turns into a confident
      // read of whatever followed it.
      if (l_seq < 0 || l_read_name == 0) return fail(BamStatus::kBadRecord);
      const std::uint64_t declared =
          32ull + l_read_name + 4ull * n_cigar +
          (static_cast<std::uint64_t>(l_seq) + 1ull) / 2ull + static_cast<std::uint64_t>(l_seq);
      if (declared > rec_len) return fail(BamStatus::kBadRecord);

      // Unplaced reads sit at the end of a coordinate-sorted BAM. Legal, and
      // not a fragment.
      if (refid == kNoTid) {
        ++skipped_unplaced_;
        continue;
      }
      if (refid < 0 || static_cast<std::size_t>(refid) >= lengths_.size()) {
        return fail(BamStatus::kBadRecord);
      }
      if (pos < 0) return fail(BamStatus::kBadRecord);

      if ((flag & filter_.exclude_flags) != 0 ||
          (flag & filter_.require_flags) != filter_.require_flags) {
        ++skipped_flags_;
        continue;
      }
      if (mapq < filter_.min_mapq) {
        ++skipped_mapq_;
        continue;
      }
      // Only the leftmost mate carries a positive TLEN, so this is also what
      // stops a pair being counted twice. Counted separately from the flag
      // rejections: a single-end BAM has TLEN 0 on every record, and folding
      // that into skipped_flags_ tells the operator to go and look at their
      // flag filter, which is not where the problem is.
      if (tlen <= 0) {
        ++skipped_unpaired_;
        continue;
      }
      if (tlen > filter_.max_span) {
        ++skipped_span_;
        continue;
      }
      const std::int64_t end = static_cast<std::int64_t>(pos) + tlen;
      if (end > lengths_[static_cast<std::size_t>(refid)]) {
        ++skipped_span_;
        continue;
      }

      out.tid = refid;
      out.start = static_cast<Coord>(pos);
      out.end = static_cast<Coord>(end);
      ++fragments_emitted_;
      return BamStatus::kOk;
    }
  }

  [[nodiscard]] std::int64_t records_seen() const noexcept { return records_seen_; }
  [[nodiscard]] std::int64_t fragments_emitted() const noexcept { return fragments_emitted_; }
  [[nodiscard]] std::int64_t skipped_flags() const noexcept { return skipped_flags_; }
  [[nodiscard]] std::int64_t skipped_mapq() const noexcept { return skipped_mapq_; }
  [[nodiscard]] std::int64_t skipped_span() const noexcept { return skipped_span_; }
  [[nodiscard]] std::int64_t skipped_unpaired() const noexcept { return skipped_unpaired_; }
  [[nodiscard]] std::int64_t skipped_unplaced() const noexcept { return skipped_unplaced_; }

  static const char* describe(BamStatus st) noexcept {
    switch (st) {
      case BamStatus::kOk: return "ok";
      case BamStatus::kEndOfFile: return "end of file";
      case BamStatus::kBadMagic: return "not a BAM file (bad magic)";
      case BamStatus::kTruncated: return "file is truncated (no BGZF end-of-file block)";
      case BamStatus::kBadBlock: return "malformed BGZF block";
      case BamStatus::kBadChecksum: return "BGZF block failed its CRC32";
      case BamStatus::kBadHeader: return "malformed BAM header";
      case BamStatus::kBadRecord: return "malformed alignment record";
      case BamStatus::kUnsupported: return "unsupported BAM feature";
    }
    return "unknown status";
  }

 private:
  BamStatus fail(BamStatus st) noexcept {
    if (st != BamStatus::kOk) failed_ = st;
    return st;
  }

  // Every path below indexes raw_ and data_ at their constructed sizes and
  // dereferences the decompressor without checking. This is the one place
  // that checks, and it covers all three ways those assumptions can be false:
  // a decompressor that failed to allocate, a moved-from object (whose
  // vectors are empty and whose data() is therefore null), and an object
  // still being constructed.
  [[nodiscard]] bool usable() const noexcept {
    if (raw_.size() != detail::kMaxBlockSize) return false;
    if (data_.size() != detail::kDataCapacity) return false;
#if defined(PEAKS_BGZF_LIBDEFLATE)
    return inflator_ != nullptr;
#elif defined(PEAKS_BGZF_ZLIB)
    return zs_ != nullptr;
#else
    return false;
#endif
  }

  // --- BGZF ---------------------------------------------------------------

  [[nodiscard]] std::size_t available_bytes() const noexcept { return data_len_ - data_pos_; }

  // Make `n` contiguous decompressed bytes available, inflating blocks as
  // needed. False means the stream ran out first.
  bool ensure(std::size_t n) {
    if (n > detail::kDataCapacity) return false;
    while (available_bytes() < n) {
      if (data_pos_ > 0) {
        // Compact: the unconsumed tail moves to the front so a record that
        // straddles blocks becomes contiguous.
        std::memmove(data_.data(), data_.data() + data_pos_, available_bytes());
        data_len_ -= data_pos_;
        data_pos_ = 0;
      }
      if (data_len_ + detail::kMaxBlockSize > detail::kDataCapacity) return false;
      if (!inflate_one_block()) return false;
    }
    return true;
  }

  const unsigned char* take(std::size_t n) noexcept {
    const unsigned char* p = data_.data() + data_pos_;
    data_pos_ += n;
    return p;
  }

  bool skip_bytes(std::int64_t n) {
    while (n > 0) {
      if (available_bytes() == 0 && !ensure(1)) return false;
      const std::size_t chunk =
          static_cast<std::size_t>(std::min<std::int64_t>(n, static_cast<std::int64_t>(available_bytes())));
      data_pos_ += chunk;
      n -= static_cast<std::int64_t>(chunk);
    }
    return true;
  }

  // Reads and inflates one BGZF block onto the end of data_. Sets block_error_
  // when the block itself is malformed, which ensure()'s callers turn into the
  // status they report.
  bool inflate_one_block() {
    if (stream_done_) return false;
    unsigned char head[18];
    in_->read(reinterpret_cast<char*>(head), 18);
    if (in_->gcount() == 0) {
      stream_done_ = true;
      return false;
    }
    if (in_->gcount() != 18) {
      block_error_ = BamStatus::kTruncated;
      stream_done_ = true;
      return false;
    }
    // FLG must be FEXTRA and nothing else. FHCRC, FNAME and FCOMMENT each
    // insert bytes between the extra field and the deflate payload, so
    // tolerating them means handing those bytes to the inflater and blaming
    // the result on the compressed data.
    if (head[0] != 0x1f || head[1] != 0x8b || head[2] != 0x08 ||
        (head[3] & 0x04) == 0 || (head[3] & 0x1a) != 0) {
      block_error_ = BamStatus::kBadBlock;
      stream_done_ = true;
      return false;
    }
    const std::uint16_t xlen = detail::le16(head + 10);
    // The BC extra subfield must be first, which is what the format
    // guarantees and what every writer emits.
    if (xlen < 6 || head[12] != 'B' || head[13] != 'C' || detail::le16(head + 14) != 2) {
      block_error_ = BamStatus::kBadBlock;
      stream_done_ = true;
      return false;
    }
    const std::size_t bsize = static_cast<std::size_t>(detail::le16(head + 16)) + 1;
    if (bsize < 18u + xlen - 6u + 8u || bsize > detail::kMaxBlockSize) {
      block_error_ = BamStatus::kBadBlock;
      stream_done_ = true;
      return false;
    }
    // Any extra subfields past BC, then the deflate payload and the trailer.
    const std::size_t rest = bsize - 18;
    in_->read(reinterpret_cast<char*>(raw_.data()), static_cast<std::streamsize>(rest));
    if (static_cast<std::size_t>(in_->gcount()) != rest) {
      block_error_ = BamStatus::kTruncated;
      stream_done_ = true;
      return false;
    }
    const std::size_t skip_extra = static_cast<std::size_t>(xlen) - 6u;
    if (rest < skip_extra + 8) {
      block_error_ = BamStatus::kBadBlock;
      stream_done_ = true;
      return false;
    }
    const unsigned char* cdata = raw_.data() + skip_extra;
    const std::size_t clen = rest - skip_extra - 8;
    const std::uint32_t want_crc = detail::le32(raw_.data() + rest - 8);
    const std::uint32_t isize = detail::le32(raw_.data() + rest - 4);
    if (isize > detail::kMaxBlockSize) {
      block_error_ = BamStatus::kBadBlock;
      stream_done_ = true;
      return false;
    }

    unsigned char* dst = data_.data() + data_len_;
    std::size_t produced = 0;
    if (isize > 0) {
      if (!raw_inflate(cdata, clen, dst, isize, produced) || produced != isize) {
        block_error_ = BamStatus::kBadBlock;
        stream_done_ = true;
        return false;
      }
      if (crc32_of(dst, produced) != want_crc) {
        block_error_ = BamStatus::kBadChecksum;
        stream_done_ = true;
        return false;
      }
    }
    data_len_ += produced;
    // The end-of-file marker is a block with an empty payload. Its presence is
    // the only evidence a BAM was not cut short at a block boundary, which is
    // what `samtools quickcheck` looks for and what a failed download looks
    // like: every record decodes, and the last chromosome is missing.
    //
    // Two things here are deliberate and were each broken once by a change
    // that looked like a hardening fix:
    //
    //   * an empty block does NOT end the stream. Empty blocks are legal
    //     anywhere -- the marker itself is one, and concatenating BGZF streams
    //     puts them in the middle -- so stopping at the first one drops every
    //     record after it AND reports success. Reading on costs one pass over
    //     bytes that are already in the file: 200k empty blocks decode in 6 ms,
    //     because each iteration consumes 18 bytes of input and the loop ends
    //     at EOF like any other. There is no unbounded loop to guard against.
    //   * the flag is ASSIGNED, not latched. It answers "was the LAST block
    //     empty", which is the question truncation detection asks. Making it
    //     sticky means a file with an empty block somewhere in the middle
    //     reports a clean end however abruptly it actually stops.
    //
    // `test_an_empty_block_mid_stream_is_transparent` pins both.
    saw_eof_marker_ = (isize == 0);
    return produced > 0 || isize == 0;
  }

  bool raw_inflate(const unsigned char* in, std::size_t in_len, unsigned char* out,
                   std::size_t out_cap, std::size_t& produced) {
#if defined(PEAKS_BGZF_LIBDEFLATE)
    if (!inflator_) return false;
    std::size_t actual = 0;
    const auto rc = libdeflate_deflate_decompress(inflator_.get(), in, in_len, out,
                                                  out_cap, &actual);
    if (rc != LIBDEFLATE_SUCCESS) return false;
    produced = actual;
    return true;
#elif defined(PEAKS_BGZF_ZLIB)
    if (!zs_) return false;
    inflateReset(zs_.get());
    zs_->next_in = const_cast<Bytef*>(in);
    zs_->avail_in = static_cast<uInt>(in_len);
    zs_->next_out = out;
    zs_->avail_out = static_cast<uInt>(out_cap);
    const int rc = inflate(zs_.get(), Z_FINISH);
    if (rc != Z_STREAM_END) return false;
    produced = out_cap - zs_->avail_out;
    return true;
#else
    (void)in; (void)in_len; (void)out; (void)out_cap; (void)produced;
    return false;
#endif
  }

  static std::uint32_t crc32_of(const unsigned char* p, std::size_t n) noexcept {
#if defined(PEAKS_BGZF_LIBDEFLATE)
    return static_cast<std::uint32_t>(libdeflate_crc32(0, p, n));
#elif defined(PEAKS_BGZF_ZLIB)
    return static_cast<std::uint32_t>(crc32(0, p, static_cast<uInt>(n)));
#else
    (void)p; (void)n;
    return 0;
#endif
  }

  // Whatever ensure() ran out on: a malformed block if there was one, a clean
  // end if the file ended tidily, truncation otherwise.
  [[nodiscard]] BamStatus exhausted_status() const noexcept {
    if (block_error_ != BamStatus::kOk) return block_error_;
    if (!saw_eof_marker_) return BamStatus::kTruncated;
    return BamStatus::kEndOfFile;
  }

  // --- BAM ----------------------------------------------------------------

  // A header read that ran out of stream. If the stream itself ended tidily
  // the header is what is wrong -- its own lengths do not fit the file it is
  // in. If it did not, the file is what is wrong, and answering "malformed
  // BAM header" to a failed download sends the operator to the wrong problem.
  [[nodiscard]] BamStatus header_exhausted() const noexcept {
    const BamStatus st = exhausted_status();
    return st == BamStatus::kEndOfFile ? BamStatus::kBadHeader : st;
  }

  BamStatus read_header() {
    if (!ensure(4)) {
      const BamStatus st = exhausted_status();
      return fail(st == BamStatus::kEndOfFile ? BamStatus::kTruncated : st);
    }
    const unsigned char* magic = take(4);
    if (magic[0] != 'B' || magic[1] != 'A' || magic[2] != 'M' || magic[3] != '\1') {
      return fail(BamStatus::kBadMagic);
    }

    if (!ensure(4)) return fail(header_exhausted());
    const std::int32_t l_text = detail::lei32(take(4));
    if (l_text < 0) return fail(BamStatus::kBadHeader);
    // Skipped rather than held: a BAM with thousands of contigs has megabytes
    // of header text, and none of it is needed to place a fragment.
    if (!skip_bytes(l_text)) return fail(header_exhausted());

    if (!ensure(4)) return fail(header_exhausted());
    const std::int32_t n_ref = detail::lei32(take(4));
    if (n_ref < 0 || n_ref > detail::kMaxContigs) return fail(BamStatus::kBadHeader);

    names_.clear();
    lengths_.clear();
    names_.reserve(static_cast<std::size_t>(n_ref));
    lengths_.reserve(static_cast<std::size_t>(n_ref));
    for (std::int32_t i = 0; i < n_ref; ++i) {
      if (!ensure(4)) return fail(header_exhausted());
      const std::int32_t l_name = detail::lei32(take(4));
      // l_name counts the terminating NUL, so 1 is an EMPTY name. It used to be
      // accepted and written as an empty column 1 in the narrowPeak -- a row no
      // browser or bedtools reads. The BED front end (ContigTable::load) has
      // always refused an empty name; the two front ends now agree.
      if (l_name <= 1 || l_name > detail::kMaxContigNameLength) {
        return fail(BamStatus::kBadHeader);
      }
      if (!ensure(static_cast<std::size_t>(l_name) + 4)) return fail(header_exhausted());
      const unsigned char* name = take(static_cast<std::size_t>(l_name));
      const std::int32_t l_ref = detail::lei32(take(4));
      // SAM: @SQ LN is in [1, 2^31-1]. A zero-length contig was accepted here
      // and refused by the chrom.sizes loader, so the same genome was valid or
      // not depending on which file format described it.
      if (l_ref <= 0) return fail(BamStatus::kBadHeader);
      // These names are published through contig_names() and get written into
      // a TAB-separated track, so a name is not just a label -- a tab or a
      // newline in one forges whole rows in the caller's output, and the
      // header is the least trustworthy part of the file. Control bytes are
      // refused; everything printable is let through, because the SAM
      // grammar for rname is narrower than that and refusing a valid file is
      // the worse error. The terminating NUL is checked rather than assumed:
      // the length below trims a byte on the strength of it.
      if (name[static_cast<std::size_t>(l_name) - 1] != '\0') {
        return fail(BamStatus::kBadHeader);
      }
      for (std::int32_t j = 0; j + 1 < l_name; ++j) {
        const unsigned char ch = name[static_cast<std::size_t>(j)];
        if (ch < 0x20 || ch == 0x7f) return fail(BamStatus::kBadHeader);
      }
      // The stored length includes the terminating NUL.
      names_.emplace_back(reinterpret_cast<const char*>(name),
                          static_cast<std::size_t>(l_name) - 1);
      lengths_.push_back(static_cast<Coord>(l_ref));
    }
    // Contig names must be unique (SAM: @SQ SN). A duplicate was accepted, and
    // since the CLI writes names rather than tids, two contigs' peaks came out
    // under one name, indistinguishable -- and the BED front end refuses the
    // same duplicate in a chrom.sizes file. Checked through a sorted index
    // once the header is read, so nothing is allocated per name.
    std::vector<std::uint32_t> order(names_.size());
    std::iota(order.begin(), order.end(), 0u);
    std::sort(order.begin(), order.end(), [this](std::uint32_t a, std::uint32_t b) {
      return names_[a] < names_[b] || (names_[a] == names_[b] && a < b);
    });
    for (std::size_t i = 1; i < order.size(); ++i) {
      if (names_[order[i]] == names_[order[i - 1]]) return fail(BamStatus::kBadHeader);
    }
    return BamStatus::kOk;
  }

  // Points `rec` at a record's 32-byte core plus its variable tail, in place.
  BamStatus next_record(const unsigned char*& rec, std::size_t& rec_len) {
    if (!ensure(4)) {
      const BamStatus st = exhausted_status();
      // Bytes left over that cannot even form a record's length prefix are a
      // record cut in half, however tidily the last BGZF block ended. Without
      // this, a payload ending 1-3 bytes into a record reports a clean end of
      // file and silently drops it -- the same class of failure the
      // end-of-file marker check above exists to catch, one level down.
      if (st == BamStatus::kEndOfFile && available_bytes() != 0) {
        return BamStatus::kTruncated;
      }
      return st;
    }
    const std::int32_t block_size = detail::lei32(data_.data() + data_pos_);
    if (block_size < 32 || static_cast<std::size_t>(block_size) > detail::kMaxRecordSize) {
      return BamStatus::kBadRecord;
    }
    const std::size_t need = 4 + static_cast<std::size_t>(block_size);
    if (!ensure(need)) {
      const BamStatus st = exhausted_status();
      // A record that runs off the end is truncation, not a clean end.
      return st == BamStatus::kEndOfFile ? BamStatus::kTruncated : st;
    }
    const unsigned char* base = take(need);
    rec = base + 4;
    rec_len = static_cast<std::size_t>(block_size);
    return BamStatus::kOk;
  }

#if defined(PEAKS_BGZF_LIBDEFLATE)
  struct InflatorDeleter {
    void operator()(libdeflate_decompressor* d) const noexcept {
      if (d != nullptr) libdeflate_free_decompressor(d);
    }
  };
  std::unique_ptr<libdeflate_decompressor, InflatorDeleter> inflator_;
#elif defined(PEAKS_BGZF_ZLIB)
  // The mirror of InflatorDeleter. inflateEnd releases zlib's internal state,
  // which `delete` alone does not, and putting it in the deleter is what makes
  // the defaulted move assignment above correct -- unique_ptr::operator=
  // disposes of the old pointee through it.
  struct ZStreamDeleter {
    void operator()(z_stream* z) const noexcept {
      if (z != nullptr) {
        inflateEnd(z);
        delete z;
      }
    }
  };
  std::unique_ptr<z_stream, ZStreamDeleter> zs_;
#endif

  std::istream* in_ = nullptr;
  BamFilter filter_{};
  std::vector<unsigned char> raw_;
  std::vector<unsigned char> data_;
  std::size_t data_pos_ = 0;
  std::size_t data_len_ = 0;
  std::vector<std::string> names_;
  std::vector<Coord> lengths_;
  BamStatus failed_ = BamStatus::kOk;
  BamStatus block_error_ = BamStatus::kOk;
  bool stream_done_ = false;
  bool saw_eof_marker_ = false;
  std::int64_t records_seen_ = 0;
  std::int64_t fragments_emitted_ = 0;
  std::int64_t skipped_flags_ = 0;
  std::int64_t skipped_mapq_ = 0;
  std::int64_t skipped_span_ = 0;
  std::int64_t skipped_unpaired_ = 0;
  std::int64_t skipped_unplaced_ = 0;
};

#endif  // PEAKS_HAVE_BGZF

}  // namespace pto::peaks
