// SPDX-License-Identifier: MIT
// fastq_stream — compressed input framing.
//
// IMPORTANT DESIGN CONSTRAINT (see docs/ARCHITECTURE.md §3): libdeflate has no
// streaming interface. It decompresses one whole member at a time into a
// buffer whose size the caller must bound in advance. That is a perfect fit for
// BGZF, where every member declares its compressed size in the header and
// expands to at most 64 KiB — so members can be handed to a pool of threads and
// inflated concurrently. It is *not* usable for a single-member plain gzip
// stream, which would have to be resident in full.
//
// Therefore two paths exist and the tool reports which one it took:
//
//   kBgzf   — block-parallel. libdeflate, N inflater threads, O(1) memory.
//   kPlain  — single-member (or non-BGZF multi-member) gzip. Streaming zlib
//             inflate on one thread; parallelism comes from the QC/trim stage
//             downstream. This is the path plain Illumina bcl2fastq output
//             takes, and the benchmark must report it as such.
//   kRaw    — uncompressed FASTQ, passed through.
#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include "fastq_stream/buffer.hpp"

namespace fq {

enum class GzipFormat { kRaw, kPlain, kBgzf };

const char* format_name(GzipFormat f) noexcept;

// Largest uncompressed payload of a single BGZF member, per the SAM spec.
inline constexpr std::size_t kBgzfMaxBlock = 65536;
// Members batched per task, chosen so a batch expands to exactly kChunkBytes.
inline constexpr std::size_t kBgzfBatchMembers = kChunkBytes / kBgzfMaxBlock;

class GzipInput {
 public:
  // `path` may be "-" for stdin.
  GzipInput(const std::string& path, BufferPool& pool);
  ~GzipInput();

  GzipInput(const GzipInput&) = delete;
  GzipInput& operator=(const GzipInput&) = delete;

  GzipFormat format() const noexcept;
  bool block_parallel() const noexcept { return format() == GzipFormat::kBgzf; }

  // In kBgzf mode `out` receives whole, still-compressed BGZF members to be
  // inflated by a worker. In kPlain/kRaw mode `out` receives decompressed
  // bytes directly. Returns false at end of stream.
  bool next(Buffer& out);

  uint64_t compressed_bytes() const noexcept;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

// Per-thread libdeflate decompressor for BGZF member batches.
class BgzfInflater {
 public:
  BgzfInflater();
  ~BgzfInflater();

  BgzfInflater(const BgzfInflater&) = delete;
  BgzfInflater& operator=(const BgzfInflater&) = delete;

  // Inflates every member in `in` back-to-back into `out`. `out` must have
  // capacity for kBgzfBatchMembers * kBgzfMaxBlock bytes.
  bool inflate(const Buffer& in, Buffer& out, std::string* err);

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace fq
