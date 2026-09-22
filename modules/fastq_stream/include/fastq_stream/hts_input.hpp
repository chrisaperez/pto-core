// SPDX-License-Identifier: MIT
// fastq_stream — BAM/CRAM ingestion.
//
// WHERE THIS SITS. Everything downstream of the assembler in pipeline.hpp is
// defined over *bytes of FASTQ text*: last_record_boundary() counts newlines,
// RecordIterator splits on '\n', RecordProcessor writes '@id/seq/+/qual' back
// out. Threading an alignment-shaped record through those ranks would mean a
// second RecordProcessor and a variant Buffer payload, for no gain. So this
// replaces the READER RANK ONLY -- it hands the rest of the pipeline chunks of
// ordinary FASTQ text and nothing downstream knows the input was a BAM.
//
// ON "ZERO COPY". BAM ingestion cannot be zero-copy at the record level and it
// is worth being precise about why: bam1_t stores SEQ 4-bit packed and QUAL as
// raw Phred values, neither of which is the byte layout the SIMD trim kernels
// read. The nibble expansion and the +33 are required transforms, not copies we
// chose. What is preserved is the property that matters -- the expansion writes
// DIRECTLY into a pooled Buffer, and from there to the writer the payload moves
// by pointer only, exactly like the gzip path.
//
// Chunks always end on a record boundary, so the assembler's carry buffer stays
// permanently empty in this mode. It still runs: one memchr pass at memory
// bandwidth in exchange for a single code path downstream.
#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include "fastq_stream/buffer.hpp"

namespace fq {

// Longest read this path will emit. Reads beyond it are rejected by name rather
// than silently truncated. 32 kbp covers Illumina by three orders of magnitude
// and matches the ceiling the assembler's kPrefix carry already imposes on the
// FASTQ path, so the two entry points agree on what "too long" means.
inline constexpr std::size_t kHtsMaxReadLen = 32 * 1024;

enum class HtsLayout {
  // One FASTQ record per alignment. Mates, if present, are emitted
  // independently in file order. Works on any input, sorted or not.
  kSingle,
  // Mates emitted back to back as name/1 then name/2. Requires name-grouped
  // input; see the sort-order note on HtsOptions::layout below.
  kInterleaved,
};

struct HtsOptions {
  HtsLayout layout = HtsLayout::kSingle;

  // Local FASTA for CRAM reference resolution. Empty + CRAM input is a hard
  // error at open time, never a silent fetch: htslib's default REF_PATH
  // includes a remote registry, so decoding a CRAM without this would send the
  // file's reference checksums off-host. See the constructor's comment.
  std::string reference_fasta;

  // Threads for htslib's own BGZF decode pool. The pipeline hands its inflater
  // budget here, because BAM framing is record-level and cannot use this
  // module's block-parallel libdeflate path.
  int decode_threads = 0;

  // Substituted for every base when QUAL is absent (0xff). 'I' == Q40.
  char missing_qual = 'I';
};

struct HtsStats {
  uint64_t records_read = 0;
  uint64_t emitted = 0;
  uint64_t skipped_secondary = 0;      // FLAG 0x100
  uint64_t skipped_supplementary = 0;  // FLAG 0x800
  uint64_t skipped_no_seq = 0;         // SEQ absent ('*')
  uint64_t synthesized_qual = 0;       // QUAL absent (0xff)
  uint64_t reverse_complemented = 0;   // FLAG 0x10
  uint64_t pairs = 0;                  // kInterleaved: mates emitted together
  uint64_t singletons = 0;             // kInterleaved: mate never arrived
};

// Reader-rank source. Deliberately mirrors GzipInput's shape so pipeline.cpp
// can hold either behind one next(Buffer&) call.
class HtsInput {
 public:
  HtsInput(const std::string& path, BufferPool& pool, const HtsOptions& opts);
  ~HtsInput();

  HtsInput(const HtsInput&) = delete;
  HtsInput& operator=(const HtsInput&) = delete;

  // Fills `out` with whole FASTQ records only. Returns false at end of stream.
  bool next(Buffer& out);

  [[nodiscard]] const HtsStats& stats() const noexcept;
  [[nodiscard]] uint64_t compressed_bytes() const noexcept;
  // "bam", "cram" or "sam", from htslib's format probe.
  [[nodiscard]] const char* format_name() const noexcept;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

// Whether `path` is an alignment file, by htslib's content probe rather than by
// extension -- a BAM named .fq.gz is still a BAM, and misdetecting it produces
// a baffling "malformed FASTQ record" error several ranks downstream.
// Returns false for anything htslib cannot open, leaving the error to the
// normal path.
[[nodiscard]] bool looks_like_alignment(const std::string& path);

}  // namespace fq
