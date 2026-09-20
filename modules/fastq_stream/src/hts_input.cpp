// SPDX-License-Identifier: MIT
#include "fastq_stream/hts_input.hpp"

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>

#include <sys/stat.h>

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>

namespace fq {
namespace {

// 4-bit encoded base -> character, and its complement. Indexed by the nibble
// htslib packs SEQ into (the "=ACMGRSVTWYHKDBN" alphabet of the SAM spec).
constexpr char kNt16[17] = "=ACMGRSVTWYHKDBN";
constexpr char kNt16Comp[17] = "=TGKCYSBAWRDMHVN";

// Worst-case bytes one emitted record can occupy: '@' + qname + "/1" + '\n',
// seq + '\n', "+\n", qual + '\n'. QNAME is capped by the BAM spec at 254.
constexpr std::size_t kWorstRecord = 254 + 4 + 2 * kHtsMaxReadLen + 6;
static_assert(kWorstRecord <= kChunkBytes,
              "a chunk must hold at least one worst-case record");

}  // namespace

struct HtsInput::Impl {
  htsFile* fp = nullptr;
  sam_hdr_t* hdr = nullptr;
  bam1_t* rec = nullptr;
  bam1_t* held = nullptr;  // kInterleaved: record waiting for its mate
  bool have_held = false;
  bool eof = false;
  BufferPool* pool = nullptr;
  HtsOptions opts;
  HtsStats stats;
  const char* fmt = "bam";

  ~Impl() {
    if (held != nullptr) bam_destroy1(held);
    if (rec != nullptr) bam_destroy1(rec);
    if (hdr != nullptr) sam_hdr_destroy(hdr);
    if (fp != nullptr) hts_close(fp);
  }

  // --------------------------------------------------------------- decoding

  // SEQ is stored in ALIGNMENT orientation. A read that aligned to the minus
  // strand must be reverse-complemented, and its QUAL reversed, to recover the
  // sequence the instrument produced.
  //
  // This is the defect to be most careful about in this file, because getting
  // it wrong is silent: the FASTQ stays well-formed, the QC report stays
  // plausible, GC content barely moves -- and roughly half the reads are
  // reverse-complemented garbage for every downstream aligner.
  void decode_seq(const bam1_t* b, char* out) {
    const std::uint8_t* packed = bam_get_seq(b);
    const std::int32_t n = b->core.l_qseq;
    if ((b->core.flag & BAM_FREVERSE) != 0) {
      for (std::int32_t i = 0; i < n; ++i) {
        out[n - 1 - i] = kNt16Comp[bam_seqi(packed, i)];
      }
      ++stats.reverse_complemented;
    } else {
      for (std::int32_t i = 0; i < n; ++i) out[i] = kNt16[bam_seqi(packed, i)];
    }
  }

  void decode_qual(const bam1_t* b, char* out) {
    const std::uint8_t* q = bam_get_qual(b);
    const std::int32_t n = b->core.l_qseq;
    // 0xff in the first byte is htslib's marker for "QUAL is '*'", i.e. absent
    // for the whole record.
    if (n > 0 && q[0] == 0xff) {
      std::memset(out, opts.missing_qual, static_cast<std::size_t>(n));
      ++stats.synthesized_qual;
      return;
    }
    // SAM allows Phred 0..93 ('!'..'~'). A larger value cannot become a FASTQ
    // quality character; the pipeline's printable-byte check did refuse the
    // run, but as "malformed FASTQ record" -- about a BAM the user never
    // converted by hand. Named here, with the read and the value.
    std::uint8_t highest = 0;
    for (std::int32_t i = 0; i < n; ++i) highest = std::max(highest, q[i]);
    if (highest > 93) {
      std::int32_t at = 0;
      while (q[at] <= 93) ++at;
      throw std::runtime_error("read '" + std::string(bam_get_qname(b)) + "': QUAL value " +
                               std::to_string(q[at]) + " at base " + std::to_string(at + 1) +
                               " is outside the SAM range 0..93");
    }
    if ((b->core.flag & BAM_FREVERSE) != 0) {
      for (std::int32_t i = 0; i < n; ++i) {
        out[n - 1 - i] = static_cast<char>(q[i] + 33);
      }
    } else {
      for (std::int32_t i = 0; i < n; ++i) out[i] = static_cast<char>(q[i] + 33);
    }
  }

  // Writes one FASTQ record at `dst`, returning the bytes used. `suffix` is 0,
  // '1' or '2'. The caller guarantees kWorstRecord bytes of room.
  std::size_t format(const bam1_t* b, char suffix, char* dst) {
    const std::int32_t n = b->core.l_qseq;
    if (static_cast<std::size_t>(n) > kHtsMaxReadLen) {
      throw std::runtime_error(
          "read '" + std::string(bam_get_qname(b)) + "' is " + std::to_string(n) +
          " bp, over this build's " + std::to_string(kHtsMaxReadLen) +
          " bp limit for BAM/CRAM ingestion");
    }
    char* w = dst;
    *w++ = '@';
    const char* qname = bam_get_qname(b);
    const std::size_t qlen = std::strlen(qname);
    std::memcpy(w, qname, qlen);
    w += qlen;
    if (suffix != 0) {
      *w++ = '/';
      *w++ = suffix;
    }
    *w++ = '\n';
    decode_seq(b, w);
    w += n;
    *w++ = '\n';
    *w++ = '+';
    *w++ = '\n';
    decode_qual(b, w);
    w += n;
    *w++ = '\n';
    ++stats.emitted;
    return static_cast<std::size_t>(w - dst);
  }

  // Returns false when the alignment cannot become a FASTQ record.
  bool usable(const bam1_t* b) {
    const std::uint16_t flag = b->core.flag;
    if ((flag & BAM_FSECONDARY) != 0) {
      ++stats.skipped_secondary;
      return false;
    }
    if ((flag & BAM_FSUPPLEMENTARY) != 0) {
      ++stats.skipped_supplementary;
      return false;
    }
    if (b->core.l_qseq <= 0) {
      // SEQ '*'. Emitting a zero-length record would fail the assembler's
      // well_formed() check several ranks later with an error naming FASTQ,
      // which is a confusing way to learn the BAM had no sequence.
      ++stats.skipped_no_seq;
      return false;
    }
    return true;
  }

  static char mate_suffix(const bam1_t* b) {
    if ((b->core.flag & BAM_FREAD1) != 0) return '1';
    if ((b->core.flag & BAM_FREAD2) != 0) return '2';
    return 0;
  }

  // bam_copy1 grows the destination's data buffer and returns null if that
  // allocation fails; the result is warn_unused_result precisely because
  // ignoring it leaves `held` holding the previous record while the caller
  // believes it holds the new one -- which would emit the wrong read.
  void hold(const bam1_t* b) {
    if (bam_copy1(held, const_cast<bam1_t*>(b)) == nullptr) {
      throw std::runtime_error("out of memory buffering a read for its mate");
    }
    have_held = true;
  }

  static bool are_mates(const bam1_t* a, const bam1_t* b) {
    if (std::strcmp(bam_get_qname(a), bam_get_qname(b)) != 0) return false;
    const std::uint16_t fa = a->core.flag, fb = b->core.flag;
    if ((fa & BAM_FPAIRED) == 0 || (fb & BAM_FPAIRED) == 0) return false;
    // One must be read1 and the other read2; two copies of the same end are
    // not a pair.
    return ((fa & BAM_FREAD1) != 0 && (fb & BAM_FREAD2) != 0) ||
           ((fa & BAM_FREAD2) != 0 && (fb & BAM_FREAD1) != 0);
  }
};

// ---------------------------------------------------------------- lifetime

HtsInput::HtsInput(const std::string& path, BufferPool& pool, const HtsOptions& opts)
    : impl_(std::make_unique<Impl>()) {
  impl_->pool = &pool;
  impl_->opts = opts;

  impl_->fp = hts_open(path.c_str(), "r");
  if (impl_->fp == nullptr) {
    throw std::runtime_error("cannot open alignment file '" + path + "'");
  }

  const htsFormat* fmt = hts_get_format(impl_->fp);
  const bool is_cram = fmt != nullptr && fmt->format == cram;
  impl_->fmt = is_cram ? "cram" : (fmt != nullptr && fmt->format == sam ? "sam" : "bam");

  // CRAM stores no sequence of its own: bases are reconstructed against the
  // reference the file was compressed with. htslib resolves that reference
  // through REF_PATH, whose compiled-in default includes a REMOTE registry at
  // the EBI -- so decoding a CRAM on a stock htslib issues an outbound request
  // carrying the file's reference checksums, with no diagnostic.
  //
  // That is a direct violation of this project's premise (nothing leaves the
  // host that holds the data), and it is the kind of thing nobody notices until
  // an air-gapped node hangs. Refuse without a local reference, and pin the
  // resolver so a stale REF_PATH in the environment cannot reintroduce it.
  if (is_cram) {
    if (opts.reference_fasta.empty()) {
      throw std::runtime_error(
          "CRAM input requires --reference <genome.fa>: without it htslib "
          "resolves reference sequences over the network, which this tool does "
          "not do. Point it at the FASTA the CRAM was compressed against.");
    }
    ::setenv("REF_PATH", opts.reference_fasta.c_str(), /*overwrite=*/1);
    if (hts_set_fai_filename(impl_->fp, opts.reference_fasta.c_str()) < 0) {
      throw std::runtime_error(
          "cannot load the FASTA index for '" + opts.reference_fasta +
          "' (run: samtools faidx " + opts.reference_fasta + ")");
    }
  }

  if (opts.decode_threads > 1) {
    hts_set_opt(impl_->fp, HTS_OPT_NTHREADS, opts.decode_threads);
  }

  impl_->hdr = sam_hdr_read(impl_->fp);
  if (impl_->hdr == nullptr) {
    throw std::runtime_error("cannot read the header of '" + path + "'");
  }

  // Interleaved output needs mates adjacent. On a coordinate-sorted file they
  // are arbitrarily far apart -- up to the whole file on human WGS -- so
  // pairing would need an unbounded name-keyed map, which is a
  // memory-exhaustion shape, not a feature. Fail at open, where it costs
  // nothing, instead of at 40 GB resident.
  if (opts.layout == HtsLayout::kInterleaved) {
    kstring_t so = KS_INITIALIZE;
    const bool have_so = sam_hdr_find_tag_hd(impl_->hdr, "SO", &so) == 0;
    const std::string order = (have_so && so.s != nullptr) ? so.s : "";
    ks_free(&so);
    if (order == "coordinate") {
      throw std::runtime_error(
          "'" + path +
          "' is coordinate-sorted (@HD SO:coordinate), so mates are not "
          "adjacent and interleaved output would need unbounded memory to "
          "pair them. Regroup by name first:\n"
          "    samtools collate -u -O " + path + " | fastq_stream --input - ...");
    }
  }

  impl_->rec = bam_init1();
  impl_->held = bam_init1();
  if (impl_->rec == nullptr || impl_->held == nullptr) {
    throw std::runtime_error("out of memory allocating BAM record buffers");
  }
}

HtsInput::~HtsInput() = default;

const HtsStats& HtsInput::stats() const noexcept { return impl_->stats; }

uint64_t HtsInput::compressed_bytes() const noexcept {
  // htslib exposes no byte counter across formats; the caller reports 0 rather
  // than a number it cannot stand behind.
  return 0;
}

const char* HtsInput::format_name() const noexcept { return impl_->fmt; }

// -------------------------------------------------------------------- next

bool HtsInput::next(Buffer& out) {
  Impl& s = *impl_;
  if (s.eof && !s.have_held) return false;

  out = s.pool->acquire();
  char* base = out.chars();
  const std::size_t cap = out.capacity();
  std::size_t used = 0;

  // Reserve room for the largest thing one loop iteration can emit. Checking
  // BEFORE decoding is what guarantees a record is never formatted with
  // nowhere to put it -- and in interleaved mode one iteration can emit two.
  const std::size_t reserve =
      (s.opts.layout == HtsLayout::kInterleaved ? 2 : 1) * kWorstRecord;

  while (cap - used >= reserve) {
    if (s.eof) {
      // Stream ended with a record still waiting for a mate that never came.
      // Emitted rather than dropped: silently losing it would change the read
      // count without saying so.
      if (s.have_held) {
        used += s.format(s.held, Impl::mate_suffix(s.held), base + used);
        s.have_held = false;
        ++s.stats.singletons;
      }
      break;
    }

    const int rc = sam_read1(s.fp, s.hdr, s.rec);
    if (rc < -1) {
      throw std::runtime_error("truncated or corrupt alignment record");
    }
    if (rc == -1) {
      s.eof = true;
      continue;  // re-enter to flush a held record
    }
    ++s.stats.records_read;
    if (!s.usable(s.rec)) continue;

    if (s.opts.layout == HtsLayout::kSingle) {
      used += s.format(s.rec, 0, base + used);
      continue;
    }

    // Interleaved: a single-slot hold. Name-grouped input needs nothing more
    // than this, and it is O(1) memory by construction rather than by hoping
    // the file is well behaved.
    if (!s.have_held) {
      s.hold(s.rec);
      continue;
    }
    if (Impl::are_mates(s.held, s.rec)) {
      const bool held_is_first = (s.held->core.flag & BAM_FREAD1) != 0;
      const bam1_t* first = held_is_first ? s.held : s.rec;
      const bam1_t* second = held_is_first ? s.rec : s.held;
      used += s.format(first, '1', base + used);
      used += s.format(second, '2', base + used);
      s.have_held = false;
      ++s.stats.pairs;
    } else {
      used += s.format(s.held, Impl::mate_suffix(s.held), base + used);
      ++s.stats.singletons;
      s.hold(s.rec);
    }
  }

  if (used == 0) {
    out = Buffer{};
    return false;
  }
  out.set_size(used);
  return true;
}

// ------------------------------------------------------------------ probe

bool looks_like_alignment(const std::string& path) {
  // Content probe, not extension matching: a BAM named .fq.gz is still a BAM,
  // and letting it reach the gzip path produces a "malformed FASTQ record"
  // error several ranks downstream that says nothing about the real cause.
  if (path == "-") return false;  // stdin cannot be probed without consuming it
  // Neither can any other non-seekable input. hts_open() reads a full hfile
  // buffer to detect the format and hts_close() discards it, so on a FIFO --
  // including bash process substitution, `-i <(zcat reads.fq.gz)`, the most
  // common way to feed this tool -- the first ~32 KiB of the stream vanished
  // before the FASTQ reader opened it, and the run failed with "malformed
  // FASTQ record". Only a regular file can be probed and then reopened from
  // byte 0; a BAM arriving over a pipe must be converted upstream.
  struct stat st{};
  if (::stat(path.c_str(), &st) != 0 || !S_ISREG(st.st_mode)) return false;
  // A FASTQ that htslib cannot classify is the expected answer here, not an
  // error worth an "[E::hts_hopen] Failed to open file" line on stderr.
  const htsLogLevel saved = hts_get_log_level();
  hts_set_log_level(HTS_LOG_OFF);
  htsFile* fp = hts_open(path.c_str(), "r");
  hts_set_log_level(saved);
  if (fp == nullptr) return false;  // let the normal path report the error
  const htsFormat* fmt = hts_get_format(fp);
  const bool yes = fmt != nullptr &&
                   (fmt->format == bam || fmt->format == cram || fmt->format == sam);
  hts_close(fp);
  return yes;
}

}  // namespace fq
