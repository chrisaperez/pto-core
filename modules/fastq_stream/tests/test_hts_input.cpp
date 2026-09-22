// SPDX-License-Identifier: MIT
//
// BAM ingestion. Dependency-free in the style of test_trimmer.cpp.
//
// Fixtures are synthesised in-process with htslib rather than checked in, for
// the same reason the gzip tests build their own BGZF: a binary fixture cannot
// be reviewed, and shelling out to `samtools` assumes a tool that is not a
// build dependency.
//
// The load-bearing case is test_reverse_strand_is_restored. A missed reverse
// complement is invisible to every other check in this file -- the FASTQ stays
// well-formed, the record count is right, the QC report looks plausible -- and
// silently ruins half the reads for any downstream aligner.

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <random>
#include <string>
#include <vector>

#include "fastq_stream/hts_input.hpp"
#include "fastq_stream/pipeline.hpp"
#include "fastq_stream/record.hpp"

namespace {

int g_failures = 0;

#define CHECK(cond)                                                            \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::fprintf(stderr, "FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__); \
      ++g_failures;                                                            \
    }                                                                          \
  } while (0)

#define CHECK_EQ_STR(a, b)                                                     \
  do {                                                                         \
    const std::string lhs_ = (a), rhs_ = (b);                                  \
    if (lhs_ != rhs_) {                                                        \
      std::fprintf(stderr, "FAILED: %s == %s\n  got  '%s'\n  want '%s' (%s:%d)\n", \
                   #a, #b, lhs_.c_str(), rhs_.c_str(), __FILE__, __LINE__);    \
      ++g_failures;                                                            \
    }                                                                          \
  } while (0)

#define CHECK_THROWS(stmt)                                                     \
  do {                                                                         \
    bool threw_ = false;                                                       \
    try {                                                                      \
      stmt;                                                                    \
    } catch (const std::exception&) {                                          \
      threw_ = true;                                                           \
    }                                                                          \
    if (!threw_) {                                                             \
      std::fprintf(stderr, "FAILED (no throw): %s (%s:%d)\n", #stmt, __FILE__, \
                   __LINE__);                                                  \
      ++g_failures;                                                            \
    }                                                                          \
  } while (0)

// --------------------------------------------------------------- fixtures

struct ReadSpec {
  std::string qname;
  std::uint16_t flag = 4;  // unmapped by default: these are reads, not alignments
  std::string seq;
  std::string qual;  // empty => QUAL absent (0xff)
};

// A BAM written to a temp path, deleted on scope exit.
class TempBam {
 public:
  TempBam(const std::vector<ReadSpec>& reads, const char* sort_order = "unsorted") {
    static std::mt19937_64 rng{std::random_device{}()};
    path_ = (std::filesystem::temp_directory_path() /
             ("fq_hts_" + std::to_string(rng()) + ".bam"))
                .string();

    sam_hdr_t* hdr = sam_hdr_init();
    sam_hdr_add_line(hdr, "HD", "VN", "1.6", "SO", sort_order, nullptr);
    sam_hdr_add_line(hdr, "SQ", "SN", "chr1", "LN", "1000000", nullptr);

    htsFile* fp = hts_open(path_.c_str(), "wb");
    if (fp == nullptr || sam_hdr_write(fp, hdr) < 0) {
      std::fprintf(stderr, "cannot create test BAM at %s\n", path_.c_str());
      std::exit(2);
    }

    bam1_t* rec = bam_init1();
    for (const auto& r : reads) {
      // A quality string of all 0xff is how the BAM format spells "absent".
      const std::string q = r.qual.empty() ? std::string(r.seq.size(), '\xff')
                                           : phred_bytes(r.qual);
      // htslib refuses to build a mapped record without a CIGAR, so anything
      // whose flag lacks BAM_FUNMAP gets a placed position and a full-length
      // match. That also keeps the reverse-strand cases realistic: FLAG 0x10 on
      // an unmapped record would be a contradiction the decoder should never
      // have to see.
      const bool mapped = (r.flag & BAM_FUNMAP) == 0 && !r.seq.empty();
      const std::uint32_t cigar =
          bam_cigar_gen(static_cast<std::uint32_t>(r.seq.size()), BAM_CMATCH);
      const int rc = bam_set1(rec, r.qname.size(), r.qname.c_str(), r.flag,
                              /*tid=*/mapped ? 0 : -1, /*pos=*/mapped ? 100 : -1,
                              /*mapq=*/mapped ? 60 : 0,
                              /*n_cigar=*/mapped ? 1 : 0, mapped ? &cigar : nullptr,
                              /*mtid=*/-1, /*mpos=*/-1,
                              /*isize=*/0, r.seq.size(), r.seq.c_str(), q.c_str(), 0);
      if (rc < 0 || sam_write1(fp, hdr, rec) < 0) {
        std::fprintf(stderr, "cannot write test record %s\n", r.qname.c_str());
        std::exit(2);
      }
    }
    bam_destroy1(rec);
    hts_close(fp);
    sam_hdr_destroy(hdr);
  }

  ~TempBam() { std::filesystem::remove(path_); }
  TempBam(const TempBam&) = delete;
  TempBam& operator=(const TempBam&) = delete;

  const std::string& path() const { return path_; }

 private:
  // bam_set1 takes raw Phred values, not ASCII. Converting here keeps the
  // ReadSpec literals readable as the FASTQ they are meant to round-trip to.
  static std::string phred_bytes(const std::string& ascii) {
    std::string out(ascii.size(), '\0');
    for (std::size_t i = 0; i < ascii.size(); ++i) {
      out[i] = static_cast<char>(ascii[i] - 33);
    }
    return out;
  }

  std::string path_;
};

// Drains an HtsInput into one string of FASTQ text.
std::string drain(const std::string& path, const fq::HtsOptions& opts,
                  fq::HtsStats* stats_out = nullptr) {
  fq::BufferPool pool(fq::kChunkBytes);
  fq::HtsInput in(path, pool, opts);
  std::string all;
  fq::Buffer buf;
  while (in.next(buf)) {
    all.append(buf.chars(), buf.size());
    buf = fq::Buffer{};
  }
  if (stats_out != nullptr) *stats_out = in.stats();
  return all;
}

// ------------------------------------------------------------------ tests

void test_single_end_round_trip() {
  const TempBam bam({
      {"readA", 4, "ACGTACGTAC", "IIIIIIIIII"},
      {"readB", 4, "TTTTGGGGCC", "01234567:;"},
  });
  fq::HtsStats st;
  const std::string out = drain(bam.path(), {}, &st);

  CHECK_EQ_STR(out,
               "@readA\nACGTACGTAC\n+\nIIIIIIIIII\n"
               "@readB\nTTTTGGGGCC\n+\n01234567:;\n");
  CHECK(st.records_read == 2);
  CHECK(st.emitted == 2);
  CHECK(st.reverse_complemented == 0);
}

// THE case. FLAG 0x10 means SEQ is stored in alignment orientation, so the
// original read is its reverse complement with QUAL reversed. Asserted against
// a hand-written expectation rather than a second implementation of the same
// transform, which would agree with itself while both were wrong.
void test_reverse_strand_is_restored() {
  //  stored SEQ  ACGTTGCA        -> revcomp TGCAACGT
  //  stored QUAL ABCDEFGH        -> reversed HGFEDCBA
  const TempBam bam({{"rev", 16, "ACGTTGCA", "ABCDEFGH"}});
  fq::HtsStats st;
  const std::string out = drain(bam.path(), {}, &st);

  CHECK_EQ_STR(out, "@rev\nTGCAACGT\n+\nHGFEDCBA\n");
  CHECK(st.reverse_complemented == 1);
}

// Forward and reverse copies of the same underlying read must produce the same
// FASTQ record body.
void test_forward_and_reverse_agree() {
  const TempBam fwd({{"x", 0, "GATTACAGATTACA", "IIIIIIIIIIIIII"}});
  const TempBam rev({{"x", 16, "TGTAATCTGTAATC", "IIIIIIIIIIIIII"}});
  CHECK_EQ_STR(drain(fwd.path(), {}), drain(rev.path(), {}));
}

void test_secondary_and_supplementary_are_dropped() {
  const TempBam bam({
      {"keep", 0, "AAAA", "IIII"},
      {"sec", 0x100, "CCCC", "IIII"},
      {"sup", 0x800, "GGGG", "IIII"},
  });
  fq::HtsStats st;
  const std::string out = drain(bam.path(), {}, &st);

  // Keeping either would silently change how many reads the run contained: a
  // secondary alignment duplicates a primary, a supplementary one is a
  // hard-clipped fragment presented as a whole read.
  CHECK_EQ_STR(out, "@keep\nAAAA\n+\nIIII\n");
  CHECK(st.skipped_secondary == 1);
  CHECK(st.skipped_supplementary == 1);
  CHECK(st.emitted == 1);
}

void test_absent_qual_is_synthesized() {
  const TempBam bam({{"noqual", 4, "ACGTAC", ""}});
  fq::HtsOptions opts;
  opts.missing_qual = 'F';
  fq::HtsStats st;
  const std::string out = drain(bam.path(), opts, &st);

  CHECK_EQ_STR(out, "@noqual\nACGTAC\n+\nFFFFFF\n");
  CHECK(st.synthesized_qual == 1);
}

void test_absent_seq_is_skipped() {
  const TempBam bam({
      {"empty", 4, "", ""},
      {"real", 4, "ACGT", "IIII"},
  });
  fq::HtsStats st;
  const std::string out = drain(bam.path(), {}, &st);

  // A zero-length record would fail well_formed() several ranks downstream with
  // an error naming FASTQ, which is a confusing way to learn the BAM had no
  // sequence.
  CHECK_EQ_STR(out, "@real\nACGT\n+\nIIII\n");
  CHECK(st.skipped_no_seq == 1);
}

// Every emitted chunk must contain only whole records, so the assembler's
// carry stays empty and RecordIterator can parse the concatenation.
void test_chunk_boundaries_never_split_a_record() {
  std::vector<ReadSpec> reads;
  const std::string seq(300, 'A');
  const std::string qual(300, 'I');
  for (int i = 0; i < 4000; ++i) {
    reads.push_back({"read" + std::to_string(i), 4, seq, qual});
  }
  const TempBam bam(reads);

  fq::BufferPool pool(fq::kChunkBytes);
  fq::HtsInput in(bam.path(), pool, {});
  std::size_t chunks = 0, records = 0;
  fq::Buffer buf;
  while (in.next(buf)) {
    ++chunks;
    // Whole records only: the last complete boundary is the end of the chunk.
    CHECK(fq::last_record_boundary(buf.chars(), buf.size()) == buf.size());
    fq::RecordIterator it(buf.bytes(), buf.size());
    fq::RecordView r;
    while (it.next(r)) {
      CHECK(r.well_formed());
      ++records;
    }
    buf = fq::Buffer{};
  }
  CHECK(chunks > 3);  // otherwise the test proves nothing about boundaries
  CHECK(records == 4000);
}

// ------------------------------------------------------------ interleaved

void test_interleaved_pairs() {
  const TempBam bam({
      // 0x1 paired, 0x40 read1, 0x80 read2.
      {"pair1", 0x1 | 0x40, "AAAA", "IIII"},
      {"pair1", 0x1 | 0x80, "CCCC", "JJJJ"},
      // Mates out of order within the group: read2 first.
      {"pair2", 0x1 | 0x80, "GGGG", "KKKK"},
      {"pair2", 0x1 | 0x40, "TTTT", "LLLL"},
  });
  fq::HtsOptions opts;
  opts.layout = fq::HtsLayout::kInterleaved;
  fq::HtsStats st;
  const std::string out = drain(bam.path(), opts, &st);

  // /1 always first, regardless of the order the file happened to store them.
  CHECK_EQ_STR(out,
               "@pair1/1\nAAAA\n+\nIIII\n"
               "@pair1/2\nCCCC\n+\nJJJJ\n"
               "@pair2/1\nTTTT\n+\nLLLL\n"
               "@pair2/2\nGGGG\n+\nKKKK\n");
  CHECK(st.pairs == 2);
  CHECK(st.singletons == 0);
}

void test_interleaved_singleton_is_emitted_not_dropped() {
  const TempBam bam({
      {"lonely", 0x1 | 0x40, "AAAA", "IIII"},
      {"pair", 0x1 | 0x40, "CCCC", "JJJJ"},
      {"pair", 0x1 | 0x80, "GGGG", "KKKK"},
      {"trailing", 0x1 | 0x80, "TTTT", "LLLL"},
  });
  fq::HtsOptions opts;
  opts.layout = fq::HtsLayout::kInterleaved;
  fq::HtsStats st;
  const std::string out = drain(bam.path(), opts, &st);

  // Both singletons survive -- including the one at EOF, which is the easy one
  // to lose. Dropping them would change the read count without saying so.
  CHECK_EQ_STR(out,
               "@lonely/1\nAAAA\n+\nIIII\n"
               "@pair/1\nCCCC\n+\nJJJJ\n"
               "@pair/2\nGGGG\n+\nKKKK\n"
               "@trailing/2\nTTTT\n+\nLLLL\n");
  CHECK(st.pairs == 1);
  CHECK(st.singletons == 2);
}

void test_coordinate_sorted_rejects_interleaved() {
  const TempBam bam({{"a", 0x1 | 0x40, "AAAA", "IIII"}}, "coordinate");
  fq::BufferPool pool(fq::kChunkBytes);
  fq::HtsOptions opts;
  opts.layout = fq::HtsLayout::kInterleaved;

  // Fails at OPEN, before a single record is read. On a coordinate-sorted file
  // mates are arbitrarily far apart, so pairing would need an unbounded
  // name-keyed map -- a memory-exhaustion shape, not a feature. Failing here
  // costs nothing; failing at 40 GB resident costs a node.
  CHECK_THROWS(fq::HtsInput(bam.path(), pool, opts));

  // ...and the same file is fine in the default single-end layout.
  fq::HtsOptions single;
  CHECK_EQ_STR(drain(bam.path(), single), "@a\nAAAA\n+\nIIII\n");
}

// -------------------------------------------------------------- detection

void test_format_probe() {
  const TempBam bam({{"a", 4, "ACGT", "IIII"}});
  CHECK(fq::looks_like_alignment(bam.path()));
  CHECK(!fq::looks_like_alignment("-"));
  CHECK(!fq::looks_like_alignment("/nonexistent/path/nope.bam"));

  const auto fq_path =
      (std::filesystem::temp_directory_path() / "fq_hts_probe.fq").string();
  {
    std::FILE* f = std::fopen(fq_path.c_str(), "w");
    std::fputs("@r\nACGT\n+\nIIII\n", f);
    std::fclose(f);
  }
  CHECK(!fq::looks_like_alignment(fq_path));
  std::filesystem::remove(fq_path);
}

// Audit 2026-09-11: a QUAL value past the SAM range is refused BY NAME. The
// pipeline's printable-byte check already refused it, but as a "malformed FASTQ
// record" -- about a BAM.
void test_qual_past_93_is_refused_by_name() {
  std::string high(4, 'I');
  high[2] = static_cast<char>(94 + 33);  // TempBam stores ascii - 33, so this is 94
  const TempBam bad({{"hiq", 4, "ACGT", high}});
  fq::HtsOptions opts;
  std::string message;
  try {
    drain(bad.path(), opts);
  } catch (const std::exception& e) {
    message = e.what();
  }
  CHECK(message.find("QUAL value 94") != std::string::npos);
  CHECK(message.find("hiq") != std::string::npos);

  // 93 is the top of the range and round-trips as '~'.
  std::string top(4, 'I');
  top[2] = '~';
  const TempBam ok({{"top", 4, "ACGT", top}});
  CHECK_EQ_STR(drain(ok.path(), opts), "@top\nACGT\n+\nII~I\n");
}

// ------------------------------------------------------------- end to end

// The test that proves the seam is behaviour-preserving: the same reads
// delivered as a BAM and as FASTQ text must produce identical QC.
void test_pipeline_matches_equivalent_fastq() {
  std::vector<ReadSpec> reads;
  std::string fastq;
  std::mt19937 rng(99);
  const char* bases = "ACGT";
  for (int i = 0; i < 500; ++i) {
    std::string seq, qual;
    for (int j = 0; j < 100; ++j) {
      seq += bases[rng() % 4];
      qual += static_cast<char>(33 + 20 + static_cast<int>(rng() % 20));
    }
    const std::string name = "r" + std::to_string(i);
    reads.push_back({name, 4, seq, qual});
    fastq += "@" + name + "\n" + seq + "\n+\n" + qual + "\n";
  }
  const TempBam bam(reads);

  const auto fq_path =
      (std::filesystem::temp_directory_path() / "fq_hts_equiv.fq").string();
  {
    std::FILE* f = std::fopen(fq_path.c_str(), "w");
    std::fwrite(fastq.data(), 1, fastq.size(), f);
    std::fclose(f);
  }

  fq::Config from_bam;
  from_bam.input = bam.path();
  from_bam.hts_input = true;
  from_bam.qc_only = true;
  from_bam.threads = 4;

  fq::Config from_text;
  from_text.input = fq_path;
  from_text.qc_only = true;
  from_text.threads = 4;

  const auto a = fq::run_pipeline(from_bam);
  const auto b = fq::run_pipeline(from_text);
  std::filesystem::remove(fq_path);

  CHECK(a.stats.reads_in == b.stats.reads_in);
  CHECK(a.stats.reads_out == b.stats.reads_out);
  CHECK(a.stats.bases_in == b.stats.bases_in);
  CHECK(a.stats.bases_out == b.stats.bases_out);
  CHECK(a.from_alignment);
  CHECK(!b.from_alignment);
}

}  // namespace

int main() {
  test_single_end_round_trip();
  test_reverse_strand_is_restored();
  test_forward_and_reverse_agree();
  test_secondary_and_supplementary_are_dropped();
  test_absent_qual_is_synthesized();
  test_absent_seq_is_skipped();
  test_chunk_boundaries_never_split_a_record();
  test_interleaved_pairs();
  test_interleaved_singleton_is_emitted_not_dropped();
  test_coordinate_sorted_rejects_interleaved();
  test_format_probe();
  test_qual_past_93_is_refused_by_name();
  test_pipeline_matches_equivalent_fastq();

  if (g_failures == 0) {
    std::printf("All tests passed.\n");
    return 0;
  }
  std::fprintf(stderr, "%d test(s) failed.\n", g_failures);
  return 1;
}
