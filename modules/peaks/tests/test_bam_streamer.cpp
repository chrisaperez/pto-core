// SPDX-License-Identifier: MIT
//
// pto-peaks: native BGZF/BAM ingestion.
//
// Every fixture in this file is built byte by byte, in the test, with no
// external file and no htslib. That is possible because a BGZF block may carry
// a DEFLATE *stored* block -- three bytes of framing around uncompressed data
// -- so the test needs no compressor, only the decompressor it is checking.
// The CRC32 below is likewise a second implementation, not a call into the one
// the streamer uses.
//
// What that buys is the whole point of the file: malformed input can be
// constructed exactly. A truncated block, a BSIZE that overruns the file, an
// n_ref that overruns the header, a CRC that does not match -- these are the
// cases a real BAM corpus does not contain and a corrupted download does, and
// every one of them must be refused rather than parsed into a plausible peak.
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <algorithm>
#include <sstream>
#include <utility>
#include <string>
#include <vector>

#include "peaks_test_util.hpp"

#include "peaks/bam_streamer.hpp"

#if !defined(PEAKS_HAVE_BGZF)

int main(int argc, char** argv) {
    // There is no fixture to emit without a decompressor, and saying so on
    // stderr with a non-zero status is what lets the CLI script tell "no BAM
    // support" apart from "here is a BAM". Writing the skip text to stdout
    // instead made the script mistake the message for a fixture.
    if (argc > 1 && std::string(argv[1]) == "--emit-fixture") {
        std::fprintf(stderr, "no BGZF support in this build; no fixture to emit\n");
        return 1;
    }
    // Skips are not failures, but a suite that silently tests nothing reads as
    // coverage. Say so, loudly, and succeed.
    std::printf("== peaks/bam_streamer ==\n");
    std::printf("  SKIP built without BGZF support (-DPEAKS_WITH_BGZF=OFF, or "
                "neither libdeflate nor zlib was found)\n");
    std::printf("0 checks, 0 failure(s)\n");
    return 0;
}

#else

namespace {

using pto::peaks::BamFilter;
using pto::peaks::BamFragment;
using pto::peaks::BamStatus;
using pto::peaks::BamStreamer;
using pto::peaks::Coord;

// --- byte plumbing ---------------------------------------------------------

void put_u16(std::string& out, std::uint16_t v) {
    out.push_back(static_cast<char>(v & 0xFF));
    out.push_back(static_cast<char>((v >> 8) & 0xFF));
}

void put_u32(std::string& out, std::uint32_t v) {
    for (int i = 0; i < 4; ++i) out.push_back(static_cast<char>((v >> (8 * i)) & 0xFF));
}

void put_i32(std::string& out, std::int32_t v) {
    put_u32(out, static_cast<std::uint32_t>(v));
}

// CRC-32/IEEE, bit-at-a-time. Independent of whatever libdeflate or zlib does.
std::uint32_t crc32_ref(const std::string& data) {
    std::uint32_t crc = 0xFFFFFFFFu;
    for (char ch : data) {
        crc ^= static_cast<unsigned char>(ch);
        for (int k = 0; k < 8; ++k) {
            crc = (crc & 1u) ? ((crc >> 1) ^ 0xEDB88320u) : (crc >> 1);
        }
    }
    return ~crc;
}

// A DEFLATE stream holding one final STORED block. No compression, and every
// inflater accepts it.
std::string deflate_stored(const std::string& payload) {
    std::string out;
    out.push_back(static_cast<char>(0x01));  // BFINAL = 1, BTYPE = 00 (stored)
    put_u16(out, static_cast<std::uint16_t>(payload.size()));
    put_u16(out, static_cast<std::uint16_t>(~payload.size() & 0xFFFF));
    out += payload;
    return out;
}

struct BlockTweak {
    bool bad_magic = false;
    bool drop_extra = false;      // FLG without FEXTRA
    std::uint8_t extra_flg = 0;   // FLG bits BGZF does not allow (FHCRC, FNAME, ...)
    int bsize_delta = 0;          // corrupt the advertised block size
    std::uint32_t crc_override = 0;
    bool use_crc_override = false;
    std::uint32_t isize_override = 0;
    bool use_isize_override = false;
};

std::string bgzf_block(const std::string& payload, const BlockTweak& t = {}) {
    const std::string cdata = deflate_stored(payload);
    std::string out;
    out.push_back(static_cast<char>(t.bad_magic ? 0x1e : 0x1f));
    out.push_back(static_cast<char>(0x8b));
    out.push_back(static_cast<char>(0x08));                       // CM = deflate
    out.push_back(static_cast<char>((t.drop_extra ? 0x00 : 0x04) | t.extra_flg));  // FLG
    put_u32(out, 0);                                              // MTIME
    out.push_back(static_cast<char>(0x00));                       // XFL
    out.push_back(static_cast<char>(0xFF));                       // OS
    put_u16(out, 6);                                              // XLEN
    out.push_back('B');
    out.push_back('C');
    put_u16(out, 2);                                              // SLEN
    const std::size_t total = 12 + 6 + cdata.size() + 8;
    put_u16(out, static_cast<std::uint16_t>(total - 1 + t.bsize_delta));  // BSIZE
    out += cdata;
    put_u32(out, t.use_crc_override ? t.crc_override : crc32_ref(payload));
    put_u32(out, t.use_isize_override ? t.isize_override
                                      : static_cast<std::uint32_t>(payload.size()));
    return out;
}

std::string bgzf_eof() { return bgzf_block(std::string{}); }

// Splits a byte string into BGZF blocks of at most `chunk` bytes, so a fixture
// can be made to straddle block boundaries at any offset.
std::string bgzf_stream(const std::string& payload, std::size_t chunk = 65280) {
    std::string out;
    std::size_t at = 0;
    if (payload.empty()) out += bgzf_block(std::string{});
    while (at < payload.size()) {
        const std::size_t n = std::min(chunk, payload.size() - at);
        out += bgzf_block(payload.substr(at, n));
        at += n;
    }
    out += bgzf_eof();
    return out;
}

struct Ref {
    std::string name;
    std::int32_t length;
};

struct HeaderTweak {
    bool bad_magic = false;
    int l_text_delta = 0;
    int n_ref_delta = 0;
    int l_name_delta = 0;
};

std::string bam_header(const std::vector<Ref>& refs, const std::string& text = "@HD\tVN:1.6\tSO:coordinate\n",
                       const HeaderTweak& t = {}) {
    std::string out;
    out += t.bad_magic ? "BAN\1" : "BAM\1";
    put_i32(out, static_cast<std::int32_t>(text.size()) + t.l_text_delta);
    out += text;
    put_i32(out, static_cast<std::int32_t>(refs.size()) + t.n_ref_delta);
    for (const Ref& r : refs) {
        put_i32(out, static_cast<std::int32_t>(r.name.size()) + 1 + t.l_name_delta);
        out += r.name;
        out.push_back('\0');
        put_i32(out, r.length);
    }
    return out;
}

// A minimal alignment record: fixed 32-byte core, a name, no cigar, no seq.
std::string bam_record(std::int32_t refid, std::int32_t pos, std::uint16_t flag,
                       std::uint8_t mapq, std::int32_t tlen,
                       const std::string& name = "r", int block_size_delta = 0,
                       std::size_t pad = 0, std::int32_t l_seq = 0) {
    std::string core;
    put_i32(core, refid);
    put_i32(core, pos);
    core.push_back(static_cast<char>(name.size() + 1));  // l_read_name
    core.push_back(static_cast<char>(mapq));
    put_u16(core, 0);      // bin
    put_u16(core, 0);      // n_cigar_op
    put_u16(core, flag);
    put_i32(core, l_seq);
    put_i32(core, refid);  // next_refID
    put_i32(core, pos);    // next_pos
    put_i32(core, tlen);
    core += name;
    core.push_back('\0');
    // Trailing bytes stand in for aux data, which is how a record reaches an
    // interesting size without a megabyte of sequence.
    core.append(pad, '\0');

    std::string out;
    put_i32(out, static_cast<std::int32_t>(core.size()) + block_size_delta);
    out += core;
    return out;
}

constexpr std::uint16_t kPaired = 0x1;
constexpr std::uint16_t kProper = 0x2;
constexpr std::uint16_t kUnmapped = 0x4;
constexpr std::uint16_t kRead1 = 0x40;
constexpr std::uint16_t kSecondary = 0x100;
constexpr std::uint16_t kQcFail = 0x200;
constexpr std::uint16_t kDuplicate = 0x400;
constexpr std::uint16_t kSupplementary = 0x800;
constexpr std::uint16_t kGoodPair = kPaired | kProper | kRead1;

// Runs a stream to exhaustion. Returns the terminating status.
BamStatus collect(const std::string& bytes, std::vector<BamFragment>& out,
                  const BamFilter& filter = {}, BamStreamer* keep = nullptr) {
    static BamStreamer local;
    BamStreamer& s = keep != nullptr ? *keep : local;
    s = BamStreamer{};
    std::istringstream in(bytes, std::ios::binary);
    out.clear();
    BamStatus st = s.open(in, filter);
    if (st != BamStatus::kOk) return st;
    BamFragment f{};
    while ((st = s.next(f)) == BamStatus::kOk) out.push_back(f);
    return st;
}

// ---------------------------------------------------------------------------

void test_header_only_bam_yields_contigs_and_no_fragments() {
    BamStreamer s;
    const std::string bytes =
        bgzf_stream(bam_header({{"chr1", 1000000}, {"chr2", 50000}}));
    std::istringstream in(bytes, std::ios::binary);
    CHECK(s.open(in, {}) == BamStatus::kOk);
    CHECK_EQ(s.contig_names().size(), 2u);
    if (s.contig_names().size() == 2) {
        CHECK(s.contig_names()[0] == "chr1");
        CHECK(s.contig_names()[1] == "chr2");
        CHECK_EQ(s.contig_lengths()[0], 1000000);
        CHECK_EQ(s.contig_lengths()[1], 50000);
    }
    BamFragment f{};
    CHECK(s.next(f) == BamStatus::kEndOfFile);
}

void test_proper_pairs_become_fragments_with_the_tlen_span() {
    std::string payload = bam_header({{"chr1", 1000000}});
    payload += bam_record(0, 100, kGoodPair, 30, 200);
    payload += bam_record(0, 5000, kGoodPair, 30, 150);
    std::vector<BamFragment> got;
    CHECK(collect(bgzf_stream(payload), got) == BamStatus::kEndOfFile);
    CHECK_EQ(got.size(), 2u);
    if (got.size() == 2) {
        CHECK_EQ(got[0].tid, 0);
        CHECK_EQ(got[0].start, 100);
        CHECK_EQ(got[0].end, 300);
        CHECK_EQ(got[1].start, 5000);
        CHECK_EQ(got[1].end, 5150);
    }
}

void test_only_the_leftmost_mate_is_emitted() {
    // Both mates of a pair carry the same fragment. Emitting both would double
    // every pileup, which looks like twice the sequencing depth and shifts
    // every p-value with it. The leftmost mate is the one with TLEN > 0.
    std::string payload = bam_header({{"chr1", 1000000}});
    payload += bam_record(0, 100, kGoodPair, 30, 200);            // leftmost
    payload += bam_record(0, 250, kPaired | kProper | 0x80, 30, -200);  // its mate
    std::vector<BamFragment> got;
    CHECK(collect(bgzf_stream(payload), got) == BamStatus::kEndOfFile);
    CHECK_EQ(got.size(), 1u);
    if (!got.empty()) {
        CHECK_EQ(got[0].start, 100);
        CHECK_EQ(got[0].end, 300);
    }
}

void test_tlen_of_zero_is_skipped() {
    std::string payload = bam_header({{"chr1", 1000000}});
    payload += bam_record(0, 100, kGoodPair, 30, 0);
    std::vector<BamFragment> got;
    CHECK(collect(bgzf_stream(payload), got) == BamStatus::kEndOfFile);
    CHECK_EQ(got.size(), 0u);
}

void test_flag_filters_skip_the_usual_suspects() {
    std::string payload = bam_header({{"chr1", 1000000}});
    payload += bam_record(0, 100, kGoodPair | kUnmapped, 30, 200);
    payload += bam_record(0, 200, kGoodPair | kSecondary, 30, 200);
    payload += bam_record(0, 300, kGoodPair | kSupplementary, 30, 200);
    payload += bam_record(0, 400, kGoodPair | kQcFail, 30, 200);
    payload += bam_record(0, 500, kGoodPair | kDuplicate, 30, 200);
    payload += bam_record(0, 600, kPaired | kRead1, 30, 200);  // not a proper pair
    payload += bam_record(0, 700, kGoodPair, 30, 200);         // the only survivor
    std::vector<BamFragment> got;
    BamStreamer s;
    CHECK(collect(bgzf_stream(payload), got, {}, &s) == BamStatus::kEndOfFile);
    CHECK_EQ(got.size(), 1u);
    if (!got.empty()) CHECK_EQ(got[0].start, 700);
    CHECK_EQ(s.records_seen(), 7);
    CHECK_EQ(s.skipped_flags(), 6);
}

void test_mapq_floor_is_applied() {
    std::string payload = bam_header({{"chr1", 1000000}});
    payload += bam_record(0, 100, kGoodPair, 5, 200);
    payload += bam_record(0, 200, kGoodPair, 30, 200);
    BamFilter filter;
    filter.min_mapq = 20;
    std::vector<BamFragment> got;
    BamStreamer s;
    CHECK(collect(bgzf_stream(payload), got, filter, &s) == BamStatus::kEndOfFile);
    CHECK_EQ(got.size(), 1u);
    if (!got.empty()) CHECK_EQ(got[0].start, 200);
    CHECK_EQ(s.skipped_mapq(), 1);
}

void test_a_span_longer_than_the_ring_is_skipped() {
    // The sliding window cannot represent a fragment longer than
    // kMaxFragmentSpan, and a chimeric TLEN of a megabase is not a fragment.
    // Dropping it here, with a count, beats handing the window a record it
    // will refuse anyway.
    std::string payload = bam_header({{"chr1", 10000000}});
    payload += bam_record(0, 100, kGoodPair, 30, 5000000);
    payload += bam_record(0, 200, kGoodPair, 30, 200);
    std::vector<BamFragment> got;
    BamStreamer s;
    CHECK(collect(bgzf_stream(payload), got, {}, &s) == BamStatus::kEndOfFile);
    CHECK_EQ(got.size(), 1u);
    CHECK_EQ(s.skipped_span(), 1);
}

void test_a_fragment_running_past_the_contig_is_clamped_or_skipped() {
    std::string payload = bam_header({{"chr1", 1000}});
    payload += bam_record(0, 900, kGoodPair, 30, 500);  // ends at 1400 > 1000
    std::vector<BamFragment> got;
    CHECK(collect(bgzf_stream(payload), got) == BamStatus::kEndOfFile);
    // Either dropped or clamped to the contig, never emitted past its end.
    for (const BamFragment& f : got) CHECK(f.end <= 1000);
}

void test_records_spanning_block_boundaries_are_reassembled() {
    // A BAM record is not aligned to BGZF blocks and routinely straddles one.
    // Splitting the fixture at every offset from 1 to 200 walks the boundary
    // through the header, through a record's fixed core, and through its name.
    std::string payload = bam_header({{"chr1", 1000000}});
    for (int i = 0; i < 12; ++i) {
        payload += bam_record(0, 1000 * (i + 1), kGoodPair, 30, 200, "read_name_" + std::to_string(i));
    }
    for (std::size_t chunk = 1; chunk <= 200; ++chunk) {
        std::vector<BamFragment> got;
        const BamStatus st = collect(bgzf_stream(payload, chunk), got);
        CHECK_MSG(st == BamStatus::kEndOfFile,
                  "chunk " + std::to_string(chunk) + ": " + BamStreamer::describe(st));
        CHECK_MSG(got.size() == 12u,
                  "chunk " + std::to_string(chunk) + ": " + std::to_string(got.size()) +
                      " fragments");
        for (std::size_t i = 0; i < got.size(); ++i) {
            CHECK_MSG(got[i].start == static_cast<Coord>(1000 * (i + 1)),
                      "chunk " + std::to_string(chunk) + " frag " + std::to_string(i));
        }
    }
}

// --- malformed input -------------------------------------------------------

void test_bad_bam_magic_is_refused() {
    HeaderTweak t;
    t.bad_magic = true;
    std::vector<BamFragment> got;
    CHECK(collect(bgzf_stream(bam_header({{"chr1", 1000}}, "@HD\n", t)), got) ==
          BamStatus::kBadMagic);
}

void test_a_header_that_overruns_is_refused() {
    HeaderTweak t;
    t.l_text_delta = 5000;  // l_text claims more text than exists
    std::vector<BamFragment> got;
    CHECK(collect(bgzf_stream(bam_header({{"chr1", 1000}}, "@HD\n", t)), got) ==
          BamStatus::kBadHeader);

    HeaderTweak t2;
    t2.n_ref_delta = 100;  // more contigs than bytes
    CHECK(collect(bgzf_stream(bam_header({{"chr1", 1000}}, "@HD\n", t2)), got) ==
          BamStatus::kBadHeader);

    HeaderTweak t3;
    t3.l_name_delta = 4000;  // a contig name longer than the header
    CHECK(collect(bgzf_stream(bam_header({{"chr1", 1000}}, "@HD\n", t3)), got) ==
          BamStatus::kBadHeader);

    HeaderTweak t4;
    t4.n_ref_delta = -2;  // negative n_ref
    CHECK(collect(bgzf_stream(bam_header({{"chr1", 1000}}, "@HD\n", t4)), got) ==
          BamStatus::kBadHeader);
}

void test_an_empty_stream_is_refused() {
    std::vector<BamFragment> got;
    CHECK(collect(std::string{}, got) == BamStatus::kTruncated);
}

void test_a_truncated_bgzf_header_is_refused() {
    std::string bytes = bgzf_stream(bam_header({{"chr1", 1000}}));
    bytes.resize(10);  // less than the 18-byte BGZF header
    std::vector<BamFragment> got;
    CHECK(collect(bytes, got) == BamStatus::kTruncated);
}

void test_a_truncated_bgzf_payload_is_refused() {
    std::string bytes = bgzf_stream(bam_header({{"chr1", 1000}}));
    bytes.resize(bytes.size() - 12);  // cut inside the last block
    std::vector<BamFragment> got;
    const BamStatus st = collect(bytes, got);
    CHECK_MSG(st == BamStatus::kTruncated, BamStreamer::describe(st));
}

void test_bad_bgzf_magic_is_refused() {
    BlockTweak t;
    t.bad_magic = true;
    std::vector<BamFragment> got;
    CHECK(collect(bgzf_block(bam_header({{"chr1", 1000}}), t) + bgzf_eof(), got) ==
          BamStatus::kBadBlock);
}

void test_a_block_without_the_bc_field_is_refused() {
    BlockTweak t;
    t.drop_extra = true;
    std::vector<BamFragment> got;
    CHECK(collect(bgzf_block(bam_header({{"chr1", 1000}}), t) + bgzf_eof(), got) ==
          BamStatus::kBadBlock);
}

void test_a_bsize_that_overruns_the_file_is_refused() {
    BlockTweak t;
    t.bsize_delta = 4000;
    std::vector<BamFragment> got;
    const BamStatus st = collect(bgzf_block(bam_header({{"chr1", 1000}}), t) + bgzf_eof(), got);
    CHECK_MSG(st == BamStatus::kTruncated || st == BamStatus::kBadBlock,
              BamStreamer::describe(st));
}

void test_an_isize_past_the_block_limit_is_refused() {
    // A BGZF block inflates to at most 64 KiB by definition. A larger ISIZE is
    // either corruption or an attempt to make the reader allocate.
    BlockTweak t;
    t.use_isize_override = true;
    t.isize_override = 70000;
    std::vector<BamFragment> got;
    CHECK(collect(bgzf_block(bam_header({{"chr1", 1000}}), t) + bgzf_eof(), got) ==
          BamStatus::kBadBlock);
}

void test_a_crc_mismatch_is_refused() {
    // The case that matters most and shows least: a block whose bytes decode
    // without complaint into the wrong alignments.
    BlockTweak t;
    t.use_crc_override = true;
    t.crc_override = 0xDEADBEEF;
    std::vector<BamFragment> got;
    CHECK(collect(bgzf_block(bam_header({{"chr1", 1000}}), t) + bgzf_eof(), got) ==
          BamStatus::kBadChecksum);
}

void test_a_record_that_overruns_the_stream_is_refused() {
    std::string payload = bam_header({{"chr1", 1000000}});
    payload += bam_record(0, 100, kGoodPair, 30, 200, "r", 4000);  // block_size too big
    std::vector<BamFragment> got;
    const BamStatus st = collect(bgzf_stream(payload), got);
    CHECK_MSG(st == BamStatus::kBadRecord || st == BamStatus::kTruncated,
              BamStreamer::describe(st));

    std::string neg = bam_header({{"chr1", 1000000}});
    neg += bam_record(0, 100, kGoodPair, 30, 200, "r", -100);  // negative block_size
    CHECK(collect(bgzf_stream(neg), got) == BamStatus::kBadRecord);
}

void test_a_refid_out_of_range_is_refused() {
    std::string payload = bam_header({{"chr1", 1000000}});
    payload += bam_record(7, 100, kGoodPair, 30, 200);  // only one contig exists
    std::vector<BamFragment> got;
    CHECK(collect(bgzf_stream(payload), got) == BamStatus::kBadRecord);
}

void test_unplaced_records_are_skipped_not_refused() {
    // refID -1 is "unplaced", which is legal and appears at the end of a
    // coordinate-sorted BAM. It is not an error; it is simply not a fragment.
    std::string payload = bam_header({{"chr1", 1000000}});
    payload += bam_record(0, 100, kGoodPair, 30, 200);
    payload += bam_record(-1, -1, kGoodPair | kUnmapped, 0, 0);
    std::vector<BamFragment> got;
    CHECK(collect(bgzf_stream(payload), got) == BamStatus::kEndOfFile);
    CHECK_EQ(got.size(), 1u);
}

void test_a_negative_position_is_refused() {
    std::string payload = bam_header({{"chr1", 1000000}});
    payload += bam_record(0, -5, kGoodPair, 30, 200);
    std::vector<BamFragment> got;
    CHECK(collect(bgzf_stream(payload), got) == BamStatus::kBadRecord);
}

void test_a_file_without_the_eof_marker_is_reported_truncated() {
    // samtools quickcheck's whole job. A BAM cut at a block boundary decodes
    // perfectly and is missing the end of the chromosome.
    std::string payload = bam_header({{"chr1", 1000000}});
    payload += bam_record(0, 100, kGoodPair, 30, 200);
    const std::string full = bgzf_stream(payload);
    const std::string cut = full.substr(0, full.size() - bgzf_eof().size());
    std::vector<BamFragment> got;
    const BamStatus st = collect(cut, got);
    CHECK_MSG(st == BamStatus::kTruncated, BamStreamer::describe(st));
    // The records before the cut are still delivered: the diagnosis is at the
    // end, not a refusal to read anything.
    CHECK_EQ(got.size(), 1u);
}

void test_a_record_with_an_extreme_l_seq_is_refused() {
    // l_seq is a signed int32 read straight off the wire, and the `l_seq < 0`
    // guard admits INT32_MAX. Computing (l_seq + 1) / 2 in int32 then
    // overflows -- undefined behaviour, and it happens before any bounds test
    // can look at the result. Run this under UBSan: a pass with no sanitizer
    // proves only that the machine wrapped quietly.
    for (std::int32_t l_seq : {std::numeric_limits<std::int32_t>::max(),
                               std::numeric_limits<std::int32_t>::max() - 1,
                               1 << 30, -1}) {
        std::string payload = bam_header({{"chr1", 1000000}});
        payload += bam_record(0, 100, kGoodPair, 30, 200, "r", 0, 0, l_seq);
        std::vector<BamFragment> got;
        const BamStatus st = collect(bgzf_stream(payload), got);
        CHECK_MSG(st == BamStatus::kBadRecord,
                  "l_seq " + std::to_string(l_seq) + ": " + BamStreamer::describe(st));
        CHECK_EQ(got.size(), 0u);
    }
}

void test_a_record_at_the_advertised_size_limit_is_accepted() {
    // The reader publishes kMaxRecordSize and has to honour it. The buffer must
    // hold the assembled record AND the block being inflated onto the end of
    // it, so a limit computed with one block of slack rather than two rejects a
    // record of exactly the advertised size -- as "malformed", which is a
    // diagnosis of the data for a defect in the reader.
    const std::size_t limit = pto::peaks::detail::kMaxRecordSize;
    const std::size_t core = 32 + 2;  // fixed part plus "r\0"

    std::string payload = bam_header({{"chr1", 1000000}});
    payload += bam_record(0, 100, kGoodPair, 30, 200, "r", 0, limit - core);
    std::vector<BamFragment> got;
    const BamStatus st = collect(bgzf_stream(payload), got);
    CHECK_MSG(st == BamStatus::kEndOfFile, BamStreamer::describe(st));
    CHECK_EQ(got.size(), 1u);

    // One byte over is refused, and refused as a bad record rather than as a
    // truncation: the reader knows its own limit.
    std::string over = bam_header({{"chr1", 1000000}});
    over += bam_record(0, 100, kGoodPair, 30, 200, "r", 0, limit - core + 1);
    CHECK(collect(bgzf_stream(over), got) == BamStatus::kBadRecord);
}

void test_an_empty_block_mid_stream_is_transparent() {
    // An empty BGZF block is legal anywhere in a stream -- the end-of-file
    // marker IS one, and concatenating BGZF streams puts them in the middle.
    // Treating the first one as the end of the file drops every record after
    // it and reports success, which is the silent wrong answer this module is
    // built against.
    std::string head = bam_header({{"chr1", 1000000}});
    for (int i = 0; i < 5; ++i) head += bam_record(0, 1000 * (i + 1), kGoodPair, 30, 200);
    std::string tail;
    for (int i = 5; i < 10; ++i) tail += bam_record(0, 1000 * (i + 1), kGoodPair, 30, 200);

    std::vector<BamFragment> got;
    CHECK(collect(bgzf_block(head) + bgzf_block(tail) + bgzf_eof(), got) ==
          BamStatus::kEndOfFile);
    CHECK_MSG(got.size() == 10u, "baseline: " + std::to_string(got.size()));

    // One empty block in the middle changes nothing.
    CHECK(collect(bgzf_block(head) + bgzf_eof() + bgzf_block(tail) + bgzf_eof(), got) ==
          BamStatus::kEndOfFile);
    CHECK_MSG(got.size() == 10u,
              "one mid-stream empty block dropped " + std::to_string(10 - got.size()) +
                  " fragments");

    // Nor do several.
    CHECK(collect(bgzf_block(head) + bgzf_eof() + bgzf_eof() + bgzf_eof() +
                      bgzf_block(tail) + bgzf_eof(),
                  got) == BamStatus::kEndOfFile);
    CHECK_MSG(got.size() == 10u,
              "three mid-stream empty blocks dropped " +
                  std::to_string(10 - got.size()) + " fragments");

    // And truncation is still detected afterwards: the marker is evidence
    // about the LAST block, not about any block, so it must not be sticky.
    CHECK(collect(bgzf_block(head) + bgzf_eof() + bgzf_block(tail), got) ==
          BamStatus::kTruncated);
    CHECK_EQ(got.size(), 10u);
}

void test_streaming_allocates_nothing_per_record() {
    std::string payload = bam_header({{"chr1", 10000000}});
    for (int i = 0; i < 4000; ++i) {
        payload += bam_record(0, 100 * (i + 1), kGoodPair, 30, 200);
    }
    const std::string bytes = bgzf_stream(payload);

    BamStreamer s;
    std::istringstream in(bytes, std::ios::binary);
    CHECK(s.open(in, {}) == BamStatus::kOk);
    {
        testing::AllocationGuard guard;
        BamFragment f{};
        long long n = 0;
        while (s.next(f) == BamStatus::kOk) ++n;
        const long long allocs = guard.count();
        CHECK_EQ(n, 4000);
        if (testing::AllocationGuard::counting()) {
            CHECK_MSG(allocs == 0,
                      std::to_string(allocs) + " allocation(s) while streaming records");
        } else {
            testing::skip("allocation counting is off under sanitizers");
        }
    }
}

void test_a_contig_name_cannot_forge_a_track_row() {
    // The names this reader publishes are written straight into a
    // TAB-separated narrowPeak by the CLI, so a tab or a newline in one is not
    // a cosmetic problem: it forges whole rows in a file a downstream pipeline
    // parses, and a BAM header is the least trustworthy part of the file. The
    // reader is where that has to stop -- every other consumer of
    // contig_names() would otherwise have to know.
    const std::string forged = "chrA\t0\t9\tINJECTED\t1000\t.\t9\t9\t9\t0\nchrB";
    std::vector<BamFragment> got;
    for (const std::string& name : {forged, std::string("chr\t1"),
                                    std::string("chr\n1"), std::string("chr\r1"),
                                    std::string("chr\x7f" "1"),
                                    std::string("chr\x01" "1")}) {
        const BamStatus st = collect(bgzf_stream(bam_header({{name, 1000}})), got);
        CHECK_MSG(st == BamStatus::kBadHeader,
                  "name with a control byte: " + std::string(BamStreamer::describe(st)));
    }

    // A name whose declared length does not end in the NUL the format promises
    // is refused too: the stored length trims a byte on the strength of it.
    std::string unterminated;
    unterminated += "BAM\1";
    put_i32(unterminated, 0);   // l_text
    put_i32(unterminated, 1);   // n_ref
    put_i32(unterminated, 4);   // l_name, counting a terminator that is not there
    unterminated += "chr1";
    put_i32(unterminated, 1000);
    CHECK(collect(bgzf_stream(unterminated), got) == BamStatus::kBadHeader);

    // Printable names are untouched, however long: the rule is about the
    // charset, not the length, and refusing a valid file is the worse error.
    BamStreamer wide;
    const std::string longish(4000, 'S');
    std::istringstream in(bgzf_stream(bam_header({{longish, 1000}})), std::ios::binary);
    CHECK(wide.open(in, {}) == BamStatus::kOk);
    CHECK_EQ(wide.contig_names().size(), 1u);
    CHECK(wide.contig_names()[0] == longish);
}

void test_a_record_cut_short_of_its_length_prefix_is_truncated() {
    // A payload that ends 1-3 bytes into a record cannot even form the length
    // prefix, and the BGZF framing around it can still be perfect -- including
    // the end-of-file marker. Reporting a clean end there drops the record AND
    // says the file was whole, which is the failure the marker check exists to
    // catch, one level down.
    std::string base = bam_header({{"chr1", 100000}});
    base += bam_record(0, 1000, kGoodPair, 60, 200);
    for (std::size_t tail = 1; tail <= 3; ++tail) {
        std::vector<BamFragment> got;
        const BamStatus st =
            collect(bgzf_stream(base + std::string(tail, '\x7f')), got);
        CHECK_MSG(st == BamStatus::kTruncated,
                  std::to_string(tail) + " trailing byte(s): " + BamStreamer::describe(st));
        CHECK_EQ(got.size(), 1u);  // what did decode is still delivered
    }
    // Nothing left over is still a clean end.
    std::vector<BamFragment> got;
    CHECK(collect(bgzf_stream(base), got) == BamStatus::kEndOfFile);
    CHECK_EQ(got.size(), 1u);
}

void test_a_header_that_stops_with_the_file_is_truncated_not_malformed() {
    // Two different faults wear the same shape. A header whose lengths run
    // past a file that ends TIDILY is a malformed header: the bytes are all
    // there and they do not add up. A header that runs past a file which
    // simply stops is a truncated download, and answering "malformed BAM
    // header" to that sends the operator to the wrong problem.
    const std::string header = bam_header({{"chr1", 100000}, {"chr2", 5000}});
    std::vector<BamFragment> got;
    for (std::size_t cut : {4u, 8u, 30u, 40u}) {
        const std::string prefix = header.substr(0, cut);
        // ... with the end-of-file marker: the header is what is wrong.
        CHECK_MSG(collect(bgzf_stream(prefix), got) == BamStatus::kBadHeader,
                  "cut " + std::to_string(cut) + " with an EOF marker");
        // ... without it: the FILE is what is wrong.
        CHECK_MSG(collect(bgzf_block(prefix), got) == BamStatus::kTruncated,
                  "cut " + std::to_string(cut) + " with no EOF marker");
    }
}

void test_an_unopened_or_moved_from_streamer_refuses_rather_than_faults() {
    // The move operations are defaulted, which steals the two buffers -- and
    // every path in the reader indexes them at their constructed sizes, so a
    // moved-from object writes a whole BGZF block through a null pointer. The
    // class advertises movability, so it has to survive being used that way.
    BamFragment f{};
    {
        BamStreamer unopened;
        CHECK(unopened.next(f) == BamStatus::kUnsupported);
    }
    {
        BamStreamer source;
        BamStreamer sink{std::move(source)};
        std::istringstream in(bgzf_stream(bam_header({{"chr1", 1000}})), std::ios::binary);
        CHECK(source.open(in, {}) == BamStatus::kUnsupported);
        CHECK(source.next(f) == BamStatus::kUnsupported);
        // The object that was moved TO is fully working.
        std::istringstream in2(bgzf_stream(bam_header({{"chr1", 1000}})), std::ios::binary);
        CHECK(sink.open(in2, {}) == BamStatus::kOk);
    }
    {
        // Move ASSIGNMENT is what collect() above does on every case in this
        // file. It must leave a usable object -- and, on the zlib backend, must
        // run inflateEnd on the stream it replaces rather than leaking it.
        BamStreamer sink;
        sink = BamStreamer{};
        std::istringstream in(bgzf_stream(bam_header({{"chr1", 1000}})), std::ios::binary);
        CHECK(sink.open(in, {}) == BamStatus::kOk);
    }
}

void test_gzip_flag_bits_beyond_fextra_are_refused() {
    // BGZF fixes FLG at FEXTRA alone. FHCRC, FNAME and FCOMMENT each insert
    // bytes between the extra field and the deflate payload, so tolerating the
    // bit means handing those bytes to the inflater and then blaming the
    // result on the compressed data.
    for (std::uint8_t bit : {0x02, 0x08, 0x10}) {
        BlockTweak t;
        t.extra_flg = bit;
        std::string bytes = bgzf_block(bam_header({{"chr1", 1000}}), t);
        bytes += bgzf_eof();
        std::vector<BamFragment> got;
        const BamStatus st = collect(bytes, got);
        CHECK_MSG(st == BamStatus::kBadBlock,
                  "FLG bit " + std::to_string(bit) + ": " + BamStreamer::describe(st));
    }
}

void test_single_end_records_are_counted_as_unpaired_not_by_flag() {
    // Every record in a single-end BAM has TLEN 0. Folding that into
    // skipped_flags_ reports the whole file as rejected "by flag", which is the
    // one number that sends the operator to their flag filter -- where there is
    // nothing to find.
    std::string payload = bam_header({{"chr1", 100000}});
    for (int i = 0; i < 10; ++i) payload += bam_record(0, 1000 + i, kGoodPair, 60, 0);
    BamStreamer s;
    std::vector<BamFragment> got;
    CHECK(collect(bgzf_stream(payload), got, BamFilter{}, &s) == BamStatus::kEndOfFile);
    CHECK_EQ(got.size(), 0u);
    CHECK_EQ(s.skipped_unpaired(), 10);
    CHECK_EQ(s.skipped_flags(), 0);
}

// --- torture regressions, 2026-09-10 ----------------------------------------
//
// The BED front end refuses an empty contig name, a duplicate one and a length
// of zero (ContigTable::load). The BAM front end accepted all three from the
// header, which is the least trustworthy part of the file: the names are
// published through contig_names() into column 1 of a narrowPeak, where an
// empty one is a row no browser reads and a duplicate makes two contigs'
// peaks indistinguishable. The SAM spec agrees on each: @SQ SN must be
// non-empty and unique, and LN is in [1, 2^31-1].

void test_an_empty_contig_name_is_refused() {
    std::vector<BamFragment> f;
    const std::string bam = bgzf_stream(bam_header({{"chr1", 1000}, {"", 1000}}) +
                                        bam_record(0, 10, kGoodPair, 60, 100));
    CHECK(collect(bam, f) == BamStatus::kBadHeader);
    CHECK(f.empty());
}

void test_a_duplicate_contig_name_is_refused() {
    std::vector<BamFragment> f;
    const std::string bam =
        bgzf_stream(bam_header({{"chr1", 1000}, {"chr2", 500}, {"chr1", 2000}}) +
                    bam_record(0, 10, kGoodPair, 60, 100));
    CHECK(collect(bam, f) == BamStatus::kBadHeader);
    CHECK(f.empty());
    // Names that differ only by case are different contigs, as in SAM.
    CHECK(collect(bgzf_stream(bam_header({{"chrM", 100}, {"chrm", 100}})), f) ==
          BamStatus::kEndOfFile);
}

void test_a_zero_length_contig_is_refused() {
    std::vector<BamFragment> f;
    CHECK(collect(bgzf_stream(bam_header({{"chr1", 1000}, {"chrZ", 0}})), f) ==
          BamStatus::kBadHeader);
    // The longest length the format allows is still fine.
    CHECK(collect(bgzf_stream(bam_header({{"chrBig", std::numeric_limits<std::int32_t>::max()}})),
                  f) == BamStatus::kEndOfFile);
}

void test_describe_names_every_status() {
    const BamStatus all[] = {BamStatus::kOk,          BamStatus::kEndOfFile,
                             BamStatus::kBadMagic,    BamStatus::kTruncated,
                             BamStatus::kBadBlock,    BamStatus::kBadChecksum,
                             BamStatus::kBadHeader,   BamStatus::kBadRecord,
                             BamStatus::kUnsupported};
    for (BamStatus st : all) {
        const char* text = BamStreamer::describe(st);
        CHECK(text != nullptr);
        CHECK(std::strlen(text) > 0);
    }
}

}  // namespace

// The CLI integration script needs a real BAM, and the byte builders above are
// the only BAM writer in this repository. `--emit-fixture` is how it gets one:
// a background tiling plus a planted peak, written to stdout. Keeping it here
// rather than in a separate tool means the bytes the CLI is driven with are the
// same bytes these cases prove the reader accepts.
int emit_fixture() {
    std::string payload = bam_header({{"chr1", 100000}, {"chr2", 50000}});
    std::vector<std::pair<std::int32_t, std::int32_t>> frags;
    for (std::int32_t s = 0; s + 200 <= 99800; s += 200) frags.emplace_back(s, 200);
    for (int i = 0; i < 300; ++i) frags.emplace_back(50000, 200);
    std::stable_sort(frags.begin(), frags.end(),
                     [](const auto& a, const auto& b) { return a.first < b.first; });
    for (const auto& [start, span] : frags) {
        payload += bam_record(0, start, kGoodPair, 30, span);
    }
    const std::string bytes = bgzf_stream(payload, 8192);
    std::fwrite(bytes.data(), 1, bytes.size(), stdout);
    return 0;
}

int main(int argc, char** argv) {
    if (argc > 1 && std::string(argv[1]) == "--emit-fixture") return emit_fixture();
    testing::Suite suite{"peaks/bam_streamer", {}};
    suite.add("header-only BAM yields contigs", test_header_only_bam_yields_contigs_and_no_fragments);
    suite.add("proper pairs become TLEN-spanned fragments", test_proper_pairs_become_fragments_with_the_tlen_span);
    suite.add("only the leftmost mate is emitted", test_only_the_leftmost_mate_is_emitted);
    suite.add("TLEN of zero is skipped", test_tlen_of_zero_is_skipped);
    suite.add("flag filters skip the usual suspects", test_flag_filters_skip_the_usual_suspects);
    suite.add("mapq floor is applied", test_mapq_floor_is_applied);
    suite.add("a span longer than the ring is skipped", test_a_span_longer_than_the_ring_is_skipped);
    suite.add("a fragment past the contig end is not emitted", test_a_fragment_running_past_the_contig_is_clamped_or_skipped);
    suite.add("records spanning block boundaries are reassembled", test_records_spanning_block_boundaries_are_reassembled);
    suite.add("bad BAM magic is refused", test_bad_bam_magic_is_refused);
    suite.add("a header that overruns is refused", test_a_header_that_overruns_is_refused);
    suite.add("an empty stream is refused", test_an_empty_stream_is_refused);
    suite.add("a truncated BGZF header is refused", test_a_truncated_bgzf_header_is_refused);
    suite.add("a truncated BGZF payload is refused", test_a_truncated_bgzf_payload_is_refused);
    suite.add("bad BGZF magic is refused", test_bad_bgzf_magic_is_refused);
    suite.add("a block without the BC field is refused", test_a_block_without_the_bc_field_is_refused);
    suite.add("a BSIZE that overruns the file is refused", test_a_bsize_that_overruns_the_file_is_refused);
    suite.add("an ISIZE past 64 KiB is refused", test_an_isize_past_the_block_limit_is_refused);
    suite.add("a CRC mismatch is refused", test_a_crc_mismatch_is_refused);
    suite.add("a record that overruns is refused", test_a_record_that_overruns_the_stream_is_refused);
    suite.add("a refID out of range is refused", test_a_refid_out_of_range_is_refused);
    suite.add("unplaced records are skipped not refused", test_unplaced_records_are_skipped_not_refused);
    suite.add("a negative position is refused", test_a_negative_position_is_refused);
    suite.add("a file without the EOF marker is truncated", test_a_file_without_the_eof_marker_is_reported_truncated);
    suite.add("an extreme l_seq is refused", test_a_record_with_an_extreme_l_seq_is_refused);
    suite.add("a record at the size limit is accepted", test_a_record_at_the_advertised_size_limit_is_accepted);
    suite.add("an empty block mid-stream is transparent", test_an_empty_block_mid_stream_is_transparent);
    suite.add("streaming allocates nothing per record", test_streaming_allocates_nothing_per_record);
    suite.add("a contig name cannot forge a track row", test_a_contig_name_cannot_forge_a_track_row);
    suite.add("a record cut short of its length prefix is truncated", test_a_record_cut_short_of_its_length_prefix_is_truncated);
    suite.add("a header that stops with the file is truncated", test_a_header_that_stops_with_the_file_is_truncated_not_malformed);
    suite.add("an unopened or moved-from streamer refuses", test_an_unopened_or_moved_from_streamer_refuses_rather_than_faults);
    suite.add("gzip FLG bits beyond FEXTRA are refused", test_gzip_flag_bits_beyond_fextra_are_refused);
    suite.add("single-end records count as unpaired", test_single_end_records_are_counted_as_unpaired_not_by_flag);
    suite.add("an empty contig name is refused", test_an_empty_contig_name_is_refused);
    suite.add("a duplicate contig name is refused", test_a_duplicate_contig_name_is_refused);
    suite.add("a zero-length contig is refused", test_a_zero_length_contig_is_refused);
    suite.add("describe names every status", test_describe_names_every_status);
    return suite.run();
}

#endif  // PEAKS_HAVE_BGZF
