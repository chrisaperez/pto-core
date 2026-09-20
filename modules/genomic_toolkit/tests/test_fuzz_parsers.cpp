// SPDX-License-Identifier: MIT
//
// Adversarial parser suite: structured fuzzing plus the regressions it found.
//
// The other four suites in this module test what the parsers are supposed to
// accept. This one tests what they are supposed to survive. It is a
// deterministic, seeded mutation fuzzer rather than a libFuzzer target on
// purpose: libFuzzer's runtime is not shipped with Apple clang, and a harness
// that only runs where a coverage-guided engine is installed is a harness that
// does not run in CI on this project's own machines. A fixed seed also means a
// failure here reproduces exactly, which a coverage-guided run does not.
//
// The invariants asserted after every input are the ones the rest of the
// module's code assumes without checking:
//
//   * every emitted Fragment has 0 <= start <= end <= kMaxCoord
//   * every emitted Fragment has length() > 0 and satisfies the filter
//   * every tid indexes the dictionary that was passed in
//   * fragments_emitted equals the number actually handed to the callback
//   * the parser either returns or throws; it never corrupts memory
//
// Runtime: ~4 s in a normal build, ~90 s under ASan+UBSan (the cost is the
// temp file per input, not the parsing). Lines per input are batched for that
// reason -- fewer, larger files cover the same grammar for a third of the I/O.
//
// Run under ASan+UBSan to get the memory-safety half of that:
//   cmake -S modules/genomic_toolkit -B build/gtk-asan
//        -DCMAKE_BUILD_TYPE=Debug
//        -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-sanitize-recover=all"

#include <cstdint>
#include <cstdio>
#include <random>
#include <span>
#include <string>
#include <vector>

#include "test_util.hpp"
#include "toolkit/contig_dict.hpp"
#include "toolkit/fragment_stream.hpp"
#include "toolkit/types.hpp"

using namespace toolkit;
using testing::TempFile;

namespace {

// Result of streaming one input, plus everything the callback saw.
struct Parsed {
    StreamStats stats;
    std::vector<Fragment> fragments;
    bool threw = false;
    std::string error;
};

Parsed parse(const std::string& body, const std::string& suffix,
             const FragmentFilter& filter) {
    TempFile file(suffix, body);
    Parsed out;
    ContigDict dict;
    try {
        auto sink = [&](std::span<Fragment> batch) {
            for (const Fragment& f : batch) out.fragments.push_back(f);
        };
        out.stats = suffix == ".bedpe"
                        ? stream_bedpe(file.path(), filter, dict, sink)
                        : stream_fragment_bed(file.path(), filter, dict, sink);
    } catch (const std::exception& ex) {
        out.threw = true;
        out.error = ex.what();
        return out;
    }

    // Invariants. A parser is allowed to reject anything; it is not allowed to
    // emit a fragment that violates the type's contract.
    //
    // Folded into one check per input rather than seven per fragment: the fuzz
    // cases below stream millions of fragments, and a per-fragment CHECK makes
    // the summary line unreadable (and used to overflow the counter). The
    // offending record is reported by index when the fold fails, so a failure
    // is still pinpointed.
    std::size_t bad = out.fragments.size();  // sentinel: none bad
    for (std::size_t i = 0; i < out.fragments.size(); ++i) {
        const Fragment& f = out.fragments[i];
        const bool ok = f.start >= 0 && f.end >= f.start &&
                        static_cast<std::int64_t>(f.end) <= kMaxCoord &&
                        f.length() > 0 && filter.length_ok(f.length()) &&
                        f.tid >= 0 && static_cast<std::size_t>(f.tid) < dict.size();
        if (!ok) { bad = i; break; }
    }
    ::testing::report(bad == out.fragments.size(), "all emitted fragments valid",
                      __FILE__, __LINE__,
                      bad == out.fragments.size()
                          ? std::string{}
                          : "first bad fragment at index " + std::to_string(bad));
    CHECK_EQ(out.stats.fragments_emitted, out.fragments.size());
    return out;
}

// --------------------------------------------------------------------------
// Regression: BEDPE MAPQ failed open on a corrupt column.
//
// stream_bedpe called parse_int(fields[7], mapq) and discarded the result.
// Two inputs slipped a corrupt record past --min-mapq as a maximum-confidence
// alignment:
//
//   "\tXX"     parse fails, mapq keeps its 255 initialiser
//   "\t60abc"  std::from_chars writes the parsed prefix *before* reporting the
//              trailing junk, so the discarded `false` left mapq = 60
//
// The sibling `copies` column in stream_fragment_bed already checked its
// parse; this is the same check in the function that was missed. An absent
// column, and the BEDPE "." placeholder, still mean "not recorded" and pass.
// --------------------------------------------------------------------------
void test_bedpe_mapq_fails_closed() {
    FragmentFilter filter;
    filter.max_length = 0;
    filter.min_mapq = 30;

    struct Case {
        const char* mapq_column;  // nullptr = column absent entirely
        bool expect_emitted;
        const char* why;
    };
    const Case cases[] = {
        {"60", true, "a real MAPQ above the threshold passes"},
        {"10", false, "a real MAPQ below the threshold is dropped as mapq"},
        {nullptr, true, "an absent column means 'not recorded', not 'bad'"},
        {".", true, "the BEDPE placeholder means 'not recorded'"},
        {"XX", false, "a non-numeric MAPQ is a corrupt record, not Q255"},
        {"60abc", false, "a partially numeric MAPQ is corrupt, not Q60"},
        {"", false, "an empty-but-present column is corrupt, not Q255"},
        {"-5", false, "a negative MAPQ cannot pass a positive threshold"},
        {"256", false, "a MAPQ past 255 is corrupt, not a high-confidence alignment"},
    };

    for (const Case& c : cases) {
        const std::string line =
            c.mapq_column == nullptr
                ? std::string("c\t0\t100\tc\t200\t300\tname\n")
                : std::string("c\t0\t100\tc\t200\t300\tname\t") + c.mapq_column + "\n";
        const Parsed p = parse(line, ".bedpe", filter);
        CHECK(!p.threw);
        const bool emitted = p.stats.fragments_emitted == 1;
        ::testing::report(emitted == c.expect_emitted, c.why, __FILE__, __LINE__,
                          std::string("emitted=") + (emitted ? "1" : "0"));
    }

    // The two corrupt forms must be accounted as malformed, not silently
    // dropped into some other bucket, so a user reading the report can see
    // that their file is broken rather than that their data is bad.
    const Parsed corrupt = parse(
        "c\t0\t100\tc\t200\t300\tname\tXX\n"
        "c\t0\t100\tc\t200\t300\tname\t60abc\n",
        ".bedpe", filter);
    CHECK_EQ(corrupt.stats.records_read, 2u);
    CHECK_EQ(corrupt.stats.malformed_lines, 2u);
    CHECK_EQ(corrupt.stats.fragments_emitted, 0u);
}

// --------------------------------------------------------------------------
// Regression: LineReader grew its buffer without bound.
//
// One unterminated 200 MB line drove peak RSS to 518 MB before the cap. The
// input below is deliberately just over the ceiling so the test costs ~64 MiB
// and a second, rather than reproducing the original at multi-gigabyte scale.
// --------------------------------------------------------------------------
void test_oversized_line_fails_closed() {
    FragmentFilter filter;
    filter.max_length = 0;

    // Comfortably under the cap: must parse normally.
    {
        std::string body = "c\t0\t100\t";
        body.append(4u << 20, 'x');  // 4 MiB of trailing column
        body += "\n";
        const Parsed p = parse(body, ".bed", filter);
        CHECK(!p.threw);
        CHECK_EQ(p.stats.fragments_emitted, 1u);
    }

    // Past the cap: must throw, not grow.
    {
        std::string body = "c\t0\t100\t";
        body.append((68u << 20), 'x');  // 68 MiB > 64 MiB ceiling
        body += "\n";
        const Parsed p = parse(body, ".bed", filter);
        CHECK(p.threw);
        CHECK(p.error.find("longer than") != std::string::npos);
    }
}

// --------------------------------------------------------------------------
// Coordinate boundaries: the int32 ceiling must be refused, never truncated.
// A truncated coordinate silently relocates a read to the wrong part of the
// chromosome, which is worse than refusing to run.
// --------------------------------------------------------------------------
void test_coordinate_boundaries() {
    FragmentFilter filter;
    filter.max_length = 0;

    const std::string body =
        "c\t0\t100\n"                                     // ok
        "c\t2147483547\t2147483647\n"                     // ok, ends at ceiling
        "c\t0\t2147483648\n"                              // one past the ceiling
        "c\t0\t9223372036854775807\n"                     // INT64_MAX
        "c\t0\t9223372036854775808\n"                     // past INT64_MAX
        "c\t-1\t100\n"                                    // negative start
        "c\t100\t50\n"                                    // inverted
        "c\t0x10\t100\n"                                  // hex
        "c\t1e5\t100\n"                                   // scientific
        "c\t+5\t100\n"                                    // leading sign
        "c\t 5\t100\n"                                    // leading space
        "c\t100abc\t200\n"                                // trailing junk
        "c\t\t100\n";                                     // empty field

    const Parsed p = parse(body, ".bed", filter);
    CHECK(!p.threw);
    CHECK_EQ(p.stats.records_read, 13u);
    CHECK_EQ(p.stats.fragments_emitted, 2u);   // only the two valid lines
    CHECK_EQ(p.stats.malformed_lines, 11u);
    if (p.fragments.size() == 2) {
        CHECK_EQ(p.fragments[1].end, static_cast<Coord>(kMaxCoord));
    }
}

// --------------------------------------------------------------------------
// Degenerate and hostile file shapes.
// --------------------------------------------------------------------------
void test_degenerate_files() {
    FragmentFilter filter;
    filter.max_length = 0;

    struct Case { const char* label; std::string body; };
    std::vector<Case> cases = {
        {"empty", ""},
        {"newlines only", "\n\n\n"},
        {"no trailing newline", "c\t0\t100"},
        {"CRLF", "c\t0\t100\r\nc\t10\t20\r\n"},
        {"comments and track lines", "#h\ntrack name=x\nbrowser position y\n"},
        {"embedded NUL", std::string("c\t0\t1") + '\0' + "00\n"},
        {"bare tabs", "\t\t\t\t\t\n"},
        {"one huge contig name", "c" + std::string(100000, 'x') + "\t0\t100\n"},
        {"many tabs", "c\t0\t100" + std::string(4096, '\t') + "\n"},
        {"binary noise", std::string("\x7f\x45\x4c\x46\x02\x01\x01\x00\n\xff\xfe\n")},
    };

    for (const Case& c : cases) {
        for (const char* suffix : {".bed", ".bedpe"}) {
            const Parsed p = parse(c.body, suffix, filter);
            // Any outcome is acceptable except a violated invariant, which
            // parse() has already checked. Assert only that it terminated in a
            // defined way and accounted for what it read.
            if (!p.threw) {
                CHECK(p.stats.malformed_lines <= p.stats.records_read);
            }
        }
    }
}

// --------------------------------------------------------------------------
// The fuzzer proper: mutate valid records with adversarial field values.
// --------------------------------------------------------------------------
const char* kNastyFields[] = {
    "", " ", "\t", ".", "-", "+", "-1", "0", "1",
    "2147483647", "2147483648", "-2147483648",
    "9223372036854775807", "9223372036854775808", "-9223372036854775809",
    "99999999999999999999999999", "0x10", "1e5", "1E-5", "nan", "inf", "NaN",
    "010", "  7", "7  ", "7abc", "abc7", "++1", "--1", "1.5", ".5", "1,000",
    "\xff\xfe", "chr1", "*",
};

std::string nasty(std::mt19937& rng) {
    return kNastyFields[rng() % (sizeof(kNastyFields) / sizeof(kNastyFields[0]))];
}

void test_fuzz_fragment_bed() {
    std::mt19937 rng(20260817u);
    FragmentFilter filter;

    for (int trial = 0; trial < 1200; ++trial) {
        // Vary the filter too: length bounds and MAPQ gates are part of the
        // surface, and a filter that rejects everything is as interesting as
        // one that accepts everything.
        filter.min_length = static_cast<Coord>(rng() % 3 == 0 ? rng() % 500 : 0);
        filter.max_length = static_cast<Coord>(rng() % 3 == 0 ? rng() % 2000 : 0);
        if (filter.max_length > 0 && filter.min_length > filter.max_length) {
            std::swap(filter.min_length, filter.max_length);
        }
        filter.min_mapq = static_cast<int>(rng() % 4 == 0 ? rng() % 260 : 0);
        // Both readings of column 5 are part of the surface.
        filter.use_count_column = rng() % 2 == 0;

        std::string body;
        const int lines = static_cast<int>(rng() % 40);
        for (int i = 0; i < lines; ++i) {
            const int shape = static_cast<int>(rng() % 4);
            if (shape == 0) {
                // Well-formed.
                const long s = static_cast<long>(rng() % 1000000);
                body += "chr" + std::to_string(rng() % 3) + "\t" + std::to_string(s) +
                        "\t" + std::to_string(s + 1 + (rng() % 800)) + "\n";
            } else if (shape == 1) {
                // Well-formed shape, adversarial values.
                body += nasty(rng) + "\t" + nasty(rng) + "\t" + nasty(rng) + "\n";
            } else if (shape == 2) {
                // Extra columns, including the 10x count column.
                body += "chr1\t10\t200\t" + nasty(rng) + "\t" + nasty(rng) + "\n";
            } else {
                // Wrong arity.
                const int fields = static_cast<int>(rng() % 8);
                for (int f = 0; f < fields; ++f) {
                    if (f) body += "\t";
                    body += nasty(rng);
                }
                body += "\n";
            }
        }
        parse(body, ".bed", filter);  // invariants asserted inside
    }
}

void test_fuzz_bedpe() {
    std::mt19937 rng(20260818u);
    FragmentFilter filter;

    for (int trial = 0; trial < 1200; ++trial) {
        filter.max_length = static_cast<Coord>(rng() % 2 == 0 ? 0 : 1000);
        filter.min_mapq = static_cast<int>(rng() % 4 == 0 ? rng() % 260 : 0);

        std::string body;
        const int lines = static_cast<int>(rng() % 40);
        for (int i = 0; i < lines; ++i) {
            if (rng() % 3 == 0) {
                const long s1 = static_cast<long>(rng() % 1000000);
                const long s2 = static_cast<long>(rng() % 1000000);
                body += "chr1\t" + std::to_string(s1) + "\t" +
                        std::to_string(s1 + 100) + "\tchr1\t" + std::to_string(s2) +
                        "\t" + std::to_string(s2 + 100) + "\tn\t" + nasty(rng) + "\t" +
                        nasty(rng) + "\n";
            } else {
                const int fields = static_cast<int>(rng() % 11);
                for (int f = 0; f < fields; ++f) {
                    if (f) body += "\t";
                    body += nasty(rng);
                }
                body += "\n";
            }
        }
        parse(body, ".bedpe", filter);
    }
}

// --------------------------------------------------------------------------
// The 10x count column must never drive unbounded emission, and a PRESENT but
// corrupt count is a malformed record -- not one copy.
//
// This test used to pin the fail-open behaviour: "0", "-5" and "abc" each
// emitted one fragment, and an out-of-range count was clamped to a million.
// That is the T2 MAPQ defect in the sibling column (a corrupt field silently
// promoted to a plausible value), and it let one line add a million fragments
// to every statistic. Changed deliberately on 2026-09-11; the audit record is
// docs/AUDIT_2026-09-11_genomic_toolkit.md.
// --------------------------------------------------------------------------
void test_count_column_is_bounded() {
    FragmentFilter filter;
    filter.max_length = 0;
    filter.use_count_column = true;  // the rule below is --with-counts' rule

    struct Case { const char* count; std::uint64_t emitted; std::uint64_t malformed; };
    const Case cases[] = {
        {"1", 1, 0}, {"3", 3, 0},
        {"0", 0, 1},                     // a fragment line is at least one read pair
        {"-5", 0, 1},
        {"abc", 0, 1},
        {"", 0, 1},                      // present but empty
        {"3x", 0, 1},                    // partial number
        {"1000001", 0, 1},               // past the per-line ceiling: refused, not clamped
        {"99999999999999999999", 0, 1},  // past int64
    };
    for (const Case& c : cases) {
        const Parsed p = parse(std::string("c\t10\t200\tBC\t") + c.count + "\n",
                               ".bed", filter);
        CHECK(!p.threw);
        CHECK_EQ(p.stats.fragments_emitted, c.emitted);
        CHECK_EQ(p.stats.malformed_lines, c.malformed);
    }

    // Without --with-counts the column is not read at all: every case above,
    // corrupt or not, is one fragment, and nothing is malformed.
    FragmentFilter plain;
    plain.max_length = 0;
    for (const Case& c : cases) {
        const Parsed p = parse(std::string("c\t10\t200\tBC\t") + c.count + "\n",
                               ".bed", plain);
        CHECK(!p.threw);
        CHECK_EQ(p.stats.fragments_emitted, 1u);
        CHECK_EQ(p.stats.malformed_lines, 0u);
    }
}

}  // namespace

int main() {
    testing::Suite suite{"genomic_toolkit parser fuzzing", {}};
    suite.add("BEDPE MAPQ fails closed", test_bedpe_mapq_fails_closed);
    suite.add("oversized line fails closed", test_oversized_line_fails_closed);
    suite.add("coordinate boundaries", test_coordinate_boundaries);
    suite.add("degenerate files", test_degenerate_files);
    suite.add("count column is bounded", test_count_column_is_bounded);
    suite.add("fuzz: fragment BED", test_fuzz_fragment_bed);
    suite.add("fuzz: BEDPE", test_fuzz_bedpe);
    return suite.run();
}
