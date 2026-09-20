// SPDX-License-Identifier: MIT
//
// IntervalIndex. The differential test at the bottom is the one that matters:
// a descent with a pruning rule is easy to get subtly wrong in a way that
// hand-written cases miss, so 20k random queries are checked against an O(n)
// scan that is obviously correct by inspection.

#include <algorithm>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "test_util.hpp"

#include "profiler/interval_index.hpp"

namespace {

using profiler::Coord;
using profiler::IntervalIndex;
using profiler::Region;

// The shared harness has no throw-checking macro; interval validation is the
// first thing in this module that needs one at the call site rather than
// inline, so it is defined locally rather than widening test_util.hpp for a
// single consumer.
#define CHECK_THROWS(expr)                                                     \
    do {                                                                       \
        bool threw_ = false;                                                   \
        try {                                                                  \
            (void)(expr);                                                      \
        } catch (const std::invalid_argument&) {                               \
            threw_ = true;                                                     \
        } catch (...) {                                                        \
        }                                                                      \
        ::testing::report(threw_, "throws: " #expr, __FILE__, __LINE__);       \
    } while (0)

Region reg(std::string chrom, std::int64_t start, std::int64_t end,
           std::string name = {}) {
    Region r;
    r.chrom = std::move(chrom);
    r.start = start;
    r.end = end;
    r.name = std::move(name);
    return r;
}

// Payload ids of the hits, sorted, so tests compare sets without depending on
// the traversal order.
std::vector<std::uint32_t> hits(const IntervalIndex& idx, const char* chrom, Coord qs,
                               Coord qe) {
    std::vector<std::uint32_t> scratch;
    const auto span = idx.overlap(chrom, qs, qe, scratch);
    std::vector<std::uint32_t> out(span.begin(), span.end());
    std::sort(out.begin(), out.end());
    return out;
}

// CHECK_EQ formats its operands via testing::describe, which has no overload
// for a vector; compare payload sets through CHECK instead, printing the two
// lists ourselves when they differ.
bool same(const std::vector<std::uint32_t>& got,
          const std::vector<std::uint32_t>& want, const char* where) {
    if (got == want) return true;
    std::fprintf(stderr, "  %s: got {", where);
    for (auto v : got) std::fprintf(stderr, "%u,", v);
    std::fprintf(stderr, "} want {");
    for (auto v : want) std::fprintf(stderr, "%u,", v);
    std::fprintf(stderr, "}\n");
    return false;
}

// ---------------------------------------------------------------- basics

void test_empty_index() {
    const auto idx = IntervalIndex::build({});
    CHECK_EQ(idx.size(), std::size_t{0});
    CHECK_EQ(idx.contig_count(), std::size_t{0});
    CHECK(hits(idx, "chr1", 0, 1000).empty());
}

void test_single_interval() {
    const std::vector<Region> regions{reg("chr1", 100, 200)};
    const auto idx = IntervalIndex::build(regions);

    CHECK_EQ(hits(idx, "chr1", 150, 160).size(), std::size_t{1});
    CHECK_EQ(hits(idx, "chr1", 0, 100).size(), std::size_t{0});   // entirely left
    CHECK_EQ(hits(idx, "chr1", 200, 300).size(), std::size_t{0}); // entirely right
    CHECK_EQ(hits(idx, "chr1", 99, 101).size(), std::size_t{1});  // clips the start
    CHECK_EQ(hits(idx, "chr1", 199, 201).size(), std::size_t{1}); // clips the end
    CHECK_EQ(hits(idx, "chr2", 150, 160).size(), std::size_t{0}); // unknown contig
}

// Half-open means touching endpoints do NOT overlap. [0,10) and [10,20) are
// adjacent, not overlapping, and BED files are full of exactly this.
void test_adjacent_intervals_do_not_overlap() {
    const std::vector<Region> regions{reg("chr1", 0, 10), reg("chr1", 10, 20)};
    const auto idx = IntervalIndex::build(regions);

    const auto left = hits(idx, "chr1", 0, 10);
    CHECK_EQ(left.size(), std::size_t{1});
    CHECK_EQ(left[0], std::uint32_t{0});

    const auto right = hits(idx, "chr1", 10, 20);
    CHECK_EQ(right.size(), std::size_t{1});
    CHECK_EQ(right[0], std::uint32_t{1});

    CHECK_EQ(hits(idx, "chr1", 5, 15).size(), std::size_t{2});
}

// The case a "binary search the start, then scan forward" implementation gets
// wrong: a long interval that starts early and swallows several short ones.
// Querying inside the nested region must still find the container.
void test_nested_intervals() {
    const std::vector<Region> regions{
        reg("chr1", 100, 100000),  // 0: the container
        reg("chr1", 120, 130),     // 1
        reg("chr1", 500, 510),     // 2
        reg("chr1", 90000, 90010), // 3
    };
    const auto idx = IntervalIndex::build(regions);

    CHECK(same(hits(idx, "chr1", 125, 126), {0, 1}, "chr1:125-126"));
    CHECK(same(hits(idx, "chr1", 90000, 90001), {0, 3}, "chr1:90000-90001"));
    // Between the nested ones: only the container.
    CHECK(same(hits(idx, "chr1", 1000, 1001), {0}, "chr1:1000-1001"));
}

void test_duplicate_starts() {
    const std::vector<Region> regions{
        reg("chr1", 100, 200), reg("chr1", 100, 300), reg("chr1", 100, 150)};
    const auto idx = IntervalIndex::build(regions);

    CHECK_EQ(hits(idx, "chr1", 100, 101).size(), std::size_t{3});
    CHECK(same(hits(idx, "chr1", 160, 161), {0, 1}, "chr1:160-161"));
    CHECK(same(hits(idx, "chr1", 250, 251), {1}, "chr1:250-251"));
}

// Zero-length intervals are legal in BED (insertion points) and overlap
// nothing under the half-open predicate. Asserted so the behaviour is a
// decision on record rather than an accident.
void test_zero_length_interval_overlaps_nothing() {
    const std::vector<Region> regions{reg("chr1", 100, 100), reg("chr1", 50, 150)};
    const auto idx = IntervalIndex::build(regions);
    CHECK(same(hits(idx, "chr1", 100, 101), {1}, "chr1:100-101"));
    CHECK(same(hits(idx, "chr1", 99, 100), {1}, "chr1:99-100"));
    // Strictly containing the zero-length record still must not report it: it
    // covers no bases. The two-term half-open predicate gets this wrong.
    CHECK(same(hits(idx, "chr1", 99, 101), {1}, "chr1:99-101"));
    // ...and a zero-length *query* is itself empty, so it matches nothing,
    // including the interval that contains its coordinate.
    CHECK_EQ(hits(idx, "chr1", 100, 100).size(), std::size_t{0});
    CHECK_EQ(hits(idx, "chr1", 60, 60).size(), std::size_t{0});
}

void test_multiple_contigs_are_independent() {
    const std::vector<Region> regions{
        reg("chr1", 100, 200), reg("chr2", 100, 200), reg("chrX", 100, 200)};
    const auto idx = IntervalIndex::build(regions);
    CHECK_EQ(idx.contig_count(), std::size_t{3});
    CHECK(same(hits(idx, "chr1", 150, 151), {0}, "chr1:150-151"));
    CHECK(same(hits(idx, "chr2", 150, 151), {1}, "chr2:150-151"));
    CHECK(same(hits(idx, "chrX", 150, 151), {2}, "chrX:150-151"));
}

void test_inverted_query_returns_nothing() {
    const std::vector<Region> regions{reg("chr1", 100, 200)};
    const auto idx = IntervalIndex::build(regions);
    CHECK_EQ(hits(idx, "chr1", 200, 100).size(), std::size_t{0});
}

// ------------------------------------------------------------ validation

void test_bad_regions_are_rejected() {
    CHECK_THROWS(IntervalIndex::build({reg("chr1", -1, 100)}));
    CHECK_THROWS(IntervalIndex::build({reg("chr1", 200, 100)}));
    // Past the int32 ceiling. Truncating would place the annotation at a wrong
    // coordinate; the message has to name the contig so the user can act on it.
    CHECK_THROWS(IntervalIndex::build({reg("chrHuge", 0, 3'000'000'000LL)}));
}

// ------------------------------------------------------------- tid binding

void test_tid_binding() {
    const std::vector<Region> regions{reg("chr1", 100, 200), reg("chr2", 300, 400)};
    auto idx = IntervalIndex::build(regions);
    // Header order deliberately differs from the annotation's first-seen order,
    // and carries a contig with no annotations.
    idx.bind_contigs({"chr2", "chrM", "chr1"});
    CHECK_EQ(idx.bound_contig_count(), std::size_t{2});

    std::vector<std::uint32_t> scratch;
    CHECK_EQ(idx.overlap(0, 350, 351, scratch).size(), std::size_t{1});
    CHECK_EQ(idx.overlap(0, 350, 351, scratch)[0], std::uint32_t{1});
    CHECK_EQ(idx.overlap(2, 150, 151, scratch)[0], std::uint32_t{0});
    CHECK_EQ(idx.overlap(1, 0, 1'000'000, scratch).size(), std::size_t{0});  // chrM
    // Out-of-range tids answer empty rather than throwing: a BAM legitimately
    // has contigs the annotation does not cover.
    CHECK_EQ(idx.overlap(-1, 0, 100, scratch).size(), std::size_t{0});
    CHECK_EQ(idx.overlap(99, 0, 100, scratch).size(), std::size_t{0});
}

// A naming mismatch is the most common way this silently returns nothing, so
// it has to be reportable rather than merely survivable.
void test_unbound_contigs_are_reported() {
    const std::vector<Region> regions{reg("1", 100, 200), reg("chr2", 100, 200)};
    auto idx = IntervalIndex::build(regions);
    idx.bind_contigs({"chr1", "chr2"});
    const auto unbound = idx.unbound_contigs();
    CHECK_EQ(unbound.size(), std::size_t{1});
    CHECK_EQ(unbound[0], std::string("1"));
}

// -------------------------------------------------------- differential test

void test_matches_brute_force_on_random_input() {
    std::mt19937 rng(20260816);
    std::uniform_int_distribution<std::int64_t> start_dist(0, 1'000'000);
    // A deliberately skewed length distribution: mostly short features with a
    // long tail of gene-body-sized ones, so nesting and deep overlap stacks
    // occur rather than being a theoretical concern.
    std::uniform_int_distribution<int> shape(0, 99);
    std::uniform_int_distribution<std::int64_t> short_len(0, 500);
    std::uniform_int_distribution<std::int64_t> long_len(500, 200'000);

    std::vector<Region> regions;
    regions.reserve(5000);
    for (int i = 0; i < 5000; ++i) {
        const std::int64_t s = start_dist(rng);
        const std::int64_t len = shape(rng) < 90 ? short_len(rng) : long_len(rng);
        regions.push_back(reg("chr1", s, s + len));
    }
    const auto idx = IntervalIndex::build(regions);

    std::size_t total_hits = 0;
    bool mismatch = false;
    for (int q = 0; q < 20000 && !mismatch; ++q) {
        const std::int64_t qs = start_dist(rng);
        const std::int64_t qe = qs + (shape(rng) < 90 ? short_len(rng) : long_len(rng));

        std::vector<std::uint32_t> expected;
        for (std::size_t i = 0; i < regions.size(); ++i) {
            // Stated as "the intersection is non-empty", deliberately NOT as
            // the implementation's three-term predicate: the point of a
            // differential test is to compare two independent formulations.
            if (std::max<std::int64_t>(regions[i].start, qs) <
                std::min<std::int64_t>(regions[i].end, qe)) {
                expected.push_back(static_cast<std::uint32_t>(i));
            }
        }
        std::sort(expected.begin(), expected.end());

        const auto got = hits(idx, "chr1", static_cast<Coord>(qs), static_cast<Coord>(qe));
        total_hits += got.size();
        if (got != expected) {
            std::fprintf(stderr,
                         "  differential mismatch at query [%lld, %lld): %zu hits vs %zu\n",
                         static_cast<long long>(qs), static_cast<long long>(qe),
                         got.size(), expected.size());
            mismatch = true;
        }
    }
    CHECK(!mismatch);
    // Guard against the test passing vacuously because every query missed.
    CHECK(total_hits > 10000);
}

}  // namespace

int main() {
    testing::Suite suite{"interval_index", {}};
    suite.add("empty index", test_empty_index);
    suite.add("single interval", test_single_interval);
    suite.add("adjacent intervals do not overlap", test_adjacent_intervals_do_not_overlap);
    suite.add("nested intervals", test_nested_intervals);
    suite.add("duplicate starts", test_duplicate_starts);
    suite.add("zero-length interval overlaps nothing",
              test_zero_length_interval_overlaps_nothing);
    suite.add("contigs are independent", test_multiple_contigs_are_independent);
    suite.add("inverted query returns nothing", test_inverted_query_returns_nothing);
    suite.add("malformed regions are rejected", test_bad_regions_are_rejected);
    suite.add("tid binding", test_tid_binding);
    suite.add("unbound contigs are reported", test_unbound_contigs_are_reported);
    suite.add("matches brute force on random input",
              test_matches_brute_force_on_random_input);
    return suite.run();
}
