// SPDX-License-Identifier: MIT
//
// parse_bounded_int -- the one parse every numeric CLI option goes through.
//
// Each refusal below is an input the previous std::stol / std::stoi path got
// wrong, by aborting, by reading a prefix, or by narrowing into int32. The
// history is in include/toolkit/cli_args.hpp; this suite is what stops it
// coming back.

#include "toolkit/cli_args.hpp"

#include <cstdint>
#include <limits>

#include "test_util.hpp"
#include "toolkit/types.hpp"

using namespace toolkit;

namespace {

constexpr std::int64_t kI64Max = std::numeric_limits<std::int64_t>::max();
constexpr std::int64_t kSentinel = 42;

// Asserts a refusal AND that the output was left alone: a parser that
// reports failure after writing a partial value is the from_chars trap the
// text readers already document.
void check_refused(const char* text, std::int64_t lo, std::int64_t hi) {
    std::int64_t v = kSentinel;
    const bool ok = parse_bounded_int(text, lo, hi, v);
    ::testing::report(!ok, text, __FILE__, __LINE__, "accepted, but must be refused");
    CHECK_EQ(v, kSentinel);
}

}  // namespace

int main() {
    testing::Suite suite{"cli_args", {}};

    suite.add("accepts whole-field integers, bounds inclusive", [] {
        std::int64_t v = -1;
        CHECK(parse_bounded_int("0", 0, 255, v));
        CHECK_EQ(v, std::int64_t{0});
        CHECK(parse_bounded_int("255", 0, 255, v));
        CHECK_EQ(v, std::int64_t{255});
        CHECK(parse_bounded_int("2147483647", 0, kMaxCoord, v));
        CHECK_EQ(v, kMaxCoord);
        CHECK(parse_bounded_int("-1", -5, 5, v));
        CHECK_EQ(v, std::int64_t{-1});
    });

    suite.add("a partial number is refused, not read as its prefix", [] {
        // std::stol read "100x" as 100 and reported success.
        check_refused("100x", 0, kMaxCoord);
        check_refused("100 ", 0, kMaxCoord);
        check_refused(" 100", 0, kMaxCoord);
        check_refused("1e3", 0, kMaxCoord);
        check_refused("1.5", 0, kMaxCoord);
        // --genome-size 0x10 was read as 0, silently disabling enrichment.
        check_refused("0x10", 0, kI64Max);
        check_refused("+5", 0, kMaxCoord);
    });

    suite.add("a non-number is refused, not thrown", [] {
        // std::stol threw std::invalid_argument from argument parsing, outside
        // main's try: SIGABRT, exit 134.
        check_refused("abc", 0, kMaxCoord);
        check_refused("", 0, kMaxCoord);
        check_refused("-", -5, 5);
    });

    suite.add("a value past int32 is refused, not narrowed into it", [] {
        // static_cast<Coord>(std::stol("4294967396")) is 100: a --min-length
        // nobody asked for, exit 0. 4294967296 narrowed to 0, which means "off"
        // for --max-length.
        check_refused("4294967396", 0, kMaxCoord);
        check_refused("4294967296", 0, kMaxCoord);
        check_refused("2147483648", 0, kMaxCoord);
    });

    suite.add("a value past int64 is refused, not thrown", [] {
        // std::stol threw std::out_of_range here, same abort as "abc".
        std::int64_t v = kSentinel;
        CHECK(parse_bounded_int("9223372036854775807", 0, kI64Max, v));
        CHECK_EQ(v, kI64Max);
        check_refused("9223372036854775808", 0, kI64Max);
        check_refused("99999999999999999999", 0, kI64Max);
        check_refused("-9223372036854775809", std::numeric_limits<std::int64_t>::min(), 0);
    });

    suite.add("a value outside the option's own range is refused", [] {
        // --min-mapq 300 and -1 were accepted; MAPQ is 0-255.
        check_refused("300", 0, 255);
        check_refused("256", 0, 255);
        check_refused("-1", 0, 255);
        check_refused("-50", 0, kMaxCoord);
    });

    return suite.run();
}
