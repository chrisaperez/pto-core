// SPDX-License-Identifier: MIT
//
// Tests for the access-control secrets behind SECURITY_HTTP_2026-08-15 L9:
// how a token is minted, and how it is compared.
//
// Both properties are the kind that fail silently. A token that is 8 bits
// short still looks like a token; a comparison that short-circuits still
// returns the right answer. Neither shows up in a functional test of the
// server, so they are pinned here.
//
// Registered in the cloud configuration too. The primitives are compiled in
// both, and a security test that quietly stops running in the build that ships
// is worse than not having written it.
#include <cstdio>

#include <algorithm>
#include <set>
#include <string>

#include "profiler/secrets.hpp"
#include "test_util.hpp"

using profiler::constant_time_equals;
using profiler::secure_token;

int main() {
    testing::Suite suite{"secrets", {}};

    // ------------------------------------------------------------ the mint
    suite.add("a token is fixed-width lowercase hex", [] {
        // The defect being pinned: the old implementation formatted with
        // `std::hex`, which suppresses leading zeros, so roughly one token in
        // sixteen came out short. Length was therefore variable and the
        // delivered entropy was below the nominal figure. 400 samples makes a
        // missing nibble a near-certain catch rather than a coin flip.
        for (int i = 0; i < 400; ++i) {
            const std::string token = secure_token(32);
            CHECK_EQ(token.size(), std::size_t{64});
            CHECK(std::all_of(token.begin(), token.end(), [](char c) {
                return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
            }));
        }
    });

    suite.add("the requested width is honoured exactly", [] {
        CHECK_EQ(secure_token(1).size(), std::size_t{2});
        CHECK_EQ(secure_token(16).size(), std::size_t{32});
        CHECK_EQ(secure_token(64).size(), std::size_t{128});
    });

    suite.add("a zero-length token is refused rather than returned empty", [] {
        // An empty token would compare equal to an absent one, which is the
        // whole authentication check answering yes.
        bool threw = false;
        try {
            (void)secure_token(0);
        } catch (const std::exception&) {
            threw = true;
        }
        CHECK(threw);
    });

    suite.add("tokens do not repeat", [] {
        // Not a randomness test -- a real one does not belong in a unit suite.
        // This catches the failure that actually happens: a generator that is
        // reseeded identically per call, or per thread, and hands out the same
        // value twice. With 256 bits, any collision here is a bug.
        std::set<std::string> seen;
        for (int i = 0; i < 1000; ++i) seen.insert(secure_token(32));
        CHECK_EQ(seen.size(), std::size_t{1000});
    });

    suite.add("no byte position is stuck", [] {
        // A generator that fills only part of the buffer, or that renders the
        // high nibble twice, leaves a position constant across every sample.
        // The old `out << std::hex << rng() << rng()` could not fail this way,
        // but a future rewrite of fill_random could.
        constexpr int kSamples = 200;
        std::string first = secure_token(32);
        std::vector<bool> varied(first.size(), false);
        for (int i = 1; i < kSamples; ++i) {
            const std::string token = secure_token(32);
            for (std::size_t c = 0; c < first.size(); ++c) {
                if (token[c] != first[c]) varied[c] = true;
            }
        }
        CHECK(std::all_of(varied.begin(), varied.end(), [](bool v) { return v; }));
    });

    // --------------------------------------------------------- the compare
    suite.add("equal secrets compare equal", [] {
        const std::string token = secure_token(32);
        CHECK(constant_time_equals(token, token));
        CHECK(constant_time_equals(token, std::string(token)));  // distinct storage
        CHECK(constant_time_equals("", ""));
    });

    suite.add("any single differing byte is caught, at any position", [] {
        // The bug a constant-time rewrite invites is an accumulator that drops
        // a difference -- an `=` where `|=` was meant, or a loop bound that
        // stops one short. Flip each position in turn and require a rejection
        // every time, which the first-difference version would also pass but a
        // broken accumulator would not.
        const std::string token = secure_token(32);
        for (std::size_t i = 0; i < token.size(); ++i) {
            std::string tampered = token;
            tampered[i] = (tampered[i] == 'a') ? 'b' : 'a';
            CHECK(!constant_time_equals(token, tampered));
        }
    });

    suite.add("a correct prefix is not a match", [] {
        // The attack the constant-time compare exists to deny is recovering a
        // token byte by byte. Whatever the timing, the ANSWER must never be
        // yes for a prefix.
        const std::string token = secure_token(32);
        for (std::size_t n = 0; n < token.size(); ++n) {
            CHECK(!constant_time_equals(token, token.substr(0, n)));
        }
        CHECK(!constant_time_equals(token, token + "0"));
    });

    suite.add("an absent token never authenticates", [] {
        // The pre-routing handler passes an empty string when the request
        // carried neither ?t= nor a Bearer header. That must not match a real
        // token, and must not match a server that failed to mint one.
        CHECK(!constant_time_equals("", secure_token(32)));
        CHECK(!constant_time_equals(secure_token(32), ""));
    });

    suite.add("comparison is byte-exact, not case- or space-insensitive", [] {
        CHECK(!constant_time_equals("abcdef", "ABCDEF"));
        CHECK(!constant_time_equals("abcdef", "abcdef "));
        CHECK(!constant_time_equals("abcdef", " abcdef"));
        // Embedded NUL: string_view carries a length, so the compare must not
        // stop at the terminator the way strcmp would.
        CHECK(constant_time_equals(std::string("ab\0cd", 5), std::string("ab\0cd", 5)));
        CHECK(!constant_time_equals(std::string("ab\0cd", 5), std::string("ab\0ce", 5)));
    });

    return suite.run();
}
