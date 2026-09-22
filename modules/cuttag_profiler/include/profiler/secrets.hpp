// SPDX-License-Identifier: MIT
//
// Access-control secrets for the embedded HTTP server: how they are minted,
// and how they are compared.
//
// Two functions, both closing SECURITY_HTTP_2026-08-15 finding L9, and both
// here rather than in http_server.cpp's anonymous namespace for the same
// reason http_security.hpp exists: a security primitive with no test is a
// primitive that rots. See tests/test_secrets.cpp.
//
// They are also deliberately OUTSIDE the PTO_CLOUD_BUILD guard that strips the
// HTTP server. The compare is useful to any future caller that handles a
// secret, and keeping the file compiled in both configurations means the tests
// that pin its behaviour run in both configurations too.
#pragma once

#include <cstddef>
#include <string>
#include <string_view>

namespace profiler {

// `bytes` of entropy from the platform CSPRNG, rendered as fixed-width
// lowercase hex (so the returned string is `bytes * 2` characters).
//
// Replaces a `std::mt19937_64` seeded from one `std::random_device` call, which
// was wrong twice over. mt19937_64 is not a CSPRNG -- its full 19937-bit state
// is recoverable from 312 consecutive outputs -- and 64 bits of seed for that
// state is the real bound on the entropy anyway. It was being used as the sole
// access control on /api/matrix, which streams an entire computed matrix.
//
// The fixed width matters independently: the old implementation formatted with
// `std::hex`, which suppresses leading zeros, so roughly one token in sixteen
// was short by a nibble or more and the tokens were of variable length. That
// makes the delivered entropy smaller than the code appears to claim, and it
// makes the length itself carry information about the value.
//
// Throws std::runtime_error when the system cannot supply entropy. Failing
// closed is the only correct behaviour: the alternative is a server that comes
// up and serves data behind a predictable token.
[[nodiscard]] std::string secure_token(std::size_t bytes = 32);

// Compares two secrets in time that does not depend on where they first
// differ.
//
// `operator==` on std::string short-circuits at the first differing byte, so
// the time it takes to reject a guess reveals how long a correct prefix was --
// which lets an attacker recover a token one byte at a time instead of
// guessing all of it at once. Over a loopback HTTP round trip the signal is
// buried in noise and this is a low-severity finding, but the fix costs one
// loop and the endpoint it guards is reachable from a browser (H1/H3).
//
// What this does NOT hide is the length: unequal lengths return immediately.
// That is deliberate and safe here, because every secret compared through this
// function is minted by `secure_token()` at a fixed, public width.
//
// `diff` is volatile so the accumulation survives optimisation. Without it a
// compiler is entitled to notice that the loop's only observable effect is the
// final comparison and rewrite it into an early exit -- turning a
// constant-time function back into a variable-time one, silently, at whatever
// optimisation level the release build happens to use.
[[nodiscard]] inline bool constant_time_equals(std::string_view a,
                                               std::string_view b) noexcept {
    if (a.size() != b.size()) return false;

    volatile unsigned char diff = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        diff = static_cast<unsigned char>(
            diff | (static_cast<unsigned char>(a[i]) ^
                    static_cast<unsigned char>(b[i])));
    }
    return diff == 0;
}

}  // namespace profiler
