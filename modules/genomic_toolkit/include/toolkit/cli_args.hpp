// SPDX-License-Identifier: MIT
//
// Checked integer parsing for command-line options.
//
// The CLI used std::stol/std::stoi/std::stoll and then narrowed the result, and
// every property of that was wrong for a tool whose answers go into papers:
//
//   --min-length abc          std::invalid_argument, thrown from argument
//                             parsing, which runs outside main's try -- SIGABRT,
//                             exit 134, no message a user can act on.
//   --min-length 100x         read as 100: std::stol stops at the first
//                             non-digit and reports success.
//   --min-length 4294967396   parsed as a long, then static_cast to int32 --
//                             100 again, a length gate nobody asked for, exit 0.
//   --max-length 4294967296   narrowed to 0, which means "off".
//   --genome-size 0x10        read as 0: enrichment silently disabled.
//   --min-mapq 300, -1        accepted; MAPQ is 0-255.
//
// Same rule as the text parsers in fragment_stream.cpp: the whole field is a
// base-10 integer or it is refused, and a value outside the option's range is
// refused by name rather than wrapped into it.
#pragma once

#include <charconv>
#include <cstdint>
#include <string_view>
#include <system_error>

namespace toolkit {

// True, with `out` set, when `text` is entirely a base-10 integer in [lo, hi].
// No whitespace, no leading '+', no trailing characters, no "0x" prefix. On
// failure `out` is left untouched.
[[nodiscard]] inline bool parse_bounded_int(std::string_view text, std::int64_t lo,
                                            std::int64_t hi, std::int64_t& out) noexcept {
    if (text.empty()) return false;
    std::int64_t value = 0;
    const char* const last = text.data() + text.size();
    const auto res = std::from_chars(text.data(), last, value);
    if (res.ec != std::errc() || res.ptr != last) return false;
    if (value < lo || value > hi) return false;
    out = value;
    return true;
}

}  // namespace toolkit
