// SPDX-License-Identifier: MIT
//
// Byte-identical replacement for std::ostream's `fixed`/`precision(6)`
// formatting of a double, using std::to_chars. `[charconv.to.chars]` specifies
// `chars_format::fixed` with a precision as "as if by printf with `%.*f`",
// which is exactly what an ofstream in the classic locale -- the only locale
// this module ever imbues, i.e. the default -- already produces. Verified
// directly (tests/test_signal_calc.cpp) rather than assumed, because two
// cases where std::to_chars is correct but NOT what ostream does would
// otherwise reach production silently:
//
//   * a value past the buffer's capacity: to_chars reports
//     errc::value_too_large rather than writing a truncated or garbage
//     result, so a buffer sized too small for a legitimate value is a
//     silent-wrong-answer risk, not merely a wasted allocation. Sized here for
//     the longest possible finite double (DBL_MAX in fixed notation is 309
//     integer digits) with headroom, and the error path aborts rather than
//     ever emitting a truncated number into a matrix.
//   * NaN: to_chars spells a negative NaN's mantissa sign bit
//     ("-nan(ind)"); ostream never signs NaN at all, always "nan" regardless
//     of that bit. Checked first and returned literally so both formatters
//     agree on every NaN -- which matters here because `missing_value`
//     (types.hpp) is a user-supplied double and NaN is a legitimate choice
//     for it.
#pragma once

#include <charconv>
#include <cmath>
#include <cstdlib>
#include <string>
#include <system_error>

#include "profiler/types.hpp"

namespace profiler {

// Appends `v` to `out` exactly as `out_stream.setf(std::ios::fixed);
// out_stream.precision(6); out_stream << v;` would have appended to a
// std::ostream. Appends rather than returning a std::string so a caller
// building one row of many columns reuses a single buffer instead of
// allocating a small string per value.
inline void append_fixed6(std::string& out, double v) {
    if (std::isnan(v)) {
        out += "nan";
        return;
    }
    // 309 integer digits (DBL_MAX) + 1 sign + 1 point + 6 fractional digits +
    // slack. Comfortably unreachable in practice for genomic signal values,
    // which is exactly why the failure path must not be a silent truncation.
    char buf[400];
    const auto res = std::to_chars(buf, buf + sizeof buf, v, std::chars_format::fixed, 6);
    if (res.ec != std::errc()) {
        std::abort();  // buffer sizing bug, not a reachable runtime input
    }
    out.append(buf, res.ptr);
}

// Formats `matrix` as one complete string, in the exact byte layout
// write_matrix_tsv (signal_calc.cpp) writes to a file -- the two are
// independent implementations (this one builds the whole matrix as a single
// std::string for an HTTP response body, that one streams row buffers to disk
// in ~1 MiB chunks and must never hold a whole large matrix in memory at
// once), but the format is one thing, verified equal directly
// (tests/test_signal_calc.cpp) rather than by each separately restating
// ostream's behaviour.
//
// Used by GET /api/matrix (http_server.cpp), which previously built the same
// content into a std::ostringstream, copied it out via .str(), and copied it
// again into httplib's response body -- three live copies of a matrix that is
// 73 MB at the module's own README benchmark geometry. Building directly into
// the string that becomes the response body (moved into set_content, not
// copied) removes two of those three.
[[nodiscard]] inline std::string format_matrix_tsv(const SignalMatrix& matrix) {
    std::string out;
    // ~10 bytes per fixed-precision-6 field is a deliberate estimate, not a
    // bound: append_fixed6 grows `out` past this for values it does not fit,
    // exactly as std::string::append always does.
    out.reserve(matrix.rows * (matrix.cols * 10 + 24));
    out += "region";
    for (std::int64_t offset : matrix.bin_offsets) {
        out += '\t';
        out += std::to_string(offset);
    }
    out += '\n';
    for (std::size_t r = 0; r < matrix.rows; ++r) {
        out += matrix.row_names[r];
        const double* row = matrix.values.data() + r * matrix.cols;
        for (std::size_t c = 0; c < matrix.cols; ++c) {
            out += '\t';
            append_fixed6(out, row[c]);
        }
        out += '\n';
    }
    return out;
}

}  // namespace profiler
