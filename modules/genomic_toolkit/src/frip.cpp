// SPDX-License-Identifier: MIT

#include "toolkit/frip.hpp"

#include <charconv>
#include <stdexcept>
#include <string_view>

#include "line_reader.hpp"

namespace toolkit {
namespace {

bool parse_int(std::string_view s, std::int64_t& out) {
    if (s.empty()) return false;
    const auto res = std::from_chars(s.data(), s.data() + s.size(), out);
    return res.ec == std::errc() && res.ptr == s.data() + s.size();
}

}  // namespace

std::vector<Region> read_peaks_bed(const std::string& path) {
    detail::LineReader reader(path);
    std::vector<Region> peaks;
    std::string_view line;
    std::size_t line_no = 0;

    while (reader.next(line)) {
        ++line_no;
        if (line.empty() || line.front() == '#' || line.starts_with("track") ||
            line.starts_with("browser")) {
            continue;
        }

        std::string_view f[4];
        std::size_t n = 0;
        std::size_t start = 0;
        while (n < 4) {
            const std::size_t tab = line.find('\t', start);
            if (tab == std::string_view::npos) {
                f[n++] = line.substr(start);
                break;
            }
            f[n++] = line.substr(start, tab - start);
            start = tab + 1;
        }

        Region r;
        // A malformed peak file is fatal, unlike a malformed fragment line
        // which is counted and skipped. The asymmetry is deliberate: a peak
        // file is small, hand-curated and the denominator of the statistic, so
        // silently dropping part of it changes the answer with no signal. A
        // fragment file has hundreds of millions of lines and a tolerated
        // handful is normal.
        if (n < 3 || !parse_int(f[1], r.start) || !parse_int(f[2], r.end)) {
            throw std::invalid_argument("malformed BED record at " + path + ":" +
                                        std::to_string(line_no));
        }
        r.chrom = std::string(f[0]);
        if (n >= 4) r.name = std::string(f[3]);
        peaks.push_back(std::move(r));
    }

    if (peaks.empty()) {
        throw std::invalid_argument("no usable peak records in '" + path + "'");
    }
    return peaks;
}

}  // namespace toolkit
