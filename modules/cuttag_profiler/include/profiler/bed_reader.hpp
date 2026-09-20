// SPDX-License-Identifier: MIT
//
// Minimal readers for the two region formats a bench scientist is likely to
// have on hand: BED (3, 4, 5 or 6+ columns) and GTF/GFF restricted to a single
// feature type. Both are streamed, so a 100k-line promoter file costs one pass
// and one allocation per region.
#pragma once

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

#include "profiler/safe_open.hpp"
#include "profiler/types.hpp"

namespace profiler {

// The longest region name kept. Names identify a locus; anything past this is
// not an identifier, it is payload.
inline constexpr std::size_t kMaxRegionNameLength = 256;

// Bounds and de-fangs a name lifted out of a region file.
//
// SECURITY_HTTP_2026-08-15 finding M5: BED column 4 is copied verbatim into
// `SignalMatrix::row_names`, and from there into `topRegions[].name` in the
// /api/profile JSON, into the first column of the /api/matrix TSV, and into
// the TSV that `profile --out-matrix` writes. It is untrusted text on a
// straight line to three outputs, so it gets bounded once, here, at the point
// it is parsed -- not at each of the three places it is emitted, which is the
// arrangement that eventually misses one.
//
// Two properties:
//
//   * LENGTH. Capped, so a 4 MB column-4 field cannot be pulled through the
//     API one region at a time. Truncation is deliberate rather than a
//     rejection: a long name is a malformed annotation, not an attack to
//     abort a whole run over.
//
//   * CONTROL CHARACTERS, replaced with '?'. This is the half that is a bug
//     rather than a mitigation. A name containing a tab or a newline breaks
//     the framing of every TSV this tool writes, letting a crafted BED forge
//     extra rows in a file a downstream tool will parse as data. DEL and the
//     C0 range also carry terminal escape sequences, and these names are
//     printed to an operator's terminal by the batch path.
//
// What it does NOT do is HTML-escape: the dashboard must not be relying on
// its input being pre-escaped, and doing it here would corrupt legitimate
// names containing '&' or '<'.
[[nodiscard]] inline std::string sanitize_region_name(std::string_view raw) {
    const std::string_view bounded = raw.substr(0, kMaxRegionNameLength);
    std::string out;
    out.reserve(bounded.size());
    for (char c : bounded) {
        const auto byte = static_cast<unsigned char>(c);
        out.push_back(byte < 0x20 || byte == 0x7F ? '?' : c);
    }
    return out;
}

// Parses a BED file. Track/browser/comment lines are skipped. Throws
// std::runtime_error on an unreadable file or a malformed coordinate field.
std::vector<Region> read_bed(const std::string& path);

// Parses a GTF/GFF, keeping rows whose feature column equals `feature`
// (default "gene") and naming each region from the given attribute key.
std::vector<Region> read_gtf(const std::string& path,
                             const std::string& feature = "gene",
                             const std::string& name_attribute = "gene_name");

// Dispatches on file extension (.bed / .gtf / .gff / .gff3), transparently
// handling a .gz suffix on either.
//
// NOTE for callers with untrusted input: the format is chosen from the
// FILENAME, and the two parsers are not equally strict. read_bed throws on the
// first line it cannot parse; read_gtf skips such lines and keeps going. So a
// caller who lets a client name the path also lets that client choose how
// permissive the parser is, and the lenient one will harvest usable rows out
// of a file that is not an annotation at all. Confinement and authentication
// are what bound this (H1/H2); prefer calling read_bed or read_gtf directly
// when the format is known.
std::vector<Region> read_regions(const std::string& path);

// As read_regions, but reads an already-open descriptor.
//
// SECURITY_HTTP finding M6: the caller checked a path and then this function
// looked the same name up again to open it, so anything able to write in the
// data root could pass the check with a real file and be read through a
// symlink. Taking the descriptor the caller already verified means there is no
// second lookup to race -- see profiler/safe_open.hpp.
//
// `name_for_format` selects the parser exactly as read_regions does and
// appears in error messages; it is never opened. Ownership of `fd` is taken
// and it is closed before this returns, on the throwing path as well.
std::vector<Region> read_regions_from_fd(UniqueFd fd,
                                         const std::string& name_for_format);

}  // namespace profiler
