// SPDX-License-Identifier: MIT
#include "profiler/bed_reader.hpp"

#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include <algorithm>
#include <charconv>
#include <cstring>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#if !defined(_WIN32)
#include <unistd.h>  // close, for the bgzf_hdopen failure path
#else
#include <io.h>
#define close _close
#endif

#include "profiler/safe_open.hpp"

namespace profiler {
namespace {

// BGZF transparently reads plain text, gzip and bgzip files, so one code path
// covers `.bed`, `.bed.gz` and `.gtf.gz` alike.
class LineReader {
public:
    explicit LineReader(const std::string& path) : path_(path) {
        fp_ = bgzf_open(path.c_str(), "r");
        if (fp_ == nullptr) {
            throw std::runtime_error("cannot open region file: " + path);
        }
        buffer_.resize(kChunkBytes);
    }

    // Reads a descriptor the caller has already opened and verified, so the
    // name is never resolved a second time (M6). bgzf_dopen wraps the
    // descriptor in an hFILE that takes ownership of it, and bgzf_close
    // releases both -- which is why the UniqueFd is released here rather than
    // left to close it as well. A double close is not a harmless no-op in a
    // threaded process: the number is recycled immediately, and the second
    // close lands on whatever another thread has just opened.
    LineReader(UniqueFd fd, std::string display_name)
        : path_(std::move(display_name)) {
        const int raw = fd.release();
        fp_ = bgzf_dopen(raw, "r");
        if (fp_ == nullptr) {
            ::close(raw);  // ownership did not transfer; still ours to release
            throw std::runtime_error("cannot open region file: " + path_);
        }
        buffer_.resize(kChunkBytes);
    }
    ~LineReader() {
        if (fp_ != nullptr) bgzf_close(fp_);
    }
    LineReader(const LineReader&) = delete;
    LineReader& operator=(const LineReader&) = delete;

    // Returns false at end of input. The view is invalidated by the next call.
    // A trailing \r is stripped, so CRLF files parse identically -- a
    // Windows-authored BED otherwise turns every last column into a value with
    // an invisible character glued to it.
    bool next(std::string_view& line) {
        for (;;) {
            for (std::size_t i = pos_; i < filled_; ++i) {
                if (buffer_[i] == '\n') {
                    std::size_t end = i;
                    if (end > pos_ && buffer_[end - 1] == '\r') --end;
                    line = std::string_view(buffer_.data() + pos_, end - pos_);
                    pos_ = i + 1;
                    return true;
                }
            }
            if (!refill()) {
                if (pos_ >= filled_) return false;
                // Final line with no terminator.
                std::size_t end = filled_;
                if (end > pos_ && buffer_[end - 1] == '\r') --end;
                line = std::string_view(buffer_.data() + pos_, end - pos_);
                pos_ = filled_;
                return true;
            }
        }
    }

private:
    // Moves the unconsumed tail to the front and reads another chunk. Returns
    // false when the source is exhausted. The buffer grows when a single line
    // does not fit, so a long line is never silently split into two malformed
    // records -- but growth stops at kMaxLineBytes and throws rather than
    // following the input up to an OOM kill.
    bool refill() {
        if (pos_ > 0) {
            const std::size_t tail = filled_ - pos_;
            if (tail > 0) std::memmove(buffer_.data(), buffer_.data() + pos_, tail);
            filled_ = tail;
            pos_ = 0;
        }
        if (filled_ == buffer_.size()) {
            if (buffer_.size() >= kMaxLineBytes) {
                throw std::runtime_error(
                    "'" + path_ + "' contains a line longer than the " +
                    std::to_string(kMaxLineBytes / (1024 * 1024)) +
                    " MiB limit; this is not a text fragment file (check for a "
                    "binary or truncated input)");
            }
            buffer_.resize(std::min(buffer_.size() * 2, kMaxLineBytes));
        }

        const std::size_t want = buffer_.size() - filled_;
        long got = static_cast<long>(bgzf_read(fp_, buffer_.data() + filled_, want));
        if (got < 0) throw std::runtime_error("read error on '" + path_ + "'");
        filled_ += static_cast<std::size_t>(got);
        return got > 0;
    }

    static constexpr std::size_t kChunkBytes = 1 << 20;   // 1 MiB
public:
    static constexpr std::size_t kMaxLineBytes = 64ull << 20;  // 64 MiB
private:
    std::string path_;
    BGZF* fp_ = nullptr;
    std::vector<char> buffer_;
    std::size_t filled_ = 0;
    std::size_t pos_ = 0;
};

std::vector<std::string_view> split(std::string_view line, char delim) {
    std::vector<std::string_view> fields;
    std::size_t start = 0;
    while (start <= line.size()) {
        const std::size_t pos = line.find(delim, start);
        if (pos == std::string_view::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    return fields;
}

bool parse_int(std::string_view text, std::int64_t& out) {
    const auto* first = text.data();
    const auto* last = text.data() + text.size();
    const auto result = std::from_chars(first, last, out);
    return result.ec == std::errc{} && result.ptr == last;
}

Strand parse_strand(std::string_view text) {
    if (text == "+") return Strand::kForward;
    if (text == "-") return Strand::kReverse;
    return Strand::kUnknown;
}

bool is_comment(std::string_view line) {
    return line.empty() || line.front() == '#' ||
           line.starts_with("track") || line.starts_with("browser");
}

bool ends_with_ci(const std::string& text, std::string_view suffix) {
    if (text.size() < suffix.size()) return false;
    return std::equal(suffix.rbegin(), suffix.rend(), text.rbegin(),
                      [](char a, char b) {
                          return std::tolower(static_cast<unsigned char>(a)) ==
                                 std::tolower(static_cast<unsigned char>(b));
                      });
}

// The parsers below take an already-open LineReader so that one body serves
// both the path-based entry points and the descriptor-based one (M6). `path`
// is used for error messages only; it is never opened here.
std::vector<Region> parse_bed(LineReader& reader, const std::string& path) {
    std::vector<Region> regions;
    regions.reserve(4096);

    std::string_view line;
    std::size_t line_no = 0;
    while (reader.next(line)) {
        ++line_no;
        if (is_comment(line)) continue;

        const auto fields = split(line, '\t');
        if (fields.size() < 3) {
            throw std::runtime_error(path + ":" + std::to_string(line_no) +
                                     ": BED requires at least 3 columns");
        }

        Region region;
        region.chrom = std::string(fields[0]);
        if (!parse_int(fields[1], region.start) ||
            !parse_int(fields[2], region.end)) {
            throw std::runtime_error(path + ":" + std::to_string(line_no) +
                                     ": non-numeric coordinate");
        }
        // Negative coordinates are rejected rather than clamped. BED is
        // 0-based and half-open, so a negative start is not a representable
        // locus; it used to be carried through to the window arithmetic in
        // signal_calc, where it is a silently wrong answer rather than a
        // refusal. from_chars happily parses "-1".
        if (region.start < 0 || region.end < 0) {
            throw std::runtime_error(path + ":" + std::to_string(line_no) +
                                     ": negative coordinate");
        }
        if (region.end < region.start) {
            throw std::runtime_error(path + ":" + std::to_string(line_no) +
                                     ": end precedes start");
        }
        if (fields.size() >= 4 && fields[3] != ".") {
            // M5: column 4 is echoed to the client and written into every TSV
            // this tool produces. Bound it here, once, where it enters.
            region.name = sanitize_region_name(fields[3]);
        }
        if (fields.size() >= 6) region.strand = parse_strand(fields[5]);
        regions.push_back(std::move(region));
    }
    return regions;
}

std::vector<Region> parse_gtf(LineReader& reader, const std::string& path,
                              const std::string& feature,
                              const std::string& name_attribute) {
    std::vector<Region> regions;
    regions.reserve(4096);

    std::string_view line;
    std::size_t line_no = 0;
    while (reader.next(line)) {
        ++line_no;
        if (line.empty() || line.front() == '#') continue;

        const auto fields = split(line, '\t');
        if (fields.size() < 8) continue;
        if (fields[2] != feature) continue;

        Region region;
        region.chrom = std::string(fields[0]);
        std::int64_t start = 0;
        std::int64_t end = 0;
        if (!parse_int(fields[3], start) || !parse_int(fields[4], end)) {
            throw std::runtime_error(path + ":" + std::to_string(line_no) +
                                     ": non-numeric coordinate");
        }
        // GTF is 1-based inclusive, so coordinate 0 is not a legal value and
        // `start - 1` would fabricate a negative one. Checked before the
        // conversion rather than after, so the subtraction never runs on a
        // value that cannot support it.
        if (start < 1 || end < 0 || end < start - 1) {
            throw std::runtime_error(path + ":" + std::to_string(line_no) +
                                     ": coordinate out of range");
        }
        region.start = start - 1;
        region.end = end;
        region.strand = parse_strand(fields[6]);

        if (fields.size() >= 9) {
            // Handles both GTF (`key "value";`) and GFF3 (`key=value;`).
            const std::string_view attrs = fields[8];
            // The key must be a whole attribute name: at the start of the
            // field or after a separator, and followed by the value delimiter.
            // A bare attrs.find() matched inside a longer key, so asking for
            // gene_name on `havana_gene_name "X"; gene_name "Y";` named the
            // region X (reproduced 2026-09-11).
            std::size_t key_pos = std::string_view::npos;
            for (std::size_t from = 0; from < attrs.size();) {
                const std::size_t at = attrs.find(name_attribute, from);
                if (at == std::string_view::npos) break;
                const std::size_t after = at + name_attribute.size();
                const bool starts_key =
                    at == 0 || attrs[at - 1] == ' ' || attrs[at - 1] == ';' ||
                    attrs[at - 1] == '\t';
                const bool ends_key = after < attrs.size() &&
                                      (attrs[after] == ' ' || attrs[after] == '=' ||
                                       attrs[after] == '"');
                if (starts_key && ends_key) {
                    key_pos = at;
                    break;
                }
                from = at + 1;
            }
            if (key_pos != std::string_view::npos) {
                std::size_t vs = key_pos + name_attribute.size();
                while (vs < attrs.size() &&
                       (attrs[vs] == ' ' || attrs[vs] == '=' || attrs[vs] == '"')) {
                    ++vs;
                }
                std::size_t ve = vs;
                while (ve < attrs.size() && attrs[ve] != '"' && attrs[ve] != ';') {
                    ++ve;
                }
                // Same bound as read_bed's column 4, and needed more here:
                // read_gtf skips lines it cannot use rather than throwing, so
                // it will harvest an attribute-shaped field out of a file that
                // is not an annotation at all.
                region.name = sanitize_region_name(attrs.substr(vs, ve - vs));
            }
        }
        regions.push_back(std::move(region));
    }
    return regions;
}

// True when `path` names a GTF/GFF by extension, transparently ignoring a
// trailing .gz. Shared by both read_regions entry points so the descriptor
// path and the name path can never dispatch to different parsers.
bool looks_like_gtf(const std::string& path) {
    std::string stem = path;
    if (ends_with_ci(stem, ".gz")) stem.resize(stem.size() - 3);
    return ends_with_ci(stem, ".gtf") || ends_with_ci(stem, ".gff") ||
           ends_with_ci(stem, ".gff3");
}

}  // namespace

std::vector<Region> read_bed(const std::string& path) {
    LineReader reader(path);
    return parse_bed(reader, path);
}

std::vector<Region> read_gtf(const std::string& path, const std::string& feature,
                             const std::string& name_attribute) {
    LineReader reader(path);
    return parse_gtf(reader, path, feature, name_attribute);
}

std::vector<Region> read_regions(const std::string& path) {
    return looks_like_gtf(path) ? read_gtf(path) : read_bed(path);
}

std::vector<Region> read_regions_from_fd(UniqueFd fd,
                                         const std::string& name_for_format) {
    // The descriptor is consumed by the LineReader (which hands it to
    // bgzf_hdopen) before either parser runs, so the file behind it is fixed
    // at the moment the caller opened it. Nothing below resolves the name.
    LineReader reader(std::move(fd), name_for_format);
    return looks_like_gtf(name_for_format)
               ? parse_gtf(reader, name_for_format, "gene", "gene_name")
               : parse_bed(reader, name_for_format);
}

}  // namespace profiler
