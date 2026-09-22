// SPDX-License-Identifier: MIT
//
// Buffered line input, internal to the module.
//
// std::getline on an ifstream is the obvious choice and the wrong one: it
// constructs a std::string per line, which for a 200M-line fragment file is
// 200M allocations to produce data that is consumed and discarded immediately.
// This reader owns one growable buffer, refills it in kChunkBytes reads, and
// hands out std::string_view into it. Nothing is copied and nothing is
// allocated per line after the buffer reaches steady state.
//
// When htslib is present the underlying handle is BGZF, which reads plain text
// and gzip/BGZF transparently -- so a .gz fragment file costs no extra code
// path. Without htslib the reader is a plain file handle and a .gz input is
// rejected with a message saying so, rather than parsed as binary noise.
#pragma once

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#ifdef GTK_HAVE_HTSLIB
#include <htslib/bgzf.h>
#endif

namespace toolkit::detail {

class LineReader {
public:
    explicit LineReader(const std::string& path) : path_(path) {
#ifdef GTK_HAVE_HTSLIB
        bgzf_ = bgzf_open(path.c_str(), "r");
        if (bgzf_ == nullptr) {
            throw std::runtime_error("cannot open '" + path + "' for reading");
        }
#else
        if (path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0) {
            throw std::runtime_error(
                "'" + path +
                "' is gzip-compressed, but this build has no htslib. Rebuild "
                "with -DGTK_WITH_HTSLIB=ON, or decompress the input first.");
        }
        file_ = std::fopen(path.c_str(), "rb");
        if (file_ == nullptr) {
            throw std::runtime_error("cannot open '" + path + "' for reading");
        }
#endif
        buffer_.resize(kChunkBytes);
    }

    ~LineReader() {
#ifdef GTK_HAVE_HTSLIB
        if (bgzf_ != nullptr) bgzf_close(bgzf_);
#else
        if (file_ != nullptr) std::fclose(file_);
#endif
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
    // following the input up to an OOM kill. See the comment at the cap.
    bool refill() {
        if (pos_ > 0) {
            const std::size_t tail = filled_ - pos_;
            if (tail > 0) std::memmove(buffer_.data(), buffer_.data() + pos_, tail);
            filled_ = tail;
            pos_ = 0;
        }
        if (filled_ == buffer_.size()) {
            // Growth is capped.
            //
            // The buffer doubles so that a single long line is never split
            // into two malformed records, and that is right -- but with no
            // ceiling it is also an unbounded allocation driven straight by
            // input. Measured on this reader before the cap: a file that is
            // one unterminated 200 MB line drives peak RSS to 518 MB (~2.6x,
            // from the doubled buffer plus the memmove), and a crafted
            // multi-gigabyte line scales the same way until the process dies.
            //
            // fastq_stream already solved this in the same shape: it refuses a
            // record past its 64 KiB carry limit with a message naming the
            // limit, rather than growing to meet it. This is that rule, in the
            // module that was missing it.
            //
            // 64 MiB is four orders of magnitude above any real record. A
            // fragment BED line is a few hundred bytes; the longest plausible
            // BEDPE line with a full name column is a few kilobytes. Anything
            // past the cap is a binary file being read as text, a truncated
            // download, or a hostile input -- and for all three, failing with
            // the byte count is more useful than an OOM kill with none.
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
        long got = 0;
#ifdef GTK_HAVE_HTSLIB
        got = static_cast<long>(bgzf_read(bgzf_, buffer_.data() + filled_, want));
        if (got < 0) throw std::runtime_error("read error on '" + path_ + "'");
#else
        got = static_cast<long>(std::fread(buffer_.data() + filled_, 1, want, file_));
        if (got == 0 && std::ferror(file_) != 0) {
            throw std::runtime_error("read error on '" + path_ + "'");
        }
#endif
        filled_ += static_cast<std::size_t>(got);
        return got > 0;
    }

    static constexpr std::size_t kChunkBytes = 1 << 20;   // 1 MiB
    // Ceiling on a single line; see refill(). Public so the test suite can
    // assert the boundary rather than hard-coding the number twice.
public:
    static constexpr std::size_t kMaxLineBytes = 64ull << 20;  // 64 MiB
private:

    std::string path_;
    std::vector<char> buffer_;
    std::size_t filled_ = 0;
    std::size_t pos_ = 0;
#ifdef GTK_HAVE_HTSLIB
    BGZF* bgzf_ = nullptr;
#else
    std::FILE* file_ = nullptr;
#endif
};

}  // namespace toolkit::detail
