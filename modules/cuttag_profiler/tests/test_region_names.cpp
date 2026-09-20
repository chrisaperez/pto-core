// SPDX-License-Identifier: MIT
//
// SECURITY_HTTP_2026-08-15 M5 and M6, for the region readers.
//
// M5: BED column 4 was copied verbatim into SignalMatrix::row_names, which is
// returned as `topRegions[].name` in the /api/profile JSON, as the first
// column of the /api/matrix TSV, and written by `profile --out-matrix`. Two
// separate problems shared that one line -- an unbounded amount of file
// content leaving through the API, and a tab or newline in a name breaking the
// framing of every TSV this tool produces, which lets a crafted BED forge rows
// in a file a downstream tool will parse as data.
//
// M6: the readers used to be handed a path that had already been checked by
// somebody else, and would look it up a second time to open it. These pin the
// descriptor-taking entry point that removed the second lookup.
//
// Registered in the cloud configuration as well as the local one: the batch
// `profile` path writes row names into a TSV, and the batch path is exactly
// what the cloud worker runs.
#include <cstdio>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#ifndef _WIN32
#include <sys/stat.h>  // mkfifo
#endif

#include "profiler/bed_reader.hpp"
#include "profiler/safe_open.hpp"
#include "test_util.hpp"

namespace fs = std::filesystem;

using profiler::open_regular_file;
using profiler::read_regions;
using profiler::read_regions_from_fd;
using profiler::sanitize_region_name;

namespace {

// A file that removes itself. The suite must not leave anything behind in the
// temp directory when a case fails partway through.
class TempFile {
public:
    TempFile(const std::string& suffix, const std::string& contents)
        : path_(fs::temp_directory_path() /
                ("pto_regions_" + testing::unique_tag() + suffix)) {
        std::ofstream out(path_, std::ios::binary);
        out << contents;
    }
    ~TempFile() {
        std::error_code ec;
        fs::remove(path_, ec);
    }
    TempFile(const TempFile&) = delete;
    TempFile& operator=(const TempFile&) = delete;

    [[nodiscard]] const std::string& path() const { return path_str_; }

private:
    fs::path path_;
    std::string path_str_ = path_.string();
};

bool throws(const std::function<void()>& body) {
    try {
        body();
        return false;
    } catch (const std::exception&) {
        return true;
    }
}

}  // namespace

int main() {
    testing::Suite suite{"region_names", {}};

    // ------------------------------------------------------- M5: the bound
    suite.add("M5: an ordinary name is untouched", [] {
        // The mitigation must be invisible to real annotation files, or it
        // will be turned off by the first person it inconveniences.
        CHECK_EQ(sanitize_region_name("ENSG00000141510"), std::string("ENSG00000141510"));
        CHECK_EQ(sanitize_region_name("TP53_promoter-1"), std::string("TP53_promoter-1"));
        CHECK_EQ(sanitize_region_name(""), std::string(""));
        // Not HTML-escaped: the dashboard must not be relying on its input
        // arriving pre-escaped, and escaping here would corrupt real names.
        CHECK_EQ(sanitize_region_name("a&b<c"), std::string("a&b<c"));
        // Non-ASCII passes through: a gene name may legitimately be UTF-8, and
        // every byte of a multi-byte sequence is >= 0x80.
        CHECK_EQ(sanitize_region_name("gène"), std::string("gène"));
    });

    suite.add("M5: a name is bounded in length", [] {
        const std::string huge(profiler::kMaxRegionNameLength * 40, 'A');
        CHECK_EQ(sanitize_region_name(huge).size(), profiler::kMaxRegionNameLength);
        // Exactly at the bound is kept whole; one over is trimmed to it.
        const std::string exact(profiler::kMaxRegionNameLength, 'B');
        CHECK_EQ(sanitize_region_name(exact).size(), profiler::kMaxRegionNameLength);
        CHECK_EQ(sanitize_region_name(exact + "C").size(),
                 profiler::kMaxRegionNameLength);
    });

    suite.add("M5: control characters cannot survive into a name", [] {
        // Tab and newline are the TSV-framing break: without this, a region
        // named "x\tforged\trow" adds columns to /api/matrix's output, and one
        // containing a newline adds a whole row.
        CHECK_EQ(sanitize_region_name("a\tb"), std::string("a?b"));
        CHECK_EQ(sanitize_region_name("a\nb"), std::string("a?b"));
        CHECK_EQ(sanitize_region_name("a\r\nb"), std::string("a??b"));
        // NUL would truncate the name at every C string boundary downstream.
        CHECK_EQ(sanitize_region_name(std::string("a\0b", 3)), std::string("a?b"));
        // ESC opens a terminal escape sequence, and the batch path prints
        // these to the operator's terminal.
        CHECK_EQ(sanitize_region_name("\x1b[31mred"), std::string("?[31mred"));
        // DEL is not below 0x20 and is handled separately; a check that only
        // tested `< 0x20` would let it through.
        CHECK_EQ(sanitize_region_name("a\x7f" "b"), std::string("a?b"));
    });

    suite.add("M5: the bound is applied by read_bed, not just available to it", [] {
        // The unit above proves the function works. This proves it is wired
        // in -- the actual defect was a call site, not a missing helper.
        const std::string payload(4096, 'X');
        TempFile bed(".bed",
                     "chr1\t100\t200\tname\twith\ta\ttab\n"
                     "chr1\t300\t400\t" + payload + "\n");
        const auto regions = read_regions(bed.path());
        CHECK_EQ(regions.size(), std::size_t{2});
        // Column 4 is "name" here; the tabs after it are columns 5-7, so what
        // is being pinned is that the parser splits before this runs.
        CHECK_EQ(regions[0].name, std::string("name"));
        CHECK_EQ(regions[1].name.size(), profiler::kMaxRegionNameLength);
    });

    suite.add("M5: a tab-bearing name cannot reach a name at all", [] {
        // A BED whose column 4 came from a quoted field with an embedded tab
        // is not representable, so the closest reachable case is the GTF
        // attribute path -- which is also the more permissive parser, and
        // therefore the better exfiltration primitive of the two.
        TempFile gtf(".gtf",
                     "chr1\tsrc\tgene\t100\t200\t.\t+\t.\t"
                     "gene_name \"bad\x1bname\";\n");
        const auto regions = read_regions(gtf.path());
        CHECK_EQ(regions.size(), std::size_t{1});
        CHECK_EQ(regions[0].name, std::string("bad?name"));
    });

    suite.add("M5: read_gtf bounds its attribute value too", [] {
        const std::string payload(2000, 'G');
        TempFile gtf(".gtf",
                     "chr1\tsrc\tgene\t100\t200\t.\t+\t.\tgene_name \"" +
                         payload + "\";\n");
        const auto regions = read_regions(gtf.path());
        CHECK_EQ(regions.size(), std::size_t{1});
        CHECK_EQ(regions[0].name.size(), profiler::kMaxRegionNameLength);
    });

    // ------------------------------------------------ coordinate soundness
    suite.add("negative BED coordinates are refused, not carried through", [] {
        // from_chars parses "-1" happily, and a negative start used to reach
        // the window arithmetic in signal_calc -- where it is a quietly wrong
        // answer rather than a refusal.
        TempFile bed(".bed", "chr1\t-5\t200\tx\n");
        CHECK(throws([&] { (void)read_regions(bed.path()); }));

        TempFile ok(".bed", "chr1\t0\t200\tx\n");
        CHECK_EQ(read_regions(ok.path()).size(), std::size_t{1});
    });

    suite.add("GTF coordinate 0 is refused rather than becoming -1", [] {
        // GTF is 1-based inclusive, so `start - 1` on a 0 fabricates a
        // negative coordinate out of a merely malformed file.
        TempFile gtf(".gtf", "chr1\tsrc\tgene\t0\t200\t.\t+\t.\tgene_name \"x\";\n");
        CHECK(throws([&] { (void)read_regions(gtf.path()); }));

        TempFile ok(".gtf", "chr1\tsrc\tgene\t1\t200\t.\t+\t.\tgene_name \"x\";\n");
        const auto regions = read_regions(ok.path());
        CHECK_EQ(regions.size(), std::size_t{1});
        CHECK_EQ(regions[0].start, std::int64_t{0});
    });

    // ---------------------------------------------------- M6: open exactly once
    suite.add("M6: a descriptor parses to the same regions as the path", [] {
        // The fd path is the one the HTTP handler now uses, so it has to be
        // equivalent to the path one and not a second, subtly different parser.
        TempFile bed(".bed", "chr1\t100\t200\talpha\t0\t+\nchr2\t5\t9\tbeta\t0\t-\n");
        const auto by_path = read_regions(bed.path());
        auto fd = open_regular_file(bed.path());
        const auto by_fd = read_regions_from_fd(std::move(fd), bed.path());

        CHECK_EQ(by_fd.size(), by_path.size());
        for (std::size_t i = 0; i < by_fd.size(); ++i) {
            CHECK_EQ(by_fd[i].chrom, by_path[i].chrom);
            CHECK_EQ(by_fd[i].start, by_path[i].start);
            CHECK_EQ(by_fd[i].end, by_path[i].end);
            CHECK_EQ(by_fd[i].name, by_path[i].name);
            CHECK(by_fd[i].strand == by_path[i].strand);
        }
    });

    suite.add("M6: the format still follows the name on the descriptor path", [] {
        // read_regions_from_fd cannot sniff an extension off a descriptor, so
        // it takes the name for that purpose only. If the two entry points
        // dispatched differently, a GTF would parse as a BED and throw.
        TempFile gtf(".gtf",
                     "chr1\tsrc\tgene\t100\t200\t.\t+\t.\tgene_name \"g1\";\n");
        auto fd = open_regular_file(gtf.path());
        const auto regions = read_regions_from_fd(std::move(fd), gtf.path());
        CHECK_EQ(regions.size(), std::size_t{1});
        CHECK_EQ(regions[0].name, std::string("g1"));
    });

    suite.add("M6: open_regular_file refuses a directory", [] {
        // Not a regular file, so it must be rejected at the open rather than
        // discovered later inside htslib.
        CHECK(throws([] {
            (void)open_regular_file(fs::temp_directory_path().string());
        }));
    });

    suite.add("M6: open_regular_file refuses a missing file", [] {
        CHECK(throws([] {
            (void)open_regular_file("/nonexistent/pto/definitely/not/here.bed");
        }));
    });

#ifndef _WIN32
    suite.add("M6: open_regular_file refuses a symlink final component", [] {
        // The race M6 describes: a name that passed confinement is swapped for
        // a symlink before the open. O_NOFOLLOW is what makes the swap fail
        // rather than succeed against the target.
        TempFile target(".bed", "chr1\t1\t2\tx\n");
        const fs::path link = fs::temp_directory_path() /
                              ("pto_link_" + testing::unique_tag() + ".bed");
        std::error_code ec;
        fs::create_symlink(target.path(), link, ec);
        if (!ec) {
            CHECK(throws([&] { (void)open_regular_file(link.string()); }));
            // ...and the target itself still opens, so the refusal is about
            // the symlink and not about the file being unreadable.
            CHECK(!throws([&] { (void)open_regular_file(target.path()); }));
            fs::remove(link, ec);
        }
    });

    suite.add("M6: open_regular_file refuses a FIFO", [] {
        // A named pipe inside the data root passes every path check, and
        // bgzf_open on it parks the calling thread until a writer appears.
        // httplib's pool is fixed-size, so a handful of these is a denial of
        // service -- and this was reachable with confinement fully enabled.
        const fs::path fifo = fs::temp_directory_path() /
                              ("pto_fifo_" + testing::unique_tag());
        if (::mkfifo(fifo.c_str(), 0600) == 0) {
            // O_NONBLOCK is not used here: the point is that the rejection
            // happens without this test blocking, which it would if the check
            // were missing.
            CHECK(throws([&] { (void)open_regular_file(fifo.string()); }));
            std::error_code ec;
            fs::remove(fifo, ec);
        }
    });
#endif

    return suite.run();
}
