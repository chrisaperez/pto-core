// SPDX-License-Identifier: MIT
//
// A deliberately tiny assertion harness. Pulling in a framework would mean
// either a network fetch at configure time or another vendored tree, and the
// project's premise is that a build works on a machine with no route out.
#pragma once

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace testing {

inline int g_failures = 0;
inline int g_checks = 0;

inline void report(bool ok, const char* expression, const char* file, int line,
                   const std::string& detail = {}) {
    ++g_checks;
    if (ok) return;
    ++g_failures;
    std::fprintf(stderr, "  FAIL %s:%d: %s%s%s\n", file, line, expression,
                 detail.empty() ? "" : " -- ", detail.c_str());
}

#define CHECK(expr) ::testing::report((expr), #expr, __FILE__, __LINE__)

#define CHECK_EQ(a, b)                                                     \
    ::testing::report((a) == (b), #a " == " #b, __FILE__, __LINE__,        \
                      ::testing::describe(a) + " vs " + ::testing::describe(b))

#define CHECK_NEAR(a, b, tol)                                              \
    ::testing::report(std::fabs((a) - (b)) <= (tol),                       \
                      #a " ~= " #b, __FILE__, __LINE__,                    \
                      ::testing::describe(a) + " vs " + ::testing::describe(b))

template <typename T>
std::string describe(const T& value) {
    if constexpr (std::is_convertible_v<T, std::string>) {
        return std::string(value);
    } else if constexpr (std::is_floating_point_v<T>) {
        char buffer[32];
        std::snprintf(buffer, sizeof(buffer), "%.6g", static_cast<double>(value));
        return buffer;
    } else {
        return std::to_string(value);
    }
}

// Runs each named case, printing a one-line summary per case.
struct Suite {
    const char* name;
    std::vector<std::pair<const char*, std::function<void()>>> cases;

    void add(const char* case_name, std::function<void()> body) {
        cases.emplace_back(case_name, std::move(body));
    }

    int run() {
        std::printf("== %s ==\n", name);
        for (auto& [case_name, body] : cases) {
            const int before = g_failures;
            body();
            std::printf("  %s %s\n", g_failures == before ? "ok  " : "FAIL",
                        case_name);
        }
        std::printf("%d checks, %d failure(s)\n", g_checks, g_failures);
        return g_failures == 0 ? 0 : 1;
    }
};

// --- Synthetic BAM construction -------------------------------------------

// Distinguishes concurrent test runs without reaching for a POSIX-only
// getpid(); the suite must build on Windows too.
inline std::string unique_tag() {
    static std::mt19937_64 rng{std::random_device{}()};
    static int counter = 0;
    return std::to_string(rng()) + "_" + std::to_string(counter++);
}

struct ReadSpec {
    std::string qname;
    int tid = 0;
    std::int64_t pos = 0;
    std::uint16_t flag = 0;
    std::uint8_t mapq = 60;
    std::vector<std::uint32_t> cigar;  // pre-packed with bam_cigar_gen
    std::int32_t mtid = -1;
    std::int64_t mpos = -1;
    std::int64_t isize = 0;
    std::size_t seq_len = 0;  // 0 => derived from the CIGAR's query length
};

inline std::uint32_t cigar_op(int op, int len) {
    return bam_cigar_gen(static_cast<std::uint32_t>(len),
                         static_cast<std::uint32_t>(op));
}

// A reference contig. `aliases` is the comma-separated @SQ AN: value; leave it
// empty for the usual SN-only line. htslib resolves AN names through
// sam_hdr_name2tid(), so BamReader has to as well.
struct Contig {
    std::string name;
    std::int64_t length = 0;
    std::string aliases;
};

// A BAM plus its index, removed from disk when the object dies.
class TempBam {
public:
    TempBam(const std::vector<Contig>& contigs,
            const std::vector<ReadSpec>& reads) {
        path_ = (std::filesystem::temp_directory_path() /
                 ("cuttag_test_" + unique_tag() + ".bam"))
                    .string();

        sam_hdr_t* header = sam_hdr_init();
        sam_hdr_add_line(header, "HD", "VN", "1.6", "SO", "coordinate", nullptr);
        for (const auto& contig : contigs) {
            const std::string length = std::to_string(contig.length);
            if (contig.aliases.empty()) {
                sam_hdr_add_line(header, "SQ", "SN", contig.name.c_str(), "LN",
                                 length.c_str(), nullptr);
            } else {
                sam_hdr_add_line(header, "SQ", "SN", contig.name.c_str(), "LN",
                                 length.c_str(), "AN", contig.aliases.c_str(),
                                 nullptr);
            }
        }

        htsFile* fp = hts_open(path_.c_str(), "wb");
        if (fp == nullptr || sam_hdr_write(fp, header) < 0) {
            std::fprintf(stderr, "cannot create test BAM at %s\n", path_.c_str());
            std::exit(2);
        }

        bam1_t* rec = bam_init1();
        for (const auto& spec : reads) {
            std::size_t query_len = spec.seq_len;
            if (query_len == 0) {
                for (std::uint32_t op : spec.cigar) {
                    if ((bam_cigar_type(bam_cigar_op(op)) & 1) != 0) {
                        query_len += bam_cigar_oplen(op);
                    }
                }
            }
            const std::string seq(query_len, 'A');
            const std::string qual(query_len, '\x28');  // Phred 40

            const int ret = bam_set1(
                rec, spec.qname.size(), spec.qname.c_str(), spec.flag, spec.tid,
                spec.pos, spec.mapq, spec.cigar.size(), spec.cigar.data(),
                spec.mtid, spec.mpos, spec.isize, seq.size(), seq.c_str(),
                qual.c_str(), 0);
            if (ret < 0 || sam_write1(fp, header, rec) < 0) {
                std::fprintf(stderr, "cannot write test record %s\n",
                             spec.qname.c_str());
                std::exit(2);
            }
        }

        bam_destroy1(rec);
        hts_close(fp);
        sam_hdr_destroy(header);

        if (sam_index_build(path_.c_str(), 0) < 0) {
            std::fprintf(stderr, "cannot index test BAM %s\n", path_.c_str());
            std::exit(2);
        }
    }

    // Convenience overload for the common SN-only case.
    TempBam(const std::vector<std::pair<std::string, std::int64_t>>& contigs,
            const std::vector<ReadSpec>& reads)
        : TempBam(without_aliases(contigs), reads) {}

    ~TempBam() {
        std::error_code ec;
        std::filesystem::remove(path_, ec);
        std::filesystem::remove(path_ + ".bai", ec);
    }

    TempBam(const TempBam&) = delete;
    TempBam& operator=(const TempBam&) = delete;

    [[nodiscard]] const std::string& path() const { return path_; }

private:
    static std::vector<Contig> without_aliases(
        const std::vector<std::pair<std::string, std::int64_t>>& contigs) {
        std::vector<Contig> out;
        out.reserve(contigs.size());
        for (const auto& [name, length] : contigs) out.push_back({name, length, {}});
        return out;
    }

    std::string path_;
};

// A scratch text file (BED, GTF) cleaned up on destruction.
class TempFile {
public:
    TempFile(const std::string& suffix, const std::string& contents) {
        path_ = (std::filesystem::temp_directory_path() /
                 ("cuttag_test_" + unique_tag() + suffix))
                    .string();
        std::FILE* fp = std::fopen(path_.c_str(), "w");
        if (fp == nullptr) std::exit(2);
        std::fwrite(contents.data(), 1, contents.size(), fp);
        std::fclose(fp);
    }
    ~TempFile() {
        std::error_code ec;
        std::filesystem::remove(path_, ec);
    }
    TempFile(const TempFile&) = delete;
    TempFile& operator=(const TempFile&) = delete;

    [[nodiscard]] const std::string& path() const { return path_; }

private:
    std::string path_;
};

}  // namespace testing
