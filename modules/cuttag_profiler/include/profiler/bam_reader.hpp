// SPDX-License-Identifier: MIT
//
// Indexed random access over a coordinate-sorted BAM/CRAM file.
//
// Threading model: htslib's `htsFile` handle and its iterators carry mutable
// decompression state and are *not* safe to share across threads. The index
// (`hts_idx_t`) is immutable once loaded and is shared freely.
//
// The header is NOT immutable, which is the subtle half. `sam_hdr_t` carries a
// lazily-built `hrecs` index that `sam_hdr_name2tid()` populates on its first
// call -- note that function takes a non-const `sam_hdr_t*` while
// `sam_hdr_tid2name()` takes a const one. That build is unsynchronised, so
// `BamReader` resolves contig names out of its own table (built once, in the
// constructor) and serialises the single htslib call that can still reach the
// header. See `tid_for()`.
//
// `BamReader` therefore owns one shared index/header plus a pool of per-thread
// file handles handed out through `HandleLease`.
#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "profiler/types.hpp"

// Opaque htslib forward declarations; the definitions stay in the .cpp.
struct htsFile;
struct hts_idx_t;
struct sam_hdr_t;
struct bam1_t;

namespace profiler {

// Signal accumulated for a single bin. Kept as two independent tallies so a
// caller can switch between `CountMode::kReads` and `CountMode::kDepth`
// without re-reading the BAM.
struct BinTally {
    double read_count = 0.0;  // alignments/fragments overlapping the bin
    double base_count = 0.0;  // aligned bases falling inside the bin
};

// How many BGZF decompression threads a single shared BamReader should build
// a pool for, given how many region-query worker threads the caller itself
// will run. A pure function of that one input so the decision has exactly one
// place to live and one thing pinning it: `worker_threads <= 1` used to select
// `hardware_concurrency()` here, on the reasoning that a lone worker could use
// the help. Measured on 100k +/-2kb region queries against an 87 MB /
// 8M-alignment BAM (docs/AUDIT_2026-09-11_cuttag_profiler.md finding 1), that
// reasoning was backwards: 17.32s wall / 34.53s user / 37.34s sys with a pool
// of 12, against 5.08s / 4.84s / 0.22s with none -- 3.4x slower wall time and
// ~14 core-seconds of system time burned per wall second, on a "single
// thread" run. The reason is structural: an indexed region query calls
// bgzf_seek per region, which invalidates whatever block was already decoded,
// so there is no run of sequential blocks under one query for a shared pool
// to parallel-inflate ahead of -- only a queue hand-off cost paid on every
// single block, for every one of a hundred thousand tiny queries.
[[nodiscard]] constexpr int choose_bgzf_thread_count(int /*worker_threads*/) noexcept {
    return 0;
}

class BamReader {
public:
    // Opens `path` and loads its index (.bai/.csi/.crai). Throws
    // std::runtime_error if either cannot be read.
    //
    // `decompression_threads` controls BGZF inflation *within* a single query:
    //   0 or 1 (default) - no shared pool; each leased handle inflates on the
    //                      thread that owns it. This is what you want whenever
    //                      the caller parallelises across regions, which is the
    //                      profiler's normal mode -- and also, per
    //                      `choose_bgzf_thread_count` above, what a
    //                      single-threaded region-query caller wants too.
    //   > 1              - build a shared htsThreadPool of that size. NOT a
    //                      general win for "a single-threaded caller" as
    //                      such: it only helps a caller whose query holds a
    //                      run of sequential BGZF blocks to overlap inflation
    //                      across, e.g. a handful of whole-chromosome-sized
    //                      windows. An indexed region-query caller calls
    //                      bgzf_seek per region and never has such a run --
    //                      see choose_bgzf_thread_count's measurements before
    //                      reaching for this.
    explicit BamReader(std::string path, int decompression_threads = 0);
    ~BamReader();

    BamReader(const BamReader&) = delete;
    BamReader& operator=(const BamReader&) = delete;

    [[nodiscard]] const std::string& path() const noexcept { return path_; }
    [[nodiscard]] const std::vector<std::string>& chromosomes() const noexcept {
        return chrom_names_;
    }
    [[nodiscard]] const std::vector<std::int64_t>& chromosome_lengths()
        const noexcept {
        return chrom_lengths_;
    }
    // Returns -1 when the name is absent from the header.
    [[nodiscard]] int tid_for(std::string_view chrom) const;
    // Returns -1 when the name is absent from the header.
    [[nodiscard]] std::int64_t chromosome_length(std::string_view chrom) const;

    // Total mapped alignments, read from the index metadata in O(#contigs)
    // rather than by streaming the file. This is the denominator for CPM/RPKM.
    [[nodiscard]] std::uint64_t total_mapped_reads() const noexcept {
        return total_mapped_reads_;
    }

    // Accumulates signal for `region` into `out`, which must already be sized
    // to the desired bin count. Bins are laid out left-to-right in ascending
    // genomic coordinate; strand flipping is the caller's responsibility.
    // Thread-safe: may be called concurrently from any number of threads.
    void query_bins(const Region& region, std::span<BinTally> out,
                    const FilterOptions& filters) const;

    // Number of alignments passing `filters` in `region`. Thread-safe.
    [[nodiscard]] std::uint64_t count_reads(const Region& region,
                                            const FilterOptions& filters) const;

private:
    // A borrowed per-thread htsFile handle, returned to the pool on scope exit.
    class HandleLease;
    struct Handle;

    [[nodiscard]] HandleLease acquire() const;
    void release(std::unique_ptr<Handle> handle) const;

    std::string path_;
    int decompression_threads_ = 0;

    hts_idx_t* index_ = nullptr;
    sam_hdr_t* header_ = nullptr;
    void* thread_pool_ = nullptr;  // htsThreadPool*, owned

    // Contig name -> tid. Transparent hash/equality so `tid_for()` can look up a
    // std::string_view without materialising a std::string per call.
    struct NameHash {
        using is_transparent = void;
        std::size_t operator()(std::string_view s) const noexcept {
            return std::hash<std::string_view>{}(s);
        }
    };
    using NameTidMap =
        std::unordered_map<std::string, int, NameHash, std::equal_to<>>;

    std::vector<std::string> chrom_names_;
    std::vector<std::int64_t> chrom_lengths_;
    std::uint64_t total_mapped_reads_ = 0;

    // @SQ SN names, filled once in the constructor and const thereafter.
    NameTidMap tid_by_name_;

    // Guards the one remaining call into sam_hdr_name2tid(), which mutates
    // `header_` on first use. `alias_tids_` memoises its answers, negatives
    // included, so a miss costs the lock at most once per distinct name.
    mutable std::mutex header_mutex_;
    mutable NameTidMap alias_tids_;

    mutable std::mutex pool_mutex_;
    mutable std::vector<std::unique_ptr<Handle>> free_handles_;
};

}  // namespace profiler
