// SPDX-License-Identifier: MIT
#include "profiler/bam_reader.hpp"

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>

#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <thread>
#include <utility>

namespace profiler {
namespace {

// Half-open interval helper used while walking CIGAR blocks.
struct Span {
    std::int64_t begin;
    std::int64_t end;
};

// Returns the reference span an alignment should contribute, honouring the
// fragment-extension policy. `keep` is set to false when the alignment must be
// ignored entirely (e.g. the right-hand mate of a pair already counted through
// its leftmost partner).
Span alignment_span(const bam1_t* rec, const FilterOptions& filters,
                    bool& keep) {
    keep = true;
    const std::int64_t pos = rec->core.pos;
    const std::int64_t aligned_end = bam_endpos(rec);
    const std::uint16_t flag = rec->core.flag;

    if (filters.extend_to_fragment && (flag & flags::kPaired) != 0 &&
        (flag & flags::kProperPair) != 0 && rec->core.isize != 0) {
        // Count each fragment exactly once, through its leftmost mate.
        if (rec->core.isize < 0) {
            keep = false;
            return {0, 0};
        }
        const std::int64_t fragment_end = pos + rec->core.isize;
        if (filters.max_fragment_length > 0 &&
            rec->core.isize > filters.max_fragment_length) {
            keep = false;
            return {0, 0};
        }
        return {pos, std::max(fragment_end, aligned_end)};
    }

    if (filters.extend_reads_to > 0) {
        // Single-end extension in the direction the read points.
        if ((flag & flags::kReverse) != 0) {
            return {std::max<std::int64_t>(0, aligned_end - filters.extend_reads_to),
                    aligned_end};
        }
        return {pos, pos + filters.extend_reads_to};
    }

    return {pos, aligned_end};
}

// Adds `weight` per base to every bin overlapping [begin, end).
void add_bases(std::span<BinTally> bins, std::int64_t window_start,
               std::int64_t window_end, Span s) {
    const std::int64_t lo = std::max(s.begin, window_start);
    const std::int64_t hi = std::min(s.end, window_end);
    if (lo >= hi) return;

    const auto nbins = static_cast<std::int64_t>(bins.size());
    const std::int64_t width = window_end - window_start;

    // Bin i covers [window_start + width*i/nbins, window_start + width*(i+1)/nbins).
    auto bin_of = [&](std::int64_t coord) {
        return ((coord - window_start) * nbins) / width;
    };
    std::int64_t first = std::clamp<std::int64_t>(bin_of(lo), 0, nbins - 1);
    std::int64_t last = std::clamp<std::int64_t>(bin_of(hi - 1), 0, nbins - 1);

    for (std::int64_t b = first; b <= last; ++b) {
        const std::int64_t b_lo = window_start + (width * b) / nbins;
        const std::int64_t b_hi = window_start + (width * (b + 1)) / nbins;
        const std::int64_t overlap =
            std::min(hi, b_hi) - std::max(lo, b_lo);
        if (overlap > 0) {
            bins[static_cast<std::size_t>(b)].base_count +=
                static_cast<double>(overlap);
        }
    }
}

}  // namespace

// --- Handle pool -----------------------------------------------------------

struct BamReader::Handle {
    htsFile* fp = nullptr;
    bam1_t* rec = nullptr;

    ~Handle() {
        if (rec != nullptr) bam_destroy1(rec);
        if (fp != nullptr) hts_close(fp);
    }
};

class BamReader::HandleLease {
public:
    HandleLease(const BamReader* owner, std::unique_ptr<Handle> handle)
        : owner_(owner), handle_(std::move(handle)) {}
    ~HandleLease() {
        if (handle_) owner_->release(std::move(handle_));
    }
    HandleLease(HandleLease&&) noexcept = default;
    HandleLease(const HandleLease&) = delete;
    HandleLease& operator=(const HandleLease&) = delete;

    Handle* operator->() const noexcept { return handle_.get(); }

private:
    const BamReader* owner_;
    std::unique_ptr<Handle> handle_;
};

BamReader::BamReader(std::string path, int decompression_threads)
    : path_(std::move(path)), decompression_threads_(decompression_threads) {
    // A shared htsThreadPool is the wrong tool for this access pattern, and
    // measurably so. When the caller already runs one worker thread per core
    // and each worker holds its own handle, routing every BGZF block through a
    // shared pool funnels all inflation into the pool's threads and adds a
    // queue hand-off per block. Measured on 100k promoter windows over an 87 MB
    // BAM (8 cores, Apple M-series):
    //
    //   shared pool, 8 threads   5.95 s wall, 17.9 s user, 34.6 s sys, 133 MB
    //   shared pool, 1 thread    6.26 s wall,  7.6 s user,  2.9 s sys,  89 MB
    //   no pool (inline)         1.38 s wall,  5.3 s user,  0.9 s sys,  85 MB
    //
    // So the default is no pool at all: each handle inflates on whichever
    // worker thread owns it, which parallelises perfectly across regions.
    // A pool is only worth creating for a caller that queries single-threaded,
    // where it can still overlap inflation within one query.
    const bool want_pool = decompression_threads_ > 1;
    if (want_pool) {
        auto* pool = new htsThreadPool{nullptr, 0};
        pool->pool = hts_tpool_init(decompression_threads_);
        if (pool->pool == nullptr) {
            delete pool;
            throw std::runtime_error("failed to create htslib thread pool");
        }
        thread_pool_ = pool;
    }

    // A bootstrap handle loads the header; the index is shared thereafter.
    htsFile* fp = hts_open(path_.c_str(), "r");
    if (fp == nullptr) {
        throw std::runtime_error("cannot open alignment file: " + path_);
    }
    if (thread_pool_ != nullptr) {
        hts_set_opt(fp, HTS_OPT_THREAD_POOL, thread_pool_);
    }

    header_ = sam_hdr_read(fp);
    if (header_ == nullptr) {
        hts_close(fp);
        throw std::runtime_error("cannot read BAM header: " + path_);
    }

    index_ = sam_index_load(fp, path_.c_str());
    if (index_ == nullptr) {
        sam_hdr_destroy(header_);
        hts_close(fp);
        throw std::runtime_error(
            "cannot load index for " + path_ +
            " (run `samtools index` first; queries require a .bai/.csi)");
    }

    const int n = sam_hdr_nref(header_);
    chrom_names_.reserve(static_cast<std::size_t>(n));
    chrom_lengths_.reserve(static_cast<std::size_t>(n));
    tid_by_name_.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        chrom_names_.emplace_back(sam_hdr_tid2name(header_, i));
        chrom_lengths_.push_back(sam_hdr_tid2len(header_, i));

        // Resolve names from this table rather than from the header, so that
        // concurrent queries never trigger htslib's lazy header parse. Both
        // sam_hdr_tid2name and sam_hdr_tid2len read the plain target arrays and
        // leave `hrecs` alone, so building it here touches nothing lazy.
        // First entry wins on a duplicate SN, which is a malformed header.
        tid_by_name_.emplace(chrom_names_.back(), i);

        // Index metadata gives the mapped/unmapped tallies without a full pass.
        std::uint64_t mapped = 0;
        std::uint64_t unmapped = 0;
        if (hts_idx_get_stat(index_, i, &mapped, &unmapped) == 0) {
            total_mapped_reads_ += mapped;
        }
    }

    // Recycle the bootstrap handle as the first pooled one.
    auto handle = std::make_unique<Handle>();
    handle->fp = fp;
    handle->rec = bam_init1();
    free_handles_.push_back(std::move(handle));
}

BamReader::~BamReader() {
    free_handles_.clear();
    if (index_ != nullptr) hts_idx_destroy(index_);
    if (header_ != nullptr) sam_hdr_destroy(header_);
    if (thread_pool_ != nullptr) {
        auto* pool = static_cast<htsThreadPool*>(thread_pool_);
        if (pool->pool != nullptr) hts_tpool_destroy(pool->pool);
        delete pool;
    }
}

BamReader::HandleLease BamReader::acquire() const {
    {
        std::lock_guard<std::mutex> lock(pool_mutex_);
        if (!free_handles_.empty()) {
            auto handle = std::move(free_handles_.back());
            free_handles_.pop_back();
            return HandleLease(this, std::move(handle));
        }
    }

    auto handle = std::make_unique<Handle>();
    handle->fp = hts_open(path_.c_str(), "r");
    if (handle->fp == nullptr) {
        throw std::runtime_error("cannot re-open alignment file: " + path_);
    }
    if (thread_pool_ != nullptr) {
        hts_set_opt(handle->fp, HTS_OPT_THREAD_POOL, thread_pool_);
    }
    // The header must be consumed so the stream is positioned for seeking; the
    // copy is discarded because the shared header_ is authoritative.
    sam_hdr_t* scratch = sam_hdr_read(handle->fp);
    if (scratch != nullptr) sam_hdr_destroy(scratch);
    handle->rec = bam_init1();
    return HandleLease(this, std::move(handle));
}

void BamReader::release(std::unique_ptr<Handle> handle) const {
    std::lock_guard<std::mutex> lock(pool_mutex_);
    free_handles_.push_back(std::move(handle));
}

int BamReader::tid_for(std::string_view chrom) const {
    // Fast path: the SN table built in the constructor. This is a correctness
    // requirement, not an optimisation.
    //
    // sam_hdr_name2tid() takes a non-const `sam_hdr_t*` because on its first
    // call it parses the header text into an hrecs index (sam_hdr_fill_hrecs),
    // and that lazy build is unsynchronised. compute_matrix() runs one worker
    // thread per core and every worker lands here, so several threads can each
    // observe `hrecs == nullptr`, each build a complete index, and store it over
    // one another: a data race on `hdr->hrecs`, plus every loser's index
    // orphaned with no pointer left to free it. ASan reported exactly that --
    // ~160 KB in 22 allocations under sam_hdr_fill_hrecs, reached from
    // signal_calc.cpp's worker pool through chromosome_length().
    if (const auto it = tid_by_name_.find(chrom); it != tid_by_name_.end()) {
        return it->second;
    }

    // A miss is not necessarily an unknown contig. sam_hdr_name2tid() also
    // resolves @SQ AN: alternative names, which are not in the SN table
    // (verified against htslib 1.22: with `@SQ SN:chr1 AN:1,NC_000001.11`, both
    // "1" and "NC_000001.11" return chr1's tid). Dropping that would silently
    // start skipping regions for anyone whose BED uses Ensembl or RefSeq contig
    // names against an aliased header, so fall through to htslib -- but under
    // the lock, so the lazy build happens exactly once.
    //
    // Memoise the result, negatives included: a region file naming contigs this
    // BAM does not have would otherwise funnel every worker through this mutex
    // once per region. The table is bounded by the number of distinct names the
    // caller asks about, which it is already holding in memory.
    const std::string name(chrom);
    std::lock_guard<std::mutex> lock(header_mutex_);
    if (const auto it = alias_tids_.find(name); it != alias_tids_.end()) {
        return it->second;
    }
    const int tid = sam_hdr_name2tid(header_, name.c_str());
    alias_tids_.emplace(name, tid);
    return tid;
}

std::int64_t BamReader::chromosome_length(std::string_view chrom) const {
    const int tid = tid_for(chrom);
    if (tid < 0) return -1;
    return chrom_lengths_[static_cast<std::size_t>(tid)];
}

void BamReader::query_bins(const Region& region, std::span<BinTally> out,
                           const FilterOptions& filters) const {
    if (out.empty()) return;
    const int tid = tid_for(region.chrom);
    if (tid < 0) return;

    // Widen the iterator window to the left so that fragments starting before
    // the region but reaching into it are still visited. Without this, extended
    // fragments are silently truncated at the window edge.
    const std::int64_t lookback =
        filters.extend_to_fragment
            ? std::max<std::int64_t>(filters.max_fragment_length, 1000)
            : std::max<std::int64_t>(filters.extend_reads_to, 0);
    const std::int64_t query_start = std::max<std::int64_t>(0, region.start - lookback);

    auto lease = acquire();
    hts_itr_t* iter =
        sam_itr_queryi(index_, tid, query_start, region.end);
    if (iter == nullptr) return;

    bam1_t* rec = lease->rec;
    const auto nbins = static_cast<std::int64_t>(out.size());
    const std::int64_t width = region.end - region.start;

    auto bin_of = [&](std::int64_t coord) {
        return ((coord - region.start) * nbins) / width;
    };

    int ret = 0;
    while ((ret = sam_itr_next(lease->fp, iter, rec)) >= 0) {
        if (!filters.passes_flags(rec->core.flag)) continue;
        if (rec->core.qual < filters.min_mapq) continue;

        bool keep = true;
        const Span span = alignment_span(rec, filters, keep);
        if (!keep) continue;
        if (span.end <= region.start || span.begin >= region.end) continue;

        // Read-count tally: one unit to every bin the fragment touches.
        const std::int64_t lo = std::max(span.begin, region.start);
        const std::int64_t hi = std::min(span.end, region.end);
        if (lo < hi) {
            const std::int64_t first =
                std::clamp<std::int64_t>(bin_of(lo), 0, nbins - 1);
            const std::int64_t last =
                std::clamp<std::int64_t>(bin_of(hi - 1), 0, nbins - 1);
            for (std::int64_t b = first; b <= last; ++b) {
                out[static_cast<std::size_t>(b)].read_count += 1.0;
            }
        }

        // Base-depth tally: CIGAR-aware when using the raw alignment span, so
        // deletions and reference skips do not inflate coverage. Extended
        // fragments are treated as contiguous, which is the intended semantics.
        const bool extended = (span.begin != rec->core.pos ||
                               span.end != bam_endpos(rec));
        if (extended) {
            add_bases(out, region.start, region.end, span);
        } else {
            std::int64_t ref_pos = rec->core.pos;
            const std::uint32_t* cigar = bam_get_cigar(rec);
            for (std::uint32_t i = 0; i < rec->core.n_cigar; ++i) {
                const int op = bam_cigar_op(cigar[i]);
                const std::int64_t len = bam_cigar_oplen(cigar[i]);
                const int type = bam_cigar_type(op);
                const bool consumes_ref = (type & 2) != 0;
                const bool consumes_query = (type & 1) != 0;
                if (consumes_ref && consumes_query) {
                    add_bases(out, region.start, region.end,
                              {ref_pos, ref_pos + len});
                }
                if (consumes_ref) ref_pos += len;
            }
        }
    }

    hts_itr_destroy(iter);
    if (ret < -1) {
        throw std::runtime_error("truncated or corrupt BAM record in " + path_);
    }
}

std::uint64_t BamReader::count_reads(const Region& region,
                                     const FilterOptions& filters) const {
    const int tid = tid_for(region.chrom);
    if (tid < 0) return 0;

    auto lease = acquire();
    hts_itr_t* iter = sam_itr_queryi(index_, tid, region.start, region.end);
    if (iter == nullptr) return 0;

    std::uint64_t count = 0;
    bam1_t* rec = lease->rec;
    while (sam_itr_next(lease->fp, iter, rec) >= 0) {
        if (!filters.passes_flags(rec->core.flag)) continue;
        if (rec->core.qual < filters.min_mapq) continue;
        bool keep = true;
        const Span span = alignment_span(rec, filters, keep);
        if (!keep) continue;
        if (span.end <= region.start || span.begin >= region.end) continue;
        ++count;
    }
    hts_itr_destroy(iter);
    return count;
}

}  // namespace profiler
