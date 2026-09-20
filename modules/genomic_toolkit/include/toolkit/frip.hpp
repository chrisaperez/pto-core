// SPDX-License-Identifier: MIT
//
// FRiP -- fraction of reads in peaks -- and the BED peak reader it needs.
//
// The conventional recipe is
//
//     bedtools intersect -a frags.bed -b peaks.bed -u | wc -l
//
// which writes every overlapping fragment back out as text purely to have
// `wc` count the lines. Here the count is a reduction over the streaming pass:
// nothing is written, and the peaks are a merged disjoint cover (see
// peak_set.hpp) so a fragment overlapping three peaks counts once, which is
// what the statistic means and what the naive pipeline gets wrong.
//
// Two denominators are reported, because the literature uses both and papers
// rarely say which:
//   * fragment FRiP  -- fragments touching a peak / fragments passing filters
//   * base FRiP      -- fragment bases inside peaks / total fragment bases
// They differ whenever fragment length correlates with peak occupancy, which
// in ATAC-seq it does.
#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "toolkit/contig_dict.hpp"
#include "toolkit/peak_set.hpp"
#include "toolkit/types.hpp"

namespace toolkit {

struct FripResult {
    std::uint64_t fragments_total = 0;
    std::uint64_t fragments_in_peaks = 0;
    std::uint64_t duplicates_excluded = 0;
    std::int64_t bases_total = 0;
    std::int64_t bases_in_peaks = 0;

    std::size_t peaks_input = 0;
    std::size_t peaks_merged = 0;
    std::int64_t peak_bases = 0;

    // Contigs present in the fragment stream with no counterpart in the peak
    // file. A "chr1" vs "1" mismatch shows up here as every contig, and is the
    // difference between a real FRiP of 0.01 and a naming bug.
    std::vector<std::string> unmatched_contigs;

    [[nodiscard]] double fragment_frip() const noexcept {
        return fragments_total == 0 ? 0.0
                                    : static_cast<double>(fragments_in_peaks) /
                                          static_cast<double>(fragments_total);
    }
    [[nodiscard]] double base_frip() const noexcept {
        return bases_total == 0 ? 0.0
                                : static_cast<double>(bases_in_peaks) /
                                      static_cast<double>(bases_total);
    }
    // FRiP divided by the FRiP a uniformly scattered library would give. Below
    // ~2 means the peaks are barely enriched over background, which a raw FRiP
    // of 0.2 on a peak set covering 10% of the genome does not make obvious.
    [[nodiscard]] double enrichment(std::int64_t genome_bases) const noexcept {
        if (genome_bases <= 0 || peak_bases <= 0) return 0.0;
        const double expected =
            static_cast<double>(peak_bases) / static_cast<double>(genome_bases);
        return expected > 0.0 ? fragment_frip() / expected : 0.0;
    }
};

// Accumulates FRiP over batches. Held by the streaming callback; call
// `finish()` once the stream is done to attach the naming diagnostics.
class FripCounter {
public:
    // `peaks` must outlive the counter. `bind()` must already have been called
    // on it against the dictionary the stream fills.
    FripCounter(const PeakSet& peaks, bool skip_duplicates)
        : peaks_(peaks), skip_duplicates_(skip_duplicates) {
        result_.peaks_input = peaks.input_count();
        result_.peaks_merged = peaks.size();
        result_.peak_bases = peaks.covered_bases();
    }

    // Counts one batch. The peak set is immutable and each fragment's flag
    // write is to its own element, so the loop parallelises with a plain
    // reduction over the five counters -- no locking, and the result is
    // bit-identical to the serial order because every accumulator is an
    // integer.
    //
    // OpenMP is optional here (unlike scrna_matrix, where it is mandatory):
    // the queries are pointer-chasing binary searches with essentially no
    // arithmetic, so on a warm page cache this is memory-latency bound and
    // scaling flattens out around four threads. A batch is 4096 fragments,
    // which is enough work to cover the fork/join cost and small enough to
    // stay in L2.
    void add_batch(std::span<Fragment> batch) {
        std::uint64_t n_total = 0, n_in = 0, n_dup = 0;
        std::int64_t b_total = 0, b_in = 0;
        const auto count = static_cast<std::ptrdiff_t>(batch.size());
        Fragment* data = batch.data();
        const PeakSet& peaks = peaks_;
        const bool skip_dups = skip_duplicates_;

        // One OverlapCursor per OpenMP thread, allocated once for the whole
        // run (not per batch, not per fragment) and never shrunk. Because
        // every batch but the last is exactly kBatchFragments long and
        // schedule(static) hands out the same fixed index range to the same
        // thread every time, thread T sees a stable, disjoint,
        // monotonically-advancing slice of the (coordinate-sorted) stream
        // across batch boundaries -- so its cursor keeps paying off run to
        // run rather than starting cold at the top of every batch. This is a
        // performance property only: PeakSet::overlap_bases resets a cursor
        // itself whenever the resolved contig changes or a query goes
        // backwards, so even a runtime that renumbered threads between
        // batches would still answer correctly, just less cheaply.
#ifdef _OPENMP
        const std::size_t max_threads = static_cast<std::size_t>(omp_get_max_threads());
#else
        const std::size_t max_threads = 1;
#endif
        if (cursors_.size() < max_threads) cursors_.resize(max_threads);
        OverlapCursor* const cursors = cursors_.data();

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (count >= 512) \
    reduction(+ : n_total, n_in, n_dup, b_total, b_in)
#endif
        for (std::ptrdiff_t i = 0; i < count; ++i) {
            Fragment& f = data[i];
            if (skip_dups && f.duplicate()) {
                ++n_dup;
                continue;
            }
            ++n_total;
            b_total += f.length();
#ifdef _OPENMP
            OverlapCursor& cursor = cursors[static_cast<std::size_t>(omp_get_thread_num())];
#else
            OverlapCursor& cursor = cursors[0];
#endif
            const std::int64_t hit = peaks.overlap_bases(f.tid, f.start, f.end, cursor);
            if (hit > 0) {
                ++n_in;
                b_in += hit;
                f.flags |= fragflags::kInPeak;
            }
        }

        // Inside the loop a batch cannot overflow: 4096 fragments of at most
        // 2^31 - 1 bp is under 2^44. Across batches it can -- see checked_add
        // -- so the running totals are checked here, once per batch, outside
        // the parallel region (an exception may not leave one).
        std::int64_t bases_total = result_.bases_total;
        std::int64_t bases_in_peaks = result_.bases_in_peaks;
        checked_add(bases_total, b_total, "FRiP total fragment bases");
        checked_add(bases_in_peaks, b_in, "FRiP in-peak bases");
        result_.fragments_total += n_total;
        result_.fragments_in_peaks += n_in;
        result_.duplicates_excluded += n_dup;
        result_.bases_total = bases_total;
        result_.bases_in_peaks = bases_in_peaks;
    }

    // Records which stream contigs the peak file never mentioned.
    FripResult finish(const ContigDict& stream_dict) {
        result_.unmatched_contigs = stream_dict.unmatched_against(peaks_.dict());
        return result_;
    }

    [[nodiscard]] const FripResult& result() const noexcept { return result_; }

private:
    const PeakSet& peaks_;
    bool skip_duplicates_;
    FripResult result_;
    // See add_batch(): sized lazily to the OpenMP thread count on first use,
    // then held for the life of the counter.
    std::vector<OverlapCursor> cursors_;
};

// Reads a BED3+ peak file. Track/browser/comment lines are skipped; throws
// std::runtime_error on an unreadable file and std::invalid_argument on a
// malformed coordinate. Peak files are small (1e4-1e6 lines) so this one does
// materialise its result -- PeakSet::build needs the whole set anyway.
std::vector<Region> read_peaks_bed(const std::string& path);

}  // namespace toolkit
