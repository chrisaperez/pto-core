// SPDX-License-Identifier: MIT
#include "profiler/signal_calc.hpp"

#include "profiler/tsv_format.hpp"
#include "profiler/worker_abort.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <mutex>
#include <new>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <vector>

namespace profiler {
namespace {

// A heap buffer whose first element sits on a cache-line boundary. Each worker
// reduces into its own instance, so no two threads ever dirty the same line.
class AlignedBuffer {
public:
    explicit AlignedBuffer(std::size_t count) : count_(count) {
        // Round the allocation up to a whole number of cache lines so that
        // neighbouring buffers cannot share a trailing line either.
        const std::size_t bytes = ((count * sizeof(double) + kCacheLineSize - 1) /
                                   kCacheLineSize) *
                                  kCacheLineSize;
        data_ = static_cast<double*>(
            ::operator new(bytes, std::align_val_t{kCacheLineSize}));
        std::fill_n(data_, bytes / sizeof(double), 0.0);
    }
    ~AlignedBuffer() {
        ::operator delete(data_, std::align_val_t{kCacheLineSize});
    }
    AlignedBuffer(const AlignedBuffer&) = delete;
    AlignedBuffer& operator=(const AlignedBuffer&) = delete;

    [[nodiscard]] double* data() noexcept { return data_; }
    [[nodiscard]] std::size_t size() const noexcept { return count_; }

private:
    double* data_ = nullptr;
    std::size_t count_ = 0;
};

std::int64_t anchor_for(const Region& region, const ProfileOptions& opts) {
    const bool reverse =
        opts.respect_strand && region.strand == Strand::kReverse;
    switch (opts.reference_point) {
        case ReferencePoint::kTSS:
            return reverse ? region.end : region.start;
        case ReferencePoint::kTES:
            return reverse ? region.start : region.end;
        case ReferencePoint::kCenter:
            // NOT (start + end) / 2. That sum overflows int64 for coordinates
            // near the type's limit, which is undefined behaviour, and the
            // saturating arithmetic in make_window() below cannot help: the
            // overflow happens HERE, before the anchor is ever handed to it.
            //
            // Confirmed under UBSan from a plain BED line --
            // "chr1 9223372036854775000 9223372036854775007" with
            // --reference-point center:
            //
            //   signed integer overflow: 9223372036854775000 +
            //   9223372036854775007 cannot be represented in type 'int64_t'
            //
            // And the consequence was worse than the UB. The wrapped sum is
            // always negative, so for a narrow band of inputs the anchor lands
            // just below zero, make_window's off-contig test passes, and the
            // region is PROFILED at coordinate ~0 instead of being skipped --
            // a locus at 9.2e18 contributing real signal from the start of the
            // contig to the meta-profile. The same region under --reference-point
            // TSS is correctly skipped, which is how the two disagreed.
            //
            // This form cannot overflow, given what the readers guarantee:
            // bed_reader's parse_bed rejects a negative coordinate and
            // `end < start`, and parse_gtf rejects `start < 1` and
            // `end < start - 1`, so every Region reaching here has
            // 0 <= start <= end. Then `end - start` is non-negative and at most
            // `end`, and `start + (end - start) / 2` is at most `end` -- both
            // within range for any representable input. It is also the exact
            // midpoint rather than merely a safe one.
            return region.start + (region.end - region.start) / 2;
    }
    return region.start;
}

}  // namespace

bool make_window(const Region& region, const ProfileOptions& opts,
                 std::int64_t chrom_length, Region& out_window) {
    const bool reverse =
        opts.respect_strand && region.strand == Strand::kReverse;
    const std::int64_t anchor = anchor_for(region, opts);

    // On the minus strand "upstream" points at higher coordinates, so the
    // window is mirrored before it is laid out in genomic order.
    const std::int64_t left = reverse ? opts.downstream : opts.upstream;
    const std::int64_t right = reverse ? opts.upstream : opts.downstream;

    // Saturating rather than wrapping arithmetic.
    //
    // A BED line may carry a coordinate near the int64 limit -- parsing accepts
    // anything that fits the type -- and `anchor + right` then overflows, which
    // is undefined behaviour for a signed type (confirmed by UBSan:
    // "9223372036854775806 + 2000 cannot be represented in type 'int64_t'").
    // The off-contig rejection below cannot help, because the overflow happens
    // first and leaves it comparing a wrapped negative value. Saturating keeps
    // the value ordered so that check still does its job.
    // Checked arithmetic rather than a hand-rolled bound: `left` and `right`
    // may be negative (nothing stops --upstream -100), and a comparison like
    // `anchor > kMax - right` itself overflows when right < 0. The builtins are
    // sign-agnostic and compile to a flags check.
    constexpr std::int64_t kMax = std::numeric_limits<std::int64_t>::max();
    constexpr std::int64_t kMin = std::numeric_limits<std::int64_t>::min();

    std::int64_t win_start = 0;
    std::int64_t win_end = 0;
    if (__builtin_sub_overflow(anchor, left, &win_start)) {
        win_start = (left > 0) ? kMin : kMax;
    }
    if (__builtin_add_overflow(anchor, right, &win_end)) {
        win_end = (right > 0) ? kMax : kMin;
    }

    out_window.chrom = region.chrom;
    out_window.start = win_start;
    out_window.end = win_end;
    out_window.name = region.name;
    out_window.strand = region.strand;

    if (chrom_length > 0 &&
        (out_window.end <= 0 || out_window.start >= chrom_length)) {
        return false;  // window lies entirely off the contig
    }
    return out_window.end > out_window.start;
}

double normalization_factor(Normalization norm, std::uint64_t total_mapped,
                            std::int64_t bin_size, double total_signal) {
    switch (norm) {
        case Normalization::kRaw:
            return 1.0;
        case Normalization::kCPM:
            return total_mapped == 0
                       ? 1.0
                       : 1.0e6 / static_cast<double>(total_mapped);
        case Normalization::kRPKM:
            if (total_mapped == 0 || bin_size <= 0) return 1.0;
            return 1.0e9 / (static_cast<double>(total_mapped) *
                            static_cast<double>(bin_size));
        case Normalization::kBPM:
            // Transcript-per-million analogue: each bin's per-bp rate divided
            // by the summed rate over every bin in the requested region set.
            return total_signal <= 0.0 ? 1.0 : 1.0e6 / total_signal;
    }
    return 1.0;
}

SignalMatrix compute_matrix(const BamReader& reader,
                            const std::vector<Region>& regions,
                            const ProfileOptions& opts,
                            const FilterOptions& filters,
                            const std::function<void(double)>& progress) {
    const auto t0 = std::chrono::steady_clock::now();

    // Defence in depth: callers are expected to have run opts.validate(), but
    // this is the function that performs the up-front regions x bins allocation,
    // so it re-checks rather than trusting its caller (REVIEW_2026-08-15
    // finding 2). int64 throughout -- bin_count() no longer narrows to int.
    opts.validate();
    // Same defence for the filters, which had no validation at all.
    filters.validate();
    const std::int64_t nbins = opts.bin_count();
    if (nbins <= 0) {
        throw std::invalid_argument(
            "profile window is empty: check --upstream/--downstream/--bin-size");
    }

    SignalMatrix matrix;
    matrix.rows = regions.size();
    matrix.cols = static_cast<std::size_t>(nbins);
    matrix.values.assign(matrix.rows * matrix.cols, 0.0);
    matrix.row_names.resize(matrix.rows);
    matrix.total_mapped_reads = reader.total_mapped_reads();

    matrix.bin_offsets.resize(static_cast<std::size_t>(nbins));
    const std::int64_t width = opts.window_length();
    for (std::int64_t b = 0; b < nbins; ++b) {
        const std::int64_t lo = (width * b) / nbins;
        const std::int64_t hi = (width * (b + 1)) / nbins;
        matrix.bin_offsets[static_cast<std::size_t>(b)] =
            -opts.upstream + (lo + hi) / 2;
    }

    if (regions.empty()) {
        matrix.column_mean.assign(static_cast<std::size_t>(nbins), 0.0);
        return matrix;
    }

    int nthreads = opts.threads > 0
                       ? opts.threads
                       : static_cast<int>(std::thread::hardware_concurrency());
    nthreads = std::clamp(nthreads, 1, 256);
    nthreads = std::min<int>(nthreads, static_cast<int>(regions.size()));

    std::atomic<std::uint64_t> skipped{0};
    // One byte per region, written only by the worker that owns the row, so no
    // two threads touch the same element. uint8_t rather than bool so the
    // vector is not the bit-packed specialisation.
    std::vector<std::uint8_t> row_skipped(regions.size(), 0);
    std::atomic<std::size_t> completed{0};
    // Records the first exception AND tells every sibling to stop. Previously
    // the exception was captured but nothing was signalled, so the remaining
    // workers ran their stripes to completion building a matrix already
    // guaranteed to be discarded -- on a large region set, the entire run.
    WorkerAbort abort;

    // Block partitioning (rather than round-robin) keeps each worker on a
    // contiguous stripe of the matrix: sequential writes, and at most one
    // shared cache line per thread boundary instead of one per region.
    const std::size_t total = regions.size();
    const std::size_t chunk = (total + static_cast<std::size_t>(nthreads) - 1) /
                              static_cast<std::size_t>(nthreads);

    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(nthreads));

    for (int t = 0; t < nthreads; ++t) {
        const std::size_t begin = static_cast<std::size_t>(t) * chunk;
        if (begin >= total) break;
        const std::size_t end = std::min(begin + chunk, total);

        workers.emplace_back([&, begin, end] {
            // Inside the try, not above it. `tallies` is an nbins-element
            // allocation driven by caller-supplied options -- up to
            // kMaxBins == 10M BinTally -- so it can throw std::bad_alloc, and
            // an exception escaping a thread's entry point is std::terminate,
            // not a catchable error. This is the same defect class as the
            // OpenMP structured-block escape in scrna_matrix's k-NN path
            // (REVIEW_2026-08-15 finding 4): the allocation was one line
            // outside the handler that exists to catch it.
            try {
                std::vector<BinTally> tallies(static_cast<std::size_t>(nbins));
                for (std::size_t i = begin; i < end; ++i) {
                    // Polled once per region, not per bin: one relaxed load
                    // against a whole BAM query is unmeasurable, and the region
                    // is the unit of work whose result would be thrown away.
                    if (abort.aborted()) break;
                    const Region& region = regions[i];
                    double* row = matrix.values.data() + i * matrix.cols;

                    // Named before the skip test. A skipped row used to keep an
                    // empty name, so the matrix TSV carried rows with a blank
                    // first field that nobody could map back to a region.
                    matrix.row_names[i] =
                        region.name.empty()
                            ? region.chrom + ':' + std::to_string(region.start) +
                                  '-' + std::to_string(region.end)
                            : region.name;

                    const std::int64_t chrom_len =
                        reader.chromosome_length(region.chrom);
                    Region window;
                    if (chrom_len < 0 ||
                        !make_window(region, opts, chrom_len, window)) {
                        skipped.fetch_add(1, std::memory_order_relaxed);
                        row_skipped[i] = 1;
                        std::fill_n(row, matrix.cols, opts.missing_value);
                        completed.fetch_add(1, std::memory_order_relaxed);
                        continue;
                    }

                    std::fill(tallies.begin(), tallies.end(), BinTally{});
                    reader.query_bins(window, tallies, filters);

                    const std::int64_t win_width = window.end - window.start;
                    for (std::int64_t b = 0; b < nbins; ++b) {
                        const std::int64_t b_lo =
                            window.start + (win_width * b) / nbins;
                        const std::int64_t b_hi =
                            window.start + (win_width * (b + 1)) / nbins;
                        const std::int64_t bin_len = std::max<std::int64_t>(1, b_hi - b_lo);

                        double value;
                        if (b_hi <= 0 || b_lo >= chrom_len) {
                            value = opts.missing_value;  // off the contig
                        } else if (opts.count_mode == CountMode::kDepth) {
                            value = tallies[static_cast<std::size_t>(b)].base_count /
                                    static_cast<double>(bin_len);
                        } else {
                            value = tallies[static_cast<std::size_t>(b)].read_count;
                        }
                        row[b] = value;
                    }

                    // Lay the row out 5'->3' so that minus-strand features line
                    // up with plus-strand ones.
                    if (opts.respect_strand && region.strand == Strand::kReverse) {
                        std::reverse(row, row + matrix.cols);
                    }

                    completed.fetch_add(1, std::memory_order_relaxed);
                }
            } catch (...) {
                abort.fail(std::current_exception());
            }
        });
    }

    if (progress) {
        while (completed.load(std::memory_order_relaxed) < total) {
            // Checked before sleeping as well as after: once the workers are
            // unwinding, `completed` stops advancing, so a poll loop that only
            // tested the counter would sleep out its full interval on a run
            // that has already failed.
            if (abort.aborted()) break;
            progress(static_cast<double>(completed.load(std::memory_order_relaxed)) /
                     static_cast<double>(total));
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }

    for (auto& worker : workers) worker.join();
    // After the join, on the caller's thread: rethrowing from inside a worker
    // would escape the thread and terminate.
    abort.rethrow_if_failed();
    if (progress) progress(1.0);

    matrix.regions_skipped = skipped.load();

    // A SKIPPED region -- its contig is not in the BAM, or its window lies
    // wholly off the contig -- has no signal to report, so it is left out of
    // everything computed across regions: BPM's total, scaling, and the column
    // mean. Before 2026-09-11 its row of missing values was averaged in, so a
    // region set half on an absent contig produced a meta-profile exactly half
    // as tall (demo data: 5.000363 -> 2.500181), exit 0. deepTools drops such
    // regions. The row itself stays in the matrix, named, so a consumer can
    // still see which regions were skipped.
    //
    // Bins that fall off the contig inside a KEPT region are unchanged: they
    // hold the missing value and count as such, which is the parity the
    // validation suite pins with deepTools' --missingDataAsZero.
    const auto kept_row = [&row_skipped](std::size_t r) { return row_skipped[r] == 0; };

    // BPM needs the summed per-bp rate before it can scale anything.
    double total_signal = 0.0;
    if (opts.normalization == Normalization::kBPM) {
        for (std::size_t r = 0; r < matrix.rows; ++r) {
            if (!kept_row(r)) continue;
            const double* row = matrix.values.data() + r * matrix.cols;
            for (std::size_t c = 0; c < matrix.cols; ++c) {
                total_signal += row[c] / static_cast<double>(opts.bin_size);
            }
        }
    }
    matrix.scale_factor = normalization_factor(
        opts.normalization, matrix.total_mapped_reads, opts.bin_size,
        total_signal);

    if (matrix.scale_factor != 1.0) {
        for (std::size_t r = 0; r < matrix.rows; ++r) {
            if (!kept_row(r)) continue;  // a missing value is not signal to scale
            double* row = matrix.values.data() + r * matrix.cols;
            for (std::size_t c = 0; c < matrix.cols; ++c) row[c] *= matrix.scale_factor;
        }
    }

    // Column means, reduced through per-thread cache-aligned accumulators.
    matrix.column_mean.assign(matrix.cols, 0.0);
    {
        // Was a bare `std::min<int>(nthreads, 8)` with no stated reason, which
        // capped a 64-core node's reduction at 8 threads. Removing the 8
        // outright is worse, though, and the reason is not obvious: each
        // reducer owns a PRIVATE cols-wide accumulator, so scratch memory is
        // `reduce_threads * cols * sizeof(double)`. At the 10M-bin ceiling one
        // accumulator is 80 MB, so 64 reducers would allocate 5 GB of scratch
        // to sum a matrix -- the thread count has to be bounded by the matrix
        // WIDTH, in the opposite direction to the intuition that a wider matrix
        // deserves more threads.
        //
        // Three bounds, all named:
        //   * the worker count the caller asked for;
        //   * enough rows per reducer to pay for starting a thread;
        //   * a fixed scratch budget, which is what the 8 was standing in for.
        constexpr std::size_t kMinRowsPerReducer = 64;
        constexpr std::size_t kReducerScratchBudget = 64u * 1024u * 1024u;

        const std::size_t row_bound = std::max<std::size_t>(1, total / kMinRowsPerReducer);
        const std::size_t bytes_per_reducer =
            std::max<std::size_t>(1, matrix.cols * sizeof(double));
        const std::size_t memory_bound =
            std::max<std::size_t>(1, kReducerScratchBudget / bytes_per_reducer);

        const int reduce_threads = static_cast<int>(std::max<std::size_t>(
            1, std::min({static_cast<std::size_t>(nthreads), row_bound, memory_bound})));
        std::vector<std::unique_ptr<AlignedBuffer>> partials;
        partials.reserve(static_cast<std::size_t>(reduce_threads));
        for (int t = 0; t < reduce_threads; ++t) {
            partials.push_back(std::make_unique<AlignedBuffer>(matrix.cols));
        }

        std::vector<std::thread> reducers;
        const std::size_t rchunk =
            (total + static_cast<std::size_t>(reduce_threads) - 1) /
            static_cast<std::size_t>(reduce_threads);
        for (int t = 0; t < reduce_threads; ++t) {
            const std::size_t begin = static_cast<std::size_t>(t) * rchunk;
            if (begin >= total) break;
            const std::size_t end = std::min(begin + rchunk, total);
            double* acc = partials[static_cast<std::size_t>(t)]->data();
            reducers.emplace_back([&matrix, &row_skipped, acc, begin, end] {
                for (std::size_t i = begin; i < end; ++i) {
                    if (row_skipped[i] != 0) continue;
                    const double* row = matrix.values.data() + i * matrix.cols;
                    for (std::size_t c = 0; c < matrix.cols; ++c) acc[c] += row[c];
                }
            });
        }
        for (auto& reducer : reducers) reducer.join();

        for (const auto& partial : partials) {
            for (std::size_t c = 0; c < matrix.cols; ++c) {
                matrix.column_mean[c] += partial->data()[c];
            }
        }
        const std::size_t kept = total - static_cast<std::size_t>(matrix.regions_skipped);
        if (kept == 0) {
            // No region contributed. The mean is undefined, so it is reported
            // as the missing value rather than as a flat line of zero signal;
            // the CLI refuses this case outright.
            matrix.column_mean.assign(matrix.cols, opts.missing_value);
        } else {
            const double denom = static_cast<double>(kept);
            for (double& v : matrix.column_mean) v /= denom;
        }
    }

    matrix.elapsed_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
            .count();
    return matrix;
}

void write_matrix_tsv(const SignalMatrix& matrix, const std::string& path) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write matrix to " + path);

    out << "region";
    for (std::int64_t offset : matrix.bin_offsets) out << '\t' << offset;
    out << '\n';

    // The per-value formatting below (append_fixed6, tsv_format.hpp) is the
    // whole reason this function is not still `out << matrix.at(r, c)`: at
    // the README's benchmark geometry (100k regions x 80 bins) that ostream
    // insertion was measured at 166-209 ns per value and roughly two thirds
    // of this command's total wall time, entirely inside the region the
    // README's headline number is quoted from. A row is built into `line`
    // and flushed in ~1 MiB chunks rather than written value by value, so the
    // write() syscall count scales with output size, not with cell count.
    std::string line;
    line.reserve(matrix.cols * 10 + 64);  // ~10 bytes/field at fixed precision 6
    std::string chunk;
    chunk.reserve(1u << 20);
    for (std::size_t r = 0; r < matrix.rows; ++r) {
        line.clear();
        line += matrix.row_names[r];
        const double* row = matrix.values.data() + r * matrix.cols;
        for (std::size_t c = 0; c < matrix.cols; ++c) {
            line += '\t';
            append_fixed6(line, row[c]);
        }
        line += '\n';
        chunk += line;
        if (chunk.size() >= (1u << 20)) {
            out.write(chunk.data(), static_cast<std::streamsize>(chunk.size()));
            chunk.clear();
        }
    }
    if (!chunk.empty()) out.write(chunk.data(), static_cast<std::streamsize>(chunk.size()));

    // An ofstream that hits ENOSPC or EFBIG sets badbit and nothing else; the
    // run used to exit 0 over a matrix cut off mid-row. close() is where the
    // buffered tail is written, so the check comes after it.
    out.close();
    if (!out) throw std::runtime_error("writing matrix to " + path + " failed");
}

void write_profile_tsv(const SignalMatrix& matrix, const std::string& path) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write profile to " + path);
    out << "offset\tmean_signal\n";
    out.setf(std::ios::fixed);
    out.precision(6);
    for (std::size_t c = 0; c < matrix.cols; ++c) {
        out << matrix.bin_offsets[c] << '\t' << matrix.column_mean[c] << '\n';
    }
    out.close();
    if (!out) throw std::runtime_error("writing profile to " + path + " failed");
}

}  // namespace profiler
