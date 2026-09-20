// SPDX-License-Identifier: MIT
//
// Token-addressed cache of recently computed matrices, bounded by BYTES rather
// than by entry count.
//
// `/api/profile` returns a token; `/api/matrix` redeems it for the full
// regions x bins TSV. Something has to hold the matrix between those two
// requests, and the previous shape was a `std::deque` capped at four entries.
//
// Four entries is not a memory bound. A `SignalMatrix` is
// `rows * cols * sizeof(double)` of dense values, and while
// `ProfileOptions::kMaxBins` caps `cols` at 10,000,000, nothing caps `rows` --
// the region count comes from the caller's BED file and no validation limits
// it. So a single cached matrix is already unbounded, and the count-based cap
// turned that from a per-request peak into a standing, cumulative one: up to
// four of them pinned for the lifetime of the process, none released when the
// response that produced them was sent.
//
// This bounds the thing that actually matters. Entry count is still capped as
// well, because holding 10,000 tiny matrices is its own kind of leak and
// because `recall` is a linear constant-time-compare scan whose cost should
// stay predictable.
//
// NOT an LRU. Redemption is a single `/api/matrix` GET moments after the
// `/api/profile` POST that minted the token, so recency of *use* and recency of
// *insertion* are the same ordering in every real interaction. FIFO gets that
// ordering without a touch-on-read step, which would need a write lock on the
// read path.
#pragma once

#include <cstddef>
#include <cstdint>
#include <deque>
#include <memory>
#include <mutex>
#include <string>
#include <utility>

#include "profiler/secrets.hpp"
#include "profiler/types.hpp"

namespace profiler {

// Total resident bytes across all cached matrices.
//
// 512 MiB, not a proportion of system memory: this is a single-operator
// workstation dashboard that shares the machine with the aligner or the browser
// the operator is reading results in, and a cache sized from total RAM competes
// with them. For scale, the benchmarked 100k-window profile is ~89 MB peak RSS
// for the whole process, so this holds several ordinary sessions and refuses to
// hold a pathological one.
inline constexpr std::size_t kDefaultMatrixCacheBytes = 512u * 1024u * 1024u;

// Entry ceiling, retained from the original design.
inline constexpr std::size_t kDefaultMatrixCacheEntries = 4;

class MatrixCache {
public:
    explicit MatrixCache(std::size_t max_bytes = kDefaultMatrixCacheBytes,
                         std::size_t max_entries = kDefaultMatrixCacheEntries)
        : max_bytes_(max_bytes), max_entries_(max_entries) {}

    MatrixCache(const MatrixCache&) = delete;
    MatrixCache& operator=(const MatrixCache&) = delete;

    // Resident cost of one matrix.
    //
    // The dense `values` array dominates, but the row names are counted too:
    // one std::string per region, and a 100k-region BED with long names is tens
    // of megabytes that the values-only estimate would not see. Saturating
    // arithmetic throughout -- this number decides an eviction, and a wrapped
    // product would decide not to evict precisely when the matrix is enormous.
    [[nodiscard]] static std::size_t matrix_bytes(const SignalMatrix& matrix) noexcept {
        std::size_t total = sizeof(SignalMatrix);
        total = saturating_add(total,
                               saturating_mul(saturating_mul(matrix.rows, matrix.cols),
                                              sizeof(double)));
        total = saturating_add(
            total, saturating_mul(matrix.column_mean.size(), sizeof(double)));
        total = saturating_add(
            total, saturating_mul(matrix.bin_offsets.size(), sizeof(std::int64_t)));
        for (const std::string& name : matrix.row_names) {
            total = saturating_add(total, sizeof(std::string));
            total = saturating_add(total, name.capacity());
        }
        return total;
    }

    // Stores `matrix` and returns the token that redeems it.
    [[nodiscard]] std::string remember(std::shared_ptr<SignalMatrix> matrix) {
        const std::size_t bytes = matrix ? matrix_bytes(*matrix) : 0;
        std::string token = secure_token(kTokenBytes);

        std::lock_guard<std::mutex> lock(mutex_);
        entries_.push_back(Entry{token, std::move(matrix), bytes});
        bytes_ = saturating_add(bytes_, bytes);

        // `size() > 1` on both loops, never `!entries_.empty()`: the entry just
        // inserted is the one whose token has already been handed to the
        // caller, so evicting it would return a token that 404s immediately.
        // A single matrix larger than the whole budget is therefore kept --
        // deliberately. The cache cannot bound a peak that has already been
        // allocated by compute_matrix; what it can do is refuse to let that
        // peak become permanent, and dropping every *other* entry is exactly
        // that.
        while (bytes_ > max_bytes_ && entries_.size() > 1) evict_oldest();
        while (entries_.size() > max_entries_ && entries_.size() > 1) evict_oldest();
        return token;
    }

    // Redeems a token. Returns nullptr when it is unknown or has been evicted.
    [[nodiscard]] std::shared_ptr<SignalMatrix> recall(const std::string& token) const {
        std::lock_guard<std::mutex> lock(mutex_);
        // An empty token must never match, including against an empty cache or
        // an entry that somehow lost its key.
        if (token.empty()) return nullptr;
        // Constant-time, and every entry is examined rather than returning at
        // the first hit: `==` leaks the length of a correct prefix through
        // response timing, and an early return leaks the matching entry's
        // position in the cache. /api/matrix streams the whole matrix, so this
        // token is the entire access control on it (SECURITY_HTTP L9).
        std::shared_ptr<SignalMatrix> found;
        for (const Entry& entry : entries_) {
            if (constant_time_equals(entry.token, token)) found = entry.matrix;
        }
        return found;
    }

    [[nodiscard]] std::size_t size() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return entries_.size();
    }

    // Resident bytes currently accounted for. Exposed so the eviction policy is
    // observable in tests rather than only inferable from behaviour.
    [[nodiscard]] std::size_t bytes() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return bytes_;
    }

    [[nodiscard]] std::size_t max_bytes() const noexcept { return max_bytes_; }

private:
    // 256 bits, 64 hex characters. Matches the session token width.
    static constexpr std::size_t kTokenBytes = 32;

    struct Entry {
        std::string token;
        std::shared_ptr<SignalMatrix> matrix;
        std::size_t bytes = 0;
    };

    static std::size_t saturating_add(std::size_t a, std::size_t b) noexcept {
        return (a > SIZE_MAX - b) ? SIZE_MAX : a + b;
    }
    static std::size_t saturating_mul(std::size_t a, std::size_t b) noexcept {
        if (a == 0 || b == 0) return 0;
        return (a > SIZE_MAX / b) ? SIZE_MAX : a * b;
    }

    // Caller holds mutex_.
    void evict_oldest() noexcept {
        // Subtract the size recorded at insertion, never a freshly recomputed
        // one: `row_names` capacities can differ between the two calls, and an
        // accounting drift here would accumulate silently until the budget
        // stopped meaning anything.
        bytes_ = (bytes_ >= entries_.front().bytes) ? bytes_ - entries_.front().bytes : 0;
        entries_.pop_front();
    }

    mutable std::mutex mutex_;
    std::deque<Entry> entries_;
    std::size_t bytes_ = 0;
    const std::size_t max_bytes_;
    const std::size_t max_entries_;
};

}  // namespace profiler
