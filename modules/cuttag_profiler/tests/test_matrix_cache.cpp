// SPDX-License-Identifier: MIT
//
// Regressions for the matrix cache's memory bound.
//
// The cache used to be a deque capped at FOUR ENTRIES, which is not a memory
// bound: `ProfileOptions::kMaxBins` caps a matrix's width, but nothing caps the
// number of regions a caller's BED file declares, so one cached matrix is
// already unbounded and four of them were pinned for the process's lifetime.
// These tests pin the property that replaced it -- eviction is decided by
// resident bytes -- and the two edge cases where a naive byte cap goes wrong.
#include "test_util.hpp"

#include <memory>
#include <string>
#include <vector>

#include "profiler/matrix_cache.hpp"

using namespace profiler;

namespace {

// A matrix of a requested shape, filled so `values` is genuinely resident.
std::shared_ptr<SignalMatrix> make_matrix(std::size_t rows, std::size_t cols) {
    auto matrix = std::make_shared<SignalMatrix>();
    matrix->rows = rows;
    matrix->cols = cols;
    matrix->values.assign(rows * cols, 1.0);
    matrix->row_names.assign(rows, "region");
    matrix->column_mean.assign(cols, 0.0);
    matrix->bin_offsets.assign(cols, 0);
    return matrix;
}

}  // namespace

int main() {
    testing::Suite suite{"matrix_cache", {}};

    suite.add("eviction is driven by bytes, not by entry count", [] {
        // Room for two of these and no more. Under the old count-based policy
        // all four would be resident, because four is four regardless of size.
        const std::size_t one = MatrixCache::matrix_bytes(*make_matrix(500, 100));
        MatrixCache cache(one * 2 + one / 2, /*max_entries=*/16);

        std::vector<std::string> tokens;
        for (int i = 0; i < 4; ++i) {
            tokens.push_back(cache.remember(make_matrix(500, 100)));
        }

        CHECK(cache.size() <= 2);
        CHECK(cache.bytes() <= cache.max_bytes());

        // FIFO: the most recent survives, the oldest is gone.
        CHECK(cache.recall(tokens.back()) != nullptr);
        CHECK(cache.recall(tokens.front()) == nullptr);
    });

    suite.add("a matrix larger than the whole budget is still redeemable", [] {
        // The token has already been handed to the caller by the time eviction
        // runs, so evicting the entry that just arrived would return a token
        // that 404s immediately -- a worse bug than the one being fixed. The
        // cache cannot un-allocate a peak compute_matrix already reached; what
        // it must do is refuse to make that peak permanent.
        MatrixCache cache(/*max_bytes=*/1024, /*max_entries=*/4);

        const std::string small = cache.remember(make_matrix(2, 2));
        const std::string huge = cache.remember(make_matrix(2000, 50));

        CHECK(cache.recall(huge) != nullptr);
        CHECK_EQ(cache.size(), std::size_t{1});
        // ...and it displaced everything else rather than being dropped itself.
        CHECK(cache.recall(small) == nullptr);
    });

    suite.add("the entry ceiling still applies under the byte budget", [] {
        // Many tiny matrices never approach the byte budget, so bytes alone
        // would let the cache grow without bound -- and `recall` is a linear
        // constant-time-compare scan whose cost must stay predictable.
        MatrixCache cache(/*max_bytes=*/SIZE_MAX, /*max_entries=*/3);

        std::vector<std::string> tokens;
        for (int i = 0; i < 10; ++i) tokens.push_back(cache.remember(make_matrix(1, 1)));

        CHECK_EQ(cache.size(), std::size_t{3});
        CHECK(cache.recall(tokens.back()) != nullptr);
        CHECK(cache.recall(tokens[0]) == nullptr);
    });

    suite.add("byte accounting returns to zero as entries are displaced", [] {
        // Eviction subtracts the size recorded at INSERTION. Recomputing it
        // instead would drift as row_name capacities changed, and the drift
        // accumulates in one direction until the budget stops meaning anything.
        MatrixCache cache(/*max_bytes=*/SIZE_MAX, /*max_entries=*/1);

        const std::size_t expected = MatrixCache::matrix_bytes(*make_matrix(10, 10));
        for (int i = 0; i < 25; ++i) {
            (void)cache.remember(make_matrix(10, 10));
            CHECK_EQ(cache.size(), std::size_t{1});
            CHECK_EQ(cache.bytes(), expected);
        }
    });

    suite.add("matrix_bytes saturates rather than wrapping", [] {
        // This number decides an eviction. A wrapped product would report a
        // pathological matrix as tiny and decline to evict anything, precisely
        // when eviction matters most. No allocation here -- the shape is
        // declared, not backed.
        SignalMatrix absurd;
        absurd.rows = SIZE_MAX / 2;
        absurd.cols = 4;
        CHECK_EQ(MatrixCache::matrix_bytes(absurd), SIZE_MAX);

        SignalMatrix empty;
        CHECK(MatrixCache::matrix_bytes(empty) > 0);  // the struct itself
    });

    suite.add("token redemption keeps its L9 properties", [] {
        // Carried over from the cache's previous home in http_server.cpp: the
        // matrix token is the entire access control on /api/matrix, so an empty
        // or wrong token must never redeem, including against an empty cache.
        MatrixCache cache;
        CHECK(cache.recall("") == nullptr);
        CHECK(cache.recall("deadbeef") == nullptr);

        const std::string token = cache.remember(make_matrix(4, 4));
        CHECK_EQ(token.size(), std::size_t{64});  // 256 bits of hex
        CHECK(cache.recall(token) != nullptr);
        CHECK(cache.recall("") == nullptr);
        // A correct prefix must not redeem.
        CHECK(cache.recall(token.substr(0, token.size() - 1)) == nullptr);
    });

    return suite.run();
}
