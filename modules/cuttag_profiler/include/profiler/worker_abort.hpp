// SPDX-License-Identifier: MIT
//
// First-error capture plus a cancellation signal for a fan-out of worker
// threads.
//
// The pattern this replaces stored the first exception under a mutex and then
// joined every worker unconditionally, so a failure was *recorded* promptly and
// *acted on* only after every sibling had finished its entire share of the
// work. On a 100k-region profile that is the whole run: the caller waits for a
// result that is already guaranteed to be thrown away, and every worker keeps
// holding its BAM handle and its per-thread tally buffer while doing it.
//
// Two pieces, deliberately separate:
//
//   * `aborted()` reads a relaxed atomic and is what the hot loop polls. It
//     must be cheap enough to check once per work item without measuring.
//   * `fail()` / `rethrow_if_failed()` carry the exception itself under a
//     mutex, which costs nothing because it happens at most once per worker.
//
// The flag is stored BEFORE the lock is taken in fail(), so a sibling spinning
// on aborted() is not waiting behind the mutex to learn it should stop.
//
// Relaxed ordering is sufficient throughout. The flag carries no data -- it
// gates a `break`, and the exception it refers to is published through the
// mutex, which supplies the happens-before edge for the exception_ptr itself.
// A worker that observes the store one iteration late does one extra iteration,
// which is the same outcome as having polled one iteration earlier.
#pragma once

#include <atomic>
#include <exception>
#include <mutex>
#include <utility>

namespace profiler {

class WorkerAbort {
public:
    WorkerAbort() = default;
    WorkerAbort(const WorkerAbort&) = delete;
    WorkerAbort& operator=(const WorkerAbort&) = delete;

    // Polled by workers between work items. Cheap by construction.
    [[nodiscard]] bool aborted() const noexcept {
        return aborted_.load(std::memory_order_relaxed);
    }

    // Records `error` if it is the first, and signals every sibling to stop.
    //
    // Safe to call from several threads at once and safe to call more than
    // once: the first exception wins, which is the one the user is shown. A
    // later failure is almost always a consequence of the first (a cancelled
    // half-finished run producing its own errors), so reporting the earliest is
    // the more useful of the two.
    void fail(std::exception_ptr error) noexcept {
        // Signal first: siblings poll the flag, not the mutex.
        aborted_.store(true, std::memory_order_relaxed);
        std::lock_guard<std::mutex> lock(mutex_);
        if (!first_error_) first_error_ = std::move(error);
    }

    // Cancels without recording an error -- a caller-requested stop rather
    // than a failure. rethrow_if_failed() stays a no-op afterwards, so a
    // cancelled run is distinguishable from a failed one.
    void cancel() noexcept { aborted_.store(true, std::memory_order_relaxed); }

    [[nodiscard]] bool has_error() const noexcept {
        std::lock_guard<std::mutex> lock(mutex_);
        return static_cast<bool>(first_error_);
    }

    // Rethrows the first recorded exception, if any. Call it after joining, on
    // the thread the caller expects to catch on -- never from inside a worker,
    // where the throw would escape the thread and reach std::terminate.
    void rethrow_if_failed() const {
        std::exception_ptr error;
        {
            std::lock_guard<std::mutex> lock(mutex_);
            error = first_error_;
        }
        if (error) std::rethrow_exception(error);
    }

private:
    mutable std::mutex mutex_;
    std::atomic<bool> aborted_{false};
    std::exception_ptr first_error_;
};

}  // namespace profiler
