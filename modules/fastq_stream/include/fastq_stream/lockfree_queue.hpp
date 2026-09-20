// SPDX-License-Identifier: MIT
// fastq_stream — bounded lock-free SPSC ring buffer.
//
// One producer thread, one consumer thread, no mutexes and no CAS on the fast
// path. Head and tail live on separate cache lines, and each side keeps a
// private cached copy of the opposite index so that a steady-state push/pop
// touches the other core's line only when the ring is near-empty or near-full.
#pragma once

#include <atomic>
#include <cstddef>
#include <new>
#include <optional>
#include <type_traits>
#include <utility>

namespace fq {

// Padding used to keep the producer's and consumer's cursors off one another's
// cache line. Getting this wrong does not break correctness, only throughput --
// two atomics sharing a line turn every push/pop into a coherence round trip.
//
// std::hardware_destructive_interference_size is deliberately NOT used, even
// where the feature-test macro advertises it. GCC warns (-Winterference-size)
// because its value is baked in per translation unit from -mtune, so two TUs
// compiled with different flags disagree about the size of any type padded with
// it -- an ODR violation and a silent ABI break in a header-only queue that is
// included from several targets here. The standard's own guidance is to use it
// only in contexts where the ABI cannot leak; a member of a shared template is
// exactly the case where it can.
//
// So: a fixed constant per architecture, chosen to match the real destructive
// interference range rather than the line size alone. Apple silicon and Neoverse
// prefetch in 128-byte pairs; x86-64 is 64 with the same caveat for Sandy Bridge
// onward, where the adjacent-line prefetcher makes 128 the safer figure.
#if defined(__aarch64__) || defined(__arm64__)
inline constexpr std::size_t kCacheLine = 128;  // Apple silicon / Neoverse
#elif defined(__x86_64__) || defined(_M_X64)
inline constexpr std::size_t kCacheLine = 128;  // adjacent-line prefetch pairs
#else
inline constexpr std::size_t kCacheLine = 64;
#endif

// Capacity must be a power of two; one slot is left unused so that
// head == tail unambiguously means "empty".
template <typename T, std::size_t Capacity>
class SpscRing {
  static_assert(Capacity >= 2, "capacity must be at least 2");
  static_assert((Capacity & (Capacity - 1)) == 0, "capacity must be a power of two");
  static_assert(std::is_nothrow_move_constructible_v<T>,
                "T must be nothrow-movable to keep push/pop non-throwing");

 public:
  SpscRing() = default;
  SpscRing(const SpscRing&) = delete;
  SpscRing& operator=(const SpscRing&) = delete;

  // Producer side. Returns false if the ring is full (caller should back off).
  bool try_push(T&& value) noexcept {
    const std::size_t tail = tail_.load(std::memory_order_relaxed);
    const std::size_t next = (tail + 1) & kMask;
    if (next == cached_head_) {
      cached_head_ = head_.load(std::memory_order_acquire);
      if (next == cached_head_) return false;
    }
    slots_[tail] = std::move(value);
    tail_.store(next, std::memory_order_release);
    return true;
  }

  // Consumer side. Returns nullopt if the ring is empty.
  std::optional<T> try_pop() noexcept {
    const std::size_t head = head_.load(std::memory_order_relaxed);
    if (head == cached_tail_) {
      cached_tail_ = tail_.load(std::memory_order_acquire);
      if (head == cached_tail_) return std::nullopt;
    }
    std::optional<T> out(std::move(slots_[head]));
    slots_[head] = T{};
    head_.store((head + 1) & kMask, std::memory_order_release);
    return out;
  }

  bool empty() const noexcept {
    return head_.load(std::memory_order_acquire) == tail_.load(std::memory_order_acquire);
  }

  std::size_t size_approx() const noexcept {
    const std::size_t t = tail_.load(std::memory_order_acquire);
    const std::size_t h = head_.load(std::memory_order_acquire);
    return (t - h) & kMask;
  }

  static constexpr std::size_t capacity() noexcept { return Capacity - 1; }

 private:
  static constexpr std::size_t kMask = Capacity - 1;

  alignas(kCacheLine) std::atomic<std::size_t> head_{0};
  alignas(kCacheLine) std::size_t cached_tail_{0};  // consumer-private
  alignas(kCacheLine) std::atomic<std::size_t> tail_{0};
  alignas(kCacheLine) std::size_t cached_head_{0};  // producer-private
  alignas(kCacheLine) T slots_[Capacity]{};
};

}  // namespace fq
