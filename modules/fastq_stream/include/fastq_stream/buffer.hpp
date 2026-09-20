// SPDX-License-Identifier: MIT
// fastq_stream — cache-aligned recycled buffers.
//
// Chunks move between pipeline stages by pointer only; the payload bytes are
// never copied. Each allocation carries kPrefix bytes of headroom in front of
// data() so the record assembler can prepend a partial record from the
// previous chunk in place instead of memmoving the whole payload.
//
// Buffers are recycled through a free list behind a spinlock. That is not on
// the hot path: one acquire/release pair per kChunkBytes (256 KiB) of stream,
// i.e. ~200k lock operations for a 50 GB dataset, against hundreds of millions
// of lock-free ring transfers.
#pragma once

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <mutex>
#include <utility>
#include <vector>

#include "fastq_stream/lockfree_queue.hpp"

namespace fq {

inline constexpr std::size_t kChunkBytes = 256 * 1024;
// Headroom in front of data() for splicing a record that straddled the
// previous chunk. This also bounds the largest single FASTQ record the
// assembler can carry: 64 KiB, i.e. reads up to ~32 kbp. Ample for Illumina;
// ONT/PacBio inputs beyond that are rejected with an explicit error rather
// than silently corrupted.
inline constexpr std::size_t kPrefix = 64 * 1024;
inline constexpr std::size_t kAlign = 64;

class BufferPool;

// Owning, movable byte buffer. Returns itself to its pool on destruction.
class Buffer {
 public:
  Buffer() = default;
  Buffer(BufferPool* pool, std::byte* base, std::size_t cap) noexcept
      : pool_(pool), base_(base), cap_(cap) {}

  Buffer(Buffer&& other) noexcept { swap(other); }
  Buffer& operator=(Buffer&& other) noexcept {
    if (this != &other) {
      Buffer tmp(std::move(other));
      swap(tmp);
    }
    return *this;
  }
  Buffer(const Buffer&) = delete;
  Buffer& operator=(const Buffer&) = delete;
  ~Buffer();

  explicit operator bool() const noexcept { return base_ != nullptr; }

  std::byte* data() noexcept { return base_ + kPrefix + front_; }
  const std::byte* data() const noexcept { return base_ + kPrefix + front_; }
  char* chars() noexcept { return reinterpret_cast<char*>(data()); }
  const char* chars() const noexcept { return reinterpret_cast<const char*>(data()); }
  const uint8_t* bytes() const noexcept {
    return reinterpret_cast<const uint8_t*>(data());
  }

  std::size_t size() const noexcept { return size_; }
  void set_size(std::size_t n) noexcept { size_ = n; }
  std::size_t capacity() const noexcept { return cap_ - kPrefix - front_; }

  // Move the logical start back by n bytes, exposing headroom to be filled.
  // Used to splice a carried partial record onto the front of a chunk.
  std::byte* extend_front(std::size_t n) noexcept {
    front_ -= static_cast<std::ptrdiff_t>(n);
    size_ += n;
    return data();
  }

  // Sequence number: restores stream order after parallel stages.
  uint64_t seq = 0;
  // Set on the final chunk so downstream stages know the stream ended.
  bool last = false;

 private:
  void swap(Buffer& o) noexcept {
    std::swap(pool_, o.pool_);
    std::swap(base_, o.base_);
    std::swap(cap_, o.cap_);
    std::swap(size_, o.size_);
    std::swap(front_, o.front_);
    std::swap(seq, o.seq);
    std::swap(last, o.last);
  }

  BufferPool* pool_ = nullptr;
  std::byte* base_ = nullptr;
  std::size_t cap_ = 0;
  std::size_t size_ = 0;
  std::ptrdiff_t front_ = 0;
};

class BufferPool {
 public:
  explicit BufferPool(std::size_t chunk_bytes = kChunkBytes)
      : cap_(kPrefix + chunk_bytes) {}

  ~BufferPool() {
    for (std::byte* p : free_) std::free(p);
  }

  BufferPool(const BufferPool&) = delete;
  BufferPool& operator=(const BufferPool&) = delete;

  Buffer acquire() {
    std::byte* p = nullptr;
    {
      std::lock_guard<std::mutex> lk(mu_);
      if (!free_.empty()) {
        p = free_.back();
        free_.pop_back();
      }
    }
    if (p == nullptr) {
      p = static_cast<std::byte*>(std::aligned_alloc(kAlign, round_up(cap_, kAlign)));
      if (p == nullptr) throw std::bad_alloc();
      allocated_.fetch_add(1, std::memory_order_relaxed);
    }
    return Buffer(this, p, cap_);
  }

  void release(std::byte* p) noexcept {
    std::lock_guard<std::mutex> lk(mu_);
    free_.push_back(p);
  }

  // Number of distinct blocks ever allocated — the pool's true high-water
  // mark, which is what the RSS ceiling is argued from.
  std::size_t blocks_allocated() const noexcept {
    return allocated_.load(std::memory_order_relaxed);
  }
  std::size_t bytes_allocated() const noexcept { return blocks_allocated() * cap_; }

 private:
  static std::size_t round_up(std::size_t v, std::size_t a) noexcept {
    return (v + a - 1) / a * a;
  }

  const std::size_t cap_;
  std::mutex mu_;
  std::vector<std::byte*> free_;
  std::atomic<std::size_t> allocated_{0};
};

inline Buffer::~Buffer() {
  if (pool_ != nullptr && base_ != nullptr) pool_->release(base_);
}

}  // namespace fq
