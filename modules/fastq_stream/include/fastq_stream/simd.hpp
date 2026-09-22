// SPDX-License-Identifier: MIT
// fastq_stream — vectorised primitives for Phred/base kernels.
//
// Runtime dispatch
// ----------------
// These kernels used to be selected at COMPILE time, on __AVX2__/__AVX512BW__.
// The wheel builds with FQ_NATIVE=OFF, so on x86-64
// neither macro was ever defined and every production binary ran the scalar
// loops on hardware with AVX2 to spare (docs/TORTURE_2026-09-10.md F14).
//
// So each ISA's kernels now carry a per-function `target` attribute and are
// compiled into every x86-64 binary, and which set runs is decided once, at the
// first call, from CPUID -- the arrangement modules/scrna_matrix's simd_math.hpp
// and modules/peaks' poisson_model.hpp already use, for the same reason. The
// attribute, not a -mavx2 flag on the target, is what matters: a flag would let
// the compiler emit AVX2 anywhere, including in code no runtime check guards,
// and the portable binary would SIGILL on a CPU without it.
//
// AArch64 needs no probe: Advanced SIMD is mandatory in ARMv8-A.
//
// AVX-512BW kernels are compiled in only with -DFQ_ENABLE_AVX512=ON, off by
// default -- the same default, for the same reason, as PEAKS_ENABLE_AVX512.
// They have not been executed on hardware by this project yet (Rosetta has no
// AVX-512); test_trimmer runs them on any CPU that reports AVX-512BW, and that
// run is what should precede turning the default on.
//
// Exactness
// ---------
// Every kernel here is exact: every ISA must produce bit-identical results to
// the scalar reference. The one intentional exception is `mismatches` past its
// early exit, which only promises "greater than cap" -- vector kernels check
// the cap once per stride. test_trimmer drives every ISA this CPU can execute
// through the *_isa entry points against a scalar reference.
//
// FQ_FORCE_ISA=scalar|neon|avx2|avx512 pins the choice, so one machine can
// exercise every compiled path. A path that was not compiled in is ignored; one
// the hardware cannot execute is honoured and will fault. That is the caller's
// explicit choice, which is why it is opt-in -- the same contract as
// SCRNA_FORCE_ISA and PTO_PEAKS_FORCE_ISA.
#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <initializer_list>

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#  include <immintrin.h>
#  define FQ_SIMD_X86 1
#  if defined(FQ_ENABLE_AVX512)
#    define FQ_SIMD_X86_AVX512 1
#  endif
#elif defined(__aarch64__) || defined(__ARM_NEON)
#  include <arm_neon.h>
#  define FQ_SIMD_NEON 1
#endif

namespace fq::simd {

enum class Isa { kScalar, kNeon, kAvx2, kAvx512 };

// Returned by find_seed4 when no offset in the range matches.
inline constexpr size_t kNotFound = static_cast<size_t>(-1);

// ------------------------------------------------------------------ scalar
// The reference every other kernel is tested against, and the fallback.
namespace scalar {

inline uint64_t sum_u8(const uint8_t* p, size_t n) noexcept {
  uint64_t total = 0;
  for (size_t i = 0; i < n; ++i) total += p[i];
  return total;
}

inline size_t count_ge(const uint8_t* p, size_t n, uint8_t thr) noexcept {
  size_t count = 0;
  for (size_t i = 0; i < n; ++i) count += (p[i] >= thr);
  return count;
}

inline void add_widen_u8_to_u32(const uint8_t* src, size_t n, uint32_t* dst) noexcept {
  for (size_t i = 0; i < n; ++i) dst[i] += src[i];
}

inline size_t mismatches(const uint8_t* a, const uint8_t* b, size_t n, size_t cap) noexcept {
  size_t diff = 0;
  for (size_t i = 0; i < n; ++i) {
    diff += (a[i] != b[i]);
    if (diff > cap) return diff;
  }
  return diff;
}

// Subtracting 33 wraps both out-of-range sides above 93, so the test is one
// unsigned compare per byte.
inline bool all_printable(const uint8_t* p, size_t n) noexcept {
  bool bad = false;
  for (size_t i = 0; i < n; ++i) bad |= static_cast<uint8_t>(p[i] - 33u) > 93u;
  return !bad;
}

// Smallest offset s in [from, to_inclusive] where the 4 bytes at p+s equal
// `seed` (a little-endian word, i.e. what `std::memcpy(&seed, adapter, 4)`
// produces -- the same load AdapterTrimmer::trim already does), or kNotFound.
// The reference every vector kernel below is tested against, and the
// fallback: an all-offsets scan, exactly what the seed loop in
// adapter_trimmer.hpp did per-offset before this kernel existed.
//
// Precondition, unchanged from that loop: reading 4 bytes at `to_inclusive`
// must be in bounds, i.e. `to_inclusive + 4 <= (the buffer's length)`. Every
// vector kernel below relies on it to size its own last full block.
inline size_t find_seed4(const uint8_t* p, size_t from, size_t to_inclusive,
                         uint32_t seed) noexcept {
  if (from > to_inclusive) return kNotFound;
  for (size_t s = from; s <= to_inclusive; ++s) {
    uint32_t w;
    std::memcpy(&w, p + s, 4);
    if (w == seed) return s;
  }
  return kNotFound;
}

}  // namespace scalar

// -------------------------------------------------------------------- NEON
#if defined(FQ_SIMD_NEON)
namespace neon {

inline uint64_t sum_u8(const uint8_t* p, size_t n) noexcept {
  size_t i = 0;
  uint32x4_t acc = vdupq_n_u32(0);
  for (; i + 16 <= n; i += 16) acc = vpadalq_u16(acc, vpaddlq_u8(vld1q_u8(p + i)));
  return vaddvq_u32(acc) + scalar::sum_u8(p + i, n - i);
}

inline size_t count_ge(const uint8_t* p, size_t n, uint8_t thr) noexcept {
  size_t i = 0;
  size_t count = 0;
  const uint8x16_t t = vdupq_n_u8(thr);
  const uint8x16_t ones = vdupq_n_u8(1);
  for (; i + 16 <= n; i += 16) {
    count += vaddvq_u8(vandq_u8(vcgeq_u8(vld1q_u8(p + i), t), ones));  // <=16, no overflow
  }
  return count + scalar::count_ge(p + i, n - i, thr);
}

inline void add_widen_u8_to_u32(const uint8_t* src, size_t n, uint32_t* dst) noexcept {
  size_t i = 0;
  for (; i + 16 <= n; i += 16) {
    const uint8x16_t v = vld1q_u8(src + i);
    const uint16x8_t lo = vmovl_u8(vget_low_u8(v));
    const uint16x8_t hi = vmovl_u8(vget_high_u8(v));
    vst1q_u32(dst + i + 0, vaddq_u32(vld1q_u32(dst + i + 0), vmovl_u16(vget_low_u16(lo))));
    vst1q_u32(dst + i + 4, vaddq_u32(vld1q_u32(dst + i + 4), vmovl_u16(vget_high_u16(lo))));
    vst1q_u32(dst + i + 8, vaddq_u32(vld1q_u32(dst + i + 8), vmovl_u16(vget_low_u16(hi))));
    vst1q_u32(dst + i + 12, vaddq_u32(vld1q_u32(dst + i + 12), vmovl_u16(vget_high_u16(hi))));
  }
  scalar::add_widen_u8_to_u32(src + i, n - i, dst + i);
}

inline size_t mismatches(const uint8_t* a, const uint8_t* b, size_t n, size_t cap) noexcept {
  size_t i = 0;
  size_t diff = 0;
  const uint8x16_t ones = vdupq_n_u8(1);
  for (; i + 16 <= n; i += 16) {
    diff += vaddvq_u8(vbicq_u8(ones, vceqq_u8(vld1q_u8(a + i), vld1q_u8(b + i))));
    if (diff > cap) return diff;
  }
  for (; i < n; ++i) {
    diff += (a[i] != b[i]);
    if (diff > cap) return diff;
  }
  return diff;
}

inline bool all_printable(const uint8_t* p, size_t n) noexcept {
  size_t i = 0;
  const uint8x16_t off = vdupq_n_u8(33);
  const uint8x16_t lim = vdupq_n_u8(93);
  uint8x16_t acc = vdupq_n_u8(0);
  for (; i + 16 <= n; i += 16) acc = vorrq_u8(acc, vcgtq_u8(vsubq_u8(vld1q_u8(p + i), off), lim));
  return vmaxvq_u8(acc) == 0 && scalar::all_printable(p + i, n - i);
}

// Tests 16 candidate offsets at once with four overlapping 16-byte loads, one
// per byte of the seed, ANDed together into a per-lane 0xFF/0x00 verdict.
// NEON has no cheap "extract a bit per lane" instruction the way x86's
// movemask does, so this stays two-tier rather than resolving the exact
// offset from the vector compare directly: vmaxvq_u8 answers "did ANY of
// these 16 match" in one reduction, and only on the rare block where it did
// -- on random DNA, a 4-byte seed matches about 1 in 256 offsets, so ~94% of
// 16-wide blocks answer no here -- does the scalar loop below run, and it
// runs the exact same per-byte check scalar::find_seed4 does, over only the
// 16 offsets this block covers, so the result is identical either way.
inline size_t find_seed4(const uint8_t* p, size_t from, size_t to_inclusive,
                         uint32_t seed) noexcept {
  if (from > to_inclusive) return kNotFound;
  uint8_t sb[4];
  std::memcpy(sb, &seed, 4);
  const uint8x16_t b0 = vdupq_n_u8(sb[0]);
  const uint8x16_t b1 = vdupq_n_u8(sb[1]);
  const uint8x16_t b2 = vdupq_n_u8(sb[2]);
  const uint8x16_t b3 = vdupq_n_u8(sb[3]);
  size_t s = from;
  // The last offset a 16-wide block covers is s+15; it reads its own 4 bytes
  // at p[s+15 .. s+19), so the block is only taken while s+15 <= to_inclusive
  // -- the same bound scalar::find_seed4 applies to a single offset, just
  // checked once for the whole block instead of once per offset in it.
  for (; s + 15 <= to_inclusive; s += 16) {
    const uint8x16_t v0 = vld1q_u8(p + s);
    const uint8x16_t v1 = vld1q_u8(p + s + 1);
    const uint8x16_t v2 = vld1q_u8(p + s + 2);
    const uint8x16_t v3 = vld1q_u8(p + s + 3);
    const uint8x16_t m = vandq_u8(vandq_u8(vceqq_u8(v0, b0), vceqq_u8(v1, b1)),
                                  vandq_u8(vceqq_u8(v2, b2), vceqq_u8(v3, b3)));
    if (vmaxvq_u8(m) == 0) continue;  // usually never entered: no match here
    for (size_t j = 0; j < 16; ++j) {
      uint32_t w;
      std::memcpy(&w, p + s + j, 4);
      if (w == seed) return s + j;
    }
  }
  return scalar::find_seed4(p, s, to_inclusive, seed);
}

}  // namespace neon
#endif  // FQ_SIMD_NEON

// -------------------------------------------------------------------- AVX2
// Integer-only kernels: AVX2 alone, no FMA, so a CPU (or an emulator) that
// reports AVX2 without FMA still gets them.
#if defined(FQ_SIMD_X86)
namespace avx2 {

__attribute__((target("avx2"))) inline uint64_t sum_u8(const uint8_t* p, size_t n) noexcept {
  size_t i = 0;
  const __m256i zero = _mm256_setzero_si256();
  __m256i acc = zero;
  for (; i + 32 <= n; i += 32) {
    const __m256i v = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p + i));
    acc = _mm256_add_epi64(acc, _mm256_sad_epu8(v, zero));
  }
  alignas(32) uint64_t lanes[4];
  _mm256_store_si256(reinterpret_cast<__m256i*>(lanes), acc);
  uint64_t total = lanes[0] + lanes[1] + lanes[2] + lanes[3];
  for (; i < n; ++i) total += p[i];
  return total;
}

__attribute__((target("avx2")))
inline size_t count_ge(const uint8_t* p, size_t n, uint8_t thr) noexcept {
  size_t i = 0;
  const __m256i zero = _mm256_setzero_si256();
  const __m256i ones = _mm256_set1_epi8(1);
  const __m256i t = _mm256_set1_epi8(static_cast<char>(thr));
  __m256i acc = zero;
  for (; i + 32 <= n; i += 32) {
    const __m256i v = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p + i));
    // max(v,t)==v  <=>  v >= t, without the signed-compare trap.
    const __m256i ge = _mm256_cmpeq_epi8(_mm256_max_epu8(v, t), v);
    acc = _mm256_add_epi64(acc, _mm256_sad_epu8(_mm256_and_si256(ge, ones), zero));
  }
  alignas(32) uint64_t lanes[4];
  _mm256_store_si256(reinterpret_cast<__m256i*>(lanes), acc);
  size_t count = static_cast<size_t>(lanes[0] + lanes[1] + lanes[2] + lanes[3]);
  for (; i < n; ++i) count += (p[i] >= thr);
  return count;
}

__attribute__((target("avx2")))
inline void add_widen_u8_to_u32(const uint8_t* src, size_t n, uint32_t* dst) noexcept {
  size_t i = 0;
  for (; i + 8 <= n; i += 8) {
    const __m128i v = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src + i));
    const __m256i w = _mm256_cvtepu8_epi32(v);
    const __m256i d = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(dst + i));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst + i), _mm256_add_epi32(d, w));
  }
  for (; i < n; ++i) dst[i] += src[i];
}

__attribute__((target("avx2")))
inline size_t mismatches(const uint8_t* a, const uint8_t* b, size_t n, size_t cap) noexcept {
  size_t i = 0;
  size_t diff = 0;
  const __m256i zero = _mm256_setzero_si256();
  const __m256i ones = _mm256_set1_epi8(1);
  for (; i + 32 <= n; i += 32) {
    const __m256i x = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(a + i));
    const __m256i y = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(b + i));
    const __m256i ne = _mm256_andnot_si256(_mm256_cmpeq_epi8(x, y), ones);
    alignas(32) uint64_t lanes[4];
    _mm256_store_si256(reinterpret_cast<__m256i*>(lanes), _mm256_sad_epu8(ne, zero));
    diff += static_cast<size_t>(lanes[0] + lanes[1] + lanes[2] + lanes[3]);
    if (diff > cap) return diff;
  }
  for (; i < n; ++i) {
    diff += (a[i] != b[i]);
    if (diff > cap) return diff;
  }
  return diff;
}

__attribute__((target("avx2"))) inline bool all_printable(const uint8_t* p, size_t n) noexcept {
  size_t i = 0;
  const __m256i off = _mm256_set1_epi8(33);
  const __m256i lim = _mm256_set1_epi8(93);
  __m256i acc = _mm256_setzero_si256();
  for (; i + 32 <= n; i += 32) {
    const __m256i v =
        _mm256_sub_epi8(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(p + i)), off);
    // Unsigned v > 93  <=>  max(v, 93) != 93, without the signed-compare trap.
    acc = _mm256_or_si256(acc, _mm256_xor_si256(_mm256_max_epu8(v, lim), lim));
  }
  bool bad = !_mm256_testz_si256(acc, acc);
  for (; i < n; ++i) bad |= static_cast<uint8_t>(p[i] - 33u) > 93u;
  return !bad;
}

// Tests 32 candidate offsets at once with four overlapping 32-byte loads, one
// per byte of the seed. Unlike NEON, x86 gets the exact matching offset
// directly and cheaply: movemask_epi8 turns the per-lane 0xFF/0x00 verdict
// into a 32-bit mask (bit j set iff offset s+j matched all four byte planes),
// and ctz on that mask is the 5'-most match in the block -- exactly what the
// scalar reference's smallest-offset-first scan would have found first,
// without a second, block-local scalar pass to locate it.
__attribute__((target("avx2")))
inline size_t find_seed4(const uint8_t* p, size_t from, size_t to_inclusive,
                         uint32_t seed) noexcept {
  if (from > to_inclusive) return kNotFound;
  uint8_t sb[4];
  std::memcpy(sb, &seed, 4);
  const __m256i b0 = _mm256_set1_epi8(static_cast<char>(sb[0]));
  const __m256i b1 = _mm256_set1_epi8(static_cast<char>(sb[1]));
  const __m256i b2 = _mm256_set1_epi8(static_cast<char>(sb[2]));
  const __m256i b3 = _mm256_set1_epi8(static_cast<char>(sb[3]));
  size_t s = from;
  // The last offset a 32-wide block covers is s+31, reading p[s+31 .. s+35);
  // taken only while that is in bounds, the same precondition
  // scalar::find_seed4 applies per offset.
  for (; s + 31 <= to_inclusive; s += 32) {
    const __m256i v0 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p + s));
    const __m256i v1 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p + s + 1));
    const __m256i v2 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p + s + 2));
    const __m256i v3 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p + s + 3));
    const __m256i m = _mm256_and_si256(_mm256_and_si256(_mm256_cmpeq_epi8(v0, b0),
                                                         _mm256_cmpeq_epi8(v1, b1)),
                                       _mm256_and_si256(_mm256_cmpeq_epi8(v2, b2),
                                                         _mm256_cmpeq_epi8(v3, b3)));
    const unsigned mask = static_cast<unsigned>(_mm256_movemask_epi8(m));
    if (mask != 0) return s + static_cast<size_t>(__builtin_ctz(mask));
  }
  return scalar::find_seed4(p, s, to_inclusive, seed);
}

}  // namespace avx2
#endif  // FQ_SIMD_X86

// ----------------------------------------------------------------- AVX-512
// 64-byte strides with MASKED tails: a masked load does not fault on the bytes
// it does not load, so the remainder is one more vector step, not a scalar
// loop. Every compare on a masked tail is itself masked, because the zeros a
// masked load fills in would otherwise be counted ("0 >= thr" holds for thr 0,
// and 0 - 33 is not printable).
#if defined(FQ_SIMD_X86_AVX512)
namespace avx512 {

inline __mmask64 tail_mask(size_t rem) noexcept {  // rem in [1, 63]
  return static_cast<__mmask64>(~0ULL >> (64 - rem));
}

__attribute__((target("avx512f,avx512bw")))
inline uint64_t sum_u8(const uint8_t* p, size_t n) noexcept {
  size_t i = 0;
  const __m512i zero = _mm512_setzero_si512();
  __m512i acc = zero;
  for (; i + 64 <= n; i += 64) {
    acc = _mm512_add_epi64(acc, _mm512_sad_epu8(_mm512_loadu_si512(p + i), zero));
  }
  if (i < n) {
    const __m512i v = _mm512_maskz_loadu_epi8(tail_mask(n - i), p + i);
    acc = _mm512_add_epi64(acc, _mm512_sad_epu8(v, zero));
  }
  return static_cast<uint64_t>(_mm512_reduce_add_epi64(acc));
}

__attribute__((target("avx512f,avx512bw")))
inline size_t count_ge(const uint8_t* p, size_t n, uint8_t thr) noexcept {
  size_t i = 0;
  size_t count = 0;
  const __m512i t = _mm512_set1_epi8(static_cast<char>(thr));
  for (; i + 64 <= n; i += 64) {
    count += static_cast<size_t>(
        __builtin_popcountll(_mm512_cmpge_epu8_mask(_mm512_loadu_si512(p + i), t)));
  }
  if (i < n) {
    const __mmask64 m = tail_mask(n - i);
    count += static_cast<size_t>(__builtin_popcountll(
        _mm512_mask_cmpge_epu8_mask(m, _mm512_maskz_loadu_epi8(m, p + i), t)));
  }
  return count;
}

__attribute__((target("avx512f,avx512bw")))
inline void add_widen_u8_to_u32(const uint8_t* src, size_t n, uint32_t* dst) noexcept {
  size_t i = 0;
  for (; i + 16 <= n; i += 16) {
    const __m512i w =
        _mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src + i)));
    _mm512_storeu_si512(dst + i, _mm512_add_epi32(_mm512_loadu_si512(dst + i), w));
  }
  for (; i < n; ++i) dst[i] += src[i];
}

__attribute__((target("avx512f,avx512bw")))
inline size_t mismatches(const uint8_t* a, const uint8_t* b, size_t n, size_t cap) noexcept {
  size_t i = 0;
  size_t diff = 0;
  for (; i + 64 <= n; i += 64) {
    diff += static_cast<size_t>(__builtin_popcountll(
        _mm512_cmpneq_epu8_mask(_mm512_loadu_si512(a + i), _mm512_loadu_si512(b + i))));
    if (diff > cap) return diff;
  }
  if (i < n) {
    const __mmask64 m = tail_mask(n - i);
    diff += static_cast<size_t>(__builtin_popcountll(_mm512_mask_cmpneq_epu8_mask(
        m, _mm512_maskz_loadu_epi8(m, a + i), _mm512_maskz_loadu_epi8(m, b + i))));
  }
  return diff;
}

__attribute__((target("avx512f,avx512bw")))
inline bool all_printable(const uint8_t* p, size_t n) noexcept {
  size_t i = 0;
  const __m512i off = _mm512_set1_epi8(33);
  const __m512i lim = _mm512_set1_epi8(93);
  __mmask64 bad = 0;
  for (; i + 64 <= n; i += 64) {
    bad |= _mm512_cmpgt_epu8_mask(_mm512_sub_epi8(_mm512_loadu_si512(p + i), off), lim);
  }
  if (i < n) {
    const __mmask64 m = tail_mask(n - i);
    bad |= _mm512_mask_cmpgt_epu8_mask(
        m, _mm512_sub_epi8(_mm512_maskz_loadu_epi8(m, p + i), off), lim);
  }
  return bad == 0;
}

// 64 candidate offsets at once. _mm512_cmpeq_epi8_mask already returns a
// per-lane mask directly (no movemask step needed), so the AND of four such
// masks is the same "matched all four byte planes" verdict AVX2's kernel
// builds from vectors, and ctzll on it is the 5'-most match in the block.
__attribute__((target("avx512f,avx512bw")))
inline size_t find_seed4(const uint8_t* p, size_t from, size_t to_inclusive,
                         uint32_t seed) noexcept {
  if (from > to_inclusive) return kNotFound;
  uint8_t sb[4];
  std::memcpy(sb, &seed, 4);
  const __m512i b0 = _mm512_set1_epi8(static_cast<char>(sb[0]));
  const __m512i b1 = _mm512_set1_epi8(static_cast<char>(sb[1]));
  const __m512i b2 = _mm512_set1_epi8(static_cast<char>(sb[2]));
  const __m512i b3 = _mm512_set1_epi8(static_cast<char>(sb[3]));
  size_t s = from;
  // The last offset a 64-wide block covers is s+63, reading p[s+63 .. s+67);
  // taken only while that is in bounds, the same precondition
  // scalar::find_seed4 applies per offset.
  for (; s + 63 <= to_inclusive; s += 64) {
    const __m512i v0 = _mm512_loadu_si512(p + s);
    const __m512i v1 = _mm512_loadu_si512(p + s + 1);
    const __m512i v2 = _mm512_loadu_si512(p + s + 2);
    const __m512i v3 = _mm512_loadu_si512(p + s + 3);
    const __mmask64 m = _mm512_cmpeq_epi8_mask(v0, b0) & _mm512_cmpeq_epi8_mask(v1, b1) &
                        _mm512_cmpeq_epi8_mask(v2, b2) & _mm512_cmpeq_epi8_mask(v3, b3);
    if (m != 0) return s + static_cast<size_t>(__builtin_ctzll(m));
  }
  return scalar::find_seed4(p, s, to_inclusive, seed);
}

}  // namespace avx512
#endif  // FQ_SIMD_X86_AVX512

// ---------------------------------------------------------------- dispatch

// Compiled into this binary.
inline bool isa_compiled(Isa isa) noexcept {
  switch (isa) {
    case Isa::kScalar:
      return true;
    case Isa::kNeon:
#if defined(FQ_SIMD_NEON)
      return true;
#else
      return false;
#endif
    case Isa::kAvx2:
#if defined(FQ_SIMD_X86)
      return true;
#else
      return false;
#endif
    case Isa::kAvx512:
#if defined(FQ_SIMD_X86_AVX512)
      return true;
#else
      return false;
#endif
  }
  return false;
}

// Compiled in AND executable on this CPU. `__builtin_cpu_supports` also checks
// that the OS saves the wider registers (XGETBV), not only the CPUID bit.
inline bool isa_supported(Isa isa) noexcept {
  if (!isa_compiled(isa)) return false;
  switch (isa) {
    case Isa::kScalar:
    case Isa::kNeon:
      return true;
    case Isa::kAvx2:
#if defined(FQ_SIMD_X86)
      return __builtin_cpu_supports("avx2");
#else
      return false;
#endif
    case Isa::kAvx512:
#if defined(FQ_SIMD_X86_AVX512)
      return __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512bw");
#else
      return false;
#endif
  }
  return false;
}

inline const char* isa_name(Isa isa) noexcept {
  // The strings the JSON report's "simd_backend" field has always carried.
  switch (isa) {
    case Isa::kAvx512: return "AVX-512BW";
    case Isa::kAvx2:   return "AVX2";
    case Isa::kNeon:   return "NEON";
    case Isa::kScalar: break;
  }
  return "scalar";
}

inline Isa detect_isa() noexcept {
  if (const char* forced = std::getenv("FQ_FORCE_ISA"); forced != nullptr) {
    const struct { const char* name; Isa isa; } kNames[] = {
        {"scalar", Isa::kScalar}, {"neon", Isa::kNeon}, {"avx2", Isa::kAvx2}, {"avx512", Isa::kAvx512}};
    for (const auto& entry : kNames) {
      if (std::strcmp(forced, entry.name) == 0 && isa_compiled(entry.isa)) return entry.isa;
    }
  }
  for (const Isa isa : {Isa::kAvx512, Isa::kAvx2, Isa::kNeon}) {
    if (isa_supported(isa)) return isa;
  }
  return Isa::kScalar;
}

// Resolved once, on first use. A magic static: thread safe, and evaluated after
// static initialisation, so there is no __builtin_cpu_init ordering to manage.
inline Isa isa_level() noexcept {
  static const Isa level = detect_isa();
  return level;
}

inline const char* backend() noexcept { return isa_name(isa_level()); }

// Each kernel twice: an explicit-ISA entry point, which is what the tests drive
// for every ISA the CPU can execute, and the dispatching one production calls.
// An ISA that is not compiled in falls through to the scalar reference.
inline uint64_t sum_u8_isa(Isa isa, const uint8_t* p, size_t n) noexcept {
  switch (isa) {
#if defined(FQ_SIMD_X86_AVX512)
    case Isa::kAvx512: return avx512::sum_u8(p, n);
#endif
#if defined(FQ_SIMD_X86)
    case Isa::kAvx2: return avx2::sum_u8(p, n);
#endif
#if defined(FQ_SIMD_NEON)
    case Isa::kNeon: return neon::sum_u8(p, n);
#endif
    default: return scalar::sum_u8(p, n);
  }
}

inline size_t count_ge_isa(Isa isa, const uint8_t* p, size_t n, uint8_t thr) noexcept {
  switch (isa) {
#if defined(FQ_SIMD_X86_AVX512)
    case Isa::kAvx512: return avx512::count_ge(p, n, thr);
#endif
#if defined(FQ_SIMD_X86)
    case Isa::kAvx2: return avx2::count_ge(p, n, thr);
#endif
#if defined(FQ_SIMD_NEON)
    case Isa::kNeon: return neon::count_ge(p, n, thr);
#endif
    default: return scalar::count_ge(p, n, thr);
  }
}

inline void add_widen_u8_to_u32_isa(Isa isa, const uint8_t* src, size_t n,
                                    uint32_t* dst) noexcept {
  switch (isa) {
#if defined(FQ_SIMD_X86_AVX512)
    case Isa::kAvx512: avx512::add_widen_u8_to_u32(src, n, dst); return;
#endif
#if defined(FQ_SIMD_X86)
    case Isa::kAvx2: avx2::add_widen_u8_to_u32(src, n, dst); return;
#endif
#if defined(FQ_SIMD_NEON)
    case Isa::kNeon: neon::add_widen_u8_to_u32(src, n, dst); return;
#endif
    default: scalar::add_widen_u8_to_u32(src, n, dst); return;
  }
}

inline size_t mismatches_isa(Isa isa, const uint8_t* a, const uint8_t* b, size_t n,
                             size_t cap) noexcept {
  switch (isa) {
#if defined(FQ_SIMD_X86_AVX512)
    case Isa::kAvx512: return avx512::mismatches(a, b, n, cap);
#endif
#if defined(FQ_SIMD_X86)
    case Isa::kAvx2: return avx2::mismatches(a, b, n, cap);
#endif
#if defined(FQ_SIMD_NEON)
    case Isa::kNeon: return neon::mismatches(a, b, n, cap);
#endif
    default: return scalar::mismatches(a, b, n, cap);
  }
}

inline bool all_printable_isa(Isa isa, const uint8_t* p, size_t n) noexcept {
  switch (isa) {
#if defined(FQ_SIMD_X86_AVX512)
    case Isa::kAvx512: return avx512::all_printable(p, n);
#endif
#if defined(FQ_SIMD_X86)
    case Isa::kAvx2: return avx2::all_printable(p, n);
#endif
#if defined(FQ_SIMD_NEON)
    case Isa::kNeon: return neon::all_printable(p, n);
#endif
    default: return scalar::all_printable(p, n);
  }
}

inline size_t find_seed4_isa(Isa isa, const uint8_t* p, size_t from, size_t to_inclusive,
                             uint32_t seed) noexcept {
  switch (isa) {
#if defined(FQ_SIMD_X86_AVX512)
    case Isa::kAvx512: return avx512::find_seed4(p, from, to_inclusive, seed);
#endif
#if defined(FQ_SIMD_X86)
    case Isa::kAvx2: return avx2::find_seed4(p, from, to_inclusive, seed);
#endif
#if defined(FQ_SIMD_NEON)
    case Isa::kNeon: return neon::find_seed4(p, from, to_inclusive, seed);
#endif
    default: return scalar::find_seed4(p, from, to_inclusive, seed);
  }
}

// Sum of raw ASCII quality bytes.
//
// NOT a Phred sum, and it must not be turned into one by subtracting
// len*kPhredOffset from the result: that is an unsigned subtraction, and a
// single quality byte below '!' wraps it to ~1.8e19. That defect has now been
// fixed twice (AUDIT.md S3 in accumulate_read, REVIEW_2026-08-15 finding 3 in
// mean_phred). Any production path needing Phred values must clamp per byte,
// which is why this kernel currently has no production caller and is exercised
// only by tests/benchmark.cpp and tests/test_trimmer.cpp. It is kept for the
// kernel-exactness tests; do not reintroduce the subtraction.
inline uint64_t sum_u8(const uint8_t* p, size_t n) noexcept {
  return sum_u8_isa(isa_level(), p, n);
}

// Number of bytes >= thr (unsigned compare).
inline size_t count_ge(const uint8_t* p, size_t n, uint8_t thr) noexcept {
  return count_ge_isa(isa_level(), p, n, thr);
}

// dst[i] += src[i] for i in [0,n). Used for per-cycle quality accumulation.
inline void add_widen_u8_to_u32(const uint8_t* src, size_t n, uint32_t* dst) noexcept {
  add_widen_u8_to_u32_isa(isa_level(), src, n, dst);
}

// Hamming distance over n bytes, saturating early once `cap` is exceeded.
// Early exit keeps the adapter scan O(read) in the common no-match case.
inline size_t mismatches(const uint8_t* a, const uint8_t* b, size_t n, size_t cap) noexcept {
  return mismatches_isa(isa_level(), a, b, n, cap);
}

// True when every byte is in '!'..'~' (33..126), the range FASTQ allows in
// SEQ and QUAL.
inline bool all_printable(const uint8_t* p, size_t n) noexcept {
  return all_printable_isa(isa_level(), p, n);
}

// Smallest offset in [from, to_inclusive] whose 4-byte word equals `seed`, or
// kNotFound. See scalar::find_seed4 for the exact contract (including the
// `to_inclusive + 4 <= buffer length` precondition every ISA relies on).
// AdapterTrimmer::trim uses this to skip past non-matching offsets in bulk
// instead of testing one 4-byte word at a time.
inline size_t find_seed4(const uint8_t* p, size_t from, size_t to_inclusive,
                         uint32_t seed) noexcept {
  return find_seed4_isa(isa_level(), p, from, to_inclusive, seed);
}

}  // namespace fq::simd
