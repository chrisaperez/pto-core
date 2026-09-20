cmake_minimum_required(VERSION 3.20)

project(peaks
        VERSION 0.1.0
        DESCRIPTION "Streaming epigenomic peak caller: single-pass pileup, multi-scale background, SIMD Poisson"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(PEAKS_BUILD_TESTS   "Build the test suite" ON)
option(PEAKS_ENABLE_AVX512 "Enable AVX-512 code paths" OFF)

# ---------------------------------------------------------------------------
# BGZF/BAM input: optional, like genomic_toolkit's htslib and for the same
# reason.
#
# This module's premise is that it builds on a cluster with a compiler and
# nothing else, and the fragment-BED path needs no decompressor at all. So BAM
# support is discovered, not required, and its absence is a message(STATUS)
# rather than an error -- the whole test suite still runs, because the BGZF
# suite skips itself loudly and successfully.
#
# libdeflate is preferred (it is what fastq_stream uses and it is roughly twice
# zlib's inflate throughput on BGZF blocks); zlib backs it because zlib is on
# every machine that has a compiler, and raw DEFLATE is what a BGZF block
# contains -- inflateInit2 with a negative windowBits reads it directly.
#
# Deliberately NO FetchContent fallback, unlike fastq_stream: fetching sources
# at configure time needs a network, and an offline build that silently loses
# BAM support is better than one that fails.
set(PEAKS_WITH_BGZF "AUTO" CACHE STRING "BAM/BGZF input: ON, OFF or AUTO")
set_property(CACHE PEAKS_WITH_BGZF PROPERTY STRINGS ON OFF AUTO)

# ---------------------------------------------------------------------------
# Header-only core
# ---------------------------------------------------------------------------
add_library(peaks_core INTERFACE)
target_include_directories(peaks_core INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_compile_features(peaks_core INTERFACE cxx_std_20)

set(PEAKS_BGZF_BACKEND "none")
if(NOT PEAKS_WITH_BGZF STREQUAL "OFF")
  find_path(PEAKS_LIBDEFLATE_INCLUDE_DIR libdeflate.h)
  find_library(PEAKS_LIBDEFLATE_LIB NAMES deflate libdeflate libdeflate.a)
  if(PEAKS_LIBDEFLATE_INCLUDE_DIR AND PEAKS_LIBDEFLATE_LIB)
    target_include_directories(peaks_core SYSTEM INTERFACE ${PEAKS_LIBDEFLATE_INCLUDE_DIR})
    target_link_libraries(peaks_core INTERFACE ${PEAKS_LIBDEFLATE_LIB})
    target_compile_definitions(peaks_core INTERFACE PEAKS_HAVE_BGZF=1 PEAKS_BGZF_LIBDEFLATE=1)
    set(PEAKS_BGZF_BACKEND "libdeflate")
  else()
    find_package(ZLIB QUIET)
    if(ZLIB_FOUND)
      target_link_libraries(peaks_core INTERFACE ZLIB::ZLIB)
      target_compile_definitions(peaks_core INTERFACE PEAKS_HAVE_BGZF=1 PEAKS_BGZF_ZLIB=1)
      set(PEAKS_BGZF_BACKEND "zlib")
    endif()
  endif()

  if(PEAKS_BGZF_BACKEND STREQUAL "none" AND PEAKS_WITH_BGZF STREQUAL "ON")
    message(FATAL_ERROR
      "PEAKS_WITH_BGZF=ON but neither libdeflate nor zlib was found. Install "
      "one (brew install libdeflate / apt install libdeflate-dev zlib1g-dev), "
      "or configure with -DPEAKS_WITH_BGZF=OFF to build the fragment-BED path "
      "only.")
  endif()
endif()

if(PEAKS_BGZF_BACKEND STREQUAL "none")
  message(STATUS "peaks: BAM/BGZF input DISABLED; fragment BED still available")
else()
  message(STATUS "peaks: BAM/BGZF input via ${PEAKS_BGZF_BACKEND}")
endif()

# ---------------------------------------------------------------------------
# SIMD capability detection
#
# Same arrangement, and the same reasoning, as modules/scrna_matrix: what is
# probed for is not "does this compiler accept -mavx2" but "can it build a
# function carrying a target attribute", because that -- not a flag on the
# target -- is how poisson_model.hpp reaches AVX2. Putting the ISA flag on the
# interface would let the compiler emit AVX2 anywhere it liked, including in
# code no runtime check guards, and the binary would SIGILL on a CPU without
# it. NEON needs no probe: Advanced SIMD is mandatory in ARMv8-A, so the header
# gates it on __aarch64__ alone.
# ---------------------------------------------------------------------------
include(CheckCXXSourceCompiles)

check_cxx_source_compiles("
#include <immintrin.h>
__attribute__((target(\"avx2,fma\")))
static double probe(const double* p, const int* q) {
  __m256d v = _mm256_loadu_pd(p);
  __m128i i = _mm_loadu_si128(reinterpret_cast<const __m128i*>(q));
  __m256d g = _mm256_i32gather_pd(p, i, 8);
  return _mm256_cvtsd_f64(_mm256_fmadd_pd(v, g, v));
}
int main() { double a[4]={0}; int b[4]={0}; return (int)probe(a,b); }
" PEAKS_HAS_AVX2_TARGET_ATTR)

check_cxx_source_compiles("
#include <immintrin.h>
__attribute__((target(\"avx512f,avx512dq\")))
static double probe(const double* p, const int* q) {
  __m512d v = _mm512_loadu_pd(p);
  __m256i i = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(q));
  __m512d g = _mm512_i32gather_pd(i, p, 8);
  return _mm512_cvtsd_f64(_mm512_fmadd_pd(v, g, v));
}
int main() { double a[8]={0}; int b[8]={0}; return (int)probe(a,b); }
" PEAKS_HAS_AVX512_TARGET_ATTR)

if(PEAKS_HAS_AVX2_TARGET_ATTR)
  target_compile_definitions(peaks_core INTERFACE PEAKS_HAVE_AVX2=1)
  message(STATUS "peaks: AVX2 Poisson/lgamma kernels compiled in (runtime-dispatched).")
else()
  message(STATUS "peaks: AVX2 unavailable for this target.")
endif()

if(PEAKS_ENABLE_AVX512 AND PEAKS_HAS_AVX512_TARGET_ATTR)
  target_compile_definitions(peaks_core INTERFACE PEAKS_HAVE_AVX512=1)
  message(STATUS "peaks: AVX-512 Poisson/lgamma kernels compiled in (runtime-dispatched).")
endif()

# ---------------------------------------------------------------------------
# Executable
#
# Registered here rather than in the top-level CMakeLists, which only
# aggregates: every module in this tree owns its own targets and must stay
# buildable on its own (`cmake -S modules/peaks -B build/peaks`).
# ---------------------------------------------------------------------------
add_executable(pto-peaks src/main.cpp)
target_link_libraries(pto-peaks PRIVATE peaks_core)
target_compile_definitions(pto-peaks PRIVATE PEAKS_VERSION="${PROJECT_VERSION}")
if(NOT MSVC)
  target_compile_options(pto-peaks PRIVATE -Wall -Wextra -Wshadow
                         $<$<CONFIG:Release>:-O3>)
endif()
install(TARGETS pto-peaks RUNTIME DESTINATION bin)

# ---------------------------------------------------------------------------
# Tests
#
# Dependency-free, like every other module's: a bare compiler is enough.
# ---------------------------------------------------------------------------
if(PEAKS_BUILD_TESTS)
  enable_testing()

  add_executable(test_peaks_poisson_model tests/test_poisson_model.cpp)
  target_link_libraries(test_peaks_poisson_model PRIVATE peaks_core)
  # The last case in this binary scans the header for a scalar libm call, so it
  # needs to know where the header is. Passed as a definition rather than
  # assumed relative to the cwd: ctest runs from the build tree, and a guard
  # that silently finds no file reads as a pass.
  target_compile_definitions(test_peaks_poisson_model
                             PRIVATE PTO_PEAKS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
  if(NOT MSVC)
    target_compile_options(test_peaks_poisson_model PRIVATE -Wall -Wextra -Wshadow)
  endif()
  add_test(NAME peaks_poisson_model COMMAND test_peaks_poisson_model)

  # The streaming halves. One binary and one add_test per structure, so
  # `ctest -R peaks_sliding_window` runs a single suite.
  foreach(SUITE sliding_window caller bam_streamer file_identity)
    add_executable(test_peaks_${SUITE} tests/test_${SUITE}.cpp)
    target_link_libraries(test_peaks_${SUITE} PRIVATE peaks_core)
    target_include_directories(test_peaks_${SUITE} PRIVATE
                               ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    # The sliding-window suite scans its own header for an integer division on
    # the hot path, so it needs to know where the sources are.
    target_compile_definitions(test_peaks_${SUITE}
                               PRIVATE PTO_PEAKS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    if(NOT MSVC)
      target_compile_options(test_peaks_${SUITE} PRIVATE -Wall -Wextra -Wshadow)
    endif()
    add_test(NAME peaks_${SUITE} COMMAND test_peaks_${SUITE})
  endforeach()

  # The unit suites drive the headers directly. This one drives the BINARY,
  # which is the only way to cover argument parsing, the BED reader, the exit
  # codes pto-cloud's job contract depends on, and the narrowPeak writer.
  add_test(NAME peaks_cli
           COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/test_cli.sh
                   $<TARGET_FILE:pto-peaks>
                   $<TARGET_FILE:test_peaks_bam_streamer>)
endif()

message(STATUS "peaks ${PROJECT_VERSION} (${CMAKE_BUILD_TYPE})")
