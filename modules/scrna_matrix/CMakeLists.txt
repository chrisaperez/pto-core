cmake_minimum_required(VERSION 3.20)
project(scrna_matrix LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(SCRNA_BUILD_TESTS "Build C++ unit tests" ON)
option(SCRNA_BUILD_PYTHON "Build pybind11 Python bindings" ON)
option(SCRNA_ENABLE_AVX512 "Enable AVX-512 code paths" OFF)
option(SCRNA_ENABLE_HNSW "Enable approximate k-NN via vendored hnswlib" ON)
option(SCRNA_HNSW_NATIVE_ISA
       "Compile hnswlib's own AVX/AVX-512 kernels by putting ISA flags on the target. \
Off by default: it makes the binary require those extensions at runtime. \
Our ScrnaCosineSpace already reaches AVX-512 safely via target attributes." OFF)

# ---------------------------------------------------------------------------
# Header-only core library
# ---------------------------------------------------------------------------
add_library(scrna_matrix INTERFACE)
target_include_directories(scrna_matrix INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_compile_features(scrna_matrix INTERFACE cxx_std_20)

# ---------------------------------------------------------------------------
# SIMD / arch capability detection
#
# Deliberately NO -mavx2 / -mavx512f on the target. Applying an ISA flag to the
# whole interface lets the compiler emit those instructions anywhere it likes,
# including in code that never goes through a runtime capability check -- so the
# resulting binary SIGILLs on any CPU that lacks the extension. Instead the
# kernels in simd_math.hpp carry per-function `target` attributes and are
# selected at runtime, which keeps the rest of the binary at the baseline ISA.
#
# So what we probe for here is not "does the compiler accept -mavx2" but "can
# this compiler build a target-attributed AVX2 function", which is the mechanism
# actually relied upon. On non-x86 targets (arm64/Apple Silicon) both probes fail
# and only the scalar path is compiled.
# ---------------------------------------------------------------------------
include(CheckCXXSourceCompiles)

check_cxx_source_compiles("
#include <immintrin.h>
__attribute__((target(\"avx2,fma\")))
static float probe(const float* p, const int* q) {
  __m256i vi = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(q));
  __m256 vg = _mm256_i32gather_ps(p, vi, 4);
  return _mm256_cvtss_f32(_mm256_fmadd_ps(vg, vg, vg));
}
int main() { float a[8]={0}; int b[8]={0}; return (int)probe(a,b); }
" SCRNA_HAS_AVX2_TARGET_ATTR)

check_cxx_source_compiles("
#include <immintrin.h>
__attribute__((target(\"avx512f,avx512bw\")))
static float probe(const float* p, const int* q) {
  __m512i vi = _mm512_loadu_si512(reinterpret_cast<const void*>(q));
  __m512 vg = _mm512_i32gather_ps(vi, p, 4);
  return _mm512_reduce_add_ps(vg);
}
int main() { float a[16]={0}; int b[16]={0}; return (int)probe(a,b); }
" SCRNA_HAS_AVX512_TARGET_ATTR)

if(SCRNA_HAS_AVX2_TARGET_ATTR)
  target_compile_definitions(scrna_matrix INTERFACE SCRNA_HAVE_AVX2=1)
  message(STATUS "scrna_matrix: AVX2 gather kernel compiled in (runtime-dispatched).")
else()
  message(STATUS "scrna_matrix: AVX2 unavailable for this target; scalar path only.")
endif()

if(SCRNA_ENABLE_AVX512 AND SCRNA_HAS_AVX512_TARGET_ATTR)
  target_compile_definitions(scrna_matrix INTERFACE SCRNA_HAVE_AVX512=1)
  message(STATUS "scrna_matrix: AVX-512 gather kernel compiled in (runtime-dispatched).")
endif()

# ---------------------------------------------------------------------------
# hnswlib (vendored, header-only, Apache-2.0) for approximate k-NN
# ---------------------------------------------------------------------------
if(SCRNA_ENABLE_HNSW)
  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third_party/hnswlib/hnswlib.h)
    message(FATAL_ERROR "SCRNA_ENABLE_HNSW=ON but third_party/hnswlib/hnswlib.h is missing.")
  endif()
  # SYSTEM: hnswlib is vendored and we do not own its warnings. Without this,
  # -Wshadow and -Wunused-parameter fire inside hnswalg.h and drown the
  # first-party diagnostics that actually need acting on -- and make a
  # -Werror build impossible for reasons that have nothing to do with our code.
  target_include_directories(scrna_matrix SYSTEM INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party>)
  target_compile_definitions(scrna_matrix INTERFACE SCRNA_ENABLE_HNSW=1)

  if(SCRNA_HNSW_NATIVE_ISA)
    # hnswlib gates its AVX kernels on compile-time __AVX__/__AVX512F__ and does
    # not use per-function target attributes, so reaching them means putting the
    # ISA flag on the target -- which is exactly the arrangement the rest of this
    # build avoids, because it lets the compiler emit those instructions in code
    # that no runtime check guards.
    if(SCRNA_HAS_AVX2_TARGET_ATTR)
      target_compile_options(scrna_matrix INTERFACE -mavx2 -mfma)
    endif()
    if(SCRNA_ENABLE_AVX512 AND SCRNA_HAS_AVX512_TARGET_ATTR)
      target_compile_options(scrna_matrix INTERFACE -mavx512f -mavx512bw)
    endif()
    message(WARNING
      "SCRNA_HNSW_NATIVE_ISA=ON: the resulting binary may execute AVX2/AVX-512 "
      "instructions without a runtime guard and will fault on CPUs lacking them. "
      "Only use this when the build and target CPUs are known to match.")
  else()
    message(STATUS
      "scrna_matrix: hnswlib enabled; its distances run through ScrnaCosineSpace "
      "(runtime-dispatched AVX2/AVX-512). hnswlib's own kernels stay at baseline ISA.")
  endif()
endif()

# ---------------------------------------------------------------------------
# OpenMP -- MANDATORY.
#
# This was `find_package(OpenMP QUIET)` with a `message(WARNING)` fallback, and
# that fallback was actively harmful. Without OpenMP, `_OPENMP` is undefined and
# every `#pragma omp` in knn_graph.hpp is ignored by the preprocessor, so the
# library still builds, still passes its whole test suite, and still reports
# "100% tests passed" -- while never executing a single parallel iteration.
#
# That is exactly what had happened here: the shipped configuration on the
# development host had OpenMP_CXX_FLAGS=NOTFOUND, so the concurrency paths
# (parallel_for_checked, the per-thread ThreadState indexing, the atomics
# ordering in the error latch) were dead code that CI was silently signing off
# on. A warning in a thousand-line configure log is not a control.
#
# So it is REQUIRED, and configuration fails on a host that cannot supply it.
# The single-threaded build is not a supported configuration: this library's
# stated purpose is k-NN over 1M+ cells, where losing the parallel loop is a
# correctness-of-benchmark and time-to-result problem, not a minor degradation.
#
# There is deliberately NO opt-out variable. If you need to build without
# OpenMP, that is a decision to make explicitly by editing this line, so it
# shows up in review and in `git blame`, rather than by setting a flag that
# silently reintroduces untested single-threaded stubs.
#
# Toolchain notes (this is the common failure and the message must fix it, not
# merely report it):
#   * Apple Clang ships no OpenMP runtime.   -> brew install libomp (auto-wired below)
#   * Homebrew GCC has it built in.          -> -DCMAKE_CXX_COMPILER=g++-15
#   * Debian/Ubuntu Clang.                   -> apt install libomp-dev

# Resolve the toolchain before searching. On macOS + Apple Clang this locates
# Homebrew's keg-only libomp and hands FindOpenMP the flags it cannot derive
# itself; everywhere else it is a no-op. See cmake/PtoOpenMP.cmake -- it
# makes the requirement satisfiable without relaxing it.
include("${CMAKE_CURRENT_LIST_DIR}/../../cmake/PtoOpenMP.cmake" OPTIONAL
        RESULT_VARIABLE _scrna_omp_helper)
if(_scrna_omp_helper)
  pto_openmp_apply_hints()
endif()

find_package(OpenMP QUIET COMPONENTS CXX)
if(NOT OpenMP_CXX_FOUND)
  if(COMMAND pto_openmp_failure_message)
    pto_openmp_failure_message(_scrna_omp_help)
    message(FATAL_ERROR "\n${_scrna_omp_help}\n")
  endif()
  message(FATAL_ERROR
    "OpenMP is REQUIRED for scrna_matrix and was not found. "
    "macOS: brew install libomp. Debian/Ubuntu: apt-get install libomp-dev.")
endif()

if(NOT TARGET OpenMP::OpenMP_CXX)
  message(FATAL_ERROR
    "OpenMP was reported as found but the OpenMP::OpenMP_CXX imported target "
    "does not exist. Refusing to continue: linking would silently drop the "
    "runtime and every '#pragma omp' would compile to serial code.")
endif()

target_link_libraries(scrna_matrix INTERFACE OpenMP::OpenMP_CXX)

# `find_package(OpenMP)` succeeding is necessary but not sufficient: it proves
# the compiler accepts the flag, not that `_OPENMP` is actually defined when
# compiling this target with these flags. Prove the macro end-to-end, because
# `_OPENMP` is the thing every parallel region in this library is guarded on.
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)
cmake_push_check_state(RESET)
set(CMAKE_REQUIRED_FLAGS "${OpenMP_CXX_FLAGS}")
set(CMAKE_REQUIRED_LIBRARIES OpenMP::OpenMP_CXX)
check_cxx_source_compiles("
  #include <omp.h>
  #ifndef _OPENMP
  #error _OPENMP is not defined; parallel regions would compile to serial code
  #endif
  int main() {
    int n = 0;
    #pragma omp parallel reduction(+ : n)
    { n += 1; }
    return n > 0 ? 0 : 1;
  }" SCRNA_OPENMP_PRAGMAS_ACTIVE)
cmake_pop_check_state()

if(NOT SCRNA_OPENMP_PRAGMAS_ACTIVE)
  message(FATAL_ERROR
    "OpenMP::OpenMP_CXX was found, but a translation unit compiled with "
    "OpenMP_CXX_FLAGS ('${OpenMP_CXX_FLAGS}') either does not define _OPENMP or "
    "cannot compile a '#pragma omp parallel' region.\n"
    "Building in this state would produce a library whose k-NN construction is "
    "silently single-threaded and whose concurrency code is never executed. "
    "Refusing.\n"
    "  Apple Clang : brew install libomp\n"
    "  Homebrew GCC: cmake -DCMAKE_CXX_COMPILER=g++-15 ...\n"
    "  Debian/Ubuntu Clang: apt-get install libomp-dev")
endif()

message(STATUS "  OpenMP     : ${OpenMP_CXX_FLAGS} (pragmas verified active)")

# ---------------------------------------------------------------------------
# Eigen3 (optional; linked if present)
#
# Nothing in include/matrix currently uses Eigen -- the scalar fallbacks are
# plain C++ and have no third-party dependency. This stays as an opt-in hook for
# future dense-side work; its absence disables nothing.
# ---------------------------------------------------------------------------
find_package(Eigen3 QUIET NO_MODULE)
if(TARGET Eigen3::Eigen)
  target_link_libraries(scrna_matrix INTERFACE Eigen3::Eigen)
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(SCRNA_BUILD_TESTS)
  enable_testing()
  add_executable(test_matrix_ops tests/test_matrix_ops.cpp)
  target_link_libraries(test_matrix_ops PRIVATE scrna_matrix)
  add_test(NAME test_matrix_ops COMMAND test_matrix_ops)

  # Wide accumulation (docs/AUDIT.md S2, REVIEW finding 5 + its residual, T5).
  # The last case in this binary scans the module's own sources for a bare
  # `float` accumulator, so it needs to know where they are. Passed as a
  # definition rather than assumed relative to the cwd: ctest runs from the
  # build tree, and a guard that silently finds no files reads as a pass.
  add_executable(test_wide_accumulation tests/test_wide_accumulation.cpp)
  target_link_libraries(test_wide_accumulation PRIVATE scrna_matrix)
  target_compile_definitions(test_wide_accumulation
                             PRIVATE SCRNA_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_test(NAME test_wide_accumulation COMMAND test_wide_accumulation)

  # Index serialization. Built unconditionally; the translation unit compiles
  # to a no-op main() when SCRNA_ENABLE_HNSW is off, so the test list does not
  # change shape with the option.
  add_executable(test_hnsw_persist tests/test_hnsw_persist.cpp)
  target_link_libraries(test_hnsw_persist PRIVATE scrna_matrix)
  add_test(NAME test_hnsw_persist COMMAND test_hnsw_persist)
endif()

# ---------------------------------------------------------------------------
# Python bindings (pybind11)
# ---------------------------------------------------------------------------
if(SCRNA_BUILD_PYTHON)
  find_package(pybind11 CONFIG QUIET)
  if(pybind11_FOUND)
    pybind11_add_module(scrna_matrix_py src/python_bindings.cpp)
    target_link_libraries(scrna_matrix_py PRIVATE scrna_matrix)

    # Wheel packaging. scikit-build-core collects whatever CMake installs and
    # lays it into the wheel, so the extension goes *inside* the Python package
    # directory rather than at the wheel root -- python/scrna_matrix/__init__.py
    # re-exports it from there. Harmless for an ordinary `cmake --install`.
    install(TARGETS scrna_matrix_py
            LIBRARY DESTINATION scrna_matrix
            RUNTIME DESTINATION scrna_matrix)

    # pybind11 spells this PYTHON_EXECUTABLE; FindPython spells it Python_EXECUTABLE.
    set(_scrna_python "${PYTHON_EXECUTABLE}")
    if(NOT _scrna_python)
      set(_scrna_python "${Python_EXECUTABLE}")
    endif()

    # -------------------------------------------------------------------------
    # De-duplicate the OpenMP runtime the extension loads (macOS only).
    #
    # scrna_matrix links OpenMP::OpenMP_CXX, which on Apple Clang resolves to
    # Homebrew's keg-only libomp (PtoOpenMP.cmake, above). That is correct for
    # the plain C++ build. But when this extension is imported into a Python
    # whose numpy/scipy already loaded their OWN OpenMP runtime -- e.g. any
    # conda-forge interpreter, since conda-forge's numpy/scipy are built against
    # conda-forge's own `llvm-openmp` package -- dyld ends up with two DIFFERENT
    # libomp.dylib images resident in one process. Both are the same upstream
    # LLVM OpenMP runtime (verified: identical compatibility version 5.0.0,
    # identical __kmpc_* export set), so they are ABI-compatible, but the
    # runtime's own duplicate-initialization guard treats a second distinct
    # image as unsafe and aborts the process with "OMP: Error #15". The
    # unsafe workaround is KMP_DUPLICATE_LIB_OK=TRUE; the real fix is to never
    # load two copies in the first place.
    #
    # So: if the target Python's own environment ships a libomp.dylib distinct
    # from the one scrna_matrix compiled against, repoint the BUILT extension's
    # load command at that one instead, post-link. Only the .dylib reference
    # changes -- compilation still uses Homebrew's headers/flags, which is safe
    # because the ABI is the same runtime. This is scoped to scrna_matrix_py
    # only: the core library and its C++ test suite are untouched and keep
    # using Homebrew's libomp exactly as before.
    if(APPLE AND _scrna_python AND PTO_LIBOMP_PREFIX)
      execute_process(
        COMMAND "${_scrna_python}" -c "import sys; print(sys.prefix)"
        OUTPUT_VARIABLE _scrna_py_prefix
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
      set(_scrna_py_libomp "${_scrna_py_prefix}/lib/libomp.dylib")
      set(_scrna_built_libomp "${PTO_LIBOMP_PREFIX}/lib/libomp.dylib")
      if(_scrna_py_prefix AND EXISTS "${_scrna_py_libomp}"
         AND NOT "${_scrna_py_libomp}" STREQUAL "${_scrna_built_libomp}")
        find_program(_scrna_install_name_tool install_name_tool)
        if(_scrna_install_name_tool)
          add_custom_command(TARGET scrna_matrix_py POST_BUILD
            COMMAND "${_scrna_install_name_tool}" -change
                    "${_scrna_built_libomp}" "${_scrna_py_libomp}"
                    "$<TARGET_FILE:scrna_matrix_py>"
            COMMENT "scrna_matrix_py: repointing OpenMP runtime at ${_scrna_py_libomp} (the target Python's own copy) to avoid a duplicate-runtime abort")
          message(STATUS
            "scrna_matrix_py: will load OpenMP from ${_scrna_py_libomp} at import "
            "time instead of ${_scrna_built_libomp}, matching the runtime "
            "numpy/scipy already initialize in this interpreter.")
        endif()
      endif()
    endif()

    # Binding-level tests (adoption vs copy, buffer lifetime, GIL release).
    # Skips itself if pytest/numpy/scipy are unavailable.
    if(SCRNA_BUILD_TESTS AND _scrna_python)
      add_test(NAME test_python_bindings
               COMMAND ${_scrna_python} -m pytest
                       ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_python_bindings.py -q)
      set_tests_properties(test_python_bindings PROPERTIES
        ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:scrna_matrix_py>")
    endif()
  else()
    message(WARNING "pybind11 not found (pip install pybind11 or add as submodule); skipping Python module.")
  endif()
endif()
