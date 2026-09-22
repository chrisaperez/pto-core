# SPDX-License-Identifier: MIT
#
# Shared OpenMP toolchain resolution for pto-core.
#
# WHY THIS FILE EXISTS
# --------------------
# OpenMP is mandatory for modules/scrna_matrix (see that module's CMakeLists for
# the reasoning: without it `_OPENMP` is undefined, every `#pragma omp` compiles
# to serial code, and the test suite reports a clean run while never executing a
# parallel iteration).
#
# But "mandatory" is only a useful control if it is *satisfiable*. On the most
# common developer machine in this project -- macOS with Apple Clang -- CMake's
# FindOpenMP reports NOT FOUND even after `brew install libomp`, because Apple
# Clang needs `-Xclang -fopenmp` rather than `-fopenmp`, and because Homebrew
# installs libomp as a keg-only formula whose headers and library are not on any
# default search path. A hard error there is a broken clean clone, not a safety
# feature.
#
# pto_openmp_apply_hints() closes that gap: it locates the keg and hands
# FindOpenMP the variables it cannot derive on its own. It is NOT a fallback and
# it never relaxes anything -- callers still run find_package(OpenMP REQUIRED)
# afterwards, and that search still hard-fails if the result is unusable.
#
# Idempotent, and safe to call from several directories.

include_guard(GLOBAL)

# Populate OpenMP_* hint variables on platforms where FindOpenMP cannot manage
# by itself. No-op everywhere else (Linux GCC/Clang with libgomp/libomp-dev, and
# Homebrew GCC on macOS, all work unaided).
function(pto_openmp_apply_hints)
  if(OpenMP_CXX_FOUND)
    return()
  endif()
  if(NOT APPLE)
    return()
  endif()
  # Homebrew GCC brings its own libgomp and needs no help; only the Clang family
  # (AppleClang, and Homebrew LLVM clang) needs the keg wired up.
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    return()
  endif()

  set(_hints "")
  # Explicit override wins, so a user with a hand-built runtime can point at it:
  #   cmake -DLIBOMP_ROOT=/path/to/libomp ...
  if(DEFINED LIBOMP_ROOT)
    list(APPEND _hints "${LIBOMP_ROOT}")
  endif()

  # `brew --prefix libomp` is authoritative when Homebrew lives somewhere
  # non-standard (a surprisingly common setup on managed laptops).
  find_program(_pto_brew brew)
  if(_pto_brew)
    execute_process(
      COMMAND "${_pto_brew}" --prefix libomp
      OUTPUT_VARIABLE _brew_prefix
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
      RESULT_VARIABLE _brew_rc)
    if(_brew_rc EQUAL 0 AND _brew_prefix)
      list(APPEND _hints "${_brew_prefix}")
    endif()
  endif()

  list(APPEND _hints
       /opt/homebrew/opt/libomp   # Homebrew, Apple silicon
       /usr/local/opt/libomp)     # Homebrew, Intel
  list(REMOVE_DUPLICATES _hints)

  foreach(_prefix IN LISTS _hints)
    if(NOT EXISTS "${_prefix}/include/omp.h")
      continue()
    endif()
    unset(_omp_lib CACHE)
    find_library(_omp_lib NAMES omp iomp5 gomp HINTS "${_prefix}/lib" NO_DEFAULT_PATH)
    if(NOT _omp_lib)
      continue()
    endif()

    # The variables FindOpenMP consumes when it cannot derive them itself.
    # PARENT_SCOPE is not enough: FindOpenMP runs in whichever directory calls
    # find_package, so these go into the cache.
    #
    # OpenMP_<lang>_INCLUDE_DIR is the one FindOpenMP reads to populate
    # OpenMP_<lang>_INCLUDE_DIRS and, from that, the imported target's
    # INTERFACE_INCLUDE_DIRECTORIES (see FindOpenMP.cmake's final
    # set_property(... INTERFACE_INCLUDE_DIRECTORIES ...) step). Without it,
    # FindOpenMP falls back to an unhinted find_path(), which cannot see a
    # keg-only Homebrew install and leaves the imported target with no
    # include path at all -- '#include <omp.h>' then fails not just for
    # ordinary targets but for anything that links against
    # OpenMP::OpenMP_CXX in isolation, such as check_cxx_source_compiles().
    set(OpenMP_CXX_FLAGS "-Xclang -fopenmp" CACHE STRING "OpenMP C++ flags" FORCE)
    set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "OpenMP C++ libraries" FORCE)
    set(OpenMP_omp_LIBRARY "${_omp_lib}" CACHE FILEPATH "libomp" FORCE)
    set(OpenMP_CXX_INCLUDE_DIR "${_prefix}/include" CACHE PATH "OpenMP C++ include dir" FORCE)
    set(OpenMP_C_FLAGS "-Xclang -fopenmp" CACHE STRING "OpenMP C flags" FORCE)
    set(OpenMP_C_LIB_NAMES "omp" CACHE STRING "OpenMP C libraries" FORCE)
    set(OpenMP_C_INCLUDE_DIR "${_prefix}/include" CACHE PATH "OpenMP C include dir" FORCE)

    # Belt-and-suspenders for anything that consults directory-level state
    # instead of the imported target (e.g. plain add_executable() targets
    # that never touch OpenMP::OpenMP_CXX directly).
    include_directories(SYSTEM "${_prefix}/include")
    link_directories("${_prefix}/lib")

    set(PTO_LIBOMP_PREFIX "${_prefix}" CACHE INTERNAL "detected libomp keg")
    message(STATUS "libomp: ${_prefix} (Homebrew keg, auto-detected for Apple Clang)")
    return()
  endforeach()
endfunction()

# The remediation text shown when OpenMP is genuinely unavailable. Kept in one
# place so the root's early diagnosis and the module's hard failure cannot drift
# apart, and so it stays actionable rather than merely descriptive.
function(pto_openmp_failure_message out_var)
  set(_msg
"OpenMP is REQUIRED to build modules/scrna_matrix, and no usable OpenMP \
toolchain was found.

This is a hard error by design. Without OpenMP, `_OPENMP` is undefined, every \
`#pragma omp` in knn_graph.hpp is ignored, and the build produces a library \
that passes its entire test suite without ever running a parallel iteration -- \
a silently single-threaded engine whose concurrency code is never exercised.

Fix the toolchain:

  macOS (Apple Clang, recommended):
      brew install libomp
    The keg is then detected automatically; no extra flags are needed. If \
Homebrew is installed in a non-standard location, pass it explicitly:
      cmake -S . -B build -DLIBOMP_ROOT=$(brew --prefix libomp)

  macOS (Homebrew GCC, OpenMP built in):
      brew install gcc
      cmake -S . -B build -DCMAKE_C_COMPILER=gcc-15 -DCMAKE_CXX_COMPILER=g++-15

  Debian/Ubuntu (Clang):
      sudo apt-get install libomp-dev
  Debian/Ubuntu (GCC):
      already present (libgomp)

Or build only the modules that do not use OpenMP -- fastq_stream and \
cuttag_profiler are std::thread-based and contain no `#pragma omp`:
      cmake -S . -B build -DPTO_BUILD_SCRNA_MATRIX=OFF")
  set(${out_var} "${_msg}" PARENT_SCOPE)
endfunction()
