# SPDX-License-Identifier: MIT
#
# pto-core -- top-level aggregate build for the five modules.
#
# Each module under modules/ is a self-contained CMake project with its own
# project(), options and install rules, and each is still buildable on its own:
#
#     cmake -S modules/fastq_stream -B build/fastq_stream && cmake --build ...
#
# That is deliberate -- the modules are deployed independently and were
# developed as separate repositories. This file adds them as subdirectories so
# the whole monorepo can be configured, built and tested in one pass, without
# taking any ownership of their internals.
#
# Because each module calls project() itself, module-scoped variables
# (CMAKE_CXX_STANDARD, compile options, etc.) stay module-scoped. Do not hoist
# them here; a module's build must not depend on being built from this file.
cmake_minimum_required(VERSION 3.20)

project(pto-core
        VERSION 0.1.0
        DESCRIPTION "Cache-aware C++20 engines for single-cell and epigenomic data"
        LANGUAGES C CXX)

# Modules are individually switchable so a host missing one module's
# dependencies (htslib for cuttag_profiler, pybind11 for scrna_matrix's Python
# bindings) can still build and test the rest.
option(PTO_BUILD_SCRNA_MATRIX    "Build modules/scrna_matrix"    ON)
option(PTO_BUILD_FASTQ_STREAM    "Build modules/fastq_stream"    ON)
option(PTO_BUILD_CUTTAG_PROFILER "Build modules/cuttag_profiler" ON)
option(PTO_BUILD_GENOMIC_TOOLKIT "Build modules/genomic_toolkit" ON)
option(PTO_BUILD_PEAKS           "Build modules/peaks"           ON)

# Cloud-worker build: compile out every network listener.
#
# This is not a runtime switch and deliberately not a flag the worker passes.
# `cuttag_profiler serve` runs an HTTP dashboard whose authentication is a
# session token bound to one operator's browser (SECURITY_HTTP_2026-08-15 H1);
# that is a sound control on a workstation and a tenant-isolation hole in a
# multi-tenant service, where one isolated task per job is the isolation
# boundary and a listening port inside it is a way around that boundary.
#
# With this ON the server source is not compiled, httplib is not included, the
# dashboard assets are not embedded, and `serve` is refused by name at dispatch
# (exit 1, and it is omitted from --help): there is no listener to reach by any
# argv, config file or future code path. A defence you can grep for in the
# linker output beats one you have to trust a caller to apply.
#
#     cmake -S . -B build -DPTO_CLOUD_BUILD=ON
#
# Only cuttag_profiler ships a server today; the option is declared here rather
# than in that module so a cloud image configures the whole tree with one flag
# and any module that later grows a listener inherits the same switch.
option(PTO_CLOUD_BUILD
       "Cloud-worker build: strip every embedded HTTP server from the binaries"
       OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# ---------------------------------------------------------------------------
# Toolchain policy: OpenMP is mandatory for any module that uses it.
#
# scrna_matrix enforces this itself (see modules/scrna_matrix/CMakeLists.txt,
# which additionally proves `_OPENMP` is defined and that a parallel region
# compiles). The check is repeated here only so the failure lands at the top of
# the configure output with actionable instructions, instead of surfacing a
# find_package error partway through a subdirectory.
#
# It is scoped to the modules that actually contain `#pragma omp`: today that is
# scrna_matrix alone. fastq_stream and cuttag_profiler use std::thread and a
# hand-written lock-free queue, and have zero OpenMP code -- requiring an
# OpenMP toolchain to build them would be a fabricated dependency that makes the
# build harder to satisfy without making anything safer.
if(PTO_BUILD_SCRNA_MATRIX)
  include(PtoOpenMP)
  # Locate Homebrew's keg-only libomp on macOS/Apple Clang, so that
  # `brew install libomp` followed by a plain `cmake -S . -B build` works with
  # no extra flags. A no-op on every other platform.
  pto_openmp_apply_hints()

  find_package(OpenMP QUIET COMPONENTS CXX)
  if(NOT OpenMP_CXX_FOUND)
    # Diagnosed here, at the top of the configure output, rather than surfacing
    # partway through a subdirectory. scrna_matrix repeats the check and is the
    # authoritative enforcement point.
    pto_openmp_failure_message(_pto_omp_help)
    message(FATAL_ERROR "\n${_pto_omp_help}\n")
  endif()
endif()

# One test registry for the whole tree, so `ctest` at the top level runs every
# module's suite. Each module also calls enable_testing() for its standalone
# build; calling it here first makes the aggregate build the one that owns the
# top-level CTestTestfile.cmake.
enable_testing()

set(_pto_enabled "")

if(PTO_BUILD_SCRNA_MATRIX)
  add_subdirectory(modules/scrna_matrix)
  list(APPEND _pto_enabled scrna_matrix)
endif()

if(PTO_BUILD_FASTQ_STREAM)
  add_subdirectory(modules/fastq_stream)
  list(APPEND _pto_enabled fastq_stream)
endif()

if(PTO_BUILD_CUTTAG_PROFILER)
  add_subdirectory(modules/cuttag_profiler)
  list(APPEND _pto_enabled cuttag_profiler)
endif()

if(PTO_BUILD_GENOMIC_TOOLKIT)
  add_subdirectory(modules/genomic_toolkit)
  list(APPEND _pto_enabled genomic_toolkit)
endif()

if(PTO_BUILD_PEAKS)
  add_subdirectory(modules/peaks)
  list(APPEND _pto_enabled peaks)
endif()

if(NOT _pto_enabled)
  message(FATAL_ERROR
    "Every module is disabled; there is nothing to build. Enable at least one "
    "of PTO_BUILD_SCRNA_MATRIX / PTO_BUILD_FASTQ_STREAM / "
    "PTO_BUILD_CUTTAG_PROFILER / PTO_BUILD_GENOMIC_TOOLKIT / PTO_BUILD_PEAKS.")
endif()

message(STATUS "pto-core ${PROJECT_VERSION}")
message(STATUS "  build type : ${CMAKE_BUILD_TYPE}")
message(STATUS "  modules    : ${_pto_enabled}")
if(PTO_CLOUD_BUILD)
  message(STATUS "  http       : DISABLED (PTO_CLOUD_BUILD=ON; no listener is compiled in)")
else()
  message(STATUS "  http       : enabled (cuttag_profiler serve)")
endif()
