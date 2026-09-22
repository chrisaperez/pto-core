cmake_minimum_required(VERSION 3.20)

project(genomic_toolkit
        VERSION 0.1.0
        DESCRIPTION "Streaming fragment QC: size distributions, duplicate marking and FRiP in one pass"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(GTK_BUILD_TESTS "Build the test suite" ON)
option(GTK_NATIVE_ARCH "Add -march=native; disable for redistributable binaries" OFF)

# ---------------------------------------------------------------------------
# htslib: optional, unlike cuttag_profiler, where it is the whole premise.
#
# BAM/CRAM input needs it; BEDPE and fragment-BED input do not, and the
# fragment file is what an scATAC user has on disk anyway. Making it optional
# means this module builds on a cluster with a compiler and nothing else --
# which is the deployment constraint the repo exists for -- and the entire test
# suite still runs, because every test drives the text paths.
#
# AUTO: use htslib if pkg-config finds it, otherwise build without BAM support
# and say so. ON: hard-fail if it is missing. OFF: never look.
# ---------------------------------------------------------------------------
set(GTK_WITH_HTSLIB "AUTO" CACHE STRING "BAM/CRAM support: ON, OFF or AUTO")
set_property(CACHE GTK_WITH_HTSLIB PROPERTY STRINGS ON OFF AUTO)

set(GTK_HTSLIB_FOUND FALSE)
if(NOT GTK_WITH_HTSLIB STREQUAL "OFF")
  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND)
    if(GTK_WITH_HTSLIB STREQUAL "ON")
      pkg_check_modules(HTSLIB REQUIRED htslib)
    else()
      pkg_check_modules(HTSLIB QUIET htslib)
    endif()
  endif()
  if(HTSLIB_FOUND)
    set(GTK_HTSLIB_FOUND TRUE)
  elseif(GTK_WITH_HTSLIB STREQUAL "ON")
    message(FATAL_ERROR
      "GTK_WITH_HTSLIB=ON but htslib was not found by pkg-config. Install it "
      "(brew install htslib / apt install libhts-dev), point PKG_CONFIG_PATH at "
      "it, or configure with -DGTK_WITH_HTSLIB=OFF to build the text-input "
      "paths only.")
  endif()
endif()

# ---------------------------------------------------------------------------
# OpenMP: optional here, in contrast to scrna_matrix's hard requirement.
#
# The only parallel region is the FRiP reduction, which is memory-latency
# bound; without OpenMP the same loop runs serially and produces bit-identical
# results (every accumulator is an integer). A hard failure would be a
# fabricated dependency -- the same reasoning the top-level CMakeLists applies
# when it scopes its OpenMP check to scrna_matrix alone.
# ---------------------------------------------------------------------------
if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/../../cmake/PtoOpenMP.cmake")
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../../cmake")
  include(PtoOpenMP OPTIONAL)
  if(COMMAND pto_openmp_apply_hints)
    pto_openmp_apply_hints()
  endif()
endif()
find_package(OpenMP QUIET COMPONENTS CXX)

# ---------------------------------------------------------------------------
# Core library
# ---------------------------------------------------------------------------
add_library(toolkit_core STATIC
            src/fragment_stream.cpp
            src/frip.cpp)

target_include_directories(toolkit_core
                           PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_compile_definitions(toolkit_core PUBLIC GTK_VERSION="${PROJECT_VERSION}")

if(GTK_HTSLIB_FOUND)
  target_include_directories(toolkit_core SYSTEM PRIVATE ${HTSLIB_INCLUDE_DIRS})
  # STATIC_LIBRARY_DIRS/STATIC_LIBRARIES: see fastq_stream/CMakeLists.txt's
  # identical comment. The plain LIBRARIES/LIBRARY_DIRS set is what
  # `pkg-config --libs htslib` reports, and it omits lzma/bz2/deflate -- fine
  # when htslib is a shared library (its own dependency list resolves them),
  # a link failure with a static-only htslib.a, which is what
  # packaging/scripts/build_htslib_manylinux.sh builds. The Homebrew glob is
  # needed alongside it because libdeflate/xz live in their own per-formula
  # keg, not a directory the linker searches by default; it's a no-op on
  # manylinux, which installs everything under /usr/local/lib64 instead.
  file(GLOB GTK_HOMEBREW_OPT_LIBS "/opt/homebrew/opt/*/lib" "/usr/local/opt/*/lib")
  target_link_directories(toolkit_core PUBLIC ${HTSLIB_STATIC_LIBRARY_DIRS} ${GTK_HOMEBREW_OPT_LIBS})
  target_link_libraries(toolkit_core PUBLIC ${HTSLIB_STATIC_LIBRARIES})
  # PUBLIC: the header's `kHaveBam` and the stream_bam declaration are guarded
  # by this, so a consumer compiled without it would see a different API than
  # the library was built with -- an ODR violation that shows up as a link
  # error at best.
  target_compile_definitions(toolkit_core PUBLIC GTK_HAVE_HTSLIB)
  message(STATUS "genomic_toolkit: BAM/CRAM input enabled (htslib ${HTSLIB_VERSION})")
else()
  message(STATUS "genomic_toolkit: BAM/CRAM input DISABLED (no htslib); "
                 "BEDPE and fragment-BED input still available")
endif()

if(OpenMP_CXX_FOUND)
  target_link_libraries(toolkit_core PUBLIC OpenMP::OpenMP_CXX)
  message(STATUS "genomic_toolkit: OpenMP enabled for the FRiP reduction")
else()
  message(STATUS "genomic_toolkit: no OpenMP; FRiP runs single-threaded")
endif()

if(MSVC)
  target_compile_options(toolkit_core PRIVATE /W4 /permissive-)
else()
  target_compile_options(toolkit_core PRIVATE
                         -Wall -Wextra -Wpedantic -Wshadow
                         $<$<CONFIG:Release>:-O3>)
  if(GTK_NATIVE_ARCH)
    target_compile_options(toolkit_core PRIVATE -march=native)
  endif()
endif()

# ---------------------------------------------------------------------------
# Executable
# ---------------------------------------------------------------------------
add_executable(genomic_toolkit src/main.cpp)
target_link_libraries(genomic_toolkit PRIVATE toolkit_core)
install(TARGETS genomic_toolkit RUNTIME DESTINATION bin)

# ---------------------------------------------------------------------------
# Tests
#
# Dependency-free, like every other module's: a bare compiler is enough. Each
# suite is one binary and one add_test, so `ctest -R gtk_peak_set` runs a
# single structure's tests.
# ---------------------------------------------------------------------------
if(GTK_BUILD_TESTS)
  enable_testing()

  foreach(SUITE peak_set dup_marker fragment_stream frip fuzz_parsers file_identity cli_args)
    add_executable(test_gtk_${SUITE} tests/test_${SUITE}.cpp)
    target_link_libraries(test_gtk_${SUITE} PRIVATE toolkit_core)
    target_include_directories(test_gtk_${SUITE} PRIVATE
                               ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    if(NOT MSVC)
      target_compile_options(test_gtk_${SUITE} PRIVATE -Wall -Wextra -Wshadow)
    endif()
    add_test(NAME gtk_${SUITE} COMMAND test_gtk_${SUITE})
  endforeach()

  # The unit suites deliberately carry no htslib (tests/test_util.hpp), so the
  # one property that only an INDEXED BAM can demonstrate -- that --region
  # partitions a file rather than filtering it -- is checked by driving the
  # built binary instead. Skips itself, loudly and successfully, when samtools
  # is absent, the same contract as validate_against_bedtools.sh.
  if(GTK_HTSLIB_FOUND)
    add_test(
      NAME gtk_region_partition
      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_region_partition.sh
              $<TARGET_FILE:genomic_toolkit>)
    # Every other BAM fixture has one contig, where the header's tids and the
    # contig dictionary's ids coincide by accident. This one has two, with the
    # peak file listing them in the opposite order.
    add_test(
      NAME gtk_bam_contig_order
      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_bam_contig_order.sh
              $<TARGET_FILE:genomic_toolkit>)
  endif()

  # src/main.cpp -- argument parsing, exit codes, and the refusals frip and
  # markdup make after streaming -- is linked by no unit suite, so it is driven
  # as a binary. Registered in every configuration: the text cases need only
  # coreutils, and the BAM cases skip themselves without htslib or samtools.
  add_test(
    NAME gtk_cli
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/test_cli.sh $<TARGET_FILE:genomic_toolkit>)
endif()

message(STATUS "genomic_toolkit ${PROJECT_VERSION} (${CMAKE_BUILD_TYPE})")
