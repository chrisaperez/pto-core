cmake_minimum_required(VERSION 3.20)

project(fastq_stream
        VERSION 0.1.0
        DESCRIPTION "Low-latency in-memory FASTQ/BAM streaming and QC engine"
        LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(FQ_NATIVE       "Tune for the build host (-march=native)" ON)
option(FQ_BUILD_TESTS  "Build unit tests and the microbenchmark" ON)
option(FQ_STATIC_LINK  "Link the C++ runtime statically where supported" OFF)
option(FQ_SANITIZE     "Build with ASan+UBSan (developer builds)" OFF)
# htslib is a heavy dependency for a module whose entire current dependency set
# is zlib + libdeflate + pthreads, and a cluster streaming raw Illumina output
# does not need it. OFF by default: with it off, hts_input.cpp is not compiled
# and a BAM input fails with a build-configuration message rather than a link
# error. (cuttag_profiler requires htslib unconditionally; that is correct
# there and would be wrong here.)
option(FQ_ENABLE_HTS   "Enable BAM/CRAM ingestion via htslib" OFF)

# ---------------------------------------------------------------- dependencies
#
# libdeflate does the BGZF block inflate. Preference order:
#   1. an installed CMake package,
#   2. an installed static/shared library found by hand,
#   3. a vendored source build via FetchContent (offline-friendly if the
#      sources are already in _deps, otherwise requires network).
find_package(libdeflate CONFIG QUIET)

if(TARGET libdeflate::libdeflate_static)
  set(FQ_DEFLATE_TARGET libdeflate::libdeflate_static)
elseif(TARGET libdeflate::libdeflate_shared)
  set(FQ_DEFLATE_TARGET libdeflate::libdeflate_shared)
else()
  find_path(LIBDEFLATE_INCLUDE_DIR libdeflate.h)
  find_library(LIBDEFLATE_STATIC NAMES libdeflate.a deflatestatic)
  find_library(LIBDEFLATE_ANY NAMES deflate libdeflate)

  if(LIBDEFLATE_INCLUDE_DIR AND (LIBDEFLATE_STATIC OR LIBDEFLATE_ANY))
    add_library(fq_libdeflate INTERFACE)
    target_include_directories(fq_libdeflate INTERFACE ${LIBDEFLATE_INCLUDE_DIR})
    if(LIBDEFLATE_STATIC)
      target_link_libraries(fq_libdeflate INTERFACE ${LIBDEFLATE_STATIC})
    else()
      target_link_libraries(fq_libdeflate INTERFACE ${LIBDEFLATE_ANY})
    endif()
    set(FQ_DEFLATE_TARGET fq_libdeflate)
  else()
    message(STATUS "libdeflate not found locally; fetching sources")
    include(FetchContent)
    set(LIBDEFLATE_BUILD_SHARED_LIB OFF CACHE BOOL "" FORCE)
    set(LIBDEFLATE_BUILD_GZIP OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(libdeflate
      GIT_REPOSITORY https://github.com/ebiggers/libdeflate.git
      GIT_TAG v1.24)
    FetchContent_MakeAvailable(libdeflate)
    set(FQ_DEFLATE_TARGET libdeflate::libdeflate_static)
  endif()
endif()

# zlib backs the streaming inflate path for single-member gzip, which
# libdeflate structurally cannot handle (see include/fastq_stream/reader.hpp).
find_package(ZLIB REQUIRED)

find_package(Threads REQUIRED)

# -------------------------------------------------------------------- library
add_library(fastq_stream_core STATIC
  src/reader.cpp
  src/pipeline.cpp)

target_include_directories(fastq_stream_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)

# libdeflate/zlib are PUBLIC so the test suite can synthesise BGZF and gzip
# fixtures directly rather than shelling out to bgzip, which is not guaranteed
# to be present on a build host.
target_link_libraries(fastq_stream_core
  PUBLIC Threads::Threads ${FQ_DEFLATE_TARGET} ZLIB::ZLIB)

target_compile_options(fastq_stream_core PRIVATE
  -Wall -Wextra -Wpedantic -Wshadow)

# Vectorisation. -march=native is right for a benchmark host but wrong for a
# redistributable binary, so it is a switch rather than a default assumption.
include(CheckCXXCompilerFlag)
if(FQ_NATIVE)
  check_cxx_compiler_flag("-march=native" FQ_HAS_MARCH_NATIVE)
  if(FQ_HAS_MARCH_NATIVE)
    target_compile_options(fastq_stream_core PUBLIC -march=native)
  else()
    check_cxx_compiler_flag("-mcpu=native" FQ_HAS_MCPU_NATIVE)
    if(FQ_HAS_MCPU_NATIVE)
      target_compile_options(fastq_stream_core PUBLIC -mcpu=native)
    endif()
  endif()
endif()

# AVX-512BW kernels (include/fastq_stream/simd.hpp). Runtime-dispatched like the
# AVX2 ones, but OFF by default -- the same default as PEAKS_ENABLE_AVX512 --
# until they have been executed on AVX-512 hardware: test_trimmer drives them on
# any CPU that reports AVX-512BW, and that run is what should flip this.
option(FQ_ENABLE_AVX512 "Compile runtime-dispatched AVX-512BW kernels" OFF)
if(FQ_ENABLE_AVX512)
  target_compile_definitions(fastq_stream_core PUBLIC FQ_ENABLE_AVX512=1)
endif()

if(FQ_SANITIZE)
  target_compile_options(fastq_stream_core PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer -g)
  target_link_options(fastq_stream_core PUBLIC -fsanitize=address,undefined)
endif()

# ------------------------------------------------------------------- htslib
if(FQ_ENABLE_HTS)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(HTSLIB REQUIRED htslib)

  target_sources(fastq_stream_core PRIVATE src/hts_input.cpp)
  target_compile_definitions(fastq_stream_core PUBLIC FQ_ENABLE_HTS=1)
  # SYSTEM, and PRIVATE: htslib's headers do not survive -Wall -Wextra
  # -Wpedantic, and CI fails on any warning whose path is under modules/*/.
  # Keeping htslib out of the public include set keeps that grep honest.
  target_include_directories(fastq_stream_core SYSTEM PRIVATE ${HTSLIB_INCLUDE_DIRS})
  # Homebrew installs each formula (libdeflate, xz/liblzma) into its own keg
  # under opt/<formula>/lib rather than a directory the linker searches by
  # default -- the same reason cmake/PtoOpenMP.cmake exists for libomp. On
  # manylinux none of this is needed (build_htslib_manylinux.sh installs
  # everything under /usr/local/lib64, already a default search path), so the
  # glob is a no-op there rather than conditional on it.
  file(GLOB FQ_HOMEBREW_OPT_LIBS "/opt/homebrew/opt/*/lib" "/usr/local/opt/*/lib")
  # STATIC_LIBRARY_DIRS/STATIC_LIBRARIES, not the plain ones. FindPkgConfig
  # always populates both -- the plain set is `pkg-config --libs htslib`,
  # which is enough when htslib is a shared library, because a .dylib/.so
  # carries its own dependency list. It is NOT enough when htslib is
  # static-only, which is exactly how packaging/scripts/build_htslib_manylinux.sh
  # builds it: libhts.a has no dependency list of its own, so lzma/bz2/deflate
  # never reach the link line and every symbol they own comes back as an
  # undefined reference. The static set (`pkg-config --static --libs`, i.e.
  # htslib.pc's Libs.private) names them explicitly and is a safe superset on
  # a shared-htslib host too, since those libraries are already installed
  # there as htslib's own build dependencies.
  target_link_directories(fastq_stream_core PUBLIC ${HTSLIB_STATIC_LIBRARY_DIRS} ${FQ_HOMEBREW_OPT_LIBS})
  target_link_libraries(fastq_stream_core PUBLIC ${HTSLIB_STATIC_LIBRARIES})
  message(STATUS "  htslib     : ${HTSLIB_VERSION} (BAM/CRAM ingestion enabled)")
endif()

# ---------------------------------------------------------------- executable
add_executable(fastq_stream src/main.cpp)
target_link_libraries(fastq_stream PRIVATE fastq_stream_core)

if(FQ_STATIC_LINK AND NOT APPLE)
  # Apple does not ship crt0.o for fully static linking; on Linux this yields a
  # single relocatable binary for cluster deployment.
  target_link_options(fastq_stream PRIVATE -static-libstdc++ -static-libgcc)
endif()

# --------------------------------------------------------------------- tests
if(FQ_BUILD_TESTS)
  enable_testing()

  add_executable(test_trimmer tests/test_trimmer.cpp)
  target_link_libraries(test_trimmer PRIVATE fastq_stream_core)
  add_test(NAME trimmer COMMAND test_trimmer)

  # The binary as a black box: argument validation, report files and
  # non-regular-file inputs, none of which the library tests can reach.
  add_test(NAME fastq_stream_cli
           COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cli.sh $<TARGET_FILE:fastq_stream>)

  add_executable(benchmark tests/benchmark.cpp)
  target_link_libraries(benchmark PRIVATE fastq_stream_core)

  if(FQ_ENABLE_HTS)
    add_executable(test_hts_input tests/test_hts_input.cpp)
    target_link_libraries(test_hts_input PRIVATE fastq_stream_core)
    target_include_directories(test_hts_input SYSTEM PRIVATE ${HTSLIB_INCLUDE_DIRS})
    add_test(NAME hts_input COMMAND test_hts_input)
  endif()
endif()

install(TARGETS fastq_stream RUNTIME DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)

message(STATUS "fastq_stream ${PROJECT_VERSION}")
message(STATUS "  build type : ${CMAKE_BUILD_TYPE}")
message(STATUS "  libdeflate : ${FQ_DEFLATE_TARGET}")
message(STATUS "  zlib       : ${ZLIB_VERSION_STRING}")
