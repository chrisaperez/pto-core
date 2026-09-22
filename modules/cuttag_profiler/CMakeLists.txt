cmake_minimum_required(VERSION 3.20)

project(cuttag_profiler
        VERSION 1.0.0
        DESCRIPTION "Air-gapped CUT&Tag/CUT&RUN signal profiler with an embedded dashboard"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

option(PROFILER_BUILD_TESTS "Build the test suite" ON)
option(PROFILER_STATIC_HTSLIB
       "Link libhts.a rather than the shared library" ON)
option(PROFILER_NATIVE_ARCH
       "Add -march=native; disable for redistributable binaries" OFF)

# Cloud-worker build: strip the embedded HTTP server (SECURITY_HTTP H1).
#
# Declared here as well as at the pto-core top level because this module is
# still buildable standalone; when configured from the top level that call
# already cached the variable and this one is a no-op, so the two cannot drift.
#
# The strip is structural, not conditional-at-runtime. Below, ON removes
# src/http_server.cpp from the library, skips the dashboard-asset embedding
# entirely, and does not link the vendored httplib/nlohmann headers -- so the
# artefact contains no HTTP server to reach, whatever argv it is given.
# main.cpp's `serve` subcommand is compiled out to match (#ifndef
# PTO_CLOUD_BUILD), and ctest gains `cloud_no_server`, which runs the built
# binary and pins that it refuses.
option(PTO_CLOUD_BUILD
       "Cloud-worker build: strip the embedded HTTP server from the binary"
       OFF)

# ---------------------------------------------------------------------------
# htslib
#
# Prefer the static archive so the shipped binary has no runtime .so/.dylib
# dependency beyond the platform libc. pkg-config supplies the transitive
# compression and networking libraries htslib itself was built against.
# ---------------------------------------------------------------------------
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)

if(PROFILER_STATIC_HTSLIB)
  find_library(HTSLIB_STATIC_LIBRARY
               NAMES libhts.a hts
               HINTS ${HTSLIB_LIBRARY_DIRS} ${HTSLIB_LIBDIR}
               NO_DEFAULT_PATH)
  if(NOT HTSLIB_STATIC_LIBRARY)
    find_library(HTSLIB_STATIC_LIBRARY NAMES libhts.a)
  endif()
endif()

add_library(htslib::htslib INTERFACE IMPORTED)
target_include_directories(htslib::htslib SYSTEM INTERFACE ${HTSLIB_INCLUDE_DIRS})

if(PROFILER_STATIC_HTSLIB AND HTSLIB_STATIC_LIBRARY)
  # `--static` pulls in zlib/bzip2/lzma/libdeflate/curl/crypto, which the
  # archive needs. htslib's .pc names some of these with a bare `-lfoo` and no
  # matching `-L`, so every entry is resolved to an absolute path here; that
  # also lets us prefer each dependency's static archive, which is the whole
  # point of a single-binary deployment.
  execute_process(COMMAND ${PKG_CONFIG_EXECUTABLE} --static --libs htslib
                  OUTPUT_VARIABLE HTSLIB_STATIC_LDFLAGS_RAW
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_QUIET)
  separate_arguments(HTSLIB_STATIC_LDFLAGS_LIST NATIVE_COMMAND
                     "${HTSLIB_STATIC_LDFLAGS_RAW}")

  # Pass 1: collect every -L directory as a search hint.
  set(HTSLIB_SEARCH_DIRS ${HTSLIB_LIBRARY_DIRS})
  foreach(FLAG IN LISTS HTSLIB_STATIC_LDFLAGS_LIST)
    if(FLAG MATCHES "^-L(.+)$")
      list(APPEND HTSLIB_SEARCH_DIRS "${CMAKE_MATCH_1}")
    endif()
  endforeach()
  # Homebrew keeps each formula in its own cellar; opt/ holds the stable links.
  file(GLOB HOMEBREW_OPT_LIBS "/opt/homebrew/opt/*/lib" "/usr/local/opt/*/lib")
  list(APPEND HTSLIB_SEARCH_DIRS ${HOMEBREW_OPT_LIBS})
  list(REMOVE_DUPLICATES HTSLIB_SEARCH_DIRS)

  # Pass 2: turn each -lfoo into a concrete file where one can be found.
  set(HTSLIB_RESOLVED_DEPS "")
  foreach(FLAG IN LISTS HTSLIB_STATIC_LDFLAGS_LIST)
    if(FLAG MATCHES "^-L")
      continue()  # superseded by the absolute paths below
    elseif(FLAG STREQUAL "-lhts")
      continue()  # the archive is named explicitly
    elseif(FLAG STREQUAL "-lpthread" OR FLAG STREQUAL "-pthread")
      continue()  # handled by Threads::Threads
    elseif(FLAG MATCHES "^-l(.+)$")
      set(DEP_NAME "${CMAKE_MATCH_1}")

      # Libraries that must NEVER be resolved to a static archive.
      #
      # The previous code relied on find_library FAILING for these, with the
      # comment "probably a libc-provided library (-lm)". That assumption is
      # false on glibc: Ubuntu ships /usr/lib/x86_64-linux-gnu/libm.a, so -lm
      # resolved to the static archive and got linked into a PIE. The result is
      #
      #     undefined reference to `_dl_x86_cpu_features'
      #
      # because static libm's ifunc resolvers reference a private loader symbol
      # that is only available to the dynamic linker. Reproduced on
      # ubuntu:24.04 with the full htslib dep chain installed, and it is the
      # failure that took down every CI build job.
      #
      # These are part of the platform runtime, not of htslib's dependency
      # closure -- bundling them serves no deployment purpose and breaks the
      # link. An explicit list beats hoping a lookup fails.
      if(DEP_NAME MATCHES "^(m|c|dl|rt|resolv|nsl|util|gcc|gcc_s|stdc\\+\\+)$")
        list(APPEND HTSLIB_RESOLVED_DEPS "${FLAG}")
        continue()
      endif()

      unset(DEP_PATH CACHE)
      find_library(DEP_PATH
                   NAMES "lib${DEP_NAME}.a" "${DEP_NAME}"
                   HINTS ${HTSLIB_SEARCH_DIRS})
      if(DEP_PATH)
        list(APPEND HTSLIB_RESOLVED_DEPS "${DEP_PATH}")
      else()
        # No archive found; leave the flag for the linker to resolve dynamically.
        list(APPEND HTSLIB_RESOLVED_DEPS "${FLAG}")
      endif()
      unset(DEP_PATH CACHE)
    else()
      list(APPEND HTSLIB_RESOLVED_DEPS "${FLAG}")
    endif()
  endforeach()
  # htslib's archive always carries its hfile_libcurl / hfile_s3 plugins, but
  # brew's .pc omits libcurl from --static --libs. Those code paths are only
  # reachable through an http(s):// or s3:// filename, which this tool never
  # constructs; the symbols still have to resolve. To ship a binary that
  # provably contains no HTTP client at all, rebuild htslib with
  # scripts/build_htslib_minimal.sh and point CMAKE_PREFIX_PATH at it.
  # Static archives are tried first for the same reason as above: a resolved
  # .dylib under /opt/homebrew pins the binary to a Homebrew install and
  # defeats the copy-and-run deployment story.
  foreach(EXTRA_DEP curl crypto ssl)
    unset(EXTRA_PATH CACHE)
    if(APPLE AND EXTRA_DEP STREQUAL "curl")
      # macOS ships /usr/lib/libcurl.4.dylib, whose dependency closure is
      # already closed (see pyproject.toml's [tool.cibuildwheel.macos]
      # comment). Homebrew's curl -- pulled onto this search path as an
      # htslib dependency -- is a static archive built with HTTP/3 support,
      # needing ngtcp2/libpsl/libidn2 this list does not resolve; preferring
      # it here breaks the link with undefined ngtcp2/psl/uidna symbols.
      # Leave -lcurl for the linker, which finds the system dylib instead.
      list(APPEND HTSLIB_RESOLVED_DEPS "-lcurl")
      continue()
    endif()
    find_library(EXTRA_PATH
                 NAMES "lib${EXTRA_DEP}.a" "${EXTRA_DEP}"
                 HINTS ${HTSLIB_SEARCH_DIRS})
    if(EXTRA_PATH)
      list(APPEND HTSLIB_RESOLVED_DEPS "${EXTRA_PATH}")
    endif()
    unset(EXTRA_PATH CACHE)
  endforeach()

  list(REMOVE_DUPLICATES HTSLIB_RESOLVED_DEPS)

  target_link_libraries(htslib::htslib INTERFACE
                        ${HTSLIB_STATIC_LIBRARY}
                        ${HTSLIB_RESOLVED_DEPS})
  message(STATUS "htslib: static (${HTSLIB_STATIC_LIBRARY})")
  message(STATUS "htslib deps: ${HTSLIB_RESOLVED_DEPS}")
else()
  target_link_directories(htslib::htslib INTERFACE ${HTSLIB_LIBRARY_DIRS})
  target_link_libraries(htslib::htslib INTERFACE ${HTSLIB_LIBRARIES})
  message(STATUS "htslib: shared (${HTSLIB_LIBRARIES}) - binary will not be self-contained")
endif()

# ---------------------------------------------------------------------------
# Vendored header-only dependencies
#
# Checked into third_party/ on purpose: a build on an air-gapped machine must
# not reach out to a package registry.
# ---------------------------------------------------------------------------
add_library(vendored_headers INTERFACE)
target_include_directories(vendored_headers SYSTEM INTERFACE
                           ${CMAKE_CURRENT_SOURCE_DIR}/third_party/httplib
                           ${CMAKE_CURRENT_SOURCE_DIR}/third_party)

# ---------------------------------------------------------------------------
# Frontend assets compiled into the binary
#
# Skipped entirely under PTO_CLOUD_BUILD: with no server to serve them the
# dashboard's HTML/CSS/JS is dead weight in the image, and leaving it out is
# also the cheapest way to see at a glance which build you have (the cloud
# binary is materially smaller).
# ---------------------------------------------------------------------------
if(PTO_CLOUD_BUILD)
  set(GENERATED_ASSETS "")
else()
  set(GENERATED_ASSETS ${CMAKE_CURRENT_BINARY_DIR}/generated/embedded_assets.cpp)
  file(GLOB_RECURSE WEB_ASSET_SOURCES CONFIGURE_DEPENDS
       ${CMAKE_CURRENT_SOURCE_DIR}/web/*.html
       ${CMAKE_CURRENT_SOURCE_DIR}/web/*.css
       ${CMAKE_CURRENT_SOURCE_DIR}/web/*.js)
  list(FILTER WEB_ASSET_SOURCES EXCLUDE REGEX "/web/src/")

  add_custom_command(
    OUTPUT ${GENERATED_ASSETS}
    COMMAND ${CMAKE_COMMAND}
            -DASSET_DIR=${CMAKE_CURRENT_SOURCE_DIR}/web
            -DOUTPUT_CPP=${GENERATED_ASSETS}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/embed_assets.cmake
    DEPENDS ${WEB_ASSET_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/embed_assets.cmake
    COMMENT "Embedding dashboard assets into the binary"
    VERBATIM)
  add_custom_target(generate_assets DEPENDS ${GENERATED_ASSETS})
endif()

# ---------------------------------------------------------------------------
# Core library
# ---------------------------------------------------------------------------
# secrets.cpp and safe_open.cpp are compiled in BOTH configurations on
# purpose. The constant-time compare is a general primitive, and the
# open-once/verify-the-descriptor helper closes M6 for the batch `profile` path
# as much as for the HTTP one -- the cloud worker runs `profile` over files a
# tenant uploaded, so the FIFO and symlink cases matter there too. Keeping them
# unconditional also keeps tests/test_secrets.cpp and tests/test_region_names.cpp
# running in the cloud configuration, which is where they would otherwise stop
# being run at all.
set(PROFILER_SOURCES
    src/bam_reader.cpp
    src/bed_reader.cpp
    src/safe_open.cpp
    src/secrets.cpp
    src/signal_calc.cpp)

if(NOT PTO_CLOUD_BUILD)
  list(APPEND PROFILER_SOURCES src/http_server.cpp ${GENERATED_ASSETS})
endif()

add_library(profiler_core STATIC ${PROFILER_SOURCES})

target_include_directories(profiler_core
                           PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(profiler_core PUBLIC htslib::htslib Threads::Threads)
target_compile_definitions(profiler_core
                           PUBLIC PROFILER_VERSION="${PROJECT_VERSION}")

if(PTO_CLOUD_BUILD)
  # PUBLIC so main.cpp compiles out `serve` in step with the library, and so
  # any future consumer of profiler_core sees the same configuration. Getting
  # this wrong in the other direction -- library stripped, main.cpp still
  # calling HttpServer -- is a link error, not a silent listener, which is the
  # failure mode to prefer.
  target_compile_definitions(profiler_core PUBLIC PTO_CLOUD_BUILD)
  # httplib and nlohmann/json are reachable ONLY from http_server.cpp. Not
  # linking the vendored headers at all is what makes "there is no HTTP server
  # in this binary" a property of the build graph rather than of a #ifdef
  # someone has to keep honest.
else()
  target_link_libraries(profiler_core PRIVATE vendored_headers)
  target_compile_definitions(profiler_core PRIVATE CPPHTTPLIB_THREAD_POOL_COUNT=16)
endif()

if(MSVC)
  target_compile_options(profiler_core PRIVATE /W4 /permissive-)
else()
  target_compile_options(profiler_core PRIVATE
                         -Wall -Wextra -Wpedantic -Wshadow
                         $<$<CONFIG:Release>:-O3>)
  if(PROFILER_NATIVE_ARCH)
    target_compile_options(profiler_core PRIVATE -march=native)
  endif()
endif()

# ---------------------------------------------------------------------------
# Executable
# ---------------------------------------------------------------------------
add_executable(cuttag_profiler src/main.cpp)
target_link_libraries(cuttag_profiler PRIVATE profiler_core)

if(APPLE)
  # macOS has no fully static libSystem; everything else is archived in.
  target_link_libraries(cuttag_profiler PRIVATE "-framework CoreFoundation")
elseif(UNIX AND PROFILER_STATIC_HTSLIB)
  target_link_options(cuttag_profiler PRIVATE -static-libstdc++ -static-libgcc)
endif()

install(TARGETS cuttag_profiler RUNTIME DESTINATION bin)

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(PROFILER_BUILD_TESTS)
  enable_testing()

  add_executable(test_bam_query tests/test_bam_query.cpp)
  target_link_libraries(test_bam_query PRIVATE profiler_core)
  target_include_directories(test_bam_query PRIVATE
                             ${CMAKE_CURRENT_SOURCE_DIR}/tests)
  add_test(NAME bam_query COMMAND test_bam_query)

  add_executable(test_signal_calc tests/test_signal_calc.cpp)
  target_link_libraries(test_signal_calc PRIVATE profiler_core)
  target_include_directories(test_signal_calc PRIVATE
                             ${CMAKE_CURRENT_SOURCE_DIR}/tests)
  add_test(NAME signal_calc COMMAND test_signal_calc)

  # Matrix-cache eviction. Header-only and htslib-free. Registered in BOTH
  # configurations even though only the dashboard populates the cache: the
  # bound it enforces is a memory-safety property, and matrix_bytes' saturating
  # arithmetic is exactly the kind of thing that rots unnoticed if its test is
  # compiled out of the configuration that ships.
  add_executable(test_matrix_cache tests/test_matrix_cache.cpp)
  target_link_libraries(test_matrix_cache PRIVATE profiler_core)
  target_include_directories(test_matrix_cache PRIVATE
                             ${CMAKE_CURRENT_SOURCE_DIR}/tests)
  add_test(NAME matrix_cache COMMAND test_matrix_cache)

  # Interval overlap queries. Header-only and htslib-free -- it links
  # profiler_core only for the shared harness.
  add_executable(test_interval_index tests/test_interval_index.cpp)
  target_link_libraries(test_interval_index PRIVATE profiler_core)
  target_include_directories(test_interval_index PRIVATE
                             ${CMAKE_CURRENT_SOURCE_DIR}/tests)
  add_test(NAME interval_index COMMAND test_interval_index)

  # Token minting and constant-time comparison (SECURITY_HTTP L9). Registered
  # in BOTH configurations: the primitives are compiled into the cloud build
  # too, and a security test that silently stops running in the configuration
  # that ships is worse than no test.
  add_executable(test_secrets tests/test_secrets.cpp)
  target_link_libraries(test_secrets PRIVATE profiler_core)
  target_include_directories(test_secrets PRIVATE
                             ${CMAKE_CURRENT_SOURCE_DIR}/tests)
  add_test(NAME secrets COMMAND test_secrets)

  # Region-name sanitisation (SECURITY_HTTP M5) and the open-once path
  # confinement helper (M6). Also both configurations: row names reach a TSV in
  # the batch `profile` path, which is precisely what the cloud worker runs, so
  # the TSV-injection half of M5 is not an HTTP-only concern.
  add_executable(test_region_names tests/test_region_names.cpp)
  target_link_libraries(test_region_names PRIVATE profiler_core)
  target_include_directories(test_region_names PRIVATE
                             ${CMAKE_CURRENT_SOURCE_DIR}/tests)
  add_test(NAME region_names COMMAND test_region_names)

  add_executable(test_bed_reader tests/test_bed_reader.cpp)
  target_link_libraries(test_bed_reader PRIVATE profiler_core)
  target_include_directories(test_bed_reader PRIVATE
                             ${CMAKE_CURRENT_SOURCE_DIR}/tests)
  add_test(NAME bed_reader COMMAND test_bed_reader)

  # The output-overwrites-input guard (include/profiler/file_identity.hpp).
  add_executable(test_file_identity tests/test_file_identity.cpp)
  target_link_libraries(test_file_identity PRIVATE profiler_core)
  target_include_directories(test_file_identity PRIVATE
                             ${CMAKE_CURRENT_SOURCE_DIR}/tests)
  add_test(NAME file_identity COMMAND test_file_identity)

  # `profile` as a black box: exit codes and option refusals that live in
  # src/main.cpp, which no unit suite links. Both configurations -- the cloud
  # worker runs exactly this subcommand. Writes its own BAM with samtools
  # (the demo data is untracked) and skips loudly without it.
  add_test(NAME profiler_cli
           COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/test_cli.sh
                   $<TARGET_FILE:cuttag_profiler>)

  # The structural claim, asserted against the linked artefact rather than
  # against the source. A #ifndef says what was compiled; nm says what is in
  # the file, and it is the latter that has to be true. Registered in both
  # configurations so the assertion runs in the direction that matters for
  # each: no server in a cloud build, and a server still present in an
  # ordinary one (a stray -DPTO_CLOUD_BUILD=ON would otherwise read as a
  # runtime bug in the dashboard).
  add_test(NAME server_linkage
           COMMAND ${CMAKE_COMMAND}
                   -DBINARY=$<TARGET_FILE:cuttag_profiler>
                   -DEXPECT_SERVER=$<NOT:$<BOOL:${PTO_CLOUD_BUILD}>>
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/assert_server_linkage.cmake)

  if(PTO_CLOUD_BUILD)
    # The behavioural half of the same claim: this binary has no `serve`.
    # PASS_REGULAR_EXPRESSION alone would pass on a crash that happened to
    # print the string, so the exit code is pinned too.
    add_test(NAME cloud_no_server
             COMMAND cuttag_profiler serve --data-root .)
    set_tests_properties(cloud_no_server PROPERTIES
                         PASS_REGULAR_EXPRESSION "built without the HTTP server"
                         WILL_FAIL FALSE)
    # `serve` must exit non-zero. WILL_FAIL inverts the whole result, so the
    # exit code is checked in its own test rather than combined with the above.
    add_test(NAME cloud_no_server_exit_code
             COMMAND cuttag_profiler serve)
    set_tests_properties(cloud_no_server_exit_code PROPERTIES WILL_FAIL TRUE)

    # And the dashboard must not be advertised in a build that cannot run it.
    add_test(NAME cloud_help_omits_serve COMMAND cuttag_profiler --help)
    set_tests_properties(cloud_help_omits_serve PROPERTIES
                         FAIL_REGULAR_EXPRESSION "SERVE OPTIONS")
  else()
    # Host/Origin validation, the session-token gate and the confinement
    # configuration checks (SECURITY_HTTP H1, H2, H3). Needs ServerOptions,
    # which lives in http_server.cpp and is absent from a cloud build.
    add_executable(test_http_security tests/test_http_security.cpp)
    target_link_libraries(test_http_security PRIVATE profiler_core)
    target_include_directories(test_http_security PRIVATE
                               ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    add_test(NAME http_security COMMAND test_http_security)
  endif()
endif()

message(STATUS "cuttag_profiler ${PROJECT_VERSION} (${CMAKE_BUILD_TYPE})")
