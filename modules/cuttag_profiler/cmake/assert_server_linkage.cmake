# SPDX-License-Identifier: MIT
#
# Asserts, against the linked artefact, whether the embedded HTTP server is in
# it. Run by ctest in both configurations:
#
#   -DEXPECT_SERVER=OFF   the PTO_CLOUD_BUILD binary must contain no server
#   -DEXPECT_SERVER=ON    the ordinary binary must still contain one
#
# Both directions matter. The first is the security claim; the second is what
# stops a stray -DPTO_CLOUD_BUILD=ON in a workstation build turning the
# dashboard off silently, which would otherwise look like a runtime bug.
#
# This is the check the source-level #ifndef cannot make. A preprocessor guard
# is a statement about what was compiled; this is a statement about what is in
# the file, and it is the one that survives someone reorganising the guards.
#
# Deliberately narrow: it matches only symbols this repository is responsible
# for (cpp-httplib's namespace and profiler::HttpServer). It does NOT assert
# the absence of socket(2)/bind(2), tempting as that is -- htslib is linked
# statically and carries hfile_libcurl/hfile_s3, so whether those imports
# appear depends on how htslib itself was built. A security gate that fails on
# an unrelated packaging choice gets disabled, and then it protects nothing.

if(NOT DEFINED BINARY OR NOT DEFINED EXPECT_SERVER)
  message(FATAL_ERROR "usage: -DBINARY=<path> -DEXPECT_SERVER=ON|OFF")
endif()

find_program(NM_EXECUTABLE NAMES nm llvm-nm)
if(NOT NM_EXECUTABLE)
  # Reported as a skip, not a pass. CTest maps this through
  # SKIP_RETURN_CODE so a toolchain without nm does not look green.
  message(STATUS "nm not found; cannot inspect ${BINARY}")
  return()
endif()

execute_process(COMMAND ${NM_EXECUTABLE} "${BINARY}"
                OUTPUT_VARIABLE SYMBOLS
                ERROR_VARIABLE NM_ERROR
                RESULT_VARIABLE NM_STATUS)
if(NOT NM_STATUS EQUAL 0)
  message(FATAL_ERROR "nm failed on ${BINARY}: ${NM_ERROR}")
endif()

# A stripped release binary keeps no local symbol names, which would make the
# ON direction fail for a reason that has nothing to do with the server.
string(LENGTH "${SYMBOLS}" SYMBOL_TEXT_LENGTH)
if(SYMBOL_TEXT_LENGTH LESS 100)
  message(STATUS "${BINARY} appears stripped; nothing to inspect")
  return()
endif()

set(FOUND FALSE)
foreach(PATTERN "httplib" "HttpServer")
  string(FIND "${SYMBOLS}" "${PATTERN}" HIT)
  if(NOT HIT EQUAL -1)
    set(FOUND TRUE)
    set(FOUND_PATTERN "${PATTERN}")
  endif()
endforeach()

if(EXPECT_SERVER AND NOT FOUND)
  message(FATAL_ERROR
    "${BINARY} contains no HTTP server symbols, but this build did not ask to "
    "strip it. Something set PTO_CLOUD_BUILD=ON, and `serve` will not work.")
elseif(NOT EXPECT_SERVER AND FOUND)
  message(FATAL_ERROR
    "PTO_CLOUD_BUILD=ON, but ${BINARY} still links the HTTP server "
    "(matched '${FOUND_PATTERN}').\n"
    "The point of that option is that the artefact cannot expose a port by "
    "any argv or config path. A #ifndef that no longer guards the whole of "
    "http_server.cpp, or a source added back to the library's list, would "
    "produce exactly this.")
endif()

if(EXPECT_SERVER)
  message(STATUS "${BINARY}: HTTP server present, as expected")
else()
  message(STATUS "${BINARY}: no HTTP server symbols -- strip verified")
endif()
