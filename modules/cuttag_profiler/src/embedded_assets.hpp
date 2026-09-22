// SPDX-License-Identifier: MIT
//
// Declarations for the byte arrays produced by cmake/embed_assets.cmake. The
// definitions live in the generated embedded_assets.cpp.
#pragma once

namespace profiler::assets {

struct Asset {
    const char* path;           // request path, e.g. "/index.html"
    const unsigned char* data;  // not NUL-terminated
    unsigned long size;
    const char* mime;
};

extern const Asset kAssets[];
extern const unsigned long kAssetCount;

}  // namespace profiler::assets
