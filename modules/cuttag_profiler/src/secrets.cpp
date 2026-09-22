// SPDX-License-Identifier: MIT
#include "profiler/secrets.hpp"

#include <cstdio>
#include <stdexcept>
#include <vector>

#if defined(_WIN32)
#include <windows.h>
// bcrypt.h must follow windows.h.
#include <bcrypt.h>
#elif defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || \
    defined(__NetBSD__)
#include <stdlib.h>  // arc4random_buf
#else
#include <cerrno>
#include <fstream>
#include <sys/random.h>  // getrandom
#endif

namespace profiler {
namespace {

// Fills `out` with CSPRNG bytes or throws. One function per platform, because
// there is no portable C++ CSPRNG: std::random_device is permitted to be a
// deterministic PRNG (and famously was, for years, on MinGW), so it cannot be
// used for a secret.
void fill_random(unsigned char* out, std::size_t count) {
#if defined(_WIN32)
    // The system-preferred RNG; needs no algorithm handle and cannot be
    // configured into a weaker mode by the caller.
    if (::BCryptGenRandom(nullptr, out, static_cast<ULONG>(count),
                          BCRYPT_USE_SYSTEM_PREFERRED_RNG) != 0) {
        throw std::runtime_error("cannot obtain secure random bytes");
    }
#elif defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || \
    defined(__NetBSD__)
    // arc4random_buf is the platform CSPRNG here, is always seeded before
    // main() and has no failure mode -- there is no error to check, which is
    // why it is preferred over getentropy() (capped at 256 bytes per call) and
    // over opening /dev/urandom (which can fail under a restrictive sandbox or
    // an exhausted fd table).
    ::arc4random_buf(out, count);
#else
    // Linux. getrandom() reads the same pool as /dev/urandom without needing a
    // file descriptor, so it works when the fd table is full and inside a
    // container with no /dev. It can return a short read and can be
    // interrupted, so it is looped.
    std::size_t filled = 0;
    while (filled < count) {
        const ssize_t got = ::getrandom(out + filled, count - filled, 0);
        if (got < 0) {
            if (errno == EINTR) continue;
            break;  // ENOSYS on a pre-3.17 kernel: fall through to the file
        }
        filled += static_cast<std::size_t>(got);
    }
    if (filled == count) return;

    // Fallback for kernels or seccomp filters without getrandom().
    std::ifstream urandom("/dev/urandom", std::ios::binary);
    if (!urandom.read(reinterpret_cast<char*>(out + filled),
                      static_cast<std::streamsize>(count - filled))) {
        throw std::runtime_error("cannot obtain secure random bytes");
    }
#endif
}

}  // namespace

std::string secure_token(std::size_t bytes) {
    if (bytes == 0) {
        throw std::invalid_argument("secure_token: zero-length token requested");
    }

    std::vector<unsigned char> raw(bytes);
    fill_random(raw.data(), raw.size());

    // Hand-rolled rather than ostringstream << std::hex: that is the exact
    // formatting that dropped leading zeros in the implementation this
    // replaces. Two characters per byte, always.
    static constexpr char kHex[] = "0123456789abcdef";
    std::string out;
    out.reserve(raw.size() * 2);
    for (unsigned char c : raw) {
        out.push_back(kHex[c >> 4]);
        out.push_back(kHex[c & 0x0F]);
    }
    return out;
}

}  // namespace profiler
