// SPDX-License-Identifier: MIT
//
// The same deliberately tiny harness the other modules use. A framework would
// mean a network fetch at configure time or another vendored tree, and the
// project's premise is that a build works on a machine with no route out.
//
// One addition over genomic_toolkit's version: an allocation counter, because
// two of this module's guarantees are about what does NOT happen. The ring
// buffers are sized once and the streaming loop must not allocate -- a claim
// that no output value can confirm, and that a `reserve()` removed by a later
// refactor would silently break. Counting global operator new is the only
// check that actually observes it.
//
// The operator replacements below are deliberately NOT inline: a replacement
// operator new may not be. Each test binary in this module is a single
// translation unit, so including this header once per binary is well formed.
#pragma once

// Whether this translation unit may replace the global allocation functions.
//
// It may not under a sanitizer, and the two sanitizers object for different
// reasons -- both of which were live CI failures:
//
//   * TSan's libclang_rt.tsan_cxx defines operator new/delete as STRONG
//     symbols in a static archive. A second definition is a link error:
//     "multiple definition of `operator new(unsigned long)'".
//   * ASan's are weak, so ours win -- but only for the part of the family we
//     name. In a normal build that is harmless: the default
//     `operator new(size_t, nothrow_t)` is specified to CALL
//     `operator new(size_t)`, so it reaches the replacement below and the
//     malloc/free pairing holds. ASan supplies its own nothrow new instead,
//     which does not delegate, so libstdc++'s _Temporary_buffer (what
//     std::stable_sort allocates through) took its buffer from ASan's
//     allocator and returned it to the free() below --
//     "alloc-dealloc-mismatch (operator new vs free)".
//
// Replacing the whole twenty-function family instead would fix the second and
// not the first, and its next omission would be this bug again. So under a
// sanitizer the operators are left alone and counting is off; `counting()`
// reports that, and the cases that depend on it announce a skip rather than
// passing over nothing. The property stays covered by every non-sanitizer job
// in the matrix.
#if defined(__has_feature)
#if __has_feature(address_sanitizer) || __has_feature(thread_sanitizer) || \
    __has_feature(memory_sanitizer)
#define PEAKS_TEST_NO_ALLOC_COUNTING 1
#endif
#endif
#if defined(__SANITIZE_ADDRESS__) || defined(__SANITIZE_THREAD__)
#define PEAKS_TEST_NO_ALLOC_COUNTING 1
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <functional>
#include <new>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace testing {

// long long, not int: a sweep drives millions of checks and a signed int
// counter overflows, which is undefined behaviour and prints nonsense.
inline long long g_failures = 0;
inline long long g_checks = 0;

// Bumped by every global allocation. `AllocationGuard` reads it.
inline long long g_allocations = 0;
inline bool g_count_allocations = false;

template <typename T>
std::string describe(const T& value) {
    if constexpr (std::is_convertible_v<T, std::string>) {
        return std::string(value);
    } else if constexpr (std::is_floating_point_v<T>) {
        char buffer[40];
        std::snprintf(buffer, sizeof(buffer), "%.10g", static_cast<double>(value));
        return buffer;
    } else {
        return std::to_string(value);
    }
}

inline void report(bool ok, const char* expression, const char* file, int line,
                   const std::string& detail = {}) {
    ++g_checks;
    if (ok) return;
    ++g_failures;
    std::fprintf(stderr, "  FAIL %s:%d: %s%s%s\n", file, line, expression,
                 detail.empty() ? "" : " -- ", detail.c_str());
}

#define CHECK(expr) ::testing::report((expr), #expr, __FILE__, __LINE__)

#define CHECK_MSG(expr, detail) \
    ::testing::report((expr), #expr, __FILE__, __LINE__, (detail))

#define CHECK_EQ(a, b)                                              \
    ::testing::report((a) == (b), #a " == " #b, __FILE__, __LINE__, \
                      ::testing::describe(a) + " vs " + ::testing::describe(b))

#define CHECK_NEAR(a, b, tol)                                         \
    ::testing::report(std::fabs((a) - (b)) <= (tol), #a " ~= " #b,    \
                      __FILE__, __LINE__,                             \
                      ::testing::describe(a) + " vs " + ::testing::describe(b))

// Asserts that `body` throws `E`. Failure to throw is as much a bug as the
// wrong result: several of this module's guarantees are "refuses to run".
#define CHECK_THROWS(E, body)                                        \
    do {                                                             \
        bool threw_ = false;                                         \
        try {                                                        \
            body;                                                    \
        } catch (const E&) {                                         \
            threw_ = true;                                           \
        } catch (...) {                                              \
        }                                                            \
        ::testing::report(threw_, "throws " #E, __FILE__, __LINE__); \
    } while (0)

// Counts allocations over its lifetime. Nothing inside a streaming loop may
// allocate: the ring buffers are sized once, at reset, and a per-position
// allocation would put a lock and a syscall-shaped latency into the hot path.
class AllocationGuard {
public:
    AllocationGuard() : before_(g_allocations) { g_count_allocations = true; }
    ~AllocationGuard() { g_count_allocations = false; }
    [[nodiscard]] long long count() const { return g_allocations - before_; }

    // False when the global operators were not replaced, which is every
    // sanitizer build -- see the header comment. A caller MUST branch on this:
    // with no replacement the counter never moves, so `count() == 0` is true
    // whatever the code under test did, and the check would pass over nothing.
    [[nodiscard]] static bool counting() noexcept {
#if defined(PEAKS_TEST_NO_ALLOC_COUNTING)
        return false;
#else
        return true;
#endif
    }

private:
    long long before_;
};

// A case that could not run. Printed, never counted as a check -- the line is
// the point, so a skipped guarantee cannot be mistaken for a kept one.
inline void skip(const char* why) { std::printf("      SKIP %s\n", why); }

struct Suite {
    const char* name;
    std::vector<std::pair<const char*, std::function<void()>>> cases;

    void add(const char* case_name, std::function<void()> body) {
        cases.emplace_back(case_name, std::move(body));
    }

    int run() {
        std::printf("== %s ==\n", name);
        for (auto& [case_name, body] : cases) {
            const long long before = g_failures;
            body();
            std::printf("  %s %s\n", g_failures == before ? "ok  " : "FAIL", case_name);
        }
        std::printf("%lld checks, %lld failure(s)\n", g_checks, g_failures);
        return g_failures == 0 ? 0 : 1;
    }
};

}  // namespace testing

#if !defined(PEAKS_TEST_NO_ALLOC_COUNTING)

// GCC cannot see that the four functions below are a matched malloc/free pair
// -- it models the replaced global operator new as the builtin, then reports
// every library free() of a pointer that came from it as mismatched. The
// pairing is correct by construction (every allocation here comes from
// std::malloc and every deallocation goes to std::free), and the diagnostic
// fires at the INLINED call site rather than at the definition, so it has to
// be silenced for the translation unit rather than around the definitions.
// Scoped to GCC: clang's own -Wmismatched-new-delete does not fire on this.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif

void* operator new(std::size_t n) {
    if (::testing::g_count_allocations) ++::testing::g_allocations;
    void* p = std::malloc(n == 0 ? 1 : n);
    if (p == nullptr) throw std::bad_alloc();
    return p;
}

void* operator new[](std::size_t n) { return ::operator new(n); }

void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }

#endif  // !PEAKS_TEST_NO_ALLOC_COUNTING
