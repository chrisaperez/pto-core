// SPDX-License-Identifier: MIT
//
// The same deliberately tiny harness the other modules use. A framework would
// mean a network fetch at configure time or another vendored tree, and the
// project's premise is that a build works on a machine with no route out.
//
// Unlike cuttag_profiler's version this one includes no htslib: every suite
// here drives the text paths, so the tests build and run in a configuration
// with BAM support switched off.
#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <random>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace testing {

// long long, not int: a fuzz suite drives millions of checks per run, and a
// signed int counter overflows (undefined behaviour, and a nonsense total in
// the summary line -- "-1676887498 checks" was the observed symptom).
inline long long g_failures = 0;
inline long long g_checks = 0;

template <typename T>
std::string describe(const T& value) {
    if constexpr (std::is_convertible_v<T, std::string>) {
        return std::string(value);
    } else if constexpr (std::is_floating_point_v<T>) {
        char buffer[32];
        std::snprintf(buffer, sizeof(buffer), "%.6g", static_cast<double>(value));
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

#define CHECK_EQ(a, b)                                                  \
    ::testing::report((a) == (b), #a " == " #b, __FILE__, __LINE__,     \
                      ::testing::describe(a) + " vs " + ::testing::describe(b))

#define CHECK_NEAR(a, b, tol)                                           \
    ::testing::report(std::fabs((a) - (b)) <= (tol), #a " ~= " #b,      \
                      __FILE__, __LINE__,                               \
                      ::testing::describe(a) + " vs " + ::testing::describe(b))

// Asserts that `body` throws `E`. Failure to throw is as much a bug as the
// wrong result: several of this module's guarantees (coordinate ceilings,
// inverted intervals, malformed peak files) are "refuses to run", and a test
// that only checked the happy path would not see them regress.
#define CHECK_THROWS(E, body)                                           \
    do {                                                                \
        bool threw_ = false;                                            \
        try {                                                           \
            body;                                                       \
        } catch (const E&) {                                            \
            threw_ = true;                                              \
        } catch (...) {                                                 \
        }                                                               \
        ::testing::report(threw_, "throws " #E, __FILE__, __LINE__);    \
    } while (0)

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

// Distinguishes concurrent test runs without a POSIX-only getpid().
inline std::string unique_tag() {
    static std::mt19937_64 rng{std::random_device{}()};
    static int counter = 0;
    return std::to_string(rng()) + "_" + std::to_string(counter++);
}

// A scratch text file removed on destruction.
class TempFile {
public:
    TempFile(const std::string& suffix, const std::string& contents) {
        path_ = (std::filesystem::temp_directory_path() /
                 ("gtk_test_" + unique_tag() + suffix))
                    .string();
        std::FILE* fp = std::fopen(path_.c_str(), "wb");
        if (fp == nullptr) {
            std::fprintf(stderr, "cannot create %s\n", path_.c_str());
            std::exit(2);
        }
        std::fwrite(contents.data(), 1, contents.size(), fp);
        std::fclose(fp);
    }
    ~TempFile() {
        std::error_code ec;
        std::filesystem::remove(path_, ec);
    }
    TempFile(const TempFile&) = delete;
    TempFile& operator=(const TempFile&) = delete;

    [[nodiscard]] const std::string& path() const { return path_; }

private:
    std::string path_;
};

}  // namespace testing
