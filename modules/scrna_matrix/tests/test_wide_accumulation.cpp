// Wide accumulation: the numerical contract, and the defect it closes.
//
// scrna_matrix's audit trail records four independent discoveries of ONE root
// cause -- a hand-rolled `float` accumulator saturating to +/-inf --- found by
// four different techniques across three separate review passes
// (docs/AUDIT.md S2, docs/REVIEW_2026-08-15.md finding 5 and its "new
// residual", and T5). Every one of them was found by accident, in the course
// of looking for something else.
//
// This file exists so a fifth is found by a test. It pins the primitives
// (`WideAcc`, `checked_narrow`, the `_wide` reductions) and, more importantly,
// pins the OBSERVABLE consequence: two identical rows must not report a
// similarity of zero, whatever their magnitude.
//
// The last test in this file is not a numerical test at all: it greps the
// module's own sources for the pattern. See `test_no_bare_float_accumulators`
// for why a regex earns its place next to eight behavioural tests.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <regex>
#include <string>
#include <vector>

#include "matrix/block_csr.hpp"
#include "matrix/knn_graph.hpp"
#include "matrix/simd_math.hpp"

namespace {

int g_failures = 0;
int g_checks = 0;

void report(bool ok, const char* expr, int line) {
  ++g_checks;
  if (ok) return;
  ++g_failures;
  std::fprintf(stderr, "  FAIL %s:%d: %s\n", __FILE__, line, expr);
}

#define CHECK(expr) report((expr), #expr, __LINE__)

bool near(double a, double b, double tol) { return std::fabs(a - b) <= tol; }

using scrna::simd::checked_narrow;

// ---------------------------------------------------------------------------

void test_checked_narrow_is_checked_after_the_cast() {
  // The T5 defect in miniature. A quantity can be perfectly finite as a double
  // and still overflow float on the way down, so guarding the wide value alone
  // proves nothing -- 1e300 is a healthy double and +inf as a float.
  CHECK(checked_narrow<float>(1e300) == 0.0f);
  CHECK(checked_narrow<float>(-1e300) == 0.0f);
  CHECK(checked_narrow<float>(std::numeric_limits<double>::infinity()) == 0.0f);
  CHECK(checked_narrow<float>(std::numeric_limits<double>::quiet_NaN()) == 0.0f);

  // Representable values pass through unchanged, including the extremes -- the
  // guard must not be a blanket clamp on "large".
  CHECK(checked_narrow<float>(1.0) == 1.0f);
  CHECK(checked_narrow<float>(-2.5) == -2.5f);
  CHECK(checked_narrow<float>(0.0) == 0.0f);
  CHECK(near(checked_narrow<float>(1e30), 1e30f, 1e24));

  // Identity for a type that is already its own wide accumulator.
  CHECK(checked_narrow<double>(1e300) == 1e300);
}

void test_squaring_does_not_saturate() {
  // AUDIT.md S2, kept here so the primitives file covers the whole family.
  // 1e30 is an ordinary float; 1e30^2 = 1e60 is not.
  const std::vector<float> row{1e30f, 1e30f};
  const double sq = scrna::simd::sparse_norm_sq<float>(row);
  CHECK(std::isfinite(sq));
  CHECK(near(sq, 2e60, 1e56));

  // sqrt before the narrowing is what keeps the norm representable.
  CHECK(std::isfinite(scrna::simd::sparse_norm<float>(row)));
}

void test_the_dot_product_does_not_saturate() {
  // REVIEW_2026-08-15's "new residual": the norms were widened, the dot was
  // not. Both operands here are ~1e60 and neither is a float.
  const std::vector<std::int32_t> idx{0, 1, 2};
  const std::vector<float> val{1e30f, 1e30f, 1e30f};

  const double wide = scrna::simd::sparse_dot_wide<float, std::int32_t>(idx, val, idx, val);
  CHECK(std::isfinite(wide));
  CHECK(near(wide, 3e60, 1e56));

  // The narrowing entry point cannot represent that, and says so by returning
  // 0 rather than +inf -- an infinity would sort ABOVE every real similarity
  // and silently win a top-k.
  const float narrow = scrna::simd::sparse_dot<float, std::int32_t>(idx, val, idx, val);
  CHECK(narrow == 0.0f);
  CHECK(std::isfinite(narrow));
}

void test_identical_large_rows_score_one_not_zero() {
  // THE headline consequence, and the reason none of the above is academic.
  //
  // Before the wide path, a dot of +inf met cosine_from_dot's `!isfinite -> 0`
  // guard and two IDENTICAL rows reported a similarity of ZERO. Not a
  // degraded answer -- an inverted one: the row that should be every other
  // row's nearest neighbour matched nothing at all, from finite input, with no
  // diagnostic anywhere.
  const std::vector<std::int32_t> idx{0, 1, 2};
  const std::vector<float> val{1e30f, 2e30f, 3e30f};

  const float sim = scrna::simd::sparse_cosine_similarity<float, std::int32_t>(
      idx, val, idx, val);
  CHECK(near(sim, 1.0f, 1e-5f));

  // Orthogonal rows at the same magnitude still score 0 -- the fix must not
  // turn the guard into "always return something".
  const std::vector<std::int32_t> other{3, 4};
  const std::vector<float> other_val{1e30f, 1e30f};
  CHECK((scrna::simd::sparse_cosine_similarity<float, std::int32_t>(
             idx, val, other, other_val) == 0.0f));

  // And an antiparallel pair scores -1, not 0: the sign has to survive too.
  const std::vector<float> negated{-1e30f, -2e30f, -3e30f};
  CHECK(near(scrna::simd::sparse_cosine_similarity<float, std::int32_t>(
                 idx, val, idx, negated), -1.0f, 1e-5f));
}

void test_cosine_from_dot_wide_beats_the_narrow_overload() {
  // Same inputs, two entry points, deliberately different answers. This is
  // what makes "stay wide until the division" a rule rather than a preference:
  // the narrow overload is not less precise here, it is WRONG, and the
  // difference is the whole finding.
  const double dot = 3e60;
  const float norm = static_cast<float>(std::sqrt(3e60));  // ~1.73e30, representable

  CHECK(near(scrna::simd::cosine_from_dot_wide<float>(dot, norm, norm), 1.0f, 1e-5f));
  // Narrowing the dot first turns 3e60 into +inf, which the guard reports as 0.
  CHECK(scrna::simd::cosine_from_dot(static_cast<float>(dot), norm, norm) == 0.0f);
}

void test_dispatcher_falls_back_when_a_vector_kernel_saturates() {
  // The dispatcher's contract: whichever ISA it picked, a saturating result is
  // recomputed exactly rather than returned. Every compiled kernel accumulates
  // in float, so this exercises the fallback on whatever the host runs.
  const std::size_t n = 64;
  std::vector<float> dense(n, 1e30f);
  std::vector<std::int32_t> idx(n);
  std::vector<float> val(n, 1e30f);
  for (std::size_t i = 0; i < n; ++i) idx[i] = static_cast<std::int32_t>(i);

  const double wide = scrna::simd::gather_dot_wide(dense.data(), idx.data(), val.data(), n);
  CHECK(std::isfinite(wide));
  CHECK(near(wide, 64e60, 1e58));
  // Bit-identical to the reference the fallback is supposed to call.
  CHECK(wide == scrna::simd::gather_dot_scalar_wide(dense.data(), idx.data(),
                                                     val.data(), n));

  const double dense_wide = scrna::simd::dense_dot_wide(dense.data(), val.data(), n);
  CHECK(std::isfinite(dense_wide));
  CHECK(dense_wide == scrna::simd::dense_dot_scalar_wide(dense.data(), val.data(), n));

  // Ordinary magnitudes must NOT take the fallback path silently changing the
  // answer: the fast result is returned as-is when it is finite.
  std::vector<float> small_dense(n, 1.5f);
  std::vector<float> small_val(n, 2.0f);
  CHECK(near(scrna::simd::gather_dot_wide(small_dense.data(), idx.data(),
                                          small_val.data(), n),
             192.0, 1e-9));
}

void test_bruteforce_knn_survives_a_matrix_at_the_validation_ceiling() {
  // End to end through the public API, at the largest magnitude
  // Block_CSR::validate() admits (||row||^2 <= 1e34). Two identical rows plus
  // one orthogonal: the identical pair must be each other's nearest neighbour.
  const float v = 5e16f;  // ||row||^2 = 5e33 for two entries, under the ceiling
  scrna::AlignedVector<float> values = {v, v, v, v, 1.0f};
  scrna::AlignedVector<std::int32_t> col_idx = {0, 2, 0, 2, 1};
  scrna::AlignedVector<std::int32_t> row_ptr = {0, 2, 4, 5};
  scrna::Block_CSR<float> mat(3, 4, std::move(values), std::move(col_idx),
                              std::move(row_ptr));

  const auto graph = scrna::build_knn_graph_bruteforce(mat, 2);
  for (float s : graph.neighbor_sim) CHECK(std::isfinite(s));
  CHECK(graph.neighbor_idx[0] == 1);
  CHECK(near(graph.neighbor_sim[0], 1.0f, 1e-4f));
  CHECK(graph.neighbor_idx[1] == 2);
  CHECK(near(graph.neighbor_sim[1], 0.0f, 1e-4f));
}

void test_wide_reductions_agree_with_their_narrow_wrappers() {
  // For values that comfortably fit, the two entry points must not diverge --
  // otherwise the refactor moved every ordinary answer as a side effect.
  const std::size_t n = 37;  // deliberately not a multiple of any lane width
  std::vector<float> dense(256), val(n);
  std::vector<std::int32_t> idx(n);
  // The int casts are load-bearing: `i % 13 - 6` on a std::size_t underflows
  // to ~1.8e19 for the first six residues, which made this test feed the
  // kernel 4.6e18 and quietly exercise the saturation path it is meant to
  // stay clear of. Caught by the assertion below, which is the argument for
  // asserting the ordinary case as well as the extreme one.
  for (std::size_t i = 0; i < dense.size(); ++i) {
    dense[i] = static_cast<float>(static_cast<int>(i % 13) - 6) * 0.25f;
  }
  for (std::size_t i = 0; i < n; ++i) {
    idx[i] = static_cast<std::int32_t>(i * 3);
    val[i] = static_cast<float>(static_cast<int>(i % 7) - 3) * 0.5f;
  }

  const double wide = scrna::simd::gather_dot_wide(dense.data(), idx.data(), val.data(), n);
  const float narrow = scrna::simd::gather_dot(dense.data(), idx.data(), val.data(), n);
  CHECK(near(narrow, wide, 1e-4));
  CHECK(std::isfinite(narrow));
  CHECK(std::fabs(wide) < 1e4);  // the inputs are ordinary; nothing saturated
}

// ---------------------------------------------------------------------------
// The source-pattern guard
// ---------------------------------------------------------------------------
//
// Every test above pins the behaviour of a reduction that EXISTS. This one
// pins the shape of the ones written next, and it is here because of a fact
// about how the previous four instances were found: differential testing
// against a double reference (AUDIT.md S2), a boundary review (REVIEW finding
// 5), a defence-in-depth read (T5), and a regression test written for an
// unrelated finding (finding 5's "new residual"). Four techniques, three
// review passes, and not one of them was looking for this.
//
// All four are expensive and none is repeatable on demand. The pattern is
// three tokens long. So: the numerical tests cover what is here, and the regex
// covers what comes next -- which is the half that would otherwise wait for a
// fifth audit to notice.
//
// THE RULE: accumulate in WideAcc<T>, narrow once through checked_narrow().
//
// THE ESCAPE HATCH: a deliberate float accumulator marks itself.
//
//     float sum = vaddvq_f32(...);  // SCRNA_FLOAT_ACC_OK: vector kernel tail
//
// The marker is accepted on the offending line or the one above it. It is not
// a rubber stamp -- it is a grep target. `git grep SCRNA_FLOAT_ACC_OK`
// enumerates every place this module deliberately accumulates narrow, which is
// exactly the list a future audit wants and exactly the list nobody could
// produce before.
//
// NOT flagged: `float32x4_t` / `__m256` / `__m512` accumulators. Those are
// unavoidable at full lane width and are covered a level up -- the dispatchers
// recompute exactly when a vector kernel returns a non-finite result. Marking
// every SIMD kernel would teach the reader that the marker means nothing.

struct Violation {
  std::string file;
  int line;
  std::string text;
  std::string why;
};

// Strips a trailing `//` comment so prose describing the defect is not read as
// an instance of it -- this file and simd_math.hpp both discuss `float acc`
// at length. Block comments are not handled; the module does not use them for
// code-shaped text, and a false positive is a marker away from resolved.
std::string strip_line_comment(const std::string& line) {
  const auto at = line.find("//");
  return at == std::string::npos ? line : line.substr(0, at);
}

void scan_file(const std::filesystem::path& path, std::vector<Violation>& out) {
  static const std::regex kDeclZero(
      R"((^|[^A-Za-z0-9_])float[ 	]+([A-Za-z_][A-Za-z0-9_]*)[ 	]*(=[ 	]*0|\{[ 	]*0))");
  static const std::regex kDeclAny(
      R"((^|[^A-Za-z0-9_])float[ 	]+([A-Za-z_][A-Za-z0-9_]*)[ 	]*[=;])");
  static const std::regex kPlusEquals(R"(([A-Za-z_][A-Za-z0-9_]*)[ 	]*\+=)");

  std::ifstream in(path);
  if (!in) return;

  std::vector<std::string> float_vars;
  std::string raw, previous;
  int line_no = 0;

  while (std::getline(in, raw)) {
    ++line_no;
    const bool exempt = raw.find("SCRNA_FLOAT_ACC_OK") != std::string::npos ||
                        previous.find("SCRNA_FLOAT_ACC_OK") != std::string::npos;
    const std::string code = strip_line_comment(raw);
    previous = raw;

    std::smatch m;
    if (!exempt && std::regex_search(code, m, kDeclZero)) {
      out.push_back({path.string(), line_no, raw,
                     "`float` initialised to zero -- accumulate in WideAcc<T>"});
    }
    // Remember the name even when the declaration was not zero-initialised:
    // the vector kernels' scalar tails accumulate into a `float sum` seeded
    // from a horizontal reduction, which a declaration-only check misses.
    if (std::regex_search(code, m, kDeclAny)) float_vars.push_back(m[2].str());

    if (!exempt && std::regex_search(code, m, kPlusEquals)) {
      const std::string target = m[1].str();
      for (const std::string& v : float_vars) {
        if (v == target) {
          out.push_back({path.string(), line_no, raw,
                         "`+=` onto the `float` variable '" + target + "'"});
          break;
        }
      }
    }
  }
}

void test_no_bare_float_accumulators() {
  const std::filesystem::path root(SCRNA_SOURCE_DIR);
  std::vector<std::filesystem::path> files;
  for (const char* dir : {"include/matrix", "src"}) {
    std::error_code ec;
    for (const auto& entry : std::filesystem::directory_iterator(root / dir, ec)) {
      const auto ext = entry.path().extension();
      if (ext == ".hpp" || ext == ".cpp") files.push_back(entry.path());
    }
  }

  if (files.empty()) {
    // An installed or packaged build has no source tree. Say so rather than
    // reporting a pass over nothing -- a guard that silently scans zero files
    // is worse than no guard, because it reads as coverage.
    std::fprintf(stderr, "  WARNING: no sources found under %s; guard did not run\n",
                 SCRNA_SOURCE_DIR);
    ++g_failures;
    return;
  }

  std::vector<Violation> violations;
  for (const auto& f : files) scan_file(f, violations);

  for (const Violation& v : violations) {
    std::fprintf(stderr, "  FAIL %s:%d: %s\n        %s\n", v.file.c_str(), v.line,
                 v.why.c_str(), v.text.c_str());
  }
  if (!violations.empty()) {
    std::fprintf(stderr,
                 "\n  A float accumulator does not lose precision on large input --\n"
                 "  it saturates to +/-inf, and the downstream `!isfinite -> 0` guards\n"
                 "  then report two IDENTICAL rows as having zero similarity.\n"
                 "  Accumulate in scrna::simd::WideAcc<T>; narrow once through\n"
                 "  checked_narrow(). If narrow accumulation is genuinely correct\n"
                 "  here, say so on the line: // SCRNA_FLOAT_ACC_OK: <why>\n\n");
  }
  g_failures += static_cast<int>(violations.size());
  g_checks += static_cast<int>(files.size());
  std::printf("  scanned %zu source file(s), %zu violation(s)\n", files.size(),
              violations.size());
}

}  // namespace

int main() {
  std::printf("== wide_accumulation ==\n");
  test_checked_narrow_is_checked_after_the_cast();
  test_squaring_does_not_saturate();
  test_the_dot_product_does_not_saturate();
  test_identical_large_rows_score_one_not_zero();
  test_cosine_from_dot_wide_beats_the_narrow_overload();
  test_dispatcher_falls_back_when_a_vector_kernel_saturates();
  test_bruteforce_knn_survives_a_matrix_at_the_validation_ceiling();
  test_wide_reductions_agree_with_their_narrow_wrappers();
  test_no_bare_float_accumulators();

  std::printf("%d checks, %d failure(s)\n", g_checks, g_failures);
  if (g_failures == 0) std::printf("All tests passed.\n");
  return g_failures == 0 ? 0 : 1;
}
