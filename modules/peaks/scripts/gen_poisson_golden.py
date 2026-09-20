#!/usr/bin/env python3
"""Regenerate the golden tables in tests/test_poisson_model.cpp.

    python3 modules/peaks/scripts/gen_poisson_golden.py

Prints C++ initialisers to stdout; paste them over the corresponding tables.
Requires mpmath, which is NOT a dependency of anything that builds or runs:
the values are checked in, the suite needs a compiler and nothing else, and
this script exists so the next person can add a case rather than trust one.

Why mpmath and not scipy: 80 significant digits, and an implementation that
shares no code, no algorithm and no floating-point format with the header under
test. A golden table computed the same way as the thing it checks is a
tautology.

Two details that are easy to get wrong and were:

  * ln P(X >= k) is the LOWER regularised incomplete gamma, P(k, lambda) --
    the Poisson upper tail is the gamma lower tail. gammainc(k, 0, lam).
  * Where the upper tail is close to 1, log(P) is computed as log1p(-Q) from
    the UPPER regularised gamma instead. At k = 1000, lambda = 2000 the true
    answer is -6.8e-136, and 1 - 6.8e-136 is not representable even at 80
    digits: log(gammainc(...)) would report exactly 0, and the golden value
    would be less accurate than the double it is there to check.

Cases past k ~ 1e6 use the tail series directly, at the same precision;
gammainc() is impractical at k ~ 2^32 and the series is exact mathematics,
converging in a handful of terms whenever lambda << k.
"""

from mpmath import mp, mpf, log, log1p, gammainc, loggamma, inf

mp.dps = 80

SF_CASES = [
    (1, "1e-300"), (1, "1e-10"), (1, "0.5"), (1, "1"), (1, "10"), (1, "1000"),
    (2, "0.5"), (2, "1"), (2, "2"), (2, "3"), (2, "100"),
    (3, "1"), (3, "3"), (3, "4"), (5, "1"), (5, "5"), (5, "6"),
    (10, "1"), (10, "5"), (10, "9"), (10, "10"), (10, "11"), (10, "20"),
    (50, "10"), (50, "49"), (50, "50"), (50, "51"), (50, "80"),
    (100, "1"), (100, "50"), (100, "99"), (100, "100"), (100, "101"), (100, "200"),
    # 511-514 straddle the end of the log-factorial table.
    (511, "100"), (512, "100"), (513, "100"), (514, "100"),
    (513, "512"), (513, "513"), (513, "514"),
    (1000, "1"), (1000, "100"), (1000, "999"), (1000, "1000"), (1000, "1001"),
    (1000, "2000"),
    (10000, "9999"), (10000, "10000"), (10000, "10001"),
    (100000, "1000"), (100000, "99999"), (100000, "100001"),
    (1000000, "999000"), (1000000, "1001000"),
]

# Counts where gammainc() is impractical; summed from the tail series instead.
SF_SERIES_CASES = [
    (65536, "1024"), (100000, "50000"), (1000000, "1"),
    (2147483648, "1000"), (4294967295, "1e6"), (4294967295, "1e9"),
]

LGAMMA_CASES = ["0.5", "1", "1.5", "2", "2.5", "3", "3.7", "7.25", "15.999",
                "16", "16.001", "17", "100.5", "512.5", "1000.5", "1e6",
                "1000000.5", "1e15", "1e100", "1e300"]

FACTORIAL_CASES = [0, 1, 2, 3, 10, 20, 21, 170, 511, 512, 513, 1000, 100000,
                   4294967295]


def ln_sf(k, lam):
    """ln P(X >= k) for X ~ Poisson(lam), via the incomplete gamma."""
    k, lam = mpf(k), mpf(lam)
    q = gammainc(k, lam, inf, regularized=True)  # P(X <= k-1)
    if q < mpf("0.5"):
        return log1p(-q)
    return log(gammainc(k, 0, lam, regularized=True))


def ln_sf_series(k, lam):
    """Same quantity, summed from P(X>=k) = pmf(k) * sum lam^n k!/(k+n)!."""
    k, lam = mpf(k), mpf(lam)
    lead = -lam + k * log(lam) - loggamma(k + 1)
    term, total, n = mpf(1), mpf(1), 0
    while True:
        n += 1
        term *= lam / (k + n)
        total += term
        if term < total * mpf("1e-70") or n > 100000:
            break
    return lead + log(total)


def emit_sf(k, lam, value):
    """One C++ initialiser, with values a double cannot hold written as zero.

    ln P(X >= 1 | lambda = 1000) is -5.08e-435. That is not a small double, it
    is no double at all -- the smallest denormal is 4.9e-324 -- so the literal
    would underflow at compile time, with a warning, to the -0.0 written here.
    Saying so beats emitting a number the compiler silently discards.
    """
    if value != 0 and abs(value) < mpf("4.9406564584124654e-324"):
        return ("    // true value %s, below the smallest denormal\n"
                "    {%du, %s, -0.0}," % (mp.nstr(value, 20), k, repr(float(mpf(lam)))))
    return "    {%du, %s, %s}," % (k, repr(float(mpf(lam))), mp.nstr(value, 20))


def main():
    print("// k, lambda, ln P(X >= k)")
    for k, lam in SF_CASES:
        print(emit_sf(k, lam, ln_sf(k, lam)))
    for k, lam in SF_SERIES_CASES:
        print(emit_sf(k, lam, ln_sf_series(k, lam)))

    print("\n// z, lnGamma(z)")
    for z in LGAMMA_CASES:
        print("    {%s, %s}," % (repr(float(mpf(z))), mp.nstr(loggamma(mpf(z)), 20)))

    print("\n// k, ln(k!)")
    for k in FACTORIAL_CASES:
        print("    {%du, %s}," % (k, mp.nstr(loggamma(mpf(k) + 1), 20)))


if __name__ == "__main__":
    main()
