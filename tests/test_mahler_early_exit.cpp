/*
 * This file is part of "Polynomials with Small Mahler Measure" (PSMM) project.
 *
 * Copyright (C) 2026, Advanpix LLC.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 */

//
// A1 regression test.
//
// estimate_mahler_reciprocal_polynomial_d must NOT report
// "over threshold" for polynomials whose true Mahler measure is
// strictly below the threshold. Before the A1 fix the early-exit
// predicate inside MPSolve's mps_mahler_check_stop trips on
//   R > T || abs(R - r) > T
// where R is an Aberth disk centre and r its radius, plus abs(...)
// resolves to the integer <stdlib.h> abs. This is optimistic and can
// falsely mark near-threshold polynomials as over-threshold,
// silently dropping them from the search. The correct sufficient
// condition is R - r > T.
//
// This test loads the near-threshold (M ∈ [1.29, 1.30]) entries from
// the reference set and asserts none of them are reported as over
// threshold when the estimator is called with T = 1.3.
//
// When the test fails on main/pre-A1-fix, that is the bug. When it
// passes post-fix, the bug is fixed.
//

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "psmm.h"
#include "rootfinder.h"
#include "mahlerestimator.h"
#include "test_helpers.h"

#include <string>

TEST_CASE("Near-threshold polynomials are NOT falsely rejected by fast estimator")
{
    const auto refs = load_reference_polynomials(
        std::string(PSMM_TEST_DATA_DIR) + "/reference-polynomials.txt");

    REQUIRE_MESSAGE(!refs.empty(), "reference file missing or empty");

    const double T = 1.3;

    // The estimator returns 2*T when it decides the polynomial is over
    // threshold, otherwise the computed Mahler measure.
    int false_rejections = 0;
    int checked          = 0;

    for (const auto& ref : refs) {
        if (ref.M >= T) continue; // only meaningful for polys below T
        ++checked;

        const double m = estimate_mahler_reciprocal_polynomial_d(
            ref.coeffs, T, /*nthreads=*/1);

        CAPTURE(ref.N);
        CAPTURE(ref.M_string);
        CAPTURE(m);
        CAPTURE(T);

        const bool falsely_over_threshold = (m > T);
        if (falsely_over_threshold) ++false_rejections;

        CHECK_MESSAGE(!falsely_over_threshold,
                      "fast estimator falsely rejected a sub-threshold polynomial");
    }

    REQUIRE(checked > 0);
    MESSAGE("Near-threshold polynomials checked: " << checked
            << ", false rejections: " << false_rejections);
}

TEST_CASE("Over-threshold polynomials are still reported as over threshold")
{
    // Pick a polynomial whose known Mahler measure is > 1.3, then ensure
    // the estimator does NOT regress in the other direction — i.e. the
    // fix doesn't disable the early exit entirely.
    //
    // Use x^10 + 1 (cyclotomic, all roots on unit circle => M = 1) as
    // a sanity polynomial and x^4 - x - 1 which has M ~ 1.2207 below
    // threshold; then craft a polynomial whose M is clearly above 1.3.
    //
    // Reciprocal polynomial a[0] .. a[N/2]:
    //   p(x) = x^10 + 2*x^5 + 1 has M = ? > 1 since not cyclotomic.
    //   Simpler: p(x) = x^2 - 3*x + 1. That's not reciprocal though.
    //   For the estimator API we need reciprocal form.
    //
    // Use p(x) = x^4 - 3*x^3 + x^2 - 3*x + 1 (reciprocal, large coeff).
    // a[0]=1, a[1]=-3, a[2]=1 in half-form. N=4.

    std::vector<int> big = {1, -3, 1};
    const double T = 1.3;

    const double m = estimate_mahler_reciprocal_polynomial_d(big, T, 1);
    CAPTURE(m);
    CAPTURE(T);
    CHECK(m > T); // must report over threshold (estimator returns 2*T or exact)
}
