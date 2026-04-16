/*
 * This file is part of "Polynomials with Small Mahler Measure" (PSMM) project.
 *
 * Copyright (C) 2026, Advanpix LLC.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 */

//
// Cross-check: for every polynomial in the reference set, computing
// its Mahler measure through PSMM's double-precision and extended-
// precision paths must agree with the stored value.
//
// This is the ground-truth integration test for the MPSolve binding.
//

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "psmm.h"
#include "rootfinder.h"
#include "test_helpers.h"

#include <string>

TEST_CASE("Mahler measure (double) matches reference within 1e-12")
{
    const auto refs = load_reference_polynomials(
        std::string(PSMM_TEST_DATA_DIR) + "/reference-polynomials.txt");

    REQUIRE_MESSAGE(!refs.empty(), "reference file missing or empty");

    for (const auto& ref : refs) {
        CAPTURE(ref.N);
        CAPTURE(ref.M_string);

        const double M = compute_mahler_reciprocal_polynomial_d(
            half_to_double(ref.coeffs), /*nthreads=*/1);

        // Double-precision gives ~15 correct digits. The stored value is
        // accurate to 72 digits so the truncated double is exact.
        CHECK(M == doctest::Approx(ref.M).epsilon(1e-12));
    }
}

TEST_CASE("Mahler measure (extended precision) matches reference to 60 digits")
{
    const auto refs = load_reference_polynomials(
        std::string(PSMM_TEST_DATA_DIR) + "/reference-polynomials.txt");

    REQUIRE_MESSAGE(!refs.empty(), "reference file missing or empty");

    const int digits   = 72;
    const int bits     = static_cast<int>(std::ceil(digits * std::log2(10.0)));

    mpf_t expected, got, diff, tol;
    mpf_init2(expected, bits);
    mpf_init2(got,      bits);
    mpf_init2(diff,     bits);
    mpf_init2(tol,      bits);
    // Tolerance ~ 1e-60 — four bits below the 72-digit reference precision.
    mpf_set_ui(tol, 1);
    mpf_div_2exp(tol, tol, static_cast<unsigned>(std::floor(60 * std::log2(10.0))));

    for (const auto& ref : refs) {
        CAPTURE(ref.N);
        CAPTURE(ref.M_string);

        // Expand to full polynomial for mpsolve_compute_mahler_with_properties.
        const auto full = expand_reciprocal(ref.coeffs, ref.N);

        std::size_t K, U, Q, R;
        int err = mpsolve_compute_mahler_with_properties(
            got, ref.N, full.data(), K, U, Q, R,
            /*target_precision=*/bits, /*nthreads=*/1);
        REQUIRE(err == 0);

        str2mpf(expected, ref.M_string.c_str());
        mpf_sub(diff, got, expected);
        mpf_abs(diff, diff);

        const bool within_tol = mpf_cmp(diff, tol) <= 0;
        CHECK_MESSAGE(within_tol,
                      "extended-precision Mahler measure disagrees with reference");
    }

    mpf_clears(expected, got, diff, tol, NULL);
}
