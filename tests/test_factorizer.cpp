/*
 * This file is part of "Polynomials with Small Mahler Measure" (PSMM) project.
 *
 * Copyright (C) 2026, Advanpix LLC.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 */

//
// A2 regression test. factor_polynomial must correctly report
// irreducibility:
//
//   - an irreducible polynomial returns TRUE;
//   - a polynomial with a single repeated factor (multiplicity > 1)
//     returns FALSE;
//   - a polynomial with several distinct factors returns FALSE;
//   - when reducible, the `factors` vector contains each distinct
//     factor repeated according to its multiplicity.
//

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "psmm.h"
#include "factorizer.h"

#include <vector>

TEST_CASE("Irreducible polynomial reports irreducible")
{
    // x^2 + x + 1 (the cyclotomic Phi_3) — irreducible over Z.
    std::vector<int> p = {1, 1, 1};
    std::vector<std::vector<int>> factors;

    const bool irr = factor_polynomial(p, factors);
    CHECK(irr == true);
    CHECK(factors.size() == 1);
}

TEST_CASE("Polynomial with a repeated factor is NOT irreducible")
{
    // (x^2 + 1)^2 = x^4 + 2x^2 + 1. One distinct factor, multiplicity 2.
    // NTL returns length() == 1, so the old `length() == 1` test said
    // "irreducible" — that's the bug this test pins down.
    std::vector<int> p = {1, 0, 2, 0, 1};
    std::vector<std::vector<int>> factors;

    const bool irr = factor_polynomial(p, factors);
    CHECK(irr == false);

    // Factors vector must carry the repeated factor twice, so the caller
    // iterating over it can process each copy.
    CHECK(factors.size() == 2);
    if (factors.size() == 2) {
        CHECK(factors[0] == factors[1]);
        // The factor itself should be x^2 + 1 = [1, 0, 1].
        CHECK(factors[0] == std::vector<int>{1, 0, 1});
    }
}

TEST_CASE("Polynomial with several distinct factors is reducible")
{
    // x^4 - 1 = (x - 1)(x + 1)(x^2 + 1), three distinct factors, each mult 1.
    std::vector<int> p = {-1, 0, 0, 0, 1};
    std::vector<std::vector<int>> factors;

    const bool irr = factor_polynomial(p, factors);
    CHECK(irr == false);
    CHECK(factors.size() == 3);
}
