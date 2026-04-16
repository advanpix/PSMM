/*
 * This file is part of "Polynomials with Small Mahler Measure" (PSMM) project.
 *
 * Copyright (C) 2026, Advanpix LLC.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 */

//
// Sanity tests for reciprocal_polynomials_iterator:
//  - total enumerated polynomials matches compute_total_number_of_polynomials;
//  - no polynomial is emitted twice (uniqueness of enumerated patterns);
//  - every emitted polynomial has the expected nnz among a[1..N/2].
//

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "psmm.h"
#include "utilities.h"
#include "polyiterator.h"

#include <set>
#include <string>
#include <sstream>
#include <vector>

namespace {

std::string encode(const std::vector<int>& poly)
{
    std::ostringstream oss;
    for (int c : poly) oss << c << ',';
    return oss.str();
}

std::size_t count_enumerated(int degree, int nnz, const std::vector<int>& coeffs)
{
    std::vector<int> buf(degree / 2 + 1);
    reciprocal_polynomials_iterator it(degree, nnz, coeffs);
    std::size_t n = 0;
    while (it.next_polynomial(buf)) ++n;
    return n;
}

std::size_t count_unique(int degree, int nnz, const std::vector<int>& coeffs)
{
    std::set<std::string> seen;
    std::vector<int> buf(degree / 2 + 1);
    reciprocal_polynomials_iterator it(degree, nnz, coeffs);
    while (it.next_polynomial(buf)) seen.insert(encode(buf));
    return seen.size();
}

std::size_t expected_total(int degree, int nnz_val, int base)
{
    mpz_t m;
    mpz_init(m);
    std::vector<int> nnz{nnz_val};
    compute_total_number_of_polynomials(m, degree, base, nnz);
    const std::size_t r = mpz_get_ui(m);
    mpz_clear(m);
    return r;
}

} // namespace

TEST_CASE("Iterator enumerates the expected number of polynomials")
{
    const std::vector<int> coeffs_pm1  = {-1, 1};          // base = 2
    const std::vector<int> coeffs_pm10 = {-1, 0, 1};       // base = 3

    struct Case { int degree; int nnz; const std::vector<int>* alphabet; };
    const std::vector<Case> cases = {
        {10, 1, &coeffs_pm1},
        {10, 2, &coeffs_pm1},
        {12, 3, &coeffs_pm1},
        {14, 2, &coeffs_pm10},
        {16, 3, &coeffs_pm10},
    };

    for (const auto& c : cases) {
        CAPTURE(c.degree);
        CAPTURE(c.nnz);

        const std::size_t enumerated = count_enumerated(c.degree, c.nnz, *c.alphabet);
        const std::size_t expected   = expected_total(c.degree, c.nnz,
                                                     static_cast<int>(c.alphabet->size()));
        CHECK(enumerated == expected);
    }
}

TEST_CASE("Iterator emits no duplicates when 0 is not in alphabet")
{
    // When the alphabet does not contain 0, every pattern/number pair maps
    // to a unique polynomial. (If 0 is in the alphabet, the iterator
    // intentionally over-enumerates; the main search loop uses the
    // alphabet {-1, 1} and relies on this invariant.)
    const std::vector<int> alphabet = {-1, 1};
    const std::size_t enumerated = count_enumerated(12, 3, alphabet);
    const std::size_t unique     = count_unique(12, 3, alphabet);
    CHECK(enumerated == unique);
}

TEST_CASE("Every emitted polynomial has a[0] = 1 and at most nnz nonzeros")
{
    // The iterator's `nnz` parameter is the number of pattern positions
    // that are "used" — if 0 is among the possible coefficients the
    // output can still have fewer than `nnz` nonzero entries, so we
    // check an upper bound. When 0 is absent from the alphabet, the
    // bound is tight.
    const std::vector<int> alphabet_no_zero = {-1, 1};
    const std::vector<int> alphabet_zero    = {-1, 0, 1};
    const int degree = 10, nnz = 2;

    auto run = [&](const std::vector<int>& alphabet, bool tight){
        std::vector<int> buf(degree / 2 + 1);
        reciprocal_polynomials_iterator it(degree, nnz, alphabet);
        std::size_t n = 0;
        while (it.next_polynomial(buf)) {
            CHECK(buf.front() == 1);
            int count = 0;
            for (std::size_t k = 1; k < buf.size(); ++k) count += (buf[k] != 0);
            CHECK(count <= nnz);
            if (tight) CHECK(count == nnz);
            ++n;
        }
        REQUIRE(n > 0);
    };

    run(alphabet_no_zero, /*tight=*/true);
    run(alphabet_zero,    /*tight=*/false);
}
