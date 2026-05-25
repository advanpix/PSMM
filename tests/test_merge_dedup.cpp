/*
 * This file is part of "Polynomials with Small Mahler Measure" (PSMM) project.
 *
 * Copyright (C) 2026, Advanpix LLC.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 */

//
// Regression test for merge_files_with_results (polyhelpers.h).
//
// Reproduces the 2026-05-25 incident (commit bc1e760) where merge silently
// dropped distinct polynomials that happened to share their Mahler measure
// with a lower-degree entry. The bug was in same_polynomial_found_m being
// used as a dedup predicate, which collapsed (N, M)-equivalence classes
// across different N.
//
// The fix dedupes by (N, half_coefficients) and x->-x equivalence only.
// Two polynomials with different (N, coeffs) but coincident M must both
// survive a merge.
//

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "psmm.h"
#include "utilities.h"
#include "polyhelpers.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>


namespace {

// Write a small DB-format file. Each entry is one line:
//   N M NNZ H L K U Q R coeffs...
// The metadata fields are placeholders -- only N, M, and the coefficient
// vector matter for the merge-dedup behavior under test.
void write_entries(const std::string& path,
                   const std::vector<std::pair<std::size_t, std::string>>& entries,
                   const std::vector<std::vector<int>>& coeffs)
{
    REQUIRE(entries.size() == coeffs.size());
    std::ofstream ofs(path);
    REQUIRE(ofs.is_open());
    ofs << "\n# test fixture\n";
    for(std::size_t i = 0; i < entries.size(); ++i)
    {
        ofs << entries[i].first << ' ' << entries[i].second
            << " 1 1 5 1 0 0 0";   // placeholder NNZ, H, L, K, U, Q, R
        for(int c : coeffs[i]) ofs << ' ' << c;
        ofs << " \n";
    }
}

std::vector<std::tuple<std::size_t, std::string, std::vector<int>>>
read_entries(const std::string& path)
{
    std::vector<std::tuple<std::size_t, std::string, std::vector<int>>> out;
    std::ifstream ifs(path);
    REQUIRE(ifs.is_open());
    std::string line;
    while(std::getline(ifs, line))
    {
        // skip blank/header
        std::size_t first = line.find_first_not_of(" \t");
        if(first == std::string::npos) continue;
        if(line[first] == '#') continue;
        std::istringstream iss(line);
        std::size_t N;
        std::string M;
        int nnz, H, L, K, U, Q, R;
        if(!(iss >> N >> M >> nnz >> H >> L >> K >> U >> Q >> R)) continue;
        std::vector<int> coeffs;
        int c;
        while(iss >> c) coeffs.push_back(c);
        out.emplace_back(N, M, coeffs);
    }
    return out;
}

std::string tmp_path(const std::string& tag)
{
    auto dir = std::filesystem::temp_directory_path();
    auto p = dir / ("psmm_merge_test_" + tag + "_"
                    + std::to_string(::getpid()) + ".txt");
    return p.string();
}

} // namespace


TEST_CASE("merge preserves distinct polynomials with coincident Mahler measure")
{
    // Construct three entries:
    //   A: N=10, M=1.176280818259917506544070338474035050693415806564695259830  (Lehmer)
    //   B: N=20, M=1.176280818259917506544070338474035050693415806564695259830  (FAKE entry sharing M)
    //   C: N=10, M=1.230391434407224702790177938975279017566574489661756241401  (distinct M)
    //
    // The buggy merge would drop B because it scans verified for any p with
    // p.N <= 20 and same M, finds A, and skips B. The fixed merge dedupes by
    // (N, half_coefficients) only and keeps all three.

    const std::string M_LEHMER =
        "1.176280818259917506544070338474035050693415806564695259830106347029688377";
    const std::string M_OTHER =
        "1.230391434407224702790177938975279017566574489661756241401914236172813448";

    std::vector<std::pair<std::size_t, std::string>> meta = {
        {10, M_LEHMER},
        {20, M_LEHMER},
        {10, M_OTHER},
    };
    std::vector<std::vector<int>> coeffs = {
        {1, 1, 0, -1, -1, -1},                     // Lehmer half (N=10)
        {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1},        // N=20 placeholder (distinct half)
        {1, 0, 0, 1, 0, 1},                        // N=10 distinct half
    };

    const std::string in_path  = tmp_path("in");
    const std::string out_path = tmp_path("out");
    write_entries(in_path, meta, coeffs);

    merge_files_with_results(in_path, out_path,
                             /*verify_precision=*/128,
                             /*output_digits=*/72,
                             /*precision=*/256);

    auto rows = read_entries(out_path);
    CHECK(rows.size() == 3);

    bool saw_A = false, saw_B = false, saw_C = false;
    for(const auto& [N, M, c] : rows)
    {
        if(N == 10 && c == coeffs[0]) saw_A = true;
        if(N == 20 && c == coeffs[1]) saw_B = true;
        if(N == 10 && c == coeffs[2]) saw_C = true;
    }
    CHECK(saw_A);
    CHECK(saw_B);  // The regression target: must survive even though M == M(A).
    CHECK(saw_C);

    std::filesystem::remove(in_path);
    std::filesystem::remove(out_path);
}


TEST_CASE("merge collapses x->-x flips of the same polynomial")
{
    // Lehmer's classical canonical form versus its x->-x flip:
    //   classical : a[0..5] = [1, 1, 0, -1, -1, -1]
    //   x->-x     : negate odd-index entries -> [1, -1, 0, 1, -1, 1]
    // Both should be reduced to a single entry under (N, half) + x->-x dedup.

    const std::string M_LEHMER =
        "1.176280818259917506544070338474035050693415806564695259830106347029688377";

    std::vector<std::pair<std::size_t, std::string>> meta = {
        {10, M_LEHMER},
        {10, M_LEHMER},
    };
    std::vector<std::vector<int>> coeffs = {
        {1,  1, 0, -1, -1, -1},   // classical Lehmer
        {1, -1, 0,  1, -1,  1},   // x->-x flip
    };

    const std::string in_path  = tmp_path("xneg_in");
    const std::string out_path = tmp_path("xneg_out");
    write_entries(in_path, meta, coeffs);

    merge_files_with_results(in_path, out_path, 128, 72, 256);

    auto rows = read_entries(out_path);
    CHECK(rows.size() == 1);

    std::filesystem::remove(in_path);
    std::filesystem::remove(out_path);
}


TEST_CASE("xneg_flip_half negates only odd-index entries")
{
    std::vector<int> in = {1, 2, 3, 4, 5, 6};
    auto out = xneg_flip_half(in);
    REQUIRE(out.size() == in.size());
    CHECK(out[0] ==  1);
    CHECK(out[1] == -2);
    CHECK(out[2] ==  3);
    CHECK(out[3] == -4);
    CHECK(out[4] ==  5);
    CHECK(out[5] == -6);
}


namespace {
// Construct a polynomial entry with a specific (N, M, M_mpf) for the
// search-phase dedup tests. nnz / coeffs irrelevant for these tests.
reciprocal_polynomial_t make_entry(std::size_t N, double M, int bits = 128)
{
    reciprocal_polynomial_t p;
    p.N      = N;
    p.M      = M;
    p.F.set_prec(bits);
    p.F      = M;
    p.nnz    = 1;
    p.H      = 1;
    p.L      = 5;
    p.K      = 1;
    p.U      = 0;
    p.Q      = 0;
    p.R      = 0;
    p.coeffs.assign(N / 2 + 1, 0);
    p.coeffs[0] = 1;
    return p;
}
} // namespace


TEST_CASE("same_polynomial_found: heuristic skip filters ONLY lower-or-equal-degree DB entries")
{
    // This is the invariant Pavel flagged 2026-05-25: with high-degree
    // polynomials now in the DB, the search-phase reducibility-skip MUST
    // NOT falsely skip a low-degree candidate just because a higher-degree
    // DB entry happens to share its Mahler measure (e.g. when the higher-N
    // entry is the reducible extension of the candidate's irreducible
    // kernel). The condition `p.N <= n` in same_polynomial_found{,_m}
    // enforces this -- DB entries with p.N > n are skipped.

    std::vector<reciprocal_polynomial_t> db;
    db.push_back(make_entry(10, 1.176280818259917));   // Lehmer (low N)
    db.push_back(make_entry(916, 1.286084660801062)); // high-N entry sharing M with hypothetical N=458 candidate

    const double tol = 1e-9;

    SUBCASE("candidate at lower N than the only matching DB entry -> NOT skipped")
    {
        // Candidate at N=458 with M close to the N=916 DB entry.
        // 916 > 458, so the N=916 entry must be EXCLUDED from the check.
        int hit = same_polynomial_found(/*candidate_N=*/458, /*M=*/1.286084660801062, tol, db);
        CHECK(hit == 0);
    }

    SUBCASE("candidate at higher N than a matching DB entry -> SKIPPED (returns DB entry's N)")
    {
        // Candidate at N=20 with M matching Lehmer (N=10). 10 < 20, so the
        // Lehmer entry IS in the search set -> match -> return 10.
        int hit = same_polynomial_found(/*candidate_N=*/20, /*M=*/1.176280818259917, tol, db);
        CHECK(hit == 10);
    }

    SUBCASE("candidate at same N as matching DB entry -> SKIPPED")
    {
        // Same degree: 10 <= 10 -> in the search set.
        int hit = same_polynomial_found(/*candidate_N=*/10, /*M=*/1.176280818259917, tol, db);
        CHECK(hit == 10);
    }

    SUBCASE("candidate at lower N than the only matching DB entry, EQUAL N case")
    {
        // Empty DB except a single high-N entry -- candidate at exactly N=916
        // matches because 916 <= 916. Candidate at any N < 916 does NOT.
        std::vector<reciprocal_polynomial_t> db_high;
        db_high.push_back(make_entry(916, 1.286084660801062));
        CHECK(same_polynomial_found(916, 1.286084660801062, tol, db_high) == 916);
        CHECK(same_polynomial_found(915, 1.286084660801062, tol, db_high) == 0);
        CHECK(same_polynomial_found(458, 1.286084660801062, tol, db_high) == 0);
        CHECK(same_polynomial_found(10,  1.286084660801062, tol, db_high) == 0);
    }
}


TEST_CASE("same_polynomial_by_coeffs: strict (N, coeffs) + x->-x dedup, no M dependence")
{
    // This is the new helper used by the verify-phase final gatekeeper
    // (psmm.cpp:475/526) and by merge_files_with_results.  It must NOT
    // collapse two polynomials with different (N, coeffs) just because
    // they share a Mahler measure (the bug pattern the old verify path
    // had via same_polynomial_found_m).

    std::vector<reciprocal_polynomial_t> db;
    db.push_back(make_entry(10, 1.176280818259917));   // Lehmer-like
    db[0].coeffs = {1, 1, 0, -1, -1, -1};              // overwrite half-coeffs

    SUBCASE("exact (N, coeffs) match -> 1")
    {
        std::vector<int> c = {1, 1, 0, -1, -1, -1};
        CHECK(same_polynomial_by_coeffs(10, c, db) == 1);
    }

    SUBCASE("x->-x flip of an existing entry -> 1")
    {
        // x->-x of [1,1,0,-1,-1,-1] negates odd indices -> [1,-1,0,1,-1,1]
        std::vector<int> flipped = {1, -1, 0, 1, -1, 1};
        CHECK(same_polynomial_by_coeffs(10, flipped, db) == 1);
    }

    SUBCASE("different coeffs at same N -> 0 (even if M coincides)")
    {
        // Candidate with same M but different half-coeffs (fake; the point
        // is the predicate ignores M entirely).
        std::vector<int> different = {1, 0, 0, 1, 0, 1};
        CHECK(same_polynomial_by_coeffs(10, different, db) == 0);
    }

    SUBCASE("different N -> 0 (the bug we're fixing)")
    {
        // Same coefficients but a different N is a different polynomial.
        // (Synthetic: half-coeff vector length should equal N/2+1, but the
        // predicate just compares vectors as-is -- the length mismatch
        // alone makes vector equality false.)
        std::vector<int> c = {1, 1, 0, -1, -1, -1};
        CHECK(same_polynomial_by_coeffs(12, c, db) == 0);
    }

    SUBCASE("empty polynomials list -> 0")
    {
        std::vector<reciprocal_polynomial_t> empty;
        std::vector<int> c = {1, 1, 0, -1, -1, -1};
        CHECK(same_polynomial_by_coeffs(10, c, empty) == 0);
    }

    SUBCASE("regression: distinct polynomials with coincident M at different N -> 0")
    {
        // The 2026-05-25 N=470 / N=916 incident shape: an entry at N=916 in
        // the DB shares its Mahler measure with a candidate at N=470 (which
        // is the irreducible kernel of the N=916 reducible). Coefficients
        // differ, lengths differ. same_polynomial_by_coeffs must say "not
        // a duplicate", letting the verify path add the candidate.
        std::vector<reciprocal_polynomial_t> db_high;
        reciprocal_polynomial_t hi = make_entry(916, 1.286084660801062);
        hi.coeffs.assign(459, 0);
        hi.coeffs[0] = 1; hi.coeffs[35] = 1; hi.coeffs[458] = -1;  // synthetic
        db_high.push_back(hi);

        // Candidate at N=470 with completely different (smaller) coeffs.
        std::vector<int> cand = std::vector<int>(236, 0);
        cand[0] = 1; cand[5] = -1; cand[235] = -1;
        CHECK(same_polynomial_by_coeffs(470, cand, db_high) == 0);  // <- the regression target
    }
}


TEST_CASE("same_polynomial_found_m: high-precision sibling enforces the same N <= n filter")
{
    std::vector<reciprocal_polynomial_t> db;
    db.push_back(make_entry(10, 1.176280818259917));
    db.push_back(make_entry(916, 1.286084660801062));

    mpf_class target(0, /*bits=*/128);

    SUBCASE("candidate at N=458 with M matching the N=916 entry -> NOT flagged")
    {
        target = 1.286084660801062;
        int hit = same_polynomial_found_m(/*candidate_N=*/458, target.get_mpf_t(),
                                          /*comp_precision=*/64, db);
        CHECK(hit == 0);
    }

    SUBCASE("candidate at N=20 with M matching Lehmer (N=10) -> FLAGGED")
    {
        target = 1.176280818259917;
        int hit = same_polynomial_found_m(/*candidate_N=*/20, target.get_mpf_t(),
                                          /*comp_precision=*/64, db);
        CHECK(hit == 1);
    }

    SUBCASE("candidate at N=10 with same M (same degree) -> FLAGGED")
    {
        target = 1.176280818259917;
        int hit = same_polynomial_found_m(/*candidate_N=*/10, target.get_mpf_t(),
                                          /*comp_precision=*/64, db);
        CHECK(hit == 1);
    }

    SUBCASE("candidate at N=9 (one below the only matching DB entry) -> NOT flagged")
    {
        std::vector<reciprocal_polynomial_t> db_high;
        db_high.push_back(make_entry(10, 1.176280818259917));
        target = 1.176280818259917;
        int hit = same_polynomial_found_m(9, target.get_mpf_t(), 64, db_high);
        CHECK(hit == 0);
    }
}
