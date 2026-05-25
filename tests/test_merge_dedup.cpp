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
