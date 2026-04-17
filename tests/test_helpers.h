/*
 * This file is part of "Polynomials with Small Mahler Measure" (PSMM) project.
 *
 * Copyright (C) 2026, Advanpix LLC.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 */

#ifndef PSMM_TEST_HELPERS_H
#define PSMM_TEST_HELPERS_H

//
// Small shared utilities for PSMM tests.
//
// Reads the reference file tests/data/reference-polynomials.txt into a
// vector of structs. Each line in the file follows the AllKnownAdvanpix
// format:
//
//   N M NNZ H L K U Q R c_0 c_1 ... c_{N/2}
//
// where coefficients are the reciprocal half-representation (a[0] = 1
// then a[1] .. a[N/2]).
//

#include "psmm.h"
#include "utilities.h"

#include <string>
#include <vector>
#include <fstream>
#include <sstream>

struct ref_polynomial_t {
    int                N;         // degree
    std::string        M_string;  // Mahler measure as text (70+ digits)
    double             M;         // Mahler measure as double
    std::vector<int>   coeffs;    // half-representation, a[0] .. a[N/2]
};

inline std::vector<ref_polynomial_t> load_reference_polynomials(const std::string& path)
{
    std::vector<ref_polynomial_t> out;
    std::ifstream ifs(path);
    if (!ifs.is_open()) return out;

    std::string line;
    while (std::getline(ifs, line)) {
        // Strip leading whitespace; skip comments and blanks.
        std::size_t first = line.find_first_not_of(" \t");
        if (first == std::string::npos) continue;
        if (line[first] == '#') continue;

        std::istringstream iss(line.substr(first));
        ref_polynomial_t p;

        if (!(iss >> p.N >> p.M_string)) continue;
        p.M = std::stod(p.M_string);

        // Skip NNZ H L K U Q R (7 fields).
        int dummy;
        for (int k = 0; k < 7; ++k) iss >> dummy;

        int expected = p.N / 2 + 1;
        p.coeffs.reserve(expected);
        int c;
        while (iss >> c) p.coeffs.push_back(c);

        if (static_cast<int>(p.coeffs.size()) == expected)
            out.push_back(std::move(p));
    }
    return out;
}

// Convert the half-representation to a full (0..N) coefficient vector.
inline std::vector<int> expand_reciprocal(const std::vector<int>& half, int N)
{
    std::vector<int> full(N + 1, 0);
    for (int k = 0; k <= N / 2; ++k) full[k] = half[k];
    for (int k = 1; k <= N / 2; ++k) full[N / 2 + k] = full[N / 2 - k];
    return full;
}

#endif // PSMM_TEST_HELPERS_H
