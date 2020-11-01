/*
 * This file is part of "Polynomials with Small Mahler Measure" (PSMM) project.
 *
 * Copyright (C) 2020, Advanpix LLC.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Pavel Holoborodko <pavel@advanpix.com>
 */

#ifndef PSMM_MAIN_HEADER_H
#define PSMM_MAIN_HEADER_H

#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <vector>
#include <map>
#include <algorithm>
#include <gmp.h>

typedef struct{

    // Properties of the roots
    std::size_t K;              // Number of roots outside the unit circle.
    std::size_t U;              // Number of complex unity roots (go in pairs):           z = exp(i*t), z* = exp(-i*t).
    std::size_t Q;              // Number of complex non-unity roots (go in quadruplets): z = r*exp(i*t), z* = r*exp(-i*t), 1/z = (1/r)*exp(-i*t), 1/z* = (1/r)*exp(i*t).
    std::size_t R;              // Number of real non-unity roots (go in pairs):          z = r, z = 1/r.

}roots_properties_t;

typedef struct{

    // Main properties
    std::size_t N;              // Polynomial degree.
    std::size_t nnz;            // Number of non-zero coefficients among a[1]..a[N/2]. The a[0] = 1 is fixed hence not counted.
    std::size_t H;              // Polynomial height (maximum of the magnitudes of its coefficients).
    std::size_t L;              // Polynomial length (sum of the magnitudes of the coefficients).
    std::vector<double> coeffs; // Polynomial coefficients. In case of reciprocal polynomials stores only half of coefficients a[0],...,a[N/2]

    // Properties of the roots
    roots_properties_t r;

    // Mahler measure
    mpf_t  F;                   // Mahler measure in extended precision (might not be available).
    double M;                   // Mahler measure in double precision.

}polynomial_t;

#endif // PSMM_MAIN_HEADER_H