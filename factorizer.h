/*
 * This file is part of "Polynomials with Small Mahler Measure" (PSMM) project.
 *
 * Copyright (C) 2020, Advanpix LLC.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Pavel Holoborodko <pavel@advanpix.com>
 */

#ifndef PSMM_POLYNOMIAL_FACTORIZER_H
#define PSMM_POLYNOMIAL_FACTORIZER_H
#include <NTL/ZZXFactoring.h>

inline bool factor_polynomial(std::vector<double>& coeffs, std::vector<std::vector<double>>& factors)
{
   //
   // Calls NTL to factorize the polynomial over Z.
   // Coefficients stored from lowers to highest: a[0],...,a[N].
   //
   // Returns true iff the polynomial is irreducible (one factor with
   // multiplicity one). When reducible, `factors` is populated with each
   // distinct factor repeated according to its multiplicity, so that the
   // caller can iterate once per "copy" and dedup by Mahler measure.
   //
   NTL::ZZX poly;
   poly.SetLength(coeffs.size());

   for(std::size_t i = 0; i < coeffs.size(); i++)
        SetCoeff(poly, i, int(coeffs[i]));

   NTL::Vec< NTL::Pair< NTL::ZZX, long > > allfactors;
   NTL::ZZ c;

   factor(c, allfactors, poly);

   for(int i = 0; i < allfactors.length(); i++)
   {
       NTL::ZZX& p          = allfactors[i].a;
       const long multiplicity = allfactors[i].b;
       const long degree       = deg(p);

       std::vector<double> f;
       f.reserve(degree + 1);
       for(long j = 0; j <= degree; j++)
           f.push_back(to_long(coeff(p,j)));

       for(long m = 0; m < multiplicity; ++m)
           factors.push_back(f);
   }

   // Irreducible iff exactly one distinct factor with multiplicity 1.
   return (allfactors.length() == 1 && allfactors[0].b == 1);
}

inline bool factor_reciprocal_polynomial(std::vector<double>& coeffs, std::vector<std::vector<double>>& factors)
{
    //
    // Convert reciprocal polynomial to full and pass it to NTL.
    //
    int N = 2 * (coeffs.size()-1);
    std::vector<double> full_coeffs(N+1);

    for(int k = 0; k <= N/2; k++) full_coeffs[k]     = coeffs[k];
    for(int k = 1; k <= N/2; k++) full_coeffs[N/2+k] = full_coeffs[N/2-k];

    //
    // Returns factors with full set of coefficients.
    //
    return factor_polynomial(full_coeffs,factors);
}

#endif // PSMM_POLYNOMIAL_FACTORIZER_H