/*
 * This file is part of "Polynomials with Small Mahler Measure" (PSMM) project.
 *
 * Copyright (C) 2020, Advanpix LLC.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Pavel Holoborodko <pavel@advanpix.com>
 */

#ifndef __PSMM_POLYNOMIAL_PROPERTIES_H__
#define __PSMM_POLYNOMIAL_PROPERTIES_H__

inline void compute_all_properties_of_reciprocal_polynomial(reciprocal_polynomial_t& p, int bits = 256, int nthreads = 1)
{
    //
    // Polynomial degree and coefficients must be set before calling the function.
    //

    int N = 2 * (p.coeffs.size()-1); // = p.N;
    std::vector<int> a(N+1);

    // Expand coefficients to full polynomial
    for(int k = 0; k <= N/2; k++) a[k]     = p.coeffs[k];
    for(int k = 1; k <= N/2; k++) a[N/2+k] = a[N/2-k];

    //
    // Compute roots, their properties and Mahler measure in extended precision.
    //
    mpsolve_compute_mahler_with_properties(p.F, N, a.data(), p.K, p.U, p.Q, p.R, bits, nthreads);
    p.M = mpf_get_d(p.F);

    //
    // Compute NNZ for the half of the coefficients (excluding the a[0] = 1 that is constant).
    //
    p.nnz = 0;
    for(std::size_t k = 1; k < p.coeffs.size(); ++k) p.nnz += (p.coeffs[k] != 0);

    //
    // Compute H(p) and L(p) over the full polynomial.
    //
    p.L = p.H = 0;
    for(std::size_t k = 0; k < a.size(); k++)
    {
        const std::size_t c = static_cast<std::size_t>(std::abs(a[k]));
        p.H  = std::max(p.H,c);
        p.L += c;
    }
}

inline bool is_primitive_polynomial(const std::vector<int>& coeffs, std::vector<int>& divisors)
{
    divisors.clear();

    std::vector<int> degrees;
    for(std::size_t i = 1; i < coeffs.size(); i++) // skip zero-degree
    {
        if(coeffs[i]!=0)
            degrees.push_back(i);
    }

    if(degrees.size() > 0)
    {
        int min_degree = degrees[0];

        for(int i = 2; i <= min_degree; i++) // divisors can be [2, min_degree]
        {
            bool divisible = true;

            for(std::size_t j = 0; j < degrees.size() && divisible; j++)
            {
                divisible = ((degrees[j] % i) == 0);
            }

            if(divisible)
                divisors.push_back(i);
        }
    }

    return (divisors.size() == 0);
}

inline bool is_primitive_reciprocal_polynomial(const std::vector<int>& coeffs, std::vector<int>& divisors)
{
    int N = 2 * (coeffs.size()-1);
    std::vector<int> a(N+1);

    // Expand coefficients to full polynomial
    for(int k = 0; k <= N/2; k++) a[k]     = coeffs[k];
    for(int k = 1; k <= N/2; k++) a[N/2+k] = a[N/2-k];

    return is_primitive_polynomial(a,divisors);
}

inline void show_statistics_of_polynomials(const std::string& filename, int extended_prec, int extended_digits)
{
    //
    // Load polynomials from file.
    //
    std::vector<reciprocal_polynomial_t> poly;

    if(!filename.empty())
        load_polynomials(filename,poly,extended_prec);

    printf("-----------------------------------------------------------------\n");
    printf("Polynomials loaded from %s (%zu polynomials in total) have following properties:\n",filename.c_str(),poly.size());

    // Yes, I know, this is terrible and awfully inefficient, etc.
    // But I need this piece of code as an indication that I control my innner perfectionist.
    std::sort(poly.begin(),poly.end(),[](reciprocal_polynomial_t& a,reciprocal_polynomial_t &b) { return (mpf_cmp(a.F, b.F) < 0);});
    printf("\nPolynomials with Minimum Mahler measure:\n");

    for(std::size_t i = 0; i < std::min(std::size_t(50),poly.size()); i++)
        printp(poly[i],extended_digits);

    //
    // Compute and show some statistics
    //
    int min_diff_idx = -1;
    mpf_t d, t;
    mpf_init2(d, extended_prec);
    mpf_init2(t, extended_prec);
    mpf_set_si(d,1);

    int maxNNZ(0), maxH(0),maxL(0),maxK(0), maxU(0), maxQ(0), maxR(0);
    for(std::size_t i = 0; i < poly.size(); i++)
    {
        std::vector<int> divisors;
        bool is_primitive = is_primitive_reciprocal_polynomial(poly[i].coeffs,divisors);

        maxNNZ = (poly[maxNNZ].nnz < poly[i].nnz ? i : maxNNZ);
        maxH   = (poly[maxH].H     < poly[i].H   ? i : maxH  );
        maxL   = (poly[maxL].L     < poly[i].L   ? i : maxL  );
        maxK   = (poly[maxK].K     < poly[i].K   ? i : maxK  );
        maxU   = (poly[maxU].U     < poly[i].U   ? i : maxU  );
        maxQ   = (poly[maxQ].Q     < poly[i].Q   ? i : maxQ  );
        maxR   = (poly[maxR].R     < poly[i].R   ? i : maxR  );

        if(!is_primitive)
        {
            printf("NON-PRIMITIVE (%d): ",int(divisors.back()));
            printp(poly[i],extended_digits);
        }

        if(i+1 < poly.size())
        {
            mpf_sub(t,poly[i+1].F,poly[i].F);
            if(mpf_cmp(t,d) < 0) // t < d
            {
                mpf_set(d,t);
                min_diff_idx = i;
            }
        }
    }

    printf("\nPolynomials with nearest Mahler measures (diff = %.2e):\n", mpf_get_d(d));
    printp(poly[min_diff_idx]  ,extended_digits);
    printp(poly[min_diff_idx+1],extended_digits);

    printf("\nMaximum number of non-zero coefficients (NNZ = %zu):\n", poly[maxNNZ].nnz); for(std::size_t i = 0; i < poly.size(); i++) if(poly[maxNNZ].nnz == poly[i].nnz) printp(poly[i],extended_digits);
    printf("\nMaximum Height (H = %zu):\n",poly[maxH].H);                                 for(std::size_t i = 0; i < poly.size(); i++) if(poly[maxH].H     == poly[i].H  ) printp(poly[i],extended_digits);
    printf("\nMaximum Length (L = %zu):\n",poly[maxL].L);                                 for(std::size_t i = 0; i < poly.size(); i++) if(poly[maxL].L     == poly[i].L  ) printp(poly[i],extended_digits);
    printf("\nMaximum number of roots outside unit disk (K = %zu):\n",poly[maxK].K);      for(std::size_t i = 0; i < poly.size(); i++) if(poly[maxK].K     == poly[i].K  ) printp(poly[i],extended_digits);
    printf("\nMaximum number of roots of unity (U = %zu):\n",poly[maxU].U);               for(std::size_t i = 0; i < poly.size(); i++) if(poly[maxU].U     == poly[i].U  ) printp(poly[i],extended_digits);
    printf("\nMaximum number of complex non-unity roots (Q = %zu):\n",poly[maxQ].Q);      for(std::size_t i = 0; i < poly.size(); i++) if(poly[maxQ].Q     == poly[i].Q  ) printp(poly[i],extended_digits);
    printf("\nMaximum number of real non-unity roots (R = %zu):\n",poly[maxR].R);         for(std::size_t i = 0; i < poly.size(); i++) if(poly[maxR].R     == poly[i].R  ) printp(poly[i],extended_digits);
    printf("-----------------------------------------------------------------\n");

    mpf_clears(d,t,NULL);
}

#endif // __PSMM_POLYNOMIAL_PROPERTIES_H__