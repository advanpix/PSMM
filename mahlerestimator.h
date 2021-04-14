/*
 * This file is part of "Polynomials with Small Mahler Measure" (PSMM) project.
 *
 * Copyright (C) 2020-2021, Advanpix LLC.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Pavel Holoborodko <pavel@advanpix.com>
 */

#ifndef PSMM_MAHLER_ESTIMATOR_H
#define PSMM_MAHLER_ESTIMATOR_H

inline double estimate_mahler_reciprocal_polynomial_d(const std::vector<double>& coeffs, double threshold, int nthreads = 1)
{
    //
    // The function estimates Mahler measure of the reciprocal polynomial in double precision.
    // We use MPSolve customized to stop as soon as any of converged roots have magnitude > T.
    //
    // If this happens, then returned value = 2*T, simply indicating that M(p) is above threshold.
    // Otherwise (if all roots have magnitude < T) then M(p) computed properly.
    //
    // This function is faster than complete root-finder. It is supposed to be used in main search loop for speed.
    // After many tests, we see that function indeed faster but only for dense polynomials with large coefficients.
    // In case of sparse polynomials with small coeffs (e.g. -1,1) - it doesn't show speed improvement.
    // That is because such polynomials tend to have a lot of small roots (|z| < T) but overall M(p) = prod(|z|) > T.
    //
    // Following functions were implemented as part of MPSolve (otherwise linker was producing incorrect code):
    //
    //      mps_mahler_is_over_threshold
    //      mps_mahler_mpsolve
    //      mps_mahler_set_threshold
    //

    int n = 2 * (coeffs.size()-1); // Degree of the polynomial.

    mps_context* s = mps_context_new();
    mps_polynomial* poly = (mps_polynomial*)mps_monomial_poly_new(s,n);

    for(int i = 0; i <= n/2; i++) mps_monomial_poly_set_coefficient_int(s,(mps_monomial_poly*)poly,    i,coeffs[i    ],0);
    for(int i = 1; i <= n/2; i++) mps_monomial_poly_set_coefficient_int(s,(mps_monomial_poly*)poly,n/2+i,coeffs[n/2-i],0);

    //
    // We use full 64-bit limb in GMP to get "double" precision (= 53 bits for mantissa).
    // This pushes MPSolve to use GMP internally, and hence roots are stored in mpc_t format.
    // We can get them directly by mps_context_get_roots_m.
    // There is no speed merit in using the mps_context_get_roots_d (roots are stored in mpc_t format anyway).
    //
    const int working_precision = 64;

    mps_context_set_input_poly            (s, poly);
    mps_context_set_output_prec           (s, working_precision);
    mps_context_set_output_goal           (s, MPS_OUTPUT_GOAL_APPROXIMATE);
    mps_context_select_algorithm          (s, MPS_ALGORITHM_STANDARD_MPSOLVE);
    mps_thread_pool_set_concurrency_limit (s, NULL, nthreads);

    mps_mahler_set_threshold(threshold);
    mps_mahler_mpsolve(s);

    double mahler (0.0);
    if(mps_mahler_is_over_threshold())
    {
        mahler = 2*threshold; // out-of threshold value.
    }
    else
    {
        mahler = 1.0;
        if(working_precision > 53)
        {
            //
            // In this case MPSolve uses GMP anyway, so it is better to work with multiprecision roots directly.
            // As a result, computed Mahler measure has (much) higher accuracy, (only last digit is a bit off).
            //
            mpf_t m, r, t;
            mpf_init2(m,working_precision);
            mpf_init2(t,working_precision);
            mpf_init2(r,working_precision);

            mpc_t *results = (mpc_t*) std::malloc(n*sizeof(mpc_t));
            mpc_vinit2(results,n,working_precision);

            mps_context_get_roots_m(s, &results, NULL);

            mpf_set_si(m,1);
            for(int i = 0; i < n; i++)
            {
                mpf_t& x = mpc_Re(results[i]);
                mpf_t& y = mpc_Im(results[i]);

                // Fast reject using Manhattan distance (to avoid slow sqrt):
                //
                //             |x|+|y| <= sqrt(x^2+y^2)
                //
                mpf_abs(x,x);
                mpf_abs(y,y);
                mpf_add(r,x,y);

                if(mpf_cmp_si(r,1) > 0)
                {
                    mpf_mul (t, x, x); // t = x^2
                    mpf_mul (r, y, y); // r = y^2
                    mpf_add (r, r, t); // r = x^2+y^2
                    mpf_sqrt(r, r);    // r = sqrt(x^2+y^2)

                    if(mpf_cmp_si(r,1) > 0) mpf_mul(m,m,r); // compute M(p) if |z| > 1
                }
            }

            mahler = mpf_get_d(m);

            mpf_clears(m,t,r,NULL);
            mpc_vclear(results,n);
            std::free(results);
        }
        else
        {
            //
            // Work in double-precision.
            //
            cplx_t *results = (cplx_t*) std::malloc(n*sizeof(cplx_t));
            mps_context_get_roots_d(s, &results, NULL);
            for(int i = 0; i < n; i++)
            {
                double r = std::hypot(cplx_Re(results[i]),cplx_Im(results[i]));

                if(r > 1)
                    mahler *= r;
            }
            std::free(results);
        }
    }

    mps_polynomial_free (s, poly);
    mps_context_free (s);

    return mahler;
}

#endif // PSMM_MAHLER_ESTIMATOR_H
