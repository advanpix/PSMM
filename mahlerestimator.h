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

#include <cmath>  // std::nan

// Thread-local hooks exported by mpsolve patch
// mpsolve/patches/0003-advanpix-defensive-mpf-set-rdpe.patch. We
// REGISTER the polynomial being processed before each mps_*_mpsolve
// call so the patch's diagnostic in mpf_set_rdpe can print the
// offending polynomial, and we READ the failure flag after the call
// to decide whether the mpsolve result is trustworthy.
extern "C" {
    extern __thread const int *  mps_advanpix_current_poly_coeffs;
    extern __thread int          mps_advanpix_current_poly_coeffs_size;
    extern __thread int          mps_advanpix_current_poly_n;
    extern __thread const char * mps_advanpix_current_poly_tag;
    extern __thread int          mps_advanpix_computation_failed;
}

inline double estimate_mahler_reciprocal_polynomial_d(const std::vector<int>& coeffs, double threshold, int nthreads = 1, mps_context* reuse_ctx = nullptr)
{
    //
    // The function estimates Mahler measure of the reciprocal polynomial in double precision.
    // We use MPSolve customized to stop as soon as any of converged roots have magnitude > T.
    //
    // If this happens, then returned value = 2*T, simply indicating that M(p) is above threshold.
    // Otherwise (if all roots have magnitude < T) then M(p) computed properly.
    //
    // If MPSolve hits the non-finite-rdpe bug (Advanpix patch 0003
    // detects it and sets mps_advanpix_computation_failed), the
    // function returns NaN. The caller MUST check std::isnan(...)
    // and treat such polynomials as uncomputable (skip + log) rather
    // than as below-threshold finds.
    //
    // When reuse_ctx != nullptr, the caller owns the context (pre-configured with output_prec / goal /
    // algorithm / concurrency). The function only creates+frees the polynomial. The threshold is set
    // per call via TLS (mps_mahler_set_threshold).
    //

    int n = 2 * (coeffs.size()-1); // Degree of the polynomial.

    const bool own_context = (reuse_ctx == nullptr);
    mps_context* s = own_context ? mps_context_new() : reuse_ctx;

    mps_polynomial* poly = (mps_polynomial*)mps_monomial_poly_new(s,n);

    for(int i = 0; i <= n/2; i++) mps_monomial_poly_set_coefficient_int(s,(mps_monomial_poly*)poly,    i,coeffs[i    ],0);
    for(int i = 1; i <= n/2; i++) mps_monomial_poly_set_coefficient_int(s,(mps_monomial_poly*)poly,n/2+i,coeffs[n/2-i],0);

    const int working_precision = 64;

    mps_context_set_input_poly(s, poly);

    if(own_context)
    {
        mps_context_set_output_prec           (s, working_precision);
        mps_context_set_output_goal           (s, MPS_OUTPUT_GOAL_APPROXIMATE);
        mps_context_select_algorithm          (s, MPS_ALGORITHM_STANDARD_MPSOLVE);
        mps_thread_pool_set_concurrency_limit (s, NULL, nthreads);
    }

    // Register polynomial for the mpsolve diagnostic + reset failure flag.
    // We register the HALF-coefficient vector (what PSMM keeps in memory);
    // the diagnostic tag tells offline reproducer scripts to expand half->full
    // before re-running.
    mps_advanpix_current_poly_coeffs       = coeffs.data();
    mps_advanpix_current_poly_coeffs_size  = (int) coeffs.size();
    mps_advanpix_current_poly_n            = n;
    mps_advanpix_current_poly_tag          = "estimate_mahler_d:half";
    mps_advanpix_computation_failed        = 0;

    mps_mahler_set_threshold(threshold);
    mps_mahler_mpsolve(s);

    // If the non-finite-rdpe path fired during this call, the mpsolve
    // result is corrupted. Return NaN so the caller knows to skip +
    // report the polynomial rather than use a wrong M(p).
    if(mps_advanpix_computation_failed)
    {
        mps_polynomial_free (s, poly);
        if(own_context) mps_context_free (s);
        return std::nan("");
    }

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
    if(own_context) mps_context_free (s);

    return mahler;
}

#endif // PSMM_MAHLER_ESTIMATOR_H
