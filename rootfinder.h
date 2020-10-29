/*
 * This file is part of "Polynomials with Small Mahler Measure" (PSMM) project.
 *
 * Copyright (C) 2020, Advanpix LLC.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Pavel Holoborodko <pavel@advanpix.com>
 */

#ifndef PSMM_ROOT_FINDER_H
#define PSMM_ROOT_FINDER_H
#include <mps/mps.h>

//
// MPSolve doesn't provide the function in public interface, so we just forward declare it here.
//
extern "C" void mps_thread_pool_set_concurrency_limit(void * s, void * pool, unsigned int concurrency_limit);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Double-precision computations
//

//
// One-shot simple version with re-allocations on each call.
//
inline int mpsolve_find_all_roots(int n, const double* coeffs, double* wr, double* wi, double* residuals, int digits = 15, int nthreads = 1)
{
    int error = 0;
    const double log2_10 = 3.3219280948873624;

    mps_context* s = mps_context_new();
    mps_polynomial* poly = (mps_polynomial*)mps_monomial_poly_new(s,n);

    for(int i = 0; i <= n; i++)
        mps_monomial_poly_set_coefficient_d(s,(mps_monomial_poly*)poly,i,coeffs[i],0.0);

    mps_context_set_input_poly            (s, poly);
    mps_context_set_output_prec           (s, std::ceil(digits * log2_10));
    mps_context_set_output_goal           (s, MPS_OUTPUT_GOAL_APPROXIMATE);
    mps_context_select_algorithm          (s, MPS_ALGORITHM_STANDARD_MPSOLVE);
    mps_thread_pool_set_concurrency_limit (s, NULL, nthreads);

    mps_mpsolve(s);

    cplx_t *results = (cplx_t*) std::malloc(n*sizeof(cplx_t));
    mps_context_get_roots_d (s, &results, &residuals);
    for(int i = 0; i < n; i++)
    {
        wr[i] = cplx_Re(results[i]);
        wi[i] = cplx_Im(results[i]);
    }
    std::free(results);

    error = mps_context_has_errors(s);

    mps_polynomial_free (s, poly);
    mps_context_free (s);

    return error;
}

//
// Computes Mahler measure of polynomial with the accuracy specified.
// Also returns roots (real & imaginary parts) and residuals for each root.
//
inline double mpsolve_compute_mahler(int n, const double* coeffs, double* wr, double* wi, double* residuals, int digits = 15, int nthreads = 1)
{
    double mahler = 0;
    int error = mpsolve_find_all_roots(n,coeffs,wr,wi,residuals,digits,nthreads);

    if(error == 0)
    {
        mahler = 1.0;
        for(int i = 0; i < n; i++)
        {
            double r = std::hypot(wr[i],wi[i]);

            if(r > 1)
                mahler *= r;
        }
    }

    return mahler;
}

inline double compute_mahler_general_polynomial(const std::vector<double>& coeffs, int digits = 15, int nthreads = 1)
{
    int N = (coeffs.size()-1);

    double* WR  = (double*)std::malloc(N*sizeof(double));
    double* WI  = (double*)std::malloc(N*sizeof(double));
    double* RE  = (double*)std::malloc(N*sizeof(double));

    double mahler = mpsolve_compute_mahler(N,coeffs.data(), WR, WI, RE, digits, nthreads);

    std::free(WR);
    std::free(WI);
    std::free(RE);

    return mahler;
}

inline double compute_mahler_reciprocal_polynomial(const std::vector<double>& coeffs, int digits = 15, int nthreads = 1)
{
    int N = 2 * (coeffs.size()-1);
    std::vector<double> a(N+1);

    // Expand coefficients to full polynomial
    for(int k = 0; k <= N/2; k++) a[k]     = coeffs[k];
    for(int k = 1; k <= N/2; k++) a[N/2+k] = a[N/2-k];

    return compute_mahler_general_polynomial(a, digits, nthreads);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Extended-precision computations
//


//
// One-shot simple version with re-allocations on each call.
//
inline int mpsolve_compute_mahler(mpf_ptr m, int n, const double* coeffs, int bits = 256, int nthreads = 1)
{
    int error = 0;

    mps_context* s = mps_context_new();
    mps_polynomial* poly = (mps_polynomial*)mps_monomial_poly_new(s,n);

    for(int i = 0; i <= n; i++)
        mps_monomial_poly_set_coefficient_d(s,(mps_monomial_poly*)poly,i,coeffs[i],0.0);

    mps_context_set_input_poly            (s, poly);
    mps_context_set_output_prec           (s, bits);
    mps_context_set_output_goal           (s, MPS_OUTPUT_GOAL_APPROXIMATE);
    mps_context_select_algorithm          (s, MPS_ALGORITHM_STANDARD_MPSOLVE);
    mps_thread_pool_set_concurrency_limit (s, NULL, nthreads);

    mps_mpsolve(s);
    error = mps_context_has_errors(s);

    if(error == 0)
    {
        mpc_t *results = (mpc_t*) std::malloc(n*sizeof(mpc_t));
        mpc_vinit2(results,n,bits);

        mps_context_get_roots_m(s, &results, NULL);

        mpf_t mahler, rabs, temp;
        mpf_init2(mahler,bits);
        mpf_init2(rabs,  bits);
        mpf_init2(temp,  bits);

        mpf_set_si(mahler,1);
        for(int i = 0; i < n; i++)
        {
            // rabs = sqrt(x^2+y^2) - magnitude of the root.
            mpf_mul (temp, mpc_Re(results[i]), mpc_Re(results[i]));
            mpf_mul (rabs, mpc_Im(results[i]), mpc_Im(results[i]));
            mpf_add (rabs, rabs, temp);
            mpf_sqrt(rabs, rabs);

            if(mpf_cmp_si(rabs,1) > 0)
                mpf_mul(mahler,mahler,rabs);
        }

        mpf_init2(m,mpf_get_prec(mahler));
        mpf_set(m,mahler);

        mpf_clears(mahler,rabs,temp,NULL);

        mpc_vclear(results,n);
        std::free(results);
    }

    mps_polynomial_free (s, poly);
    mps_context_free(s);

    return error;
}

inline void compute_mahler_general_polynomial(mpf_ptr m, const std::vector<double>& coeffs, int bits = 256, int nthreads = 1)
{
    int N = (coeffs.size()-1);
    mpsolve_compute_mahler(m,N,coeffs.data(),bits, nthreads);
}

inline void compute_mahler_reciprocal_polynomial(mpf_ptr m, const std::vector<double>& coeffs, int bits = 256, int nthreads = 1)
{
    int N = 2 * (coeffs.size()-1);
    std::vector<double> a(N+1);

    // Expand coefficients to full polynomial
    for(int k = 0; k <= N/2; k++) a[k]     = coeffs[k];
    for(int k = 1; k <= N/2; k++) a[N/2+k] = a[N/2-k];

    compute_mahler_general_polynomial(m, a, bits, nthreads);
}

#endif // PSMM_ROOT_FINDER_H
