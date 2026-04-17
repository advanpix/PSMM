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

//
// MPSolve doesn't provide the function in public interface, so we just forward declare it here. Needed iif _MPS_PRIVATE is not defined.
//

#ifndef _MPS_PRIVATE
    extern "C" void mps_thread_pool_set_concurrency_limit(void * s, void * pool, unsigned int concurrency_limit);
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Double-precision computations
//

//
// Computes Mahler measure of the (general) polynomial with double precision accuracy.
//

inline double compute_mahler_general_polynomial_d(const std::vector<int>& coeffs, int nthreads = 1, mps_context* reuse_ctx = nullptr)
{
    //
    // Find polynomial roots, compute Mahler measure for double precision.
    //
    // When reuse_ctx != nullptr the caller owns the context and is responsible
    // for pre-setting output_prec / output_goal / algorithm / concurrency.
    // The function only creates+frees the polynomial object. This eliminates
    // the mps_context_new/free overhead in the hot search loop (C1).
    //

    int n = (coeffs.size()-1);  // Degree of the polynomial.

    const bool own_context = (reuse_ctx == nullptr);
    mps_context* s = own_context ? mps_context_new() : reuse_ctx;

    mps_polynomial* poly = (mps_polynomial*)mps_monomial_poly_new(s,n);

    for(int i = 0; i <= n; i++)
        mps_monomial_poly_set_coefficient_int(s,(mps_monomial_poly*)poly,i,coeffs[i],0);

    const int working_precision = 64;

    mps_context_set_input_poly(s, poly);

    if(own_context)
    {
        mps_context_set_output_prec           (s, working_precision);
        mps_context_set_output_goal           (s, MPS_OUTPUT_GOAL_APPROXIMATE);
        mps_context_select_algorithm          (s, MPS_ALGORITHM_STANDARD_MPSOLVE);
        mps_thread_pool_set_concurrency_limit (s, NULL, nthreads);
    }

    mps_mpsolve(s);

    double mahler = 1.0;
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

    mps_polynomial_free (s, poly);
    if(own_context) mps_context_free (s);

    return mahler;
}

inline double compute_mahler_reciprocal_polynomial_d(const std::vector<int>& coeffs, int nthreads = 1, mps_context* reuse_ctx = nullptr)
{
    //
    // The function computes Mahler measure of the reciprocal polynomial in double precision.
    // It is supposed to be fast, as it is being used in main loop over all polynomials.
    //
    // There are many ways and ideas on how to compute roots of palindromic polynomials. At the moment we rely on the basic ones:
    //
    // (A) Convert palindromic polynomial to generic and use root-finder to find all roots of generic polynomial.
    //     MPSolve is very good at finding roots of generic sparse polynomials - fast and stable.
    //
    // (B) Leverage the property that palindromic polynomial can be expressed in Chebyshev basis of N/2 degree.
    //     This requires special root finder for the Chebyshev expansion (e.g. eigenvalues of colleague matrix, see Trefethen approximation book).
    //     Then roots can be converted to the roots of original palindromic polynomial.
    //     Root finder for Chebyshev basis must be optimized for sparse expansions as well.
    //
    //     MPSolve supports Chebyshev expansions, but its stability is questionable in this mode.
    //     Only MPS_ALGORITHM_SECULAR_GA is supported, and output results are on lower accuracy (at least in double precision).
    //     But most importantly, MPSolve doesn't support sparse Chebyshev expansions out-of-the-box (we can add it by providing our custom polynomial class).
    //     This leads to x2 times slow-down compared to (A).

//#define PSMM_USE_CHEBYSHEV_ROOT_FINDER

#if !defined(PSMM_USE_CHEBYSHEV_ROOT_FINDER)

    //
    // (A) The most stable way in MPSolve 3.2.1
    //
    int N = 2 * (coeffs.size()-1);
    std::vector<int> a(N+1);

    // Expand coefficients to full polynomial
    for(int k = 0; k <= N/2; k++) a[k]     = coeffs[k];
    for(int k = 1; k <= N/2; k++) a[N/2+k] = a[N/2-k];

    return compute_mahler_general_polynomial_d(a, nthreads, reuse_ctx);
#else

    //
    // (B) Suffers from accuracy loss in MPSolve 3.2.1
    //

    // Find roots of Chebyshev expansion.
    // Convert them to the roots of original polynomial.
    // Compute Mahler measure.

    int m = 2 * (coeffs.size()-1);   // Palindromic polynomial degree.
    int n = m / 2;                   // Chebyshev expansion degree. Input vector "coeffs" includes coefficients [0..n] (n+1 in total).
    double mahler = 1.0;

    mps_context* s = mps_context_new();
    mps_chebyshev_poly* poly = mps_chebyshev_poly_new(s,n,MPS_STRUCTURE_REAL_INTEGER);

    for(int i = 0; i <= n; i++)
        mps_chebyshev_poly_set_coefficient_i(s,poly,i,coeffs[n-i],0); //  = palindromic coeffs in reversed order.

    //
    // We use full 64-bit limb in GMP to get "double" precision (= 53 bits for mantissa).
    // This pushes MPSolve to use GMP internally, and hence roots are stored in mpc_t format.
    // We can get them directly by mps_context_get_roots_m.
    // There is no speed merit in using the mps_context_get_roots_d (roots are stored in mpc_t format anyway).
    //
    const int precision = 64;

    mps_context_set_input_poly            (s, MPS_POLYNOMIAL(poly));
    mps_context_set_output_prec           (s, precision);
    mps_context_set_output_goal           (s, MPS_OUTPUT_GOAL_APPROXIMATE);
    mps_context_select_algorithm          (s, MPS_ALGORITHM_SECULAR_GA); // mps_mpsolve crashes on MPS_ALGORITHM_STANDARD_MPSOLVE
    mps_thread_pool_set_concurrency_limit (s, NULL, nthreads);

    mps_mpsolve(s);

#if 0
    mpc_t *results = (mpc_t*) std::malloc(n*sizeof(mpc_t));
    mpc_vinit2(results,n,precision);

    mps_context_get_roots_m(s, &results, NULL);

    for(int i = 0; i < n; i++)
    {
        wr[i] = mpf_get_d(mpc_Re(results[i]));
        wi[i] = mpf_get_d(mpc_Im(results[i]));
    }

    mpc_vclear(results,n);
    std::free(results);
#endif

    //
    // Double-precision version. It is natural that is has much lower accuracy, especially for high degree polynomials.
    // However this fails even on lower degrees, e.g. 10.
    //
    cplx_t *results = (cplx_t*) std::malloc(n*sizeof(cplx_t));
    mps_context_get_roots_d (s, &results, NULL);
    for(int i = 0; i < n; i++)
    {
        std::complex<double> z(cplx_Re(results[i]),cplx_Im(results[i]));
        std::complex<double> x0 = (z+sqrt(z*z-4.0))/2.0;
        std::complex<double> x1 = (z-sqrt(z*z-4.0))/2.0;

        double r0 = std::hypot(x0.real(),x0.imag());
        if(r0 > 1) mahler *= r0;

        double r1 = std::hypot(x1.real(),x1.imag());
        if(r1 > 1) mahler *= r1;
    }
    std::free(results);

    mps_polynomial_free(s, MPS_POLYNOMIAL(poly));
    mps_context_free (s);

    return mahler;
#endif

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Extended-precision computations
//


//
// One-shot simple version with re-allocations on each call.
//
inline int mpsolve_compute_mahler_with_properties(mpf_ptr mahler,                   // Mahler measure of the polynomial
                                                  int n,                            // Polynomial degree
                                                  const int* coeffs,                // Full polynomial coefficients: a[0],...,a[n]
                                                  std::size_t& K,                   // Number of roots outside the unit circle.
                                                  std::size_t& U,                   // Number of complex unity roots (go in pairs):           z = exp(i*t), z* = exp(-i*t).
                                                  std::size_t& Q,                   // Number of complex non-unity roots (go in quadruplets): z = r*exp(i*t), z* = r*exp(-i*t), 1/z = (1/r)*exp(-i*t), 1/z* = (1/r)*exp(i*t).
                                                  std::size_t& R,                   // Number of real non-unity roots (go in pairs):          z = r, z = 1/r.
                                                  int target_precision = 256,       // Precision to use for root finding and analysis.
                                                  int nthreads = 1                  // Number of threads to use in root finder.
                                                 )
{
    int error = 0;

    mps_context* s = mps_context_new();
    mps_polynomial* poly = (mps_polynomial*)mps_monomial_poly_new(s,n);

    for(int i = 0; i <= n; i++)
        mps_monomial_poly_set_coefficient_int(s,(mps_monomial_poly*)poly,i,coeffs[i],0);

    //
    // MPSolve uses GMP as extended precision engine.
    // Internally GMP does calculations with precision = number of full limbs. So that, if we request 230 bits of precision, GMP will use full 256 anyway.
    // Therefore, it is better to use precision = integer multiplier of 64 (= size of limb on x64 CPU).
    //
    // Moreover, we want to make sure roots are accurate to the full precision requested (target_precision).
    // Hence we add one more extra-limb of precision ("guard precision" so to speak). This ensures MPSolve provides high-accuracy roots.
    //
    int extra_precision = (std::ceil(target_precision/double(64)) + 1) * 64;

    mps_context_set_input_poly            (s, poly);
    mps_context_set_output_prec           (s, extra_precision);
    mps_context_set_output_goal           (s, MPS_OUTPUT_GOAL_APPROXIMATE);
    mps_context_select_algorithm          (s, MPS_ALGORITHM_STANDARD_MPSOLVE);
    mps_thread_pool_set_concurrency_limit (s, NULL, nthreads);

    mps_mpsolve(s);
    error = mps_context_has_errors(s);

    if(error == 0)
    {
        K = 0;  // Number of roots outside the unit circle.
        U = 0;  // Number of complex unity roots (go in pairs):           z = exp(i*t), z* = exp(-i*t).
        Q = 0;  // Number of complex non-unity roots (go in quadruplets): z = r*exp(i*t), z* = r*exp(-i*t), 1/z = (1/r)*exp(-i*t), 1/z* = (1/r)*exp(i*t).
        R = 0;  // Number of real non-unity roots (go in pairs):          z = r, z = 1/r.

        //
        // We do all the computations with extra_precision
        // Only final result is returned with target_precision.
        //
        mpc_t *results = (mpc_t*) std::malloc(n*sizeof(mpc_t));
        mpc_vinit2(results,n,extra_precision);

        mps_context_get_roots_m(s, &results, NULL);

        mpf_t m, r, t, d, eps;
        mpf_init2(m,  extra_precision);
        mpf_init2(t,  extra_precision);
        mpf_init2(d,  extra_precision);
        mpf_init2(r,  extra_precision);
        mpf_init2(eps,extra_precision);

        //
        // Compute machine epsilon for the requested precision, machine epsilon = 2^-(target_precision-1)
        // Machine epsilon is used as tolerance in checking the closeness of floating point numbers.
        // Please note, computed roots have higher precision (see above), but we operate with initially requested precision "target_precision".
        //
        mpf_set_ui(eps,1);
        mpf_div_2exp(eps,eps,std::max(1,(target_precision-1)));

        mpf_set_si(m,1);
        for(int i = 0; i < n; i++)
        {
            //
            // Root magnitude: r = sqrt(x^2+y^2)
            // We use naive textbook algorithm, but we will have to implement proper hypot in future.
            //
            mpf_mul (t, mpc_Re(results[i]), mpc_Re(results[i])); // t = x^2
            mpf_mul (r, mpc_Im(results[i]), mpc_Im(results[i])); // r = y^2
            mpf_add (r, r, t);                                   // r = x^2+y^2
            mpf_sqrt(r, r);                                      // r = sqrt(x^2+y^2)

            mpf_sub_ui(d,r,1);
            mpf_abs(d,d);
            if(mpf_cmp(d,eps) <= 0)  // ||z|-1| <= eps
            {
                U++; // unity root
            }
            else
            {
                if(mpf_cmp_si(r,1) > 0) // |z| > 1 && ||z|-1| > eps
                {
                    mpf_mul(m,m,r); // use for Mahler measure computation.
                    K++;
                }

                mpf_abs(d,mpc_Im(results[i]));
                if(mpf_cmp(d,eps) < 0) // |Im(z)| < eps
                {
                    R++; // real non-unity root
                }
                else
                {
                    Q++; // complex non-unity root
                }
            }
        }

        // Caller owns `mahler`'s allocation and precision (see reciprocal_polynomial_t::F /
        // mpf_class). Just copy the computed value into the caller's storage. If the caller's
        // precision is lower than target_precision, GMP rounds to the caller's precision.
        mpf_set(mahler,m);

        mpf_clears(m,r,t,d,eps,NULL);

        mpc_vclear(results,n);
        std::free(results);
    }

    mps_polynomial_free (s, poly);
    mps_context_free(s);

    return error;
}

#endif // PSMM_ROOT_FINDER_H
