/*
 * This file is part of "Polynomials with Small Mahler Measure" (PSMM) project.
 *
 * Copyright (C) 2020, Advanpix LLC.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Pavel Holoborodko <pavel@advanpix.com>
 */

#include "psmm.h"
#include "utilities.h"
#include "arguments.h"
#include "rootfinder.h"
#include "factorizer.h"
#include "polyiterator.h"
#include "polyhelpers.h"
#include "polyproperties.h"
#include "polyconvert.h"

#ifdef USE_FAST_MAHLER_ESTIMATOR
#include "mahlerestimator.h"
#endif

// References.
//
// [Br51] R. Breusch, On the distribution of the roots of a polynomial with integral coefficients, Proc. Amer. Math. Soc. 2 (1951), 939–941
// [Sm71] C.J. Smyth, On the product of the conjugates outside the unit circle of an algebraic integer, Bull. London Math. Soc. 3 (1971), 169–175.


int main(int argc, char* argv[])
{
    //
    // The list of known polynomials with Mahler measure < 1.3 ('Known180' [1]) stores Mahler measure with 13 digits (12 for fractional part).
    // Hence, if you use Known180 then 'search_tolerance' must be >= 1e-12.
    //
    // This might be too high for high-order polynomials since their Mahler measures are clustered more densely.
    //
    // EXAMPLE 1. There are two different polynomials of 138 and 180 degrees which
    //     have the same Mahler measure up to 11 digits after the comma:
    //
    //     180 1.25282688286[4]667557461368123141605672483501189167102619478483990743452812 40 1 81 20 140 40 0 1 1 1 0 0 -1 -1 -1 0 0 0 0 0 0 0 1 1 1 0 0 -1 -1 -1 0 0 0 0 0 0 0 1 1 1 1 1 0 -1 -1 -1 -1 -1 0 0 0 0 1 1 1 1 1 0 -1 -1 -1 -1 -1 0 0 0 0 1 0 0 0 1 0 0 0 0 -1 -1 0 0 0 0 1 0 0 0 1 0 0 0 0 -1 -1 0 0 0 0 1
    //     138 1.25282688286[8]347583989529787618964056215870696871794089089602385130202388 3 1 7 15 108 28 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1
    //
    //     If we would store these values with the accuracy of 12 digits only (Known180) then these two values would be considered the same for high tolerances like 1e-10 or 1e-11.
    //
    // EXAMPLE 2.
    //
    //     214 1.2857117920[6]7308994496869974333734988823588105122313807358079920158251959 2 1 5 34 146 68 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1
    //     158 1.2857117920[9]6690828765613523071678857629755985708260720994246701583641435 2 1 5 26 106 52 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1
    //
    // EXAMPLE 3.
    //
    //     196 1.285872662[7]28973405119892725557620131473130201981020148720083571083521939 2 1 5 27 142 52 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1
    //     216 1.285872662[8]03832621268592047312148080780963659492447295128617130559423450 48 1 97 36 144 72 0 1 -1 1 0 -1 1 -1 0 0 0 0 0 1 -1 1 0 -1 1 -1 0 0 0 0 0 1 -1 1 0 -1 1 -1 0 0 0 0 0 1 -1 1 0 -1 1 -1 0 0 0 0 0 1 -1 1 0 -1 1 -1 0 0 0 0 0 1 -1 1 0 -1 1 -1 0 0 0 0 0 1 0 0 1 -1 0 0 -1 0 0 0 0 1 0 0 1 -1 0 0 -1 0 0 0 0 1 0 0 1 -1 0 0 -1 0 0 0 0 1
    //
    // EXAMPLE 4.
    //
    //     286 1.28535091[1]978298816872119209775780481172687954632768373236970316904376274 2 1 5 47 192 92 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 -1
    //     286 1.28535091[2]052405489247282361244720981840042015246864037398332189073656661 2 1 5 40 206 80 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
    //
    // EXAMPLE 5.
    //
    //      62 1.2858632072[8]7390566353658372190483988074074758474638202101410392136823271 18 1 37 11 40 20 2 1 1 1 0 -1 -1 -1 0 0 0 0 0 1 1 1 0 -1 -1 -1 0 0 0 0 0 1 1 1 0 -1 -1 -1 -1
    //     270 1.2858632072[3]2419360063269405987105238227103026953231638400082996792034877 2 1 5 29 212 56 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1
    //
    // EXAMPLE 6.
    //
    //     274 1.2860926270[2]6747194069236227273054797341756306778007487778481200532426484 2 1 5 33 208 64 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1
    //     284 1.2860926270[1]2435716829466235652543676338041317848010311065680690084257741 2 1 5 45 194 88 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1
    //
    // Expectedly, there are many such polynomials probably with even smaller difference.
    // Therefore we do not recommend using Known180 because of limited number of digits stored.
    //
    // Instead, we suggest to use extended precision for Mahler measure.
    // With this purpose, we re-computed Mahler measure of all previousely known polynomials with 72 digits of accuracy (octuple precision).
    // File 'AllKnown' provides these refined values - it must be used in all calculations.
    // In this case 'search_tolerance' can be as small as 1e-14.
    //
    // Still, there is a chance of loosing some polynomial with Mahler measure difference < 1e-14, but we must balance between double precision (speed) and extended precision (accuracy).
    //
    // References.
    // 1. http://www.cecm.sfu.ca/~mjm/Lehmer/lists/
    //
    const double search_tolerance = 1e-14; // We use 64-bits of precision in MPSolve, so that our roots must have full "double" precision accuracy (53-bits).
                                           // From our experiments, 1e-14 is the lowest correct tolerance for "double" precision. It is good enough for now, but in future we need to use true extended precision for search.

    const int    extended_digits  = 72;
    const int    extended_prec    = std::ceil(extended_digits * log2(10)); // Use a bit higher precision for computations.
    const int    verify_prec      = std::floor(extended_prec-3*log2(10));  // Ignore three last decimal digits (to avoid rounding errors)

    ArgumentsParser args(argc,argv);

    if(args.argSupplied("analyze"))
    {
        //
        // Calculate and show statistics for the polynomials supplied:
        //
        //    psmm -analyze=AllKnown
        //

        show_statistics_of_polynomials(args.getArgValue("analyze"),extended_prec,extended_digits);
    }
    else if(args.argSupplied("merge") && args.argSupplied("output"))
    {
        //
        // Merge several input files into final results.
        //
        //    psmm -merge=file0,file1,file2,...,fileN -output=AllKnown
        //
        // Read all input files, sort the results by degree and M(p) and remove duplicate results.
        // This mode is very imnportant as it allows to search polynomials of different degrees in parallel (and remove duplicate results later on).
        //
        merge_files_with_results(args.getArgValue("merge"),args.getArgValue("output"),verify_prec,extended_digits,extended_prec);
    }
    else
    {
        //
        // Run search:
        //    psmm -degree=N -coeffs=-1,1 -threshold=1.3 -nnz=1,2,3,4 -threads=3 -period=5 -known=AllKnown -addto=AllKnown
        //
        // Extract & process intermediate results from PSMM log file (marked with with ***):
        //    psmm -degree=N -coeffs=-1,1 -threshold=1.3 -nnz=1,2,3,4 -threads=3 -period=5 -fromlog=358_psmm.txt -known=AllKnown -addto=AllKnown
        //
        // This is needed in case if PSMM crashed (e.g. due to lack of memory) but there are a lot of polynomials were found and stored in log file.
        // Exactly this situation happened for N=358. Search was running for 5 days, 556 polynomials were found, but PSMM crashed on the verification stage for unknown reasons.
        // I suspect the memory leaks in the MPSolve were the main reason for the crash.
        //
        // UPDATE (January 11, 2021).
        // Same crash was observed for N=406. MSVC debugger narrowed down the cause to ZXXFactoring.cpp from NTL.
        // Probably we need to write special test to find out why it is crashing (or just wait for NTL update).
        //

        if(
            !args.argSupplied("degree") ||
            !args.argSupplied("coeffs") ||
            !args.argSupplied("nnz")    ||
            !args.argSupplied("threshold")
          )
        {
            // Example:
            //
            // Run search:
            //    psmm -degree=N -coeffs=-1,1 -threshold=1.3 -nnz=1,2,3,4 -threads=3 -period=5 -known=AllKnown -addto=AllKnown
            //
            // Show statistics:
            //    psmm -analyze=AllKnown
            //
            // Merge several input files into final results.
            //
            //    psmm -merge=file0,file1,file2,...,fileN -output=AllKnown
            //

            printf("Usage:\n\tpsmm -degree=N -coeffs=c0,c1,c2,... -threshold=T -nnz=nnz1,nnz2,nnz3,... -threads=M -period=3600 -known=AllKnown -addto=AllKnown"
                         "\n\tpsmm -merge=file0,file1,file2,...,fileN -output=AllKnown"
                         "\n\tpsmm -analyze=AllKnown\n"
                         );

            return EXIT_FAILURE;
        }

        // Parse command-line arguments.
        int degree            = std::stoi(args.getArgValue("degree"));     // Polynomial degree to consider (must be even)
        double threshold      = std::stod(args.getArgValue("threshold"));  // Mahler measure upper threshold.
        std::string fknown    = args.getArgValue("known");                 // File with already known polynomials with Mahler measure < threshold.
        std::string faddto    = args.getArgValue("addto");                 // File to add the found polynomials, usually the same as "load"
        std::string fromlog   = args.getArgValue("fromlog");               // Read & process intermediate results PSMM log file. No search is performed, only full verification.
        int nthreads          = args.argSupplied("threads")    ? std::stoi(args.getArgValue("threads"))    : 1;  // CPU threads to use (1 by default).
        int period            = args.argSupplied("period")     ? std::stoi(args.getArgValue("period"))     : 5;  // Show progress report every "period" of seconds

        if((degree & 1) || degree <= 0)
        {
            puts("Polynomial degree must be even and positive.");
            return EXIT_FAILURE;
        }

        if(nthreads < 0)
        {
            puts("Number of threads must be positive.");
            return EXIT_FAILURE;
        }

        std::vector<int> coeffs;
        std::vector<int> nnz;

        // Build set of possible coefficient values (integers).
        std::vector<std::string> coeffs_values;
        f_split_string(args.getArgValue("coeffs"),',',coeffs_values);
        for(std::size_t i = 0; i < coeffs_values.size(); i++)
            if(!coeffs_values[i].empty())
                coeffs.push_back(std::stoi(coeffs_values[i]));

        // Build list of nonzeros to consider.
        std::vector<std::string> nnz_values;
        f_split_string(args.getArgValue("nnz"),',',nnz_values);
        for(std::size_t i = 0; i < nnz_values.size(); i++)
            if(!nnz_values[i].empty())
            {
                int v = std::stoi(nnz_values[i]);

                if(v > degree/2)
                {
                    printf("Incorrect argument \"nnz\". The nonzeros must be: nnz <= degree/2\n");
                    return EXIT_FAILURE;
                }
                nnz.push_back(v);
            }

        // Load known polynomials
        std::vector<reciprocal_polynomial_t> known;
        if(!fknown.empty())
            load_polynomials(fknown,known,extended_prec);

        // Initialize long integers
        mpz_t total_number_of_polynomials,polynomials_left,time_left;
        mpz_t years,days,hours,minutes,seconds, time_elapsed;
        mpz_inits(total_number_of_polynomials,polynomials_left,time_left,NULL);
        mpz_inits(years,days,hours,minutes,seconds,time_elapsed,NULL);

        // Compute total number of polynomials to check
        compute_total_number_of_polynomials(total_number_of_polynomials,degree,coeffs.size(),nnz);

        // Show preamble
        printf("-----------------------------------------------------------------\n");
        printf("PSMM - Sequential search of polynomials with small Mahler measure.\n");
        printf("Degree       = %d\n",degree);
        printf("Coefficients = %s\n",args.getArgValue("coeffs").c_str());
        printf("Nonzeros     = %s\n",args.getArgValue("nnz").c_str());
        printf("Threshold    = %s\n",args.getArgValue("threshold").c_str());
        printf("Known        = %zu\n",known.size());
        printf("Polynomials  = %s\n",mpz2string(total_number_of_polynomials).c_str());
        printf("-----------------------------------------------------------------\n");
        fflush(stdout);

        // Polynomials found in the current session.
        std::vector<reciprocal_polynomial_t> candidates;
        if(!fromlog.empty())
            load_candidates_from_log(fromlog,degree,candidates);

        std::size_t polys_per_report      = 0;
        std::size_t current_polynomial    = 0;
        std::size_t polynomials_processed = 0;

        // Main loop, iterates over all possible polynomials.
        auto main_loop_start  = std::chrono::high_resolution_clock::now();
        auto last_report_time = std::chrono::high_resolution_clock::now();

        if(fromlog.empty())
        {
            //
            // Run actual search if no input/log file was provided.
            //
            // Architecture:
            //   1. The main thread drives the iterator and fills a batch of polynomials.
            //   2. Mahler-measure computation is parallelized across `nthreads` workers
            //      using a simple atomic work-stealing loop.  Each worker calls MPSolve
            //      single-threaded (its own context), so there is no lock contention.
            //   3. Results are processed sequentially (dedup, logging, candidates list).
            //
            // The batch size balances between amortizing thread-launch overhead and
            // keeping memory footprint bounded.  256 polynomials per batch is a sweet
            // spot for typical low-nnz sparse searches.
            //

            const std::size_t BATCH_SIZE = 256;

            struct batch_item {
                std::vector<int> coeffs;
                int    nnz_value;
                bool   skip;
                double mahler;   // filled by workers
            };

            std::vector<batch_item> batch;
            batch.reserve(BATCH_SIZE);

            for(std::size_t i = 0; i < nnz.size(); i++)
            {
                std::vector<int> poly(degree/2+1);
                reciprocal_polynomials_iterator p(degree,nnz[i],coeffs);

                bool iterator_active = true;
                while(iterator_active)
                {
                    // ---- Phase 1: fill a batch from the iterator (main thread only) ----
                    batch.clear();
                    while(batch.size() < BATCH_SIZE && iterator_active)
                    {
                        bool skip = p.skip_next_polynomial();
                        iterator_active = p.next_polynomial(poly);
                        if(iterator_active)
                        {
                            batch.push_back({poly, nnz[i], skip, 0.0});
                        }
                    }

                    if(batch.empty()) break;

                    // ---- Phase 2: compute Mahler measure in parallel ----
                    //
                    // Each call creates its own mps_context (one-shot). Context reuse
                    // was investigated (C1 in the review) but MPSolve's internal state
                    // (error_state stickiness, zero_roots mismatch on degree change,
                    // and observed heap corruption on reuse even at fixed degree)
                    // makes it unsafe without deeper MPSolve patches.
                    //
                    if(nthreads > 1)
                    {
                        std::atomic<std::size_t> next_idx{0};
                        std::vector<std::thread> workers;
                        workers.reserve(nthreads);
                        for(int t = 0; t < nthreads; ++t)
                        {
                            workers.emplace_back([&](){
                                while(true)
                                {
                                    const std::size_t idx = next_idx.fetch_add(1, std::memory_order_relaxed);
                                    if(idx >= batch.size()) break;
                                    auto& item = batch[idx];
                                    if(!item.skip)
                                    {
#ifdef USE_FAST_MAHLER_ESTIMATOR
                                        item.mahler = estimate_mahler_reciprocal_polynomial_d(item.coeffs, threshold, /*mpsolve_threads=*/2);
#else
                                        item.mahler = compute_mahler_reciprocal_polynomial_d(item.coeffs, /*mpsolve_threads=*/2);
#endif
                                    }
                                }
                            });
                        }
                        for(auto& w : workers) w.join();
                    }
                    else
                    {
                        // Single-threaded path — no thread overhead.
                        for(auto& item : batch)
                        {
                            if(!item.skip)
                            {
#ifdef USE_FAST_MAHLER_ESTIMATOR
                                item.mahler = estimate_mahler_reciprocal_polynomial_d(item.coeffs, threshold, /*mpsolve_threads=*/2);
#else
                                item.mahler = compute_mahler_reciprocal_polynomial_d(item.coeffs, /*mpsolve_threads=*/2);
#endif
                            }
                        }
                    }

                    // ---- Phase 3: process results sequentially (main thread) ----
                    for(auto& item : batch)
                    {
                        if(!item.skip)
                        {
                            if(item.mahler != item.mahler) // NaN check (avoids isnan macro collision)
                            {
                                // The Mahler-estimator returned NaN because mpsolve hit
                                // the non-finite-rdpe path (mpsolve patch 0003 set the
                                // failure flag). Log the polynomial with !!! prefix so
                                // operators can pull it from the log and reproduce /
                                // report upstream. DO NOT add to candidates -- using
                                // a corrupted M(p) would silently mis-classify.
                                printf("!!! NaN-from-mpsolve\t\t");
                                printf("NNZ = %d\t[",item.nnz_value);
                                for(size_t j = 0; j < item.coeffs.size()-1; j++) printf("%2d ",item.coeffs[j]);
                                printf("%2d]\n",item.coeffs.back());
                                fflush(stdout);
                            }
                            else if(item.mahler > 1.0 && item.mahler <= threshold)
                            {
                                int known_before = same_polynomial_found(degree, item.mahler, search_tolerance, known);

                                if(known_before == 0)
                                {
                                    int candidate_found_before = same_polynomial_found(degree, item.mahler, search_tolerance, candidates);

                                    if(candidate_found_before == 0)
                                    {
                                        printf("*** %.16f\t\t",item.mahler);
                                        reciprocal_polynomial_t rp;

                                        rp.N      = degree;
                                        rp.M      = item.mahler;
                                        rp.nnz    = item.nnz_value;
                                        rp.coeffs = item.coeffs;

                                        candidates.push_back(rp);
                                    }
                                    else
                                    {
                                        printf("+++ %.16f\t\t",item.mahler);
                                    }
                                }
                                else
                                {
                                    printf("--- %.16f  (%3d)\t",item.mahler,known_before);
                                }

                                printf("NNZ = %d\t[",item.nnz_value);
                                for(size_t j = 0; j < item.coeffs.size()-1; j++) printf("%2d ",item.coeffs[j]);
                                printf("%2d]\n",item.coeffs.back());
                                fflush(stdout);
                            }

                            polynomials_processed++;
                        }

                        current_polynomial++;
                        polys_per_report++;
                    }

                    // ---- Progress reporting ----
                    auto current_time = std::chrono::high_resolution_clock::now();
                    std::chrono::seconds elapsed_since_last_report = std::chrono::duration_cast<std::chrono::seconds>(current_time-last_report_time);

                    if(elapsed_since_last_report.count() > period)
                    {
                        double pps = double(polys_per_report)/double(elapsed_since_last_report.count());

                        mpz_sub_ui(polynomials_left,total_number_of_polynomials,current_polynomial);
                        const std::size_t pps_rounded = std::max<std::size_t>(1, static_cast<std::size_t>(std::round(pps)));
                        mpz_div_ui(time_left,polynomials_left,pps_rounded);

                        std::string time_left_str = sec2yhms(time_left,years,days,hours,minutes,seconds);

                        printf("\tPPS = %.2f,\tNNZ = %3d,\tDONE = %zu,\tFOUND = %zu\tTIME LEFT = %s\n",pps,nnz[i],current_polynomial,candidates.size(),time_left_str.c_str());

                        last_report_time = current_time;
                        polys_per_report = 0;
                    }

                    fflush(stdout);
                }
            }
        }

        printf("-----------------------------------------------------------------\n");
        printf("Degree                = %d\n",degree);
        printf("Coefficients          = %s\n",args.getArgValue("coeffs").c_str());
        printf("Nonzeros              = %s\n",args.getArgValue("nnz").c_str());
        printf("Threshold             = %s\n",args.getArgValue("threshold").c_str());
        printf("Polynomials Total     = %s\n",mpz2string(total_number_of_polynomials).c_str());
        printf("Polynomials Checked   = %zu\n",polynomials_processed);
        printf("-----------------------------------------------------------------\n");
        printf("Polynomials Selected  = %zu\t(M(p) <= threshold)\n",candidates.size());
        fflush(stdout);


        // Do verification of the found results.

        std::size_t reducible_polynomials            = 0;
        std::size_t factors_checked                  = 0;
        std::size_t success_irreducible_polynomials  = 0;
        std::size_t success_factors_total            = 0;

        std::vector<reciprocal_polynomial_t> verified;

        for(std::size_t i = 0; i < candidates.size(); i++)
        {
            reciprocal_polynomial_t& candidate = candidates[i];

            //
            // Check if candidate polynomial is irreducible.
            //   If "yes", then compute its detailed properties (H(p), L(p), high-accuracy Mahler measure M(p) and root characteristics).
            //   If "no", then factor the polynomial and check each factor to find the one with small Mahler measure M(p).
            //

            //
            // Factors are stored as polynomials with full set of coefficients (not half as we do for reciprocal). This is how NTL works and we follow.
            //
            // Note. Usually we are interested in polynomials with M(p) < 1.3.
            //       In this case all factors are guaranteed to be reciprocal (since lower bound for non-reciprocal polynomials M(p) >= 1.324717..., see [Br51] and [Sm71]).
            //       Be cautious when "threshold" > 1.324717
            //
            std::vector<std::vector<int>> factors;
            bool irreducible = factor_reciprocal_polynomial(candidate.coeffs,factors);

            if(irreducible)
            {
                //
                // Compute detailed properties of the candidate and add it to the verified list.
                //
                compute_all_properties_of_reciprocal_polynomial(candidate,extended_prec,nthreads);

                //
                // Do the last verification step in extended precision.
                //
                if(candidate.F > 1)
                {
                    // Reject substitution polynomials P(x) = Q(x^d), d > 1:
                    // M is inherited from Q at smaller degree (theorem
                    // M(P(x^k)) = M(P)). The DB tracks minimal-degree
                    // representatives only -- same rule as the merge stage
                    // (polyhelpers.h, commit aa63b50). Belt-and-braces here:
                    // emitting a substitution to -addto would be filtered
                    // by the later merge, but should not be written in the
                    // first place.
                    if(is_substitution_polynomial(candidate.N, candidate.coeffs))
                        continue;

                    // Strict (N, coeffs) + x->-x dedup. Was previously
                    // same_polynomial_found_m (M-within-precision), which
                    // silently dropped distinct polynomials whose Mahler
                    // measures coincided with a lower-degree DB entry --
                    // same bug shape as the 2026-05-25 merge incident.
                    bool found = (same_polynomial_by_coeffs(candidate.N, candidate.coeffs, verified) == 1) ||
                                 (same_polynomial_by_coeffs(candidate.N, candidate.coeffs, known   ) == 1);

                    if(!found)
                    {
                        verified.push_back(candidate);
                        success_irreducible_polynomials++;
                        success_factors_total++;
                    }
                }
            }
            else
            {
                //
                // Check each factor of the polynomial otherwise.
                //
                for(std::size_t j = 0; j < factors.size(); j++)
                {
                    std::vector<int>& factor = factors[j];
                    int degree = factors[j].size()-1;

                    //
                    // Find roots, compute Mahler measure and handle the result.
                    // Factors are stored with full set of coefficients, therefore we compute Mahler for full set of coefficients.
                    //
                    double mahler = compute_mahler_general_polynomial_d(factor,nthreads);

                    if(mahler > 1.0 && mahler <= threshold)
                    {
                        // Search if any of known polynomials have the same Mahler measure.
                        int known_before = same_polynomial_found(degree, mahler, search_tolerance, known);

                        if(known_before == 0)
                        {
                            // Search if any of computed polynomials have the same Mahler measure.
                            int verified_found_before = same_polynomial_found(degree, mahler, search_tolerance, verified);

                            if(verified_found_before == 0)
                            {
                                reciprocal_polynomial_t p;

                                p.N  = degree;
                                p.coeffs.assign(&factor[0],&factor[degree/2+1]);

                                compute_all_properties_of_reciprocal_polynomial(p,extended_prec,nthreads);

                                //
                                // Do the last verification step in extended precision.
                                //
                                if(p.F > 1)
                                {
                                    // Reject substitution polynomials P(x) = Q(x^d),
                                    // d > 1 (see comment at the earlier irreducible-
                                    // path site above).
                                    if(is_substitution_polynomial(p.N, p.coeffs))
                                        continue;

                                    // Strict (N, coeffs) + x->-x dedup (see comment at the
                                    // earlier same_polynomial_by_coeffs site above).
                                    bool found = (same_polynomial_by_coeffs(p.N, p.coeffs, verified) == 1) ||
                                                 (same_polynomial_by_coeffs(p.N, p.coeffs, known   ) == 1);

                                    if(!found)
                                    {
                                        verified.push_back(p);
                                        success_factors_total++;
                                    }
                                }
                            }
                        }
                    }

                    factors_checked++;
                }

                reducible_polynomials++;
            }
        }

        // Sort list of cleaned-up results by degree, tie-break by Mahler measure.
        std::sort(verified.begin(),verified.end(),[](const reciprocal_polynomial_t& a,const reciprocal_polynomial_t& b) { return (a.N < b.N) || ((a.N == b.N) && (a.F < b.F)); });

        auto main_loop_stop = std::chrono::high_resolution_clock::now();
        mpz_set_ui(time_elapsed,std::chrono::duration_cast<std::chrono::seconds>(main_loop_stop-main_loop_start).count());
        std::string total_time_elapsed = sec2yhms(time_elapsed,years,days,hours,minutes,seconds);

        printf("--- Reducible         = %zu\n",reducible_polynomials);
        printf("--- Factors Checked   = %zu\n",factors_checked);
        printf("-----------------------------------------------------------------\n");
        printf("Polynomials Found     = %zu\t(M(p) <= threshold)\n",verified.size());
        printf("--- Target degree     = %zu\n",success_irreducible_polynomials);
        printf("--- Lower Degree      = %zu\n",success_factors_total-success_irreducible_polynomials);
        printf("-----------------------------------------------------------------\n");
        printf("Elapsed Time          = %s\n",total_time_elapsed.c_str());
        printf("-----------------------------------------------------------------\n");
        fflush(stdout);

        // Show list of verified polynomials on screen
        for(std::size_t i = 0; i < verified.size(); i++)  printp(verified[i],extended_digits);
        printf("-----------------------------------------------------------------\n");

        // Append list of verified polynomials to the file
        if(verified.size() > 0)
        {
            FILE* foutput = NULL;

            if(!faddto.empty())
                foutput = fopen(faddto.c_str(),"at");

            if(foutput != NULL)
            {
                fprintf(foutput,"\n");
                fprintf(foutput,"# coeffs = [%s] nnz = [%s] time = %s found = %zu polynomials (target degree = %zu, lower degree = %zu)\n",
                            args.getArgValue("coeffs").c_str(),
                            args.getArgValue("nnz").c_str(),
                            total_time_elapsed.c_str(),
                            verified.size(),
                            success_irreducible_polynomials,
                            success_factors_total-success_irreducible_polynomials);

                for(std::size_t i = 0; i < verified.size(); i++)
                {
                    reciprocal_polynomial_t& poly = verified[i];

                    //
                    // D M NNZ H L K U Q R Coefficients
                    //
                    fprintf(foutput, "%3zu %s %zu %zu %zu %zu %zu %zu %zu ",poly.N,mpf2string(poly.F.get_mpf_t(),extended_digits).c_str(), poly.nnz, poly.H, poly.L, poly.K, poly.U, poly.Q, poly.R);
                    for(std::size_t j = 0; j < poly.coeffs.size(); j++) fprintf(foutput, "%d ",poly.coeffs[j]);
                    fprintf(foutput,"\n");
                    fflush(foutput);
                }

                fclose(foutput);
            }
        }

        mpz_clears(total_number_of_polynomials,polynomials_left,time_left,NULL);
        mpz_clears(years,days,hours,minutes,seconds,time_elapsed,NULL);
    }
}

