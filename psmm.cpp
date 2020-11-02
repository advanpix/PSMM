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
    // EXAMPLE. There are two different polynomials of 138 and 180 degrees which
    //     have the same Mahler measure up to 11 digits after the comma:
    //
    //     180 1.25282688286[4]667557461368123141605672483501189167102619478483990743452812 40 1 1 1 0 0 -1 -1 -1 0 0 0 0 0 0 0 1 1 1 0 0 -1 -1 -1 0 0 0 0 0 0 0 1 1 1 1 1 0 -1 -1 -1 -1 -1 0 0 0 0 1 1 1 1 1 0 -1 -1 -1 -1 -1 0 0 0 0 1 0 0 0 1 0 0 0 0 -1 -1 0 0 0 0 1 0 0 0 1 0 0 0 0 -1 -1 0 0 0 0 1
    //     138 1.25282688286[8]347583989529787618964056215870696871794089089602385130202388 15 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1
    //
    //     If we would store these values with the accuracy of 12 digits only (Known180) then these two values would be considered the same for high tolerances like 1e-10 or 1e-11.
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
    const int    search_accuracy  = 15;    // Request MPSolve to compute roots with this number of correct digits during the search phase.
    const double search_tolerance = 1e-14; // Must be ~pow(10,-(search_accuracy-1)). We assume polynomials are different if |M(p1)-M(p2)| > search_tolerance.

    const int    extended_digits  = 72;
    const int    extended_prec    = std::ceil(extended_digits * log2(10)); // Use a bit higher precision for computations.

    ArgumentsParser args(argc,argv);

    if(
        !args.argSupplied("degree") ||
        !args.argSupplied("coeffs") ||
        !args.argSupplied("nnz")    ||
        !args.argSupplied("threshold")
      )
    {
        // Example:
        //
        //    psmm -degree=N -coeffs=-1,1 -threshold=1.3 -nnz=1,2,3,4 -threads=3 -period=5 -known=AllKnown -addto=AllKnown
        //
        puts("psmm -degree=N -coeffs=c0,c1,c2,... -threshold=T -nnz=nnz1,nnz2,nnz3,... -threads=M -period=3600 -known=AllKnown -addto=AllKnown");
        return EXIT_FAILURE;
    }

    // Parse command-line arguments.
    int degree         = std::stoi(args.getArgValue("degree"));     // Polynomial degree to consider (must be even)
    double threshold   = std::stod(args.getArgValue("threshold"));  // Mahler measure upper threshold.
    std::string fknown = args.getArgValue("known");                 // File with already known polynomials with Mahler measure < threshold.
    std::string faddto = args.getArgValue("addto");                 // File to add the found polynomials, usually the same as "load"
    int nthreads       = args.argSupplied("threads") ? std::stoi(args.getArgValue("threads")) : 1;  // CPU threads to use (1 by default).
    int period         = args.argSupplied("period")  ? std::stoi(args.getArgValue("period"))  : 5;  // Show progress report every "period" of seconds

    if(degree & 1)
    {
        puts("Polynomial degree must be even.");
        return EXIT_FAILURE;
    }

    if(nthreads < 0)
    {
        puts("Number of threads must be positive.");
        return EXIT_FAILURE;
    }

    std::vector<double> coeffs;
    std::vector<int>    nnz;

    // Build set of possible coefficient values.
    std::vector<std::string> coeffs_values;
    f_split_string(args.getArgValue("coeffs"),',',coeffs_values);
    for(std::size_t i = 0; i < coeffs_values.size(); i++)
        if(!coeffs_values[i].empty())
            coeffs.push_back(std::stod(coeffs_values[i]));

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
    printf("PSMM - Sequential search of polynomials with small Mahler measure.\n",degree);
    printf("Degree       = %d\n",degree);
    printf("Coefficients = %s\n",args.getArgValue("coeffs").c_str());
    printf("Nonzeros     = %s\n",args.getArgValue("nnz").c_str());
    printf("Threshold    = %s\n",args.getArgValue("threshold").c_str());
    printf("Known        = %d\n",known.size());
    printf("Polynomials  = %s\n",mpz2string(total_number_of_polynomials).c_str());
    printf("-----------------------------------------------------------------\n");
    fflush(stdout);

    // Polynomials found in the current session.
    std::vector<reciprocal_polynomial_t> candidates;

    std::size_t polys_per_report      = 0;
    std::size_t current_polynomial    = 0;
    std::size_t polynomials_processed = 0;

    // Main loop, iterates over all possible polynomials.
    auto main_loop_start  = std::chrono::high_resolution_clock::now();
    auto last_report_time = std::chrono::high_resolution_clock::now();
    for(std::size_t i = 0; i < nnz.size(); i++)
    {
        std::vector<double> poly(degree/2+1);
        reciprocal_polynomials_iterator p(degree,nnz[i],coeffs);

        bool fcontinue = true;
        while(fcontinue)
        {
            bool skip = p.skip_next_polynomial(); // Skip next polynomials in a sequence (e.g. p(-x) and p(x) have the same roots by absolute value)
                                                  // This can be extended, e.g. to detect non-primitive polynomials, or even by appliyng Graeffe's pre-screening.

            fcontinue = p.next_polynomial(poly);  // Get the next polynomial in a sequence.

            if(fcontinue)
            {
                if(!skip)
                {
                    // Find roots, compute Mahler measure and handle the result.
                    double mahler = compute_mahler_reciprocal_polynomial(poly,search_accuracy,nthreads); // Estimate Mahler measure in double precision (fast).

                    if(mahler > 1.0 && mahler <= threshold)
                    {
                        // Search if any of known polynomials have the same Mahler measure.
                        int known_before = same_polynomial_found(degree, mahler, search_tolerance, known);

                        if(known_before == 0)
                        {
                            // Search if any of recently computed polynomials have the same Mahler measure.
                            int candidate_found_before = same_polynomial_found(degree, mahler, search_tolerance, candidates);

                            if(candidate_found_before == 0)
                            {
                                // Unseen polynomial has been found, celebrate it with "***".
                                printf("*** %.16f\t\t",mahler);
                                reciprocal_polynomial_t p;

                                p.N      = degree;
                                p.M      = mahler;
                                p.nnz    = nnz[i];
                                p.coeffs = poly;

                                candidates.push_back(p);
                            }
                            else
                            {
                                // Polynomial was previousely found in the current session, mark it with "+++".
                                printf("+++ %.16f\t\t",mahler);
                            }
                        }
                        else
                        {
                            // Polynomial was known before, mark it with "---" as unimportant.
                            printf("--- %.16f  (%3d)\t",mahler,known_before);
                        }

                        printf("NNZ = %d\t[",nnz[i]);
                        for(size_t j = 0; j < poly.size()-1; j++) printf("%2d ",int(poly[j]));
                        printf("%2d]\n",int(poly.back()));
                        fflush(stdout);
                    }

                    polynomials_processed++;
                }

                // Compute & show some progress statistics
                current_polynomial++;
                polys_per_report++;

                auto current_time = std::chrono::high_resolution_clock::now();
                std::chrono::seconds elapsed_since_last_report = std::chrono::duration_cast<std::chrono::seconds>(current_time-last_report_time);

                if(elapsed_since_last_report.count() > period)
                {
                    double pps = double(polys_per_report)/double(elapsed_since_last_report.count());

                    mpz_sub_ui(polynomials_left,total_number_of_polynomials,current_polynomial);
                    mpz_div_ui(time_left,polynomials_left,std::size_t(std::round(pps)));

                    std::string time_left_str = sec2yhms(time_left,years,days,hours,minutes,seconds);

                    printf("\tPPS = %.2f,\tNNZ = %3d,\tDONE = %d,\tFOUND = %I64u\tTIME LEFT = %s\n",pps,nnz[i],current_polynomial,candidates.size(),time_left_str.c_str());

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
    printf("Polynomials Checked   = %I64u\n",polynomials_processed);
    printf("-----------------------------------------------------------------\n");
    printf("Polynomials Selected  = %I64u\t(M(p) <= threshold)\n",candidates.size());
    fflush(stdout);


    // Do the refinement of the found results. Checking the irreducibility
    // The main criteria is irreducibility F(p).
    // If polynomial of the target (new) degree is irreducible F(p) = true and 1 < M(p) <= T then it is valid and good for us.
    // However F(p) is slow and we approximate it by checking if M(p) close to the already found polynomials.

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
        std::vector<std::vector<double>> factors;
        bool irreducible = factor_reciprocal_polynomial(candidate.coeffs,factors);

        if(irreducible)
        {
            //
            // Compute detailed properties of the candidate and add it to the verified list.
            //
            compute_all_properties_of_reciprocal_polynomial(candidate,extended_prec,nthreads);

            verified.push_back(candidate);
            success_irreducible_polynomials++;
            success_factors_total++;
        }
        else
        {
            //
            // Check each factor of the polynomial otherwise.
            //
            for(std::size_t j = 0; j < factors.size(); j++)
            {
                std::vector<double>& factor = factors[j];
                int degree = factors[j].size()-1;

                //
                // Find roots, compute Mahler measure and handle the result.
                // Factors are stored with full set of coefficients, therefore we compute Mahler for full set of coefficients.
                //
                double mahler = compute_mahler_general_polynomial(factor,search_accuracy,nthreads);

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

                            verified.push_back(p);
                            success_factors_total++;
                        }
                    }
                }

                factors_checked++;
            }

            reducible_polynomials++;
        }
    }

    // Sort list of cleaned-up results by degree & Mahler measure.
    std::sort(verified.begin(),verified.end(),[](reciprocal_polynomial_t& a,reciprocal_polynomial_t &b) { return a.M < b.M;}); // by Mahler
    std::sort(verified.begin(),verified.end(),[](reciprocal_polynomial_t& a,reciprocal_polynomial_t &b) { return (a.N < b.N) || ((a.N == b.N) && (a.M < b.M)); }); // by degree

    auto main_loop_stop = std::chrono::high_resolution_clock::now();
    mpz_set_ui(time_elapsed,std::chrono::duration_cast<std::chrono::seconds>(main_loop_stop-main_loop_start).count());
    std::string total_time_elapsed = sec2yhms(time_elapsed,years,days,hours,minutes,seconds);

    printf("--- Reducible         = %I64u\n",reducible_polynomials);
    printf("--- Factors Checked   = %I64u\n",factors_checked);
    printf("-----------------------------------------------------------------\n");
    printf("Polynomials Found     = %I64u\t(M(p) <= threshold)\n",verified.size());
    printf("--- Target degree     = %I64u\n",success_irreducible_polynomials);
    printf("--- Lower Degree      = %I64u\n",success_factors_total-success_irreducible_polynomials);
    printf("-----------------------------------------------------------------\n");
    printf("Elapsed Time          = %s\n",total_time_elapsed.c_str());
    printf("-----------------------------------------------------------------\n");
    fflush(stdout);

    // Show final results on screen
    for(std::size_t i = 0; i < verified.size(); i++)
    {
        reciprocal_polynomial_t& poly = verified[i];

        //
        // D M NNZ H L K U Q R Coefficients
        //
        printf("%2d %s %d %d %d %d %d %d %d ",poly.N,mpf2string(poly.F,extended_digits).c_str(), poly.nnz, poly.H, poly.L, poly.K, poly.U, poly.Q, poly.R);
        for(std::size_t j = 0; j < poly.coeffs.size(); j++) printf("%d ",int(poly.coeffs[j]));

        printf("\n");
        fflush(stdout);
    }
    printf("-----------------------------------------------------------------\n");

    // Append final results to the file
    if(verified.size() > 0)
    {
        FILE* foutput = NULL;

        if(!faddto.empty())
            foutput = fopen(faddto.c_str(),"at");

        if(foutput != NULL)
        {
            fprintf(foutput,"\n");
            fprintf(foutput,"# coeffs = [%s] nnz = [%s] time = %s found = %I64u polynomials\n",args.getArgValue("coeffs").c_str(),args.getArgValue("nnz").c_str(),total_time_elapsed.c_str(),verified.size());
            for(std::size_t i = 0; i < verified.size(); i++)
            {
                reciprocal_polynomial_t& poly = verified[i];

                //
                // D M NNZ H L K U Q R Coefficients
                //
                fprintf(foutput, "%2d %s %d %d %d %d %d %d %d ",poly.N,mpf2string(poly.F,extended_digits).c_str(), poly.nnz, poly.H, poly.L, poly.K, poly.U, poly.Q, poly.R);
                for(std::size_t j = 0; j < poly.coeffs.size(); j++) fprintf(foutput, "%d ",int(poly.coeffs[j]));
                fprintf(foutput,"\n");
                fflush(foutput);
            }

            fclose(foutput);
        }
    }

    mpz_clears(total_number_of_polynomials,polynomials_left,time_left,NULL);
    mpz_clears(years,days,hours,minutes,seconds,time_elapsed,NULL);
}

