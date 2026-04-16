/*
 * This file is part of "Polynomials with Small Mahler Measure" (PSMM) project.
 *
 * Copyright (C) 2020, Advanpix LLC.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Pavel Holoborodko <pavel@advanpix.com>
 */

#ifndef __PSMM_POLYNOMIAL_HELPER_FUNCTIONS_H__
#define __PSMM_POLYNOMIAL_HELPER_FUNCTIONS_H__

inline std::pair<int,double> find_nearest_polynomial(double m, std::vector<reciprocal_polynomial_t>& polynomials)
{
    double min_diff = std::numeric_limits<double>::max();
    int min_diff_degree = 0;
    for(std::size_t i = 0; i < polynomials.size(); i++)
    {
        reciprocal_polynomial_t& p = polynomials[i];

        //
        // Precomputed are given with 1e-12 accuracy
        // we skip last few digits so that rounding is not counted.
        //
        double diff = abs(p.M-m);
        if(diff < min_diff)
        {
            min_diff = diff;
            min_diff_degree = p.N;
        }
    }

    return std::pair<int,double>(min_diff_degree,min_diff);
}

inline int same_polynomial_found(int n, double mahler, double tolerance, std::vector<reciprocal_polynomial_t>& polynomials)
{
    std::pair<int,double> nearest_polynomial = find_nearest_polynomial(mahler,polynomials);
    int found = (nearest_polynomial.second <= tolerance) && (nearest_polynomial.first <= n) ? nearest_polynomial.first : 0;
    return found;
}

inline int same_polynomial_found_m(int n, mpf_srcptr mahler, int comp_precision, std::vector<reciprocal_polynomial_t>& polynomials)
{
    //
    // We assume that p.F is computed
    //

    mpf_t m, d, eps;

    mpf_init2(m  ,comp_precision);
    mpf_init2(d  ,comp_precision);
    mpf_init2(eps,comp_precision);

    mpf_set_si(m,1);

    mpf_set_ui(eps,1);
    mpf_div_2exp(eps,eps,std::max(1,(comp_precision-1)));

    for(std::size_t i = 0; i < polynomials.size(); i++)
    {
        reciprocal_polynomial_t& p = polynomials[i];

        if(p.N <= n)
        {
            mpf_sub(d,mahler,p.F);
            mpf_abs(d,d);

            if(mpf_cmp(d,m) <= 0)
                mpf_set(m,d);
        }
    }

    int found = (mpf_cmp(m,eps) <= 0);

    mpf_clears(m,d,eps,NULL);

    return found;
}

inline void printp(const reciprocal_polynomial_t& poly, int digits = 72)
{
    //
    // D M NNZ H L K U Q R Coefficients
    //

    printf("%3zu %s %zu %zu %zu %zu %zu %zu %zu ",poly.N,mpf2string(poly.F,digits).c_str(), poly.nnz, poly.H, poly.L, poly.K, poly.U, poly.Q, poly.R);
    for(std::size_t j = 0; j < poly.coeffs.size(); j++) printf("%d ",poly.coeffs[j]);

    printf("\n");
    fflush(stdout);
}

inline void load_candidates_from_log(const std::string& filename, int N, std::vector<reciprocal_polynomial_t>& polynomials)
{
    std::ifstream ifs(filename);

    if(ifs.is_open())
    {
        std::string original_line, line;
        while(std::getline(ifs,original_line))
        {
            line = f_trim(original_line);
            f_remove_duplicates(line,' ');
            f_remove_duplicates(line,'\t');

            if(line[0] == '*') // consider only lines starting with ***
            {
                // Format:
                //
                // *** %.16f\tNNZ = %d\t[%2d %2d .... %2d %2d]

                std::vector<std::string> tokens;
                f_split_string(line,'\t',tokens);


                reciprocal_polynomial_t poly;

                poly.N   = N;
                poly.M   = atof(tokens[0].c_str()+4);
                poly.nnz = atoi(tokens[1].c_str()+6);

                std::vector<std::string> s_coeffs;
                std::string t_coeffs = tokens[2].substr(1,tokens[2].size()-2);
                f_split_string(f_trim(t_coeffs),' ',s_coeffs);

                std::vector<int>& coeffs = poly.coeffs;
                coeffs.resize(poly.N/2+1);

                std::size_t k = 0;
                for(std::size_t i = 0; i < s_coeffs.size(); i++)
                    coeffs[k++] = atoi(s_coeffs[i].c_str());  // k = [0..N/2]

                if(k != (poly.N/2+1))
                {
                    printf("Parsing error, N = %zu, M = %.16f k = %zu, poly.N/2+1 = %zu %s\n",poly.N,poly.M,k,poly.N/2+1,original_line.c_str());

                    for(std::size_t i = 0; i < tokens.size(); i++)
                        printf("\ttoken[%zu] = %s\n",i,tokens[i].c_str());

                    for(std::size_t i = 0; i < s_coeffs.size(); i++)
                        printf("\ts_coeffs[%zu] = %s\n",i,s_coeffs[i].c_str());

                    exit(1);
                }

                polynomials.push_back(poly);
            }
        }

        ifs.close();
    }

}

inline void load_polynomials(const std::string& filename, std::vector<reciprocal_polynomial_t>& polynomials, int bits = 256)
{
    std::ifstream ifs(filename);

    if(ifs.is_open())
    {
        std::string original_line, line;
        while(std::getline(ifs,original_line))
        {
            line = f_trim(original_line);

            std::vector<std::string> tokens;

            f_remove_duplicates(line,' ');
            f_split_string(line,' ',tokens);

            if(tokens.size() > 2 && line[0] != '#') // skip invalid lines & comments
            {
                // ** Input file format:
                //
                // 0 1  2  3 4 5 6 7 8 9...9+N/2+1
                // N M NNZ H L K U Q R Coefficients
                //
                // NNZ - Number of non-zero coefficients among a[1]..a[N/2]. The a[0] = 1 is fixed hence not counted.
                // H   - Polynomial height (maximum of the magnitudes of its coefficients).
                // L   - Polynomial length (sum of the magnitudes of the coefficients).
                // K   - Number of roots outside the unit circle.
                // U   - Number of complex unity roots (go in pairs):           z = exp(i*t), z* = exp(-i*t).
                // Q   - Number of complex non-unity roots (go in quadruplets): z = r*exp(i*t), z* = r*exp(-i*t), 1/z = (1/r)*exp(-i*t), 1/z* = (1/r)*exp(i*t).
                // R   - Number of real non-unity roots (go in pairs):          z = r, z = 1/r.
                //

                reciprocal_polynomial_t poly;
                std::vector<int>& coeffs = poly.coeffs;

                poly.N  = atoi(tokens[0].c_str());

                // Read the Mahler measure in extended precision and convert to double for fast computations later on.
                mpf_init2(poly.F,bits);
                str2mpf(poly.F,tokens[1].c_str());

                poly.M   = mpf_get_d(poly.F);
                poly.nnz = atof(tokens[2].c_str());
                poly.H   = atof(tokens[3].c_str());
                poly.L   = atof(tokens[4].c_str());
                poly.K   = atof(tokens[5].c_str());
                poly.U   = atof(tokens[6].c_str());
                poly.Q   = atof(tokens[7].c_str());
                poly.R   = atof(tokens[8].c_str());

                coeffs.resize(poly.N/2+1);

                std::size_t k = 0;
                for(std::size_t i = 9; i < tokens.size(); i++)
                    coeffs[k++] = atoi(tokens[i].c_str());  // k = [0..N/2]

                if(k != (poly.N/2+1))
                {
                    printf("Parsing error, N = %zu, M = %.16f %s\n",poly.N,poly.M,original_line.c_str());
                    for(std::size_t i = 0; i < tokens.size(); i++)
                        printf("\ttoken[%zu] = %s\n",i,tokens[i].c_str());

                    exit(1);
                }

                polynomials.push_back(poly);
            }
        }

        ifs.close();
    }
}

inline void merge_files_with_results(const std::string& input, const std::string& output, int verify_precision, int output_digits, int precision = 256)
{
    // Parse the filenames
    std::vector<std::string> filenames;
    f_split_string(input,',',filenames);

    // Read all files and merge all polynomials in one big list.
    std::vector<reciprocal_polynomial_t> polynomials;
    for(std::size_t i = 0; i < filenames.size(); i++)
        load_polynomials(filenames[i],polynomials,precision);

    // Sort list of merged polynomials by degree, tie-break by Mahler measure.
    std::sort(polynomials.begin(),polynomials.end(),[](reciprocal_polynomial_t& a,reciprocal_polynomial_t &b) { return (a.N < b.N) || ((a.N == b.N) && (mpf_cmp(a.F,b.F) < 0)); });

    std::map<std::size_t,std::size_t> nresults;

    // Add unique polynomials from lowest to highest degree to final list of results.
    // Also compute number of unique polynomials for each degree.
    std::vector<reciprocal_polynomial_t> verified;
    for(std::size_t i = 0; i < polynomials.size(); i++)
    {
        reciprocal_polynomial_t& p = polynomials[i];
        if(!same_polynomial_found_m(p.N, p.F, verify_precision, verified))
        {
            verified.push_back(p);
            nresults[p.N]++;
        }
    }

    if(verified.size() > 0)
    {
        FILE* foutput = NULL;

        if(!output.empty())
            foutput = fopen(output.c_str(),"w"); // Re-writes existing file

        if(foutput != NULL)
        {
            int previous = 0;
            for(std::size_t i = 0; i < verified.size(); i++)
            {
                reciprocal_polynomial_t& poly = verified[i];

                if(previous != poly.N)
                {
                   previous = poly.N;

                   fprintf(foutput,"\n");
                   fprintf(foutput,"# degree = %zu, %zu polynomials\n",poly.N,nresults[poly.N]);
                }

                //
                // D M NNZ H L K U Q R Coefficients
                //
                fprintf(foutput, "%3zu %s %zu %zu %zu %zu %zu %zu %zu ",poly.N,mpf2string(poly.F,output_digits).c_str(), poly.nnz, poly.H, poly.L, poly.K, poly.U, poly.Q, poly.R);
                for(std::size_t j = 0; j < poly.coeffs.size(); j++) fprintf(foutput, "%d ",poly.coeffs[j]);
                fprintf(foutput,"\n");
                fflush(foutput);
            }

            fclose(foutput);
        }
    }
 }

#endif // __PSMM_POLYNOMIAL_HELPER_FUNCTIONS_H__