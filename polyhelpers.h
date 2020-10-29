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

typedef struct{
    std::size_t N;              // Polynomial degree
    double M;                   // Mahler measure in double precision.
    mpf_t  F;                   // Mahler measure in extended precision (if any)
    std::vector<double> coeffs; // Polynomial coefficients.
    std::size_t K;
    std::size_t nnz;            // Number of non-zero coefficients
}computed_polynomial_t;

inline void load_polynomials(const std::string& filename, std::vector<computed_polynomial_t>& polynomials, int bits = 256)
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
                //
                //
                // This code is for text file from Known180.gzip (downloaded from http://www.cecm.sfu.ca/~mjm/Lehmer/lists/Known180.gz)
                // Table has spaces between coefficients.
                // Only half of coefficients are stored.
                //
                // 0 - degree
                // 1 - M(f)
                // 2 - K(?)
                // [3 ... ] - Coefficients
                //
                // Example: 16  1.224278907222   2  1 1 0-1-1 0 1 1 1
                //
                computed_polynomial_t poly;
                std::vector<double>& coeffs = poly.coeffs;

                poly.N  = atoi(tokens[0].c_str());

                // Read the Mahler measure in extended precision and convcert to double for fast computations later on.
                mpf_init2(poly.F,bits);
                str2mpf(poly.F,tokens[1].c_str());
                poly.M  = mpf_get_d(poly.F);

                poly.K  = atof(tokens[2].c_str());
                coeffs.resize(poly.N/2+1);

                std::size_t k = 0;
                for(std::size_t i = 3; i < tokens.size(); i++)
                    coeffs[k++] = atoi(tokens[i].c_str());  // k = [0..N/2]

                if(k != (poly.N/2+1))
                {
                    printf("Parsing error, N = %d, M = %.16f %s\n",poly.N,poly.M,original_line.c_str());
                    for(std::size_t i = 0; i < tokens.size(); i++)
                        printf("\ttoken[%d] = %s\n",i,tokens[i].c_str());

                    exit(1);
                }

                poly.nnz = 0;
                for(std::size_t i = 1; i <= poly.N/2; i++)
                    poly.nnz += (coeffs[i] != 0);

                polynomials.push_back(poly);
            }
        }

        ifs.close();
    }
}

inline std::pair<int,double> find_nearest_polynomial(double m, std::vector<double>& poly, std::vector<computed_polynomial_t>& polynomials)
{
    double min_diff = std::numeric_limits<double>::max();
    int min_diff_degree = 0;
    for(std::size_t i = 0; i < polynomials.size(); i++)
    {
        computed_polynomial_t& p = polynomials[i];

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

#endif // __PSMM_POLYNOMIAL_HELPER_FUNCTIONS_H__