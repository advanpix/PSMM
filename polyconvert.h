/*
 * This file is part of "Polynomials with Small Mahler Measure" (PSMM) project.
 *
 * Copyright (C) 2020, Advanpix LLC.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Pavel Holoborodko <pavel@advanpix.com>
 */

#ifndef __PSMM_POLYNOMIAL_CONVERTERS_H__
#define __PSMM_POLYNOMIAL_CONVERTERS_H__

inline void load_polynomials_old(const std::string& filename, std::vector<reciprocal_polynomial_t>& polynomials, int bits = 256)
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
                // This code is for text file in a format of Known180.gzip (downloaded from http://www.cecm.sfu.ca/~mjm/Lehmer/lists/Known180.gz)
                // Table has spaces between coefficients.
                // Only half of coefficients are stored.
                //
                // 0 - degree
                // 1 - M Mahler measure
                // 2 - K Number of roots outside the unit circle.
                // [3 ... ] - Coefficients
                //
                // Example: 16  1.224278907222   2  1 1 0-1-1 0 1 1 1
                //
                reciprocal_polynomial_t poly;
                std::vector<int>& coeffs = poly.coeffs;

                poly.N  = atoi(tokens[0].c_str());

                // Read the Mahler measure in extended precision and convert to double for fast computations later on.
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
                    printf("Parsing error, N = %zu, M = %.16f %s\n",poly.N,poly.M,original_line.c_str());
                    for(std::size_t i = 0; i < tokens.size(); i++)
                        printf("\ttoken[%zu] = %s\n",i,tokens[i].c_str());

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

inline void convert_file(const std::string& srcfile, const std::string& dstfile, int output_digits, int bits = 256, int nthreads = 1)
{
    // This function converts the "old" format file to extended one.
    // Comments and whitelines are preserved.
    //
    // ** Src file format (e.g. Known180):
    //
    // N M K Coefficients
    //
    // N - polynomial degree.
    // M - Mahler measure.
    // K - Number of roots outside the unit circle.
    //
    // ** Dst file format:
    //
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

    std::ifstream ifs(srcfile);
    std::ofstream ofs(dstfile);

    if(ifs.is_open() && ofs.is_open())
    {
        std::size_t PrevN = 0;

        std::string original_line, line;
        while(std::getline(ifs,original_line))
        {
            line = f_trim(original_line);

            if(line[0] == '#' || line.size() < 10) // skip comments and empty lines.
            {
                ofs << original_line << std::endl;
            }
            else
            {
                std::vector<std::string> tokens;

                f_remove_duplicates(line,' ');
                f_split_string(line,' ',tokens);

                if(tokens.size() > 2) // skip invalid lines (do we have them at all?)
                {
                    //
                    //
                    // This code is for text file in a format of Known180.gzip (downloaded from http://www.cecm.sfu.ca/~mjm/Lehmer/lists/Known180.gz)
                    // Table has spaces between coefficients.
                    // Only half of coefficients are stored.
                    //
                    // 0 - degree
                    // 1 - M Mahler measure
                    // 2 - K Number of roots outside the unit circle.
                    // [3 ... ] - Coefficients
                    //
                    // Example: 16  1.224278907222   2  1 1 0-1-1 0 1 1 1
                    //
                    reciprocal_polynomial_t p;

                    // Parse line from old file
                    p.N  = atoi(tokens[0].c_str());

                    p.coeffs.resize(p.N/2+1);

                    std::size_t k = 0;
                    for(std::size_t i = 3; i < tokens.size(); i++)
                        p.coeffs[k++] = atoi(tokens[i].c_str());  // k = [0..N/2]

                    if(k != (p.N/2+1))
                    {
                        printf("Parsing error, N = %zu, %s\n",p.N,original_line.c_str());
                        for(std::size_t i = 0; i < tokens.size(); i++)
                            printf("\ttoken[%zu] = %s\n",i,tokens[i].c_str());

                        exit(1);
                    }

                    // Compute all the characteristics.
                    compute_all_properties_of_reciprocal_polynomial(p, bits, nthreads);

                    // Write line to the new file
                    ofs << p.N << " " << mpf2string(p.F,output_digits)<< " " << p.nnz<< " " << p.H<< " " << p.L << " " << p.K << " " << p.U<< " " << p.Q<< " " << p.R << " ";
                    for(std::size_t i = 0; i < p.N/2; i++) ofs << p.coeffs[i] << " ";
                    ofs << p.coeffs[p.N/2] << std::endl;

                    // Historical cross-check against the K field from the source file was
                    // disabled long ago; left as a note for anyone auditing the old format.
                    //   const std::size_t K_from_file = static_cast<std::size_t>(atoi(tokens[2].c_str()));
                    //   if(K_from_file != p.K) { ... report and exit ... }

                    if(PrevN != p.N)
                    {
                        printf("N = %zu\n",p.N);
                        PrevN = p.N;
                    }
                }
            }
        }
    }
}

#endif // __PSMM_POLYNOMIAL_CONVERTERS_H__