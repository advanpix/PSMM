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

inline void compute_all_properties_of_reciprocal_polynomial(reciprocal_polynomial_t& p, int bits = 256, int nthreads = 1)
{
    //
    // Polynomial degree and coefficients must be set before calling the function.
    //

    int N = 2 * (p.coeffs.size()-1); // = p.N;
    std::vector<double> a(N+1);

    // Expand coefficients to full polynomial
    for(int k = 0; k <= N/2; k++) a[k]     = p.coeffs[k];
    for(int k = 1; k <= N/2; k++) a[N/2+k] = a[N/2-k];

    //
    // Compute roots, their properties and Mahler measure in extended precision.
    //
    mpsolve_compute_mahler_with_properties(p.F, N, a.data(), p.K, p.U, p.Q, p.R, bits, nthreads);
    p.M = mpf_get_d(p.F);

    //
    // Compute NNZ, for the half of the coefficients.
    //
    p.nnz = std::accumulate(p.coeffs.begin(),p.coeffs.end(),0.0,[](double& a,double& b)->double {return a += (b!=0);});
    p.nnz-=1;  // ignore the a[0] = 1, which is constant.

    //
    // Compute NNZ, H(p) and L(p)
    //
    // Well, we can do this using std::max_element and std::accumulate, but here we do this in one pass over the coeffcients.
    p.L = p.H = 0;
    for(std::size_t k = 0; k < a.size(); k++)
    {
        std::size_t c = std::size_t(abs(a[k]));
        p.H  = std::max(p.H,c);
        p.L += c;
    }
}

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
                std::vector<double>& coeffs = poly.coeffs;

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
                std::vector<double>& coeffs = poly.coeffs;

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
                    printf("Parsing error, N = %d, M = %.16f %s\n",poly.N,poly.M,original_line.c_str());
                    for(std::size_t i = 0; i < tokens.size(); i++)
                        printf("\ttoken[%d] = %s\n",i,tokens[i].c_str());

                    exit(1);
                }

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
                        printf("Parsing error, N = %d, %s\n",p.N,original_line.c_str());
                        for(std::size_t i = 0; i < tokens.size(); i++)
                            printf("\ttoken[%d] = %s\n",i,tokens[i].c_str());

                        exit(1);
                    }

                    std::size_t K = atoi(tokens[2].c_str());

                    // Compute all the characteristics.
                    compute_all_properties_of_reciprocal_polynomial(p, bits, nthreads);

                    // Write line to the new file
                    ofs << p.N << " " << mpf2string(p.F,output_digits)<< " " << p.nnz<< " " << p.H<< " " << p.L << " " << p.K << " " << p.U<< " " << p.Q<< " " << p.R << " ";
                    for(std::size_t i = 0; i < p.N/2; i++) ofs << p.coeffs[i] << " ";
                    ofs << p.coeffs[p.N/2] << std::endl;

                    //if(K != p.K)
                    //{
                    //    printf("Parsing error, N = %d, K(%d) != p.K(%d), %s\n", p.N, K, p.K,original_line.c_str());
                    //    exit(1);
                    //}

                    if(PrevN != p.N)
                    {
                        printf("N = %d\n",p.N);
                        PrevN = p.N;
                    }
                }
            }
        }
    }
}

#endif // __PSMM_POLYNOMIAL_HELPER_FUNCTIONS_H__