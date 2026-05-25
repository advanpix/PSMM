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

#include <set>
#include <utility>

// x -> -x map on a reciprocal polynomial's half-coefficient vector
// a[0]..a[N/2]: negate the entries at odd indices. P(-x) and P(x) have
// identical Mahler measure but differ in canonical form; the DB stores
// exactly one representative per {P(x), P(-x)} equivalence class.
inline std::vector<int> xneg_flip_half(const std::vector<int>& half)
{
    std::vector<int> out(half);
    for(std::size_t j = 0; j < out.size(); ++j)
        if(j & 1) out[j] = -out[j];
    return out;
}

// STRICT dedup predicate: returns 1 if any entry in `polynomials` has
// the same degree as `n` AND coefficients equal to `half_coeffs` or to
// its x -> -x flip; 0 otherwise.
//
// Distinguishing from same_polynomial_found{,_m}: that pair is the
// FAST search-phase reducibility-skip heuristic, comparing only by N
// and Mahler measure within tolerance/precision. The heuristic
// silently drops distinct polynomials that share their Mahler measure
// (e.g. an irreducible Salem at degree N and a higher-N polynomial
// whose factorisation includes that Salem), which is fine for skipping
// expensive reducibility checks during enumeration but WRONG as a
// final-verify gatekeeper before adding to the verified list.
//
// Use this helper at every site where we have already paid for the
// reducibility/property computation and just need a strict
// "is this polynomial already in our set?" check -- merge_files_with_results,
// the verify path at psmm.cpp:475-476 / 526-527, and any other place
// that needs (N, coeffs)-level uniqueness.
inline int same_polynomial_by_coeffs(int n,
                                     const std::vector<int>& half_coeffs,
                                     std::vector<reciprocal_polynomial_t>& polynomials)
{
    std::vector<int> flip = xneg_flip_half(half_coeffs);
    for(std::size_t i = 0; i < polynomials.size(); ++i)
    {
        const reciprocal_polynomial_t& p = polynomials[i];
        if(p.N == static_cast<std::size_t>(n))
        {
            if(p.coeffs == half_coeffs) return 1;
            if(p.coeffs == flip)        return 1;
        }
    }
    return 0;
}

inline int same_polynomial_found(int n, double mahler, double tolerance, std::vector<reciprocal_polynomial_t>& polynomials)
{
    //
    // Search for ANY polynomial with degree <= n whose Mahler measure
    // is within tolerance of `mahler`. Returns the degree of the first
    // match, or 0 if none found.
    //
    // The old implementation found the globally nearest polynomial by M
    // and then checked its degree, which missed matches when a closer
    // polynomial at higher degree shadowed a valid match at lower degree.
    //
    // INTENDED USE (and its high-precision sibling same_polynomial_found_m):
    //   Search-phase reducibility skip. A candidate of degree N whose M
    //   coincides with an existing entry at lower degree is almost
    //   certainly reducible (q = p * c with cyclotomic c gives M(q) = M(p)),
    //   so we skip the expensive factorisation check. The heuristic accepts
    //   the theoretical possibility of missing a genuinely distinct
    //   irreducible polynomial with coincident M -- not a concern in
    //   practice for the search space we care about.
    //
    // DO NOT USE FOR MERGE / DEDUP / DB CONSOLIDATION:
    //   Two polynomials with different (N, half_coefficients) are different
    //   mathematical objects and the DB tracks both, even when M coincides.
    //   For those code paths, dedup by (N, coeffs) + x->-x equivalence -- see
    //   merge_files_with_results below and work/tools/safe_merge.py.
    //   Commit e0340eb (2026-05-25) fixed merge_files_with_results that
    //   had used this primitive in error.
    //
    for(std::size_t i = 0; i < polynomials.size(); i++)
    {
        const reciprocal_polynomial_t& p = polynomials[i];

        if(p.N <= static_cast<std::size_t>(n))
        {
            if(std::fabs(p.M - mahler) <= tolerance)
                return p.N;
        }
    }

    return 0;
}

inline int same_polynomial_found_m(int n, mpf_srcptr mahler, int comp_precision, std::vector<reciprocal_polynomial_t>& polynomials)
{
    //
    // We assume that p.F is computed.
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
            mpf_sub(d,mahler,p.F.get_mpf_t());
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

    printf("%3zu %s %zu %zu %zu %zu %zu %zu %zu ",poly.N,mpf2string(poly.F.get_mpf_t(),digits).c_str(), poly.nnz, poly.H, poly.L, poly.K, poly.U, poly.Q, poly.R);
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
                poly.F.set_prec(bits);
                str2mpf(poly.F.get_mpf_t(),tokens[1].c_str());

                poly.M   = poly.F.get_d();
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

    // Sort list of merged polynomials by degree, tie-break by Mahler
    // measure, final tie-break by coefficient vector (lexicographic).
    // The coefficient tie-break is what makes the canonical-form choice
    // for x->-x pairs deterministic and matches work/tools/safe_merge.py.
    std::sort(polynomials.begin(),polynomials.end(),
        [](const reciprocal_polynomial_t& a,const reciprocal_polynomial_t& b) {
            if(a.N != b.N) return a.N < b.N;
            if(a.F != b.F) return a.F < b.F;
            return a.coeffs < b.coeffs;
        });

    std::map<std::size_t,std::size_t> nresults;

    // Dedup by (degree, half-coefficients) -- distinct polynomials sharing
    // only their Mahler measure (e.g. an irreducible Salem and a higher-N
    // polynomial whose factorisation includes the Salem) are NOT duplicates
    // and must both be preserved, regardless of which degree is smaller.
    //
    // x -> -x flips of an existing entry are also collapsed here per Rule 9
    // of the scan workflow: P(x) and P(-x) share an equivalence class, and
    // the DB stores exactly one representative. With the input sorted by
    // (N asc, M asc), whichever variant appears first becomes the canonical
    // form -- callers wanting a specific canonical form (e.g. Lehmer's
    // classical "1 1 0 -1 -1 -1") must place the DB file with that form
    // first in the merge input list so it loads before competing flips.
    //
    // The previous implementation used same_polynomial_found_m which deduped
    // by (p.N <= n, M within verify_precision) and silently dropped distinct
    // higher-degree entries whenever a lower-degree entry with coincident M
    // was merged. See commit bc1e760 for the 2026-05-25 incident.
    std::vector<reciprocal_polynomial_t> verified;
    std::set<std::pair<std::size_t, std::vector<int>>> seen_keys;
    for(std::size_t i = 0; i < polynomials.size(); i++)
    {
        reciprocal_polynomial_t& p = polynomials[i];
        std::pair<std::size_t, std::vector<int>> key(p.N, p.coeffs);
        std::pair<std::size_t, std::vector<int>> xneg_key(p.N, xneg_flip_half(p.coeffs));
        if(seen_keys.count(key)      != 0) continue;  // exact (N, coeffs) duplicate
        if(seen_keys.count(xneg_key) != 0) continue;  // x->-x flip of existing entry
        seen_keys.insert(key);
        verified.push_back(p);
        nresults[p.N]++;
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
                fprintf(foutput, "%3zu %s %zu %zu %zu %zu %zu %zu %zu ",poly.N,mpf2string(poly.F.get_mpf_t(),output_digits).c_str(), poly.nnz, poly.H, poly.L, poly.K, poly.U, poly.Q, poly.R);
                for(std::size_t j = 0; j < poly.coeffs.size(); j++) fprintf(foutput, "%d ",poly.coeffs[j]);
                fprintf(foutput,"\n");
                fflush(foutput);
            }

            fclose(foutput);
        }
    }
 }

#endif // __PSMM_POLYNOMIAL_HELPER_FUNCTIONS_H__