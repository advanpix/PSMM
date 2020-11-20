/*
 * This file is part of "Polynomials with Small Mahler Measure" (PSMM) project.
 *
 * Copyright (C) 2020, Advanpix LLC.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Pavel Holoborodko <pavel@advanpix.com>
 */

#ifndef __PSMM_RECIPROCAL_POLYNOMIALS_ITERATOR_H__
#define __PSMM_RECIPROCAL_POLYNOMIALS_ITERATOR_H__

//
// Computes total number of possible polynomials of a given degree, number of non-zero coefficients and number of possible values of coefficients (base).
//
inline void compute_total_number_of_polynomials(mpz_t m, int degree, int base, std::vector<int>& nnz)
{
    mpz_t a,b;
    mpz_init(a);
    mpz_init(b);

    mpz_set_ui(m,0);
    for(std::size_t i = 0; i < nnz.size(); i++)
    {
        mpz_bin_uiui(a,degree/2,nnz[i]); // we consider only [1..n/2] coefficients of the reciprocal polynomial.
        mpz_ui_pow_ui(b,base,nnz[i]);

        mpz_mul(a,a,b);
        mpz_add(m,m,a);
    }

    mpz_clear(a);
    mpz_clear(b);
}

//
// Allows iteration over full set of reciprocal polynomials of given degree, number of nonzero coefficients and list of possible values of coefficients.
//
//
// Reciprocal polynomial has n/2+1 unique coefficients: a[0], a[1], a[2],...,a[n/2]
// However since a[0] = 1 is fixed, we consider only a[1]...a[n/2].
//
// We represent polynomial as a sparse number with base = cardinality of set of possible coefficient values.
// Then total number of possible polynomials are:
//
//                Q = binomial(n/2,nonzeros) * base^nonzeros
//
// The class is iterator over the full set of possible polynomials.
//
// We assume coefficients have non-templated "double" type for simplicity and to make it uniform with the rest of the seartch (root finders, etc.).
//
class reciprocal_polynomials_iterator
{
public:
    reciprocal_polynomials_iterator()
    { }

    reciprocal_polynomials_iterator(
                                    int degree,                          // Degree of a polynomial
                                    int nonzeros,                        // Number of non-zero coefficients, excluding the a[0] = a[n] = 1
                                    const std::vector<double>& coeffs    // List of possible values of coefficients
                                   )
    {
        setup(degree,nonzeros,coeffs);
    }

    void setup(
                int degree,                          // Degree of a polynomial, must be even.
                int nonzeros,                        // Number of non-zero coefficients, excluding the a[0] = a[n] = 1
                const std::vector<double>& coeffs    // List of possible values of coefficients
              )
    {
        m_Degree         = degree;
        m_Nonzeros       = nonzeros;
        m_PossibleCoeffs = coeffs;
        m_Base           = m_PossibleCoeffs.size();

        // Check if all coefficients have negative pair {-c,c}.
        m_MirroredCoeffs = true;
        for(std::size_t i = 0; i < m_PossibleCoeffs.size() && m_MirroredCoeffs; i++)
            m_MirroredCoeffs = (std::find(m_PossibleCoeffs.begin(),m_PossibleCoeffs.end(),-m_PossibleCoeffs[i]) != m_PossibleCoeffs.end());

        // Build reverse look-up table to map coefficients to their indices.
        for(std::size_t i = 0; i < m_PossibleCoeffs.size(); i++)
            m_CoeffsIndices[int(m_PossibleCoeffs[i])] = i;

        m_Pattern.resize(m_Degree/2);
        m_Number.resize(m_Nonzeros);
        m_TransformedNumber.resize(m_Nonzeros);

        // By default, setup lowest possible state as initial.
        // Later on we can load the initial state from file, or else.
        set_lowest_state();
    }

    // Returns n/2+1 coefficients of the next polynomial in a search set.
    bool next_polynomial(std::vector<double>& polynomial)
    {
        bool status = m_Status;

        if(status)
        {
            // Prepare and return polynomial based on current pattern & number.
            int j = 0;
            polynomial[0] = 1; // a[0] = 1
            for(int i = 1; i <= m_Degree/2; i++)
            {
                if(m_Pattern[i-1]!=0) polynomial[i] = m_PossibleCoeffs[m_Number[j++]];
                else                  polynomial[i] = 0;
            }

            // Move on to the next number by adding 1 to the current number (avoid using expensive % and / operations).
            int carry = 1;
            for(int i = 0; i < m_Nonzeros && (carry == 1); i++)
            {
                m_Number[i] += carry;

                if(carry = (m_Number[i] >= m_Base))
                    m_Number[i] = 0;
            }

            // Switch to the next pattern in case of overflow.
            // Overflow means that all possible numbers were considered and we need to move to another pattern.
            m_Status = (carry == 0) ? true : next_pattern(); // +1 leads to overflow, move to next pattern.
        }

        return status;
    }

    bool skip_next_polynomial()
    {
        //
        // Skip non-primitive polynomials. Primitiveness is independent from coefficients and can be deduced directly from pattern.
        // So that primitiveness can be checked in "next_pattern".
        // However we do it here because we want to measure total iteration speed over the search space.
        //

        bool skip = !m_PrimitivePattern;

        if(!skip && m_MirroredCoeffs)
        {
            //
            // Skip polynomials which have the same M(p) as previousely checked.
            //
            if(m_EvenPattern) // We keep this for historical reasons. Since m_EvenPattern is automatically equivalent to non-primitive pattern (which is already taken care of above).
            {
                // Build number for polynomial P(ix)
                std::size_t j = 0;
                for(std::size_t i = 0; i < m_Pattern.size(); i++)
                {
                    if(m_Pattern[i] != 0)
                    {
                        int a = m_PossibleCoeffs[m_Number[j]];

                        if((((i+1)>>1)&1) !=0 ) m_TransformedNumber[j] = m_CoeffsIndices[-a]; // c[i] = -a[i], if i = 2*k, k is odd
                        else                    m_TransformedNumber[j] = m_Number[j];         // c[i] =  a[i], if i = 2*k, k is even

                        j++;
                    }
                }
            }
            else
            {
                // At least one odd-degree coefficient is non-zero, then investigate P(-x).
                // Build number for polynomial P(-x)
                std::size_t j = 0;
                for(std::size_t i = 0; i < m_Pattern.size(); i++)
                {
                    if(m_Pattern[i] != 0)
                    {
                        int a = m_PossibleCoeffs[m_Number[j]];

                        if(((i+1)&1)!=0) m_TransformedNumber[j] = m_CoeffsIndices[-a]; // c[i] = -a[i], if i is odd
                        else             m_TransformedNumber[j] = m_Number[j];         // c[i] =  a[i], if i is even

                        j++;
                    }
                }
            }

            // Skip polynomial if m_Number > m_TransformedNumber. This means that we already processed the polynomial.
            // Compute decimal representation of the numbers to compare them.
            std::size_t a = m_Number.back();
            std::size_t b = m_TransformedNumber.back();
            for(int i = m_Nonzeros-2; i >= 0; i--)
            {
                a = a * m_Base + m_Number[i];
                b = b * m_Base + m_TransformedNumber[i];
            }

            skip = (a > b);

#if 0
            printf(":::\n");
            for(int i = 0; i < m_Number.size(); i++) printf("%d ",m_Number[i]);
            printf(" (%I64u): ",a);
            for(int i = 0; i < m_Number.size(); i++) printf("%d ",int(m_PossibleCoeffs[m_Number[i]]));
            printf("\n");

            for(int i = 0; i < m_Number.size(); i++) printf("%d ",m_TransformedNumber[i]);
            printf(" (%I64u): ",b);
            for(int i = 0; i < m_Number.size(); i++) printf("%d ",int(m_PossibleCoeffs[m_TransformedNumber[i]]));
            printf("\n%d\n",skip);
#endif
        }

        return skip;
    }
private:

    bool next_pattern()
    {
        bool status = std::next_permutation(m_Pattern.begin(),m_Pattern.end());

        if(status)
        {
            std::fill(m_Number.begin(),m_Number.end(),0); // reset number to all zeros
            analyze_pattern();
        }

        return status;
    }

    void set_lowest_state()
    {
        // All-blank setup of the pattern & number.
        std::fill(m_Pattern.begin(),m_Pattern.end(),0);
        for(int i = 0; i < m_Nonzeros; i++)
            m_Pattern[m_Pattern.size()-(i+1)] = 1; // put 1's at the end of the pattern, as required by std::next_permutation

        std::fill(m_Number.begin(),m_Number.end(),0);
        analyze_pattern();

        m_Status = true;
    }

    void analyze_pattern()
    {
        std::vector<int> degrees, divisors;

        // Checks if pattern/polynomial includes only even degree coefficients.
        m_EvenPattern = true;
        for(int i = 0; i < m_Pattern.size(); i++)
        {
            if(m_Pattern[i]!=0)
            {
                m_EvenPattern = m_EvenPattern && (((i+1)&1) == 0);
                degrees.push_back(i+1);
            }
        }

        degrees.push_back(m_Degree); // The highest degree is not included in the m_Pattern, we must add it manually.

        // Checks if pattern/polynomial is primitive.
        for(int i = 2; i <= degrees[0]; i++) // non-trivial divisors can be [2, min_degree=degrees[0]]
        {
            bool divisible = true;

            for(std::size_t j = 0; j < degrees.size() && divisible; j++)
            {
                divisible = ((degrees[j] % i) == 0);
            }

            if(divisible)
                divisors.push_back(i);
        }

        m_PrimitivePattern = (divisors.size() == 0);
    }

private:
    int m_Degree;    // Degree of a polynomial
    int m_Nonzeros;  // Number of nonzero coefficients

    std::vector<double> m_PossibleCoeffs;
    std::map<int,int>   m_CoeffsIndices; // maps coefficients into their indices in m_PossibleCoeffs

    // Internal state of the iterator.
    std::vector<int> m_Pattern;  // Sparse pattern, =1 indicates the locations of non-zero polynomial coefficients. Length = n/2.
                                 // The 0-degree coefficient is always = 1, hence we ignore it in m_Pattern.
                                 // So that m_Pattern[i] correspond to (i+1) degree coefficient.

    std::vector<int> m_Number;   // Current combination of coefficients from the list of possible values. Length = number of nonzeros.
    int m_Base;
    bool m_Status;

    // Polynomial coefficients analysis
    bool m_MirroredCoeffs;
    bool m_EvenPattern;
    bool m_PrimitivePattern;
    std::vector<int> m_TransformedNumber;

};

#endif // __PSMM_RECIPROCAL_POLYNOMIALS_ITERATOR_H__