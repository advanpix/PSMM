/*
 * This file is part of "Polynomials with Small Mahler Measure" (PSMM) project.
 *
 * Copyright (C) 2020, Advanpix LLC.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Pavel Holoborodko <pavel@advanpix.com>
 */

#ifndef __PSMM_UTILITIES_H__
#define __PSMM_UTILITIES_H__

// Converts mpz number into std::string
inline std::string mpz2string(mpz_srcptr a)
{
    std::size_t length = mpz_sizeinbase(a,10) + 2;
    char* s = (char*) std::malloc(length);

    mpz_get_str(s,10,a);

    std::string result(s);
    std::free(s);

    return result;
}

inline std::string mpf2string(mpf_srcptr x, int digits = 16)
{
    char f[64];
    char s[4096];

    sprintf(f,"%%.%dFf",digits);
    gmp_sprintf(s,f,x);

    return std::string(s);
}

inline void str2mpf(mpf_ptr x, const char* s)
{
    gmp_sscanf(s,"%Ff",x);
}

// Converts seconds into formatted string of years,months,days,hours,minutes and seconds
inline std::string sec2yhms(mpz_t total_seconds, mpz_t years, mpz_t days, mpz_t hours, mpz_t minutes, mpz_t seconds)
{
    mpz_divmod_ui(years,  total_seconds,total_seconds,std::size_t(std::size_t(24)*3600*365));
    mpz_divmod_ui(days,   total_seconds,total_seconds,std::size_t(24)*3600);
    mpz_divmod_ui(hours,  total_seconds,total_seconds,std::size_t(3600));
    mpz_divmod_ui(minutes,seconds,total_seconds,std::size_t(60));

    std::string s;

    s += mpz2string(years  )+"Y:";
    s += mpz2string(days   )+"D:";
    s += mpz2string(hours  )+"H:";
    s += mpz2string(minutes)+"M:";
    s += mpz2string(seconds)+"S";

    return s;
}

// Split string to tokens by delimiter
inline void f_split_string(const std::string& input,char delimiter,std::vector<std::string>& tokens)
{
    std::string token;
    std::stringstream s(input);

    while(std::getline(s,token,delimiter))
        tokens.push_back(token);
}

// Trim from left
inline std::string& ltrim(std::string& s,const char* t = " \t\n\r\f\v")
{
    s.erase(0,s.find_first_not_of(t));
    return s;
}

// Trim from right
inline std::string& rtrim(std::string& s,const char* t = " \t\n\r\f\v")
{
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

// Replace sequence of the same characters to just one
inline void f_remove_duplicates(std::string& str,const char c)
{
    std::string::iterator new_end = std::unique(str.begin(),str.end(),[&c](char lhs,char rhs){return (lhs == rhs) && (lhs == c);});
    str.erase(new_end,str.end());
}

// Trim from left & right
inline std::string& f_trim(std::string& s,const char* t = " \t\n\r\f\v")
{
    return ltrim(rtrim(s,t),t);
}

#endif // __PSMM_UTILITIES_H__