/*
 * This file is part of "Polynomials with Small Mahler Measure" (PSMM) project.
 *
 * Copyright (C) 2020, Advanpix LLC.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 *
 * Authors:
 *   Pavel Holoborodko <pavel@advanpix.com>
 */

#ifndef __COMMAND_LINE_ARGUMENTS_PARSER_H__
#define __COMMAND_LINE_ARGUMENTS_PARSER_H__

//
// Basic command line arguments parser.
//
// Class can parse arguments in the following format:
//
//   -key0=value0 -key1=value1 ... -keyN=valueN
//
// For example if following arguments are supplied:
//
//  -degree=10 -coeffs=-2,-1,1,2 -threshold=1.3 -threads=16
//
// then we can get them as following:
//
//     std::string degree = p.getArgValue("degree")
//     std::string coeffs = p.getArgValue("coeffs")
//     ...
//
// Please note, there is no error detection, so use properly.
//

class ArgumentsParser{

public:
    ArgumentsParser() {}

    ArgumentsParser(int argc, char **argv)
    {
        parse(argc,argv);
    }

    void parse(int argc,char **argv)
    {
        for(std::size_t i = 1; i < argc; i++)
        {
            std::vector<std::string> tokens;
            std::string s(argv[i]);

            f_split_string(f_trim(s),'=',tokens);

            if(tokens.size() > 1 && tokens[0].size() > 1)
                m_Arguments[tokens[0].substr(1)] = tokens[1];
        }
    }

    const std::string& getArgValue(const std::string &option, const std::string &defaultValue = "") const
    {
        auto it = m_Arguments.find(option);
        if(it != m_Arguments.end()) return it->second;
        else return defaultValue;
    }

    bool argSupplied(const std::string &option) const
    {
        return (m_Arguments.find(option) != m_Arguments.end());
    }

    void showAll()
    {
        for(auto it = m_Arguments.begin(); it != m_Arguments.end(); ++it)
        {
            printf("%s\t\t-->\t\t%s\n",it->first.c_str(),it->second.c_str());
        }
    }

    std::size_t size()
    {
        return m_Arguments.size();
    }
    
private:
    std::map<std::string,std::string> m_Arguments;
};

#endif // __COMMAND_LINE_ARGUMENTS_PARSER_H__