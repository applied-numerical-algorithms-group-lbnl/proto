#pragma once
#ifndef _PROTO_INPUT_PARSER_
#define _PROTO_INPUT_PARSER_
#include "Proto.H"
#include <map>
#include <istream>
#include <iomanip>

class InputArgs
{
    public:

    void parse(int a_argc, char** a_args)
    {
        int argNum = 0;
        while (argNum < a_argc)
        {
            char* arg = a_args[argNum];
            if (arg[0] == '-')
            {
                double argVal = atof(a_args[argNum+1]);
                pout() << arg+1 << " = " << argVal << std::endl;
                args[std::string(arg+1)] = argVal;
                argNum += 2;
            } else {
                argNum ++;
            }
        }
    }

    void parse(std::string a_inputsFile = "inputs")
    {
        std::ifstream inputs;
        inputs.open(a_inputsFile);

        if (inputs.fail())
        {
            pout() << "InputArgs::parse | Warning: Could not find specified inputs file: " << a_inputsFile << std::endl;
        }

        std::string s;
        double value;
        while(!inputs.eof())
        {
            inputs >> s;
            inputs >> (value);
            args[s] = value;
        }
    }

    void print()
    {
        pout() << setw(50) << setfill('=') << "=" << std::endl;
        pout() << "Input Arguments:" << std::endl;
        pout() << setw(50) << setfill('-') << "-" << std::endl;
        for (auto iter = args.begin(); iter != args.end(); ++iter)
        {
            pout() << setw(30) << left << setfill(' ') << iter->first;
            pout() << setw(20) << left << setfill(' ') << iter->second << std::endl;
        }
        pout() << setw(50) << setfill('=') << "=" << std::endl;
    }
    
    double get(std::string a_varname)
    {
        PROTO_ASSERT(args.count(a_varname) > 0,
            "InputArgs::get | Error: Cannot find variable %s.", a_varname.c_str());
        return args[a_varname];
    }
    
    template<typename T>
    void set(std::string a_varname, T* a_var)
    {
        if (args.count(a_varname) > 0)
        {
            *a_var = (T)args[a_varname];
        } 
    }
    private:

    std::map<std::string, double> args;
};
#endif //end include guard
