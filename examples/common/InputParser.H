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

    template<typename T>
    void add(std::string a_argName, T& a_defaultArg)
    {
        getArgs<T>()[a_argName] = &a_defaultArg;
    }

    void parse(int a_argc, char** a_args)
    {
        int argNum = 0;
        while (argNum < a_argc)
        {
            char* arg = a_args[argNum];
            if (arg[0] == '-')
            {
                std::string argName(arg+1);
                if (doubleArgs.count(argName) > 0)
                {
                    double argVal = atof(a_args[argNum+1]);
                    *doubleArgs[argName] = argVal;
                    argNum += 2;
                } else if (intArgs.count(argName) > 0)
                {
                    int argVal = atoi(a_args[argNum+1]);
                    *intArgs[argName] = argVal;
                    argNum += 2;
                } else if (boolArgs.count(argName) > 0)
                {
                    bool argVal = atoi(a_args[argNum+1]);
                    *boolArgs[argName] = argVal;
                    argNum += 2;
                } else if (stringArgs.count(argName) > 0)
                {
                    std::string argVal(a_args[argNum+1]);
                    *stringArgs[argName] = argVal;
                    argNum += 2;
                } else {
                    pout() << "InputArgs::parse | Warning: Unrecognized argument: " << argName << std::endl;
                    argNum += 2;
                }
            } else {
                argNum ++;
            }
        }

    }
    
    /*
    void parse(int a_argc, char** a_args)
    {
        parse("inputs");
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
        pout() << "Finished parsing command line" << std::endl;
    }
    
    void parse(std::string a_inputsFile = "inputs")
    {
        std::ifstream inputs;
        inputs.open(a_inputsFile);

        if (inputs.fail())
        {
            pout() << "InputArgs::parse | Warning: Could not find specified inputs file: " << a_inputsFile << std::endl;
            return; 
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
    */

    void print()
    {
        pout() << setw(50) << setfill('=') << "=" << std::endl;
        pout() << "Input Arguments:" << std::endl;
        pout() << setw(50) << setfill('-') << "-" << std::endl;
        for (auto iter = stringArgs.begin(); iter != stringArgs.end(); ++iter)
        {
            pout() << setw(30) << left << setfill(' ') << iter->first;
            pout() << setw(20) << left << setfill(' ') << *iter->second << std::endl;
        }
        for (auto iter = boolArgs.begin(); iter != boolArgs.end(); ++iter)
        {
            pout() << setw(30) << left << setfill(' ') << iter->first;
            pout() << setw(20) << left << setfill(' ') << *iter->second << std::endl;
        }
        for (auto iter = intArgs.begin(); iter != intArgs.end(); ++iter)
        {
            pout() << setw(30) << left << setfill(' ') << iter->first;
            pout() << setw(20) << left << setfill(' ') << *iter->second << std::endl;
        }
        for (auto iter = doubleArgs.begin(); iter != doubleArgs.end(); ++iter)
        {
            pout() << setw(30) << left << setfill(' ') << iter->first;
            pout() << setw(20) << left << setfill(' ') << *iter->second << std::endl;
        }
        pout() << setw(50) << setfill('=') << "=" << std::endl;

    }
    /*
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
    */
    
    template<typename T>
    std::map<std::string, T*>& getArgs(){}
    
    private:

    //std::map<std::string, double> args;

    std::map<std::string, double*> doubleArgs;
    std::map<std::string, int*> intArgs;
    std::map<std::string, bool*> boolArgs;
    std::map<std::string, std::string*> stringArgs;

};

template<>
std::map<std::string, double*>& InputArgs::getArgs(){return doubleArgs;}

template<>
std::map<std::string, int*>& InputArgs::getArgs(){return intArgs;}

template<>
std::map<std::string, bool*>& InputArgs::getArgs(){return boolArgs;}

template<>
std::map<std::string, std::string*>& InputArgs::getArgs(){return stringArgs;}

#endif //end include guard
