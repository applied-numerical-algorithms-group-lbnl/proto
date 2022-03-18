#pragma once
#ifndef _PROTO_INPUT_PARSER_
#define _PROTO_INPUT_PARSER_
#include "Proto.H"
#include <map>
#include <istream>
#include <iomanip>
#include <streambuf>
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
#ifdef PR_MPI
        // also print to the command line for better user feedback
        if (procID() == 0)
        {
            std::cout << setw(50) << setfill('=') << "=" << std::endl;
            std::cout << "Input Arguments:" << std::endl;
            std::cout << setw(50) << setfill('-') << "-" << std::endl;
            for (auto iter = stringArgs.begin(); iter != stringArgs.end(); ++iter)
            {
                std::cout << setw(30) << left << setfill(' ') << iter->first;
                std::cout << setw(20) << left << setfill(' ') << *iter->second << std::endl;
            }
            for (auto iter = boolArgs.begin(); iter != boolArgs.end(); ++iter)
            {
                std::cout << setw(30) << left << setfill(' ') << iter->first;
                std::cout << setw(20) << left << setfill(' ') << *iter->second << std::endl;
            }
            for (auto iter = intArgs.begin(); iter != intArgs.end(); ++iter)
            {
                std::cout << setw(30) << left << setfill(' ') << iter->first;
                std::cout << setw(20) << left << setfill(' ') << *iter->second << std::endl;
            }
            for (auto iter = doubleArgs.begin(); iter != doubleArgs.end(); ++iter)
            {
                std::cout << setw(30) << left << setfill(' ') << iter->first;
                std::cout << setw(20) << left << setfill(' ') << *iter->second << std::endl;
            }
            std::cout << setw(50) << setfill('-') << "-" << std::endl;
            std::cout << "Set above variables on the command line using: " << std::endl;
            std::cout << ">> ./<executable> -<var_name> value" << std::endl;
            std::cout << setw(50) << setfill('=') << "=" << std::endl;
        }
#endif
    }
    
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
