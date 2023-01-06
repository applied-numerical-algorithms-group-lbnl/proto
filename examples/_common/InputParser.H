#pragma once
#ifndef _PROTO_INPUT_PARSER_
#define _PROTO_INPUT_PARSER_
#include "Proto.H"
#include <map>
#include <istream>
#include <iomanip>
#include <streambuf>
#include <fstream>
#include <regex>

using Proto::pout; using Proto::procID;

class InputArgs
{
    public:

    inline std::string trimFront(std::string input)
    {
        return std::regex_replace(input, std::regex("^\\s+"), string(""));
    }

    inline std::string trimBack(std::string input)
    {
        return std::regex_replace(input, std::regex("\\s+$"), string(""));
    }

    inline bool isWhitespace(std::string s)
    {
        bool blank = false;
        if (s.empty() || std::all_of(s.begin(), s.end(), [](char c){return std::isspace(c);})) {
            blank = true;
        }
        return blank;
    }

    template<typename T>
    inline void add(std::string a_argName, T& a_defaultArg)
    {
        getArgs<T>()[a_argName] = &a_defaultArg;
    }

    inline void parse(std::vector<std::string> a_vals)
    {
        int argNum = 0;
        while (argNum < a_vals.size())
        {
            auto arg = a_vals[argNum];
            if (arg[0] == '-')
            {
                auto argName = arg.substr(1);
                if (doubleArgs.count(argName) > 0)
                {
                    double argVal = stof(a_vals[argNum+1]);
                    *doubleArgs[argName] = argVal;
                    argNum += 2;
                } else if (intArgs.count(argName) > 0)
                {
                    int argVal = stoi(a_vals[argNum+1]);
                    *intArgs[argName] = argVal;
                    argNum += 2;
                } else if (boolArgs.count(argName) > 0)
                {
                    bool argVal = stoi(a_vals[argNum+1]);
                    *boolArgs[argName] = argVal;
                    argNum += 2;
                } else if (stringArgs.count(argName) > 0)
                {
                    std::string argVal(a_vals[argNum+1]);
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

    inline void parse(int a_argc, char** a_args)
    {
        std::vector<std::string> vals;
        if (a_argc == 2)
        {
            //args have been specified using a file
            std::ifstream file;
            file.open(a_args[1]);
            while(!file.eof())
            {
                std::string var;
                std::string val;
                std::string line;
                std::getline(file, line);
                if (!isWhitespace(line))
                {
                    line = trimFront(line);
                    line = trimBack(line);
                    std::size_t pos = line.find_first_of(' ');
                    var = line.substr(0,pos);
                    line.erase(0,pos);
                    line = trimFront(line);
                    val = line;
                    vals.push_back(var);
                    vals.push_back(val);
                }
            }
        } else {
            for (int ii = 0; ii < a_argc; ii++)
            {
                vals.push_back(std::string(a_args[ii]));
            }
        }
        parse(vals);
    }
    
    inline void print()
    {
        using Proto::pout;
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
    inline std::map<std::string, T*>& getArgs(){}
    
    private:

    //std::map<std::string, double> args;

    std::map<std::string, double*> doubleArgs;
    std::map<std::string, int*> intArgs;
    std::map<std::string, bool*> boolArgs;
    std::map<std::string, std::string*> stringArgs;

};

template<>
inline std::map<std::string, double*>& InputArgs::getArgs(){return doubleArgs;}

template<>
inline std::map<std::string, int*>& InputArgs::getArgs(){return intArgs;}

template<>
inline std::map<std::string, bool*>& InputArgs::getArgs(){return boolArgs;}

template<>
inline std::map<std::string, std::string*>& InputArgs::getArgs(){return stringArgs;}

#endif //end include guard
