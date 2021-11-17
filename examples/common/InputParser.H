#include <map>
#include <istream>

class InputArgs
{
    public:

    void parse(std::string a_inputsFile = "inputs")
    {
        std::ifstream inputs;
        inputs.open(a_inputsFile);

        if (inputs.fail())
        {
            std::cout << "InputArgs::parse | Error: Invalid inputs file: " << a_inputsFile << std::endl;
            std::abort();
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
