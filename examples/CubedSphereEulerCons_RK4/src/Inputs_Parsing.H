#pragma once

#include "Proto.H"
#include "InputParser.H"
#include <algorithm> // For std::remove
#include <iostream>  // For std::cout, etc.
#include <iomanip>   // For std::setw, etc.

using namespace Proto;
using namespace std;

/// @brief ParseInputs Singleton class
class ParseInputs {
private:
    int convTestType; // 0: No convergence test, 1: Space Convergence test, 2: Space and Time Convergence test
    int init_condition_type; // 0 for radial pulse, 1 for nonradial pulse, 2 for radial outflow
    int domainSize; // Number of cells in non radial direction
    int thickness; // Number of cells in radial direction
    int max_iter; 
    double CFL;
    int write_cadence; 

    // Private Constructor
    ParseInputs() {}

    // Private Copy Constructor and Assignment Operator
    ParseInputs(const ParseInputs&) = delete;
    ParseInputs& operator=(const ParseInputs&) = delete;

public:
    // Static method to get the instance of the class
    static ParseInputs& getInstance() {
        static ParseInputs instance; // Guaranteed to be destroyed and instantiated on first use
        return instance;
    }

    // Static Getters
    static int get_convTestType() {return getInstance().convTestType;}
    static int get_init_condition_type() {return getInstance().init_condition_type;}
    static int get_domainSize() {return getInstance().domainSize;}
    static int get_thickness() {return getInstance().thickness;}
    static int get_max_iter() {return getInstance().max_iter;}
    static double get_CFL() {return getInstance().CFL;}
    static int get_write_cadence() {return getInstance().write_cadence;}

    // Parses now with command-line arguments
    void parsenow(int argc, char* argv[]) {
        InputArgs args;
        args.add("convTestType", convTestType);
        args.add("init_condition_type", init_condition_type);
        args.add("domainSize", domainSize);
        args.add("thickness", thickness);
        args.add("max_iter", max_iter);
        args.add("CFL", CFL);
        args.add("write_cadence", write_cadence);
        args.parse(argc, argv);

        int pid = procID();
        if (pid == 0) {
            cout << "Running with the following inputs (these can be changed in the file 'inputs')" << endl;
            cout << setw(75) << setfill('-') << "-" << endl;
            cout << "convTestType: " << convTestType << endl;
            cout << "init_condition_type: " << init_condition_type << endl;
            cout << "domainSize: " << domainSize << endl;
            cout << "thickness: " << thickness << endl;
            cout << "max_iter: " << max_iter << endl;
            cout << "CFL: " << CFL << endl;
            cout << "write_cadence: " << write_cadence << endl;
            cout << setw(75) << setfill('-') << "-" << endl;
        }
    }

    // Utility function to remove spaces from a string
    string removeSpaces(string str) {
        str.erase(remove(str.begin(), str.end(), ' '), str.end());
        return str;
    }

    // Additional getters and setters as needed
};

