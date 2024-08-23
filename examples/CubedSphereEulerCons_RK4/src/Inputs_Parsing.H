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
    int boxSize_nonrad; // Size of the box for domain decomposition
    int boxSize_rad; // Size of the box for domain decomposition
    int max_iter; 
    double CFL;
  double strength; // Strength of perturbation for test problems.
    int write_cadence; 
    int checkpoint_cadence;
    string data_file_prefix;
    string checkpoint_file_prefix;
    string restart_file;
    double gamma; // Ratio of specific heats
    int temporal_order; // Temporal order of accuracy
    int divergence_cleaning; // 0: No divergence cleaning, 1: Divergence cleaning
    int linear_viscosity; // 0: No linear viscosity, 1: Linear viscosity
    int boundary_reflux; // 0: don't reflux at block boundaries, 1: reflux
    int radial_refinement; // 0: refine in all directions by 2 for convergence studies;
                           // 1 refine only in radial direction.
    int operator_switchiv;   // no inviscid (1) no viscous (2) both on (0).
    int operator_switchln;   // no linear (1) no nonlinear (2) both on (0)
    string BC_file;
    string Probe_trajectory_file;
    string Probe_data_file;
    double Probe_cadence; // In seconds
    double BC_start_time;
    double BC_cadence; // In hours

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
    static int get_boxSize_nonrad() {return getInstance().boxSize_nonrad;}
    static int get_boxSize_rad() {return getInstance().boxSize_rad;}
    static int get_max_iter() {return getInstance().max_iter;}
    static double get_CFL() {return getInstance().CFL;}
    static double get_strength() {return getInstance().strength;}
    static int get_write_cadence() {return getInstance().write_cadence;}
    static int get_checkpoint_cadence() {return getInstance().checkpoint_cadence;}
    static string get_data_file_prefix() {return getInstance().data_file_prefix;}
    static string get_checkpoint_file_prefix() {return getInstance().checkpoint_file_prefix;}
    static string get_restart_file() {return getInstance().restart_file;}
    static double get_gamma() {return getInstance().gamma;}
    static int get_temporal_order() {return getInstance().temporal_order;}
    static int get_divergence_cleaning() {return getInstance().divergence_cleaning;}
    static int get_linear_viscosity() {return getInstance().linear_viscosity;}
    static int get_boundary_reflux() {return getInstance().boundary_reflux;}
    static int get_radial_refinement() {return getInstance().radial_refinement;}
    static int get_operator_switchiv() {return getInstance().operator_switchiv;}
    static int get_operator_switchln() {return getInstance().operator_switchln;}
    static string get_BC_file() {return getInstance().BC_file;}
    static string get_Probe_trajectory_file() {return getInstance().Probe_trajectory_file;}
    static string get_Probe_data_file() {return getInstance().Probe_data_file;}
    static double get_Probe_cadence() {return getInstance().Probe_cadence;}
    static double get_BC_start_time() {return getInstance().BC_start_time;}
    static double get_BC_cadence() {return getInstance().BC_cadence;}

    // Parses now with command-line arguments
    void parsenow(int argc, char* argv[]) {
        InputArgs args;
        args.add("convTestType", convTestType);
        args.add("init_condition_type", init_condition_type);
        args.add("domainSize", domainSize);
        args.add("thickness", thickness);
        args.add("boxSize_nonrad", boxSize_nonrad);
        args.add("boxSize_rad", boxSize_rad);
        args.add("max_iter", max_iter);
        args.add("CFL", CFL);
        args.add("strength", strength);
        args.add("write_cadence", write_cadence);
        args.add("checkpoint_cadence", checkpoint_cadence);
        args.add("data_file_prefix", data_file_prefix);
        args.add("checkpoint_file_prefix", checkpoint_file_prefix);
        args.add("restart_file", restart_file);
        args.add("gamma", gamma);
        args.add("temporal_order", temporal_order);
        args.add("divergence_cleaning", divergence_cleaning);
        args.add("linear_viscosity", linear_viscosity);
        args.add("boundary_reflux", boundary_reflux);
        args.add("radial_refinement", radial_refinement);
        args.add("operator_switchiv", operator_switchiv);
        args.add("operator_switchln", operator_switchln);
        args.add("BC_file", BC_file);
        args.add("Probe_trajectory_file", Probe_trajectory_file);
        args.add("Probe_data_file", Probe_data_file);
        args.add("Probe_cadence", Probe_cadence);
        args.add("BC_start_time", BC_start_time);
        args.add("BC_cadence", BC_cadence);
        args.parse(argc, argv);

        int pid = procID();
        if (pid == 0) {
            cout << "Running with the following inputs (these can be changed in the file 'inputs')" << endl;
            cout << setw(75) << setfill('-') << "-" << endl;
            cout << "convTestType: " << convTestType << endl;
            cout << "init_condition_type: " << init_condition_type << endl;
            cout << "domainSize: " << domainSize << endl;
            cout << "thickness: " << thickness << endl;
            cout << "boxSize_nonrad: " << boxSize_nonrad << endl;
            cout << "boxSize_rad: " << boxSize_rad << endl;
            cout << "max_iter: " << max_iter << endl;
            cout << "CFL: " << CFL << endl;
            cout << "strength: " << strength << endl;
            cout << "write_cadence: " << write_cadence << endl;
            cout << "checkpoint_cadence: " << checkpoint_cadence << endl;
            cout << "data_file_prefix: " << data_file_prefix << endl;
            cout << "checkpoint_file_prefix: " << checkpoint_file_prefix << endl;
            cout << "restart_file: " << restart_file << endl;
            cout << "gamma: " << gamma << endl;
            cout << "temporal_order: " << temporal_order << endl;
            cout << "divergence_cleaning: " << divergence_cleaning << endl;
            cout << "linear_viscosity: " << linear_viscosity << endl;
            cout << "boundary_reflux: " << boundary_reflux << endl;
            cout << "radial_refinement: " << radial_refinement << endl;
            cout << "operator_switchiv: " << operator_switchiv << endl;
            cout << "operator_switchln: " << operator_switchln << endl;
            cout << "BC_file: " << BC_file << endl;
            cout << "Probe_trajectory_file: " << Probe_trajectory_file << endl;
            cout << "Probe_data_file: " << Probe_data_file << endl;
            cout << "Probe_cadence: " << Probe_cadence << endl;
            cout << "BC_start_time: " << BC_start_time << endl;
            cout << "BC_cadence: " << BC_cadence << endl;
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

