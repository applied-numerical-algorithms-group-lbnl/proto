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
    int init_condition_type; 
    // 0 for radial shell, 
    // 1 for nonradial pulse, 
    // 2 for radial outflow
    // 3 for H5 BCs
    // 4 for radial shell with small constant cartesian velocity
    // 5 for radial Shock
    // 6 for radial shell with small constant cartesian magnetic field
    int domainSize; // Number of cells in non radial direction
    int thickness; // Number of cells in radial direction
    int boxSize_nonrad; // Size of the box for domain decomposition
    int boxSize_rad; // Size of the box for domain decomposition
    int max_iter; 
    double CFL;
    int limiter; // 1:4th order, 2: Minmod
    double strength; // Strength of perturbation for test problems.
    int write_cadence; 
    int checkpoint_cadence;
    string data_file_prefix;
    string checkpoint_file_prefix;
    int max_checkpoint_files;
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
    int MBInterp_define;   // old define (0) new define (1)
    string BC_file;
    string Probe_trajectory_file;
    string Probe_data_file;
    double Probe_cadence; // In seconds
    double BC_start_time;
    double BC_cadence; // In hours
    int BC_frame_rotate; // 0: Frames don't rotate (SWiG), 1: Frames rotate
    int BC_force_single_frame; // 1: Use only one frame for BCs, 0: Use all frames for BCs (if present)
    int CME_type; // 0: No CME, 1: CME
    double CME_Enter_Time;
    double CME_r1;
    double CME_r0;
    double CME_halfangle;
    double CME_pol_flux_control;
    double CME_FRiED_n;
    double CME_lat;
    double CME_lon;
    double CME_tilt;
    double CME_Vandas_R0;
    double CME_apex_speed;
    double CME_energy_control;
    double CME_density;
    double CME_helicity_sign;
    double CME_a1;
    double CME_shift;
    
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
    static int get_limiter() {return getInstance().limiter;}
    static double get_strength() {return getInstance().strength;}
    static int get_write_cadence() {return getInstance().write_cadence;}
    static int get_checkpoint_cadence() {return getInstance().checkpoint_cadence;}
    static string get_data_file_prefix() {return getInstance().data_file_prefix;}
    static string get_checkpoint_file_prefix() {return getInstance().checkpoint_file_prefix;}
    static int get_max_checkpoint_files() {return getInstance().max_checkpoint_files;}
    static string get_restart_file() {return getInstance().restart_file;}
    static double get_gamma() {return getInstance().gamma;}
    static int get_temporal_order() {return getInstance().temporal_order;}
    static int get_divergence_cleaning() {return getInstance().divergence_cleaning;}
    static int get_linear_viscosity() {return getInstance().linear_viscosity;}
    static int get_boundary_reflux() {return getInstance().boundary_reflux;}
    static int get_radial_refinement() {return getInstance().radial_refinement;}
    static int get_operator_switchiv() {return getInstance().operator_switchiv;}
    static int get_operator_switchln() {return getInstance().operator_switchln;}
    static int get_MBInterp_define() {return getInstance().MBInterp_define;}
    static string get_BC_file() {return getInstance().BC_file;}
    static string get_Probe_trajectory_file() {return getInstance().Probe_trajectory_file;}
    static string get_Probe_data_file() {return getInstance().Probe_data_file;}
    static double get_Probe_cadence() {return getInstance().Probe_cadence;}
    static double get_BC_start_time() {return getInstance().BC_start_time;}
    static double get_BC_cadence() {return getInstance().BC_cadence;}
    static int get_BC_frame_rotate() {return getInstance().BC_frame_rotate;}
    static int get_BC_force_single_frame() {return getInstance().BC_force_single_frame;}
    static int get_CME_type() {return getInstance().CME_type;}
    static double get_CME_Enter_Time() {return getInstance().CME_Enter_Time;}
    static double get_CME_r1() {return getInstance().CME_r1;}
    static double get_CME_r0() {return getInstance().CME_r0;}
    static double get_CME_halfangle() {return getInstance().CME_halfangle;}
    static double get_CME_pol_flux_control() {return getInstance().CME_pol_flux_control;}
    static double get_CME_FRiED_n() {return getInstance().CME_FRiED_n;}
    static double get_CME_lat() {return getInstance().CME_lat;}
    static double get_CME_lon() {return getInstance().CME_lon;}
    static double get_CME_tilt() {return getInstance().CME_tilt;}
    static double get_CME_Vandas_R0() {return getInstance().CME_Vandas_R0;}
    static double get_CME_apex_speed() {return getInstance().CME_apex_speed;}
    static double get_CME_energy_control() {return getInstance().CME_energy_control;}
    static double get_CME_density() {return getInstance().CME_density;}
    static double get_CME_helicity_sign() {return getInstance().CME_helicity_sign;}
    static double get_CME_a1() {return getInstance().CME_a1;}
    static double get_CME_shift() {return getInstance().CME_shift;}

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
        args.add("limiter", limiter);
        args.add("strength", strength);
        args.add("write_cadence", write_cadence);
        args.add("checkpoint_cadence", checkpoint_cadence);
        args.add("data_file_prefix", data_file_prefix);
        args.add("checkpoint_file_prefix", checkpoint_file_prefix);
        args.add("max_checkpoint_files", max_checkpoint_files);
        args.add("restart_file", restart_file);
        args.add("gamma", gamma);
        args.add("temporal_order", temporal_order);
        args.add("divergence_cleaning", divergence_cleaning);
        args.add("linear_viscosity", linear_viscosity);
        args.add("boundary_reflux", boundary_reflux);
        args.add("radial_refinement", radial_refinement);
        args.add("operator_switchiv", operator_switchiv);
        args.add("operator_switchln", operator_switchln);
        args.add("MBInterp_define", MBInterp_define);
        args.add("BC_file", BC_file);
        args.add("Probe_trajectory_file", Probe_trajectory_file);
        args.add("Probe_data_file", Probe_data_file);
        args.add("Probe_cadence", Probe_cadence);
        args.add("BC_start_time", BC_start_time);
        args.add("BC_cadence", BC_cadence);
        args.add("BC_frame_rotate", BC_frame_rotate);
        args.add("BC_force_single_frame", BC_force_single_frame);
        args.add("CME_Enter_Time", CME_Enter_Time);
        args.add("CME_type", CME_type);
        args.add("CME_r1", CME_r1);
        args.add("CME_r0", CME_r0);
        args.add("CME_halfangle", CME_halfangle);
        args.add("CME_pol_flux_control", CME_pol_flux_control);
        args.add("CME_FRiED_n", CME_FRiED_n);
        args.add("CME_lat", CME_lat);
        args.add("CME_lon", CME_lon);
        args.add("CME_tilt", CME_tilt);
        args.add("CME_Vandas_R0", CME_Vandas_R0);
        args.add("CME_apex_speed", CME_apex_speed);
        args.add("CME_energy_control", CME_energy_control);
        args.add("CME_density", CME_density);
        args.add("CME_helicity_sign", CME_helicity_sign);
        args.add("CME_a1", CME_a1);
        args.add("CME_shift", CME_shift);
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
            cout << "limiter: " << limiter << endl;
            cout << "strength: " << strength << endl;
            cout << "write_cadence: " << write_cadence << endl;
            cout << "checkpoint_cadence: " << checkpoint_cadence << endl;
            cout << "data_file_prefix: " << data_file_prefix << endl;
            cout << "checkpoint_file_prefix: " << checkpoint_file_prefix << endl;
            cout << "max_checkpoint_files: " << max_checkpoint_files << endl;
            cout << "restart_file: " << restart_file << endl;
            cout << "gamma: " << gamma << endl;
            cout << "temporal_order: " << temporal_order << endl;
            cout << "divergence_cleaning: " << divergence_cleaning << endl;
            cout << "linear_viscosity: " << linear_viscosity << endl;
            cout << "boundary_reflux: " << boundary_reflux << endl;
            cout << "radial_refinement: " << radial_refinement << endl;
            cout << "operator_switchiv: " << operator_switchiv << endl;
            cout << "operator_switchln: " << operator_switchln << endl;
            cout << "MBInterp_define: " << MBInterp_define << endl;
            cout << "BC_file: " << BC_file << endl;
            cout << "Probe_trajectory_file: " << Probe_trajectory_file << endl;
            cout << "Probe_data_file: " << Probe_data_file << endl;
            cout << "Probe_cadence: " << Probe_cadence << endl;
            cout << "BC_start_time: " << BC_start_time << endl;
            cout << "BC_cadence: " << BC_cadence << endl;
            cout << "BC_frame_rotate: " << BC_frame_rotate << endl;
            cout << "BC_force_single_frame: " << BC_force_single_frame << endl;
            cout << "CME_type: " << CME_type << endl;
            cout << "CME_Enter_Time: " << CME_Enter_Time << endl;
            cout << "CME_r1: " << CME_r1 << endl;
            cout << "CME_r0: " << CME_r0 << endl;
            cout << "CME_halfangle: " << CME_halfangle << endl;
            cout << "CME_pol_flux_control: " << CME_pol_flux_control << endl;
            cout << "CME_FRiED_n: " << CME_FRiED_n << endl;
            cout << "CME_lat: " << CME_lat << endl;
            cout << "CME_lon: " << CME_lon << endl;
            cout << "CME_tilt: " << CME_tilt << endl;
            cout << "CME_Vandas_R0: " << CME_Vandas_R0 << endl;
            cout << "CME_apex_speed: " << CME_apex_speed << endl;
            cout << "CME_energy_control: " << CME_energy_control << endl;
            cout << "CME_density: " << CME_density << endl;
            cout << "CME_helicity_sign: " << CME_helicity_sign << endl;
            cout << "CME_a1: " << CME_a1 << endl;
            cout << "CME_shift: " << CME_shift << endl;
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

