#pragma once

#include "Proto.H"
#include "InputParser.H"

using namespace Proto;

class Parsefrominputs {
  public:
    double tstop;
	double CFL;
	double domsizex;
	double domsizey;
	double domsizez;
	double gamma;
	int domainSizex;
	int domainSizey;
	int domainSizez;
	int maxStep;
	int restartStep;
	int outputInterval;
	int CheckpointInterval;
    int BoxSize;
    int grid_type_global; // 0: 2D-Rectangular/3D-Rectangular;  1: 2D-Wavy/3D-Not implemented;  2: 2D-Polar/3D-Spherical
    int LowBoundType;  // 0 for periodic, 1 for Dirichlet, 2 for open. This is for dir==0 only
    int HighBoundType;  // 0 for periodic, 1 for Dirichlet, 2 for open. This is for dir==0 only
    
    int limiter_apply;
    int slope_flattening_apply;
    int non_linear_visc_apply;
    int linear_visc_apply;
    int pole_correction;
    int takedivBstep;
    int timeIntegratorType; // 1: Euler Step    4: RK4 
    int init_condition_type;
    /*
	   0. constant solution
	   1. 2D current sheet problem
	   2. flow from 1 side
	   3. 2D Orszag Tang problem
	   4. Alfven wave problem
	   5. Acoustic pulse problem
	   6. Acoustic pulse problem in polar grid
	   7. Cylindrical pulse problem in polar grid
	   8. Shifted Acoustic pulse problem
	   9. Euler problem
	   10. Radially out flow in polar grid
	   11. Cylindrical pulse problem in cartesian grid
	   12. Acoustic pulse problem in 3D cartesian grid
	   13. 3D MHD blast wave
	   14. 2D MHD blast wave
	   15. Acoustic pulse problem with Bx
	   16. Shell acoustic pulse problem in spherical grid
	   17. Sphere Shifted Acoustic pulse problem in Spherical grid

	   19. Velocity pulse problem in polar grid
	   21. Radially out flow in spherical grid
	 */
    int Riemann_solver_type;
    int initialize_in_spherical_coords;
    double C1_fix; // A constant in wavy grid definition.
    double r_in;
    double r_out;
    double C_rad; // A constant in exponential dr in spherical grid.
    int convTestType; // 0 for no convergence test, 1 for space convergence, 2 for space and time convergence
	int saveConvTestData;
  public:
	void parsenow(){
		InputArgs args;
		args.parse(); //assumes a file "inputs" exists.
		args.set("tstop", &tstop);
		args.set("CFL", &CFL);
		args.set("domsizex", &domsizex);
		args.set("domsizey", &domsizey);
		args.set("domsizez", &domsizez);
		args.set("gamma", &gamma);
		args.set("domainSizex", &domainSizex);
		args.set("domainSizey", &domainSizey);
		args.set("domainSizez", &domainSizez);
		args.set("maxStep", &maxStep);
		args.set("restartStep", &restartStep);
		args.set("outputInterval", &outputInterval);
		args.set("CheckpointInterval", &CheckpointInterval);
		args.set("BoxSize", &BoxSize);
		args.set("grid_type_global", &grid_type_global);
		args.set("LowBoundType", &LowBoundType);
		args.set("HighBoundType", &HighBoundType);
		args.set("limiter_apply", &limiter_apply);
		args.set("slope_flattening_apply", &slope_flattening_apply);
		args.set("non_linear_visc_apply", &non_linear_visc_apply);
		args.set("linear_visc_apply", &linear_visc_apply);
		args.set("pole_correction", &pole_correction);
		args.set("takedivBstep", &takedivBstep);
		args.set("timeIntegratorType", &timeIntegratorType);
		args.set("init_condition_type", &init_condition_type);
		args.set("Riemann_solver_type", &Riemann_solver_type);
		args.set("initialize_in_spherical_coords", &initialize_in_spherical_coords);
		args.set("C1_fix", &C1_fix);
		args.set("r_in", &r_in);
		args.set("r_out", &r_out);
		args.set("C_rad", &C_rad);
		args.set("convTestType", &convTestType);
		args.set("saveConvTestData", &saveConvTestData);

		int pid = procID();

		if (pid == 0){
		cout << "Running with the following inputs (these can be changed in the file 'inputs')" << endl;
		cout << setw(75) << setfill('-') << "-" << endl;
		cout << "tstop: " << tstop << endl;
		cout << "CFL: " << CFL << endl;
		cout << "domsizex: " << domsizex << endl;
		cout << "domsizey: " << domsizey << endl;
		cout << "domsizez: " << domsizez << endl;
		cout << "gamma: " << gamma << endl;
		cout << "domainSizex: "  << domainSizex << endl;
		cout << "domainSizey: "  << domainSizey << endl;
		cout << "domainSizez: "  << domainSizez << endl;
		cout << "maxStep: "    << maxStep << endl;
		cout << "restartStep: "    << restartStep << endl;
		cout << "outputInterval: "    << outputInterval << endl;	
		cout << "CheckpointInterval: "    << CheckpointInterval << endl;	
		cout << "BoxSize: "    << BoxSize << endl;
		cout << "grid_type_global: "    << grid_type_global << endl;
		cout << "LowBoundType: "    << LowBoundType << endl;
		cout << "HighBoundType: "    << HighBoundType << endl;
		cout << "limiter_apply: "    << limiter_apply << endl;
		cout << "slope_flattening_apply: "    << slope_flattening_apply << endl;
		cout << "non_linear_visc_apply: "    << non_linear_visc_apply << endl;
		cout << "linear_visc_apply: "    << linear_visc_apply << endl;
		cout << "pole_correction: "    << pole_correction << endl;
		cout << "takedivBstep: "    << takedivBstep << endl;
		cout << "timeIntegratorType: "    << takedivBstep << endl;
		cout << "init_condition_type: "    << init_condition_type << endl;
		cout << "Riemann_solver_type: "    << Riemann_solver_type << endl;
		cout << "initialize_in_spherical_coords: "    << initialize_in_spherical_coords << endl;
		cout << "C1_fix: "    << C1_fix << endl;
		cout << "r_in: "    << r_in << endl;
		cout << "r_out: "    << r_out << endl;
		cout << "C_rad: "    << C_rad << endl;
		cout << "convTestType: "    << convTestType << endl;
		cout << "saveConvTestData: "    << saveConvTestData << endl;
		cout << setw(75) << setfill('-') << "-" << endl;
		}

	}
    //Use the following format if encapsulation is needed (if variables needed to be made private) 
    //Getter
	int getBoxSize() {
      return BoxSize;
    }
	// Setter
    void setBoxSize(int s) {
      BoxSize = s;
    }

};
