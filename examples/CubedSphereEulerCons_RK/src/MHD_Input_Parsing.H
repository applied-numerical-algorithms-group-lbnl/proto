#pragma once

#include "Proto.H"
#include "InputParser.H"

using namespace Proto;
/// @brief Parsefrominputs class
class Parsefrominputs
{
public:
	double tstop;
	double CFL;
	double gamma;
	int domainSizex;
	int domainSizey;
	int domainSizez;
	int maxStep;
	int restartStep;
	int outputInterval;
	int CheckpointInterval;
	int MaxCheckpointFiles;
	int BoxSize;
	int grid_type_global;  // 0: 2D-Rectangular/3D-Rectangular;  1: 2D-Wavy/3D-Not implemented;  2: 2D-Polar/3D-Spherical
	int LowBoundType;	   // 0 for periodic, 1 for Dirichlet, 2 for open. This is for dir==0 only
	int HighBoundType;	   // 0 for periodic, 1 for Dirichlet, 2 for open. This is for dir==0 only
	int sph_inner_BC_hdf5; // 1 is h5 file is used to provide inner BCs
	int limiter_apply;
	int slope_flattening_apply;
	int non_linear_visc_apply;
	int linear_visc_apply;
	int takedivBstep;
	int timeIntegratorType; // 1: Euler Step    4: RK4
	int init_condition_type;
	int Riemann_solver_type; // 1: Rusanov Solver, 2: Roe Solver
	double entropy_fix_coeff;
	double r_in;
	double r_out;

public:
	string removeSpaces(string str)
	{
		str.erase(remove(str.begin(), str.end(), ' '), str.end());
		return str;
	}

	void parsenow(int argc, char *argv[])
	{
		InputArgs args;
		args.add("tstop", tstop);
		args.add("CFL", CFL);
		args.add("gamma", gamma);
		args.add("domainSizex", domainSizex);
		args.add("domainSizey", domainSizey);
		args.add("domainSizez", domainSizez);
		args.add("maxStep", maxStep);
		args.add("restartStep", restartStep);
		args.add("outputInterval", outputInterval);
		args.add("CheckpointInterval", CheckpointInterval);
		args.add("MaxCheckpointFiles", MaxCheckpointFiles);
		args.add("BoxSize", BoxSize);
		args.add("grid_type_global", grid_type_global);
		args.add("LowBoundType", LowBoundType);
		args.add("HighBoundType", HighBoundType);
		args.add("sph_inner_BC_hdf5", sph_inner_BC_hdf5);
		args.add("limiter_apply", limiter_apply);
		args.add("slope_flattening_apply", slope_flattening_apply);
		args.add("non_linear_visc_apply", non_linear_visc_apply);
		args.add("linear_visc_apply", linear_visc_apply);
		args.add("takedivBstep", takedivBstep);
		args.add("timeIntegratorType", timeIntegratorType);
		args.add("init_condition_type", init_condition_type);
		args.add("Riemann_solver_type", Riemann_solver_type);
		args.add("entropy_fix_coeff", entropy_fix_coeff);
		args.add("r_in", r_in);
		args.add("r_out", r_out);
		args.parse(argc, argv);
		// args.print(); //It prints, but not in the order I want.

		int pid = procID();

		if (pid == 0)
		{
			cout << "Running with the following inputs (these can be changed in the file 'inputs')" << endl;
			cout << setw(75) << setfill('-') << "-" << endl;
			cout << "tstop: " << tstop << endl;
			cout << "CFL: " << CFL << endl;
			cout << "gamma: " << gamma << endl;
			cout << "domainSizex: " << domainSizex << endl;
			cout << "domainSizey: " << domainSizey << endl;
			cout << "domainSizez: " << domainSizez << endl;
			cout << "maxStep: " << maxStep << endl;
			cout << "restartStep: " << restartStep << endl;
			cout << "outputInterval: " << outputInterval << endl;
			cout << "CheckpointInterval: " << CheckpointInterval << endl;
			cout << "MaxCheckpointFiles: " << MaxCheckpointFiles << endl;
			cout << "BoxSize: " << BoxSize << endl;
			cout << "grid_type_global: " << grid_type_global << endl;
			cout << "LowBoundType: " << LowBoundType << endl;
			cout << "HighBoundType: " << HighBoundType << endl;
			cout << "sph_inner_BC_hdf5: " << sph_inner_BC_hdf5 << endl;
			cout << "limiter_apply: " << limiter_apply << endl;
			cout << "slope_flattening_apply: " << slope_flattening_apply << endl;
			cout << "non_linear_visc_apply: " << non_linear_visc_apply << endl;
			cout << "linear_visc_apply: " << linear_visc_apply << endl;
			cout << "takedivBstep: " << takedivBstep << endl;
			cout << "timeIntegratorType: " << takedivBstep << endl;
			cout << "init_condition_type: " << init_condition_type << endl;
			cout << "Riemann_solver_type: " << Riemann_solver_type << endl;
			cout << "entropy_fix_coeff: " << entropy_fix_coeff << endl;
			cout << "r_in: " << r_in << endl;
			cout << "r_out: " << r_out << endl;
			cout << setw(75) << setfill('-') << "-" << endl;
		}
	}
	// Use the following format if encapsulation is needed (if variables needed to be made private)
	// Getter
	int getBoxSize()
	{
		return BoxSize;
	}
	// Setter
	void setBoxSize(int s)
	{
		BoxSize = s;
	}
};