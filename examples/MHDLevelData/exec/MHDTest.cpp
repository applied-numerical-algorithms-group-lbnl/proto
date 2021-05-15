#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "Proto.H"
#include "MHDLevelDataRK4.H"
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"
#include "MHD_Initialize.H"
#include "MHDOp.H"
#include "MHD_Mapping.H"

// For Chrono Timer (Talwinder)
#include <chrono>
#include <iostream>
//////////////////////////////

using namespace std;
using namespace Proto;

double time_globalll;
int grid_type_global;
int LowBoundType;
int HighBoundType;
bool limiter_apply;
bool slope_flattening_apply;
bool non_linear_visc_apply;
bool linear_visc_apply;
int init_condition_type;
int Riemann_solver_type;

void InitializeMHDLevelDataState(MHDLevelDataState& state)
{
	(state.m_U).setToZero();
	for(DataIterator dit=state.m_U.begin(); *dit!=dit.end(); ++dit)
		MHD_Initialize::initializeState((state.m_U)[*dit],state.m_dx,state.m_gamma);
}

/**
 * @brief Write out all components of data to filename.vtk. This routine is a no-op for any process
 * that is not process 0.
 * @param data LevelBoxData that is defined on a single box assigned to process 0
 */


void WriteSinglePatchLevelData(LevelBoxData<double,DIM+NUMCOMPS>& out_data,
                               const double& dx,
                               const string& filename_data)
{
	if(procID()==0)
	{
		DataIterator dit=out_data.begin();

#if DIM == 1
		const char* varnames[5];
		varnames[0] = "X";
		varnames[1] = "density";
		varnames[2] = "Vx";
		varnames[3] = "p";
		varnames[4] = "Bx";
#endif

#if DIM == 2
		const char* varnames[8];
		varnames[0] = "X";
		varnames[1] = "Y";
		varnames[2] = "density";
		varnames[3] = "Vx";
		varnames[4] = "Vy";
		varnames[5] = "p";
		varnames[6] = "Bx";
		varnames[7] = "By";
#endif

#if DIM == 3
		const char* varnames[11];
		varnames[0] = "X";
		varnames[1] = "Y";
		varnames[2] = "Z";
		varnames[3] = "density";
		varnames[4] = "Vx";
		varnames[5] = "Vy";
		varnames[6] = "Vz";
		varnames[7] = "p";
		varnames[8] = "Bx";
		varnames[9] = "By";
		varnames[10] = "Bz";
#endif

		double origin[DIM];
		for (int ii = 0; ii < DIM; ii++)
		{
			origin[ii] = 0.0;
		}
		WriteBoxData(filename_data.c_str(),out_data[*dit],varnames,origin,dx);
	}


}



int main(int argc, char* argv[])
{

#ifdef PR_MPI
	MPI_Init(&argc,&argv);
#endif

	int pid = procID();
	bool convTest = false;
	int maxLev;

	// Defining inputs here until parmparse (or something similar becomes available)
	grid_type_global = 1;  // 0: 2D-Rectangular;  1: 2D-Wavy;  2: 2D-Polar
	double tstop = 0.5, CFL  = 0.1, domsize = 1.0, gamma = 1.66666666666666666666666667;
	int size1D = 64, maxStep = 10000, outputInterval = 10;
	limiter_apply = true;
	slope_flattening_apply = true;
	linear_visc_apply = true;
	non_linear_visc_apply = true;
	bool takedivBstep = true;
	Riemann_solver_type = 2; // 1:Rusanov; 2:Roe8Wave
	init_condition_type = 3;
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
	   10. Radially out flow
	   11. Cylindrical pulse problem in cartesian grid
	   12. Acoustic pulse problem in 3D cartesian grid
	   13. 3D MHD blast wave
	   14. 2D MHD blast wave
	 */
	if (grid_type_global == 2) {
		LowBoundType = 2;  // 0 for periodic, 1 for Dirichlet, 2 for open. This is for dir==0 only
		HighBoundType = 2;
	} else {
		LowBoundType = 0;
		HighBoundType = 0;
	}
	bool takeviscositystep = false;
	if (non_linear_visc_apply || linear_visc_apply) takeviscositystep = true;
	int domainSize=size1D;
	//int sizeDomain=domainSize;
	int sizeDomain=64;



	LevelBoxData<double,NUMCOMPS> U[3];
	LevelBoxData<double,DIM+NUMCOMPS> OUT[3];

	if (convTest) {
		maxLev = 3;
	} else {
		maxLev = 1;
	}

	for (int lev=0; lev<maxLev; lev++)
	{
		Box domain(Point::Zeros(),Point::Ones()*(domainSize -1));
		array<bool,DIM> per;
		for(int idir = 0; idir < DIM; idir++) per[idir]=true;
		ProblemDomain pd(domain,per);
		double dt, dx;
		dx = domsize/domainSize;
		//dx = params.domsize/domainSize;
		if (!convTest) dt = CFL*(1.0/domainSize);
		//if (convTest) dt = CFL*(1.0/1024.);
		if (convTest) dt = CFL*(1.0/domainSize);
		if(pid==0) std::cout << "domainSize: " << domainSize << std::endl;
		if(pid==0) std::cout << "dt: " << dt << std::endl;

		RK4<MHDLevelDataState,MHDLevelDataRK4Op,MHDLevelDataDX> rk4;
		EulerStep<MHDLevelDataState,MHDLevelDataRK4Op,MHDLevelDataDX> eulerstep;
		EulerStep<MHDLevelDataState, MHDLevelDatadivBOp, MHDLevelDataDX> divBstep;
		EulerStep<MHDLevelDataState, MHDLevelDataViscosityOp, MHDLevelDataDX> viscositystep;


		MHDLevelDataState state(pd,sizeDomain*Point::Ones(),dx,gamma);
		InitializeMHDLevelDataState(state);

		int count=0;
		for(DataIterator dit=state.m_U.begin(); *dit!=dit.end(); ++dit)
		{
			count++;
		}
		std::cout << "proc_id, num boxes " << pid << ", " << count << std::endl;

		double time = 0.;
		if(pid==0) cout << "starting time loop, maxStep = "<< maxStep << endl;
		for (int k = 0; (k < maxStep) && (time < tstop); k++)
		{
			if (!convTest)
			{
				if((outputInterval > 0) && (k == 0))
				{
					LevelBoxData<double,NUMCOMPS> new_state(state.m_dbl,Point::Ones(NGHOST));
					LevelBoxData<double,DIM> phys_coords(state.m_dbl,Point::Ones(NGHOST));
					LevelBoxData<double,NUMCOMPS+DIM> out_data(state.m_dbl,Point::Ones(NGHOST));
					for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
						//W_bar itself is not 4th order W. But it is calculated from 4th order accurate U for output.
						MHD_Mapping::JU_to_W_bar_calc(new_state[*dit],state.m_U[*dit],dx,gamma);
						MHD_Mapping::phys_coords_calc(phys_coords[*dit],state.m_U[*dit].box(),dx);
						MHD_Mapping::out_data_calc(out_data[*dit],phys_coords[*dit],new_state[*dit]);
						//MHD_Mapping::out_data_calc(out_data[*dit],phys_coords[*dit],state.m_U[*dit]);
					}

					//Solution on a single patch
					U[lev].define(DisjointBoxLayout(pd,domainSize*Point::Ones()),Point::Zeros());
					(new_state).copyTo(U[lev]);


					if (grid_type_global == 2) {
						OUT[lev].define(DisjointBoxLayout(pd,domainSize*Point::Ones()),{{0,1}});
					} else {
						OUT[lev].define(DisjointBoxLayout(pd,domainSize*Point::Ones()),Point::Zeros());
					}
					(out_data).copyTo(OUT[lev]);
					OUT[lev].exchange();
					std::string filename_data="Output_"+std::to_string(k);
					WriteSinglePatchLevelData(OUT[lev], dx, filename_data);
					if(pid==0) cout << "Written .vtk file after step "<< k << endl;
				}
			}

			auto start = chrono::steady_clock::now();
			if (convTest) {
				rk4.advance(time,dt,state);
			} else {
				rk4.advance(time,dt,state);
			}
			if (takeviscositystep) {
				// Take step for artificial viscosity
				viscositystep.advance(time,dt,state);
			}
			if (takedivBstep) {
				// Take step for divB term
				divBstep.advance(time,dt,state);
			}
			auto end = chrono::steady_clock::now();
			time += dt;
			time_globalll = time;
			if(pid==0) cout <<"nstep = " << k+1 << " time = " << time << " time step = " << dt << " Time taken: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms"  << endl;
			if (!convTest)
			{
				//if(pid==0) cout <<"nstep = " << k << " time = " << time << " time step = " << dt << endl;
				if((outputInterval > 0) && ((k+1)%outputInterval == 0))
				{
					LevelBoxData<double,NUMCOMPS> new_state(state.m_dbl,Point::Ones(NGHOST));
					LevelBoxData<double,DIM> phys_coords(state.m_dbl,Point::Ones(NGHOST));
					LevelBoxData<double,NUMCOMPS+DIM> out_data(state.m_dbl,Point::Ones(NGHOST));
					for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
						//W_bar itself is not 4th order W. But it is calculated from 4th order accurate U for output.
						//JU_to_W_calc is not suitable here as m_U doesn't have ghost cells, and deconvolve doesn't work at boundaries.
						MHD_Mapping::JU_to_W_bar_calc(new_state[*dit],state.m_U[*dit],dx,gamma);
						MHD_Mapping::phys_coords_calc(phys_coords[*dit],state.m_U[*dit].box(),dx);
						MHD_Mapping::out_data_calc(out_data[*dit],phys_coords[*dit],new_state[*dit]);
						//MHD_Mapping::out_data_calc(out_data[*dit],phys_coords[*dit],state.m_U[*dit]);
					}
					//Solution on a single patch
					U[lev].define(DisjointBoxLayout(pd,domainSize*Point::Ones()),Point::Zeros());
					(new_state).copyTo(U[lev]);

					if (grid_type_global == 2) {
						OUT[lev].define(DisjointBoxLayout(pd,domainSize*Point::Ones()),{{0,1}});
					} else {
						OUT[lev].define(DisjointBoxLayout(pd,domainSize*Point::Ones()),Point::Zeros());
					}
					(out_data).copyTo(OUT[lev]);
					OUT[lev].exchange();
					std::string filename_data="Output_"+std::to_string(k+1);
					WriteSinglePatchLevelData(OUT[lev], dx, filename_data);

					if(pid==0) cout << "Written .vtk file after step "<< k+1 << endl;
				}
			}
		}



		if (convTest) {
			//Solution on a single patch
			U[lev].define(DisjointBoxLayout(pd,domainSize*Point::Ones()),Point::Zeros());
			(state.m_U).copyTo(U[lev]);
			domainSize *= 2;
			sizeDomain *= 2; //For debugging: if you want to keep the number of boxes the same
			maxStep *= 2;
		}


	}

	//Until we have a coarsening operation for LevelBoxData, we perform the error calculations on a single patch.
	if(pid==0 && convTest)
	{
		for (int varr = 0; varr < 6; varr++) {
			Reduction<double> rxn;
			double ErrMax[2];
			for(int ilev=0; ilev<2; ilev++)
			{
				rxn.reset();
				DataIterator dit_lev=U[ilev].begin();
				DataIterator dit_levp1=U[ilev+1].begin();

				BoxData<double,1> err=slice(U[ilev][*dit_lev],varr);
				err-=Stencil<double>::AvgDown(2)(slice(U[ilev+1][*dit_levp1],varr));
				err.absMax(rxn);
				ErrMax[ilev]=rxn.fetch();
				//ErrMax[ilev]=err.absMax();
				std::string filename="rho_err"+std::to_string(ilev);
				if (varr == 0) filename="rho_err"+std::to_string(ilev);
				if (varr == 1) filename="Vx_err"+std::to_string(ilev);
				if (varr == 2) filename="Vy_err"+std::to_string(ilev);
				if (varr == 3) filename="p_err"+std::to_string(ilev);
				if (varr == 4) filename="Bx_err"+std::to_string(ilev);
				if (varr == 5) filename="By_err"+std::to_string(ilev);
				//NOTE: this assumes that the domain length is 1.0, which is assumed throughout this code. May cause errors if this changes.
				double dx=1./(err.box().size(0));
				WriteBoxData(filename.c_str(),err,dx);
				std::cout << "Lev: " << ilev << " , " << ErrMax[ilev] << std::endl;

			}
			double rate = log(abs(ErrMax[0]/ErrMax[1]))/log(2.0);
			std::cout << "order of accuracy = " << rate << std::endl;
		}
	}

#ifdef PR_MPI
	MPI_Finalize();
#endif
}
