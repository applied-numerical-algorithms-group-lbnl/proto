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
#include <iomanip>
#include "Proto.H"
#include "MHDLevelDataRK4.H"
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"
#include "MHD_Initialize.H"
#include "MHDOp.H"
#include "MHD_Mapping.H"
#include "MHD_Output_Writer.H"
#include "MHD_Input_Parsing.H"
#include "MHD_CFL.H"
#include <chrono> // Used by timer
using namespace std;
using namespace Proto;
Parsefrominputs inputs;
int main(int argc, char* argv[])
{
#ifdef PR_MPI
	MPI_Init(&argc,&argv);
#endif
	//have to do this to get a time table
	PR_TIMER_SETFILE("proto.time.table");
	//PR_TIME("main");
	int pid = procID();
	//Reading inputs file
	inputs.parsenow();
	int maxLev;	
	// When using mapping, computational domain is always from 0 to 1. The physical grid is mapped brom this cube.
	if (inputs.grid_type_global > 1){
		inputs.domsizex = 1.0;
		inputs.domsizey = 1.0;
		inputs.domsizez = 1.0;
	}
	bool takeviscositystep = false;
	if ((inputs.non_linear_visc_apply == 1) || (inputs.linear_visc_apply == 1)) takeviscositystep = true;

	LevelBoxData<double,NUMCOMPS> U[3];  // Size 3 is needed for the convergence rate tests (If indicated in inputs file)
	if (inputs.convTestType != 0) {
		maxLev = 3;
	} else {
		maxLev = 1;
	}
	for (int lev=0; lev<maxLev; lev++)
	{
		// Creating a box for our full domain 
		#if DIM == 1
			Box domain(Point::Zeros(),Point(inputs.domainSizex-1));
		#endif
		#if DIM == 2
			Box domain(Point::Zeros(),Point(inputs.domainSizex-1, inputs.domainSizey-1));
		#endif
		#if DIM == 3
			Box domain(Point::Zeros(),Point(inputs.domainSizex-1, inputs.domainSizey-1, inputs.domainSizez-1));
		#endif
		array<bool,DIM> per;
		// All outer boundaries are set to periodic by default
		for(int idir = 0; idir < DIM; idir++){
			per[idir]=true;
	 	}
		// Creating problem domain 
		ProblemDomain pd(domain,per);
		double dx = inputs.domsizex/inputs.domainSizex, dy = inputs.domsizey/inputs.domainSizey, dz = inputs.domsizez/inputs.domainSizez;
		double dt;
		// Following is done for required dt control in convergence tests
		if (inputs.convTestType == 1)
		{
			dt = inputs.CFL*(1.0/1024.);
		} else {
			#if DIM == 1
			dt = inputs.CFL*dx;
			#endif
			#if DIM == 2
			dt = inputs.CFL*std::min({dx,dy});
			#endif
			#if DIM == 3
			dt = inputs.CFL*std::min({dx,dy,dz});
			#endif
		}

		// Create an object state. state.m_U has all the consereved variables (multiplied by Jacobian for mapped grids)
		// All the mapping variables, which are functions of mapping geometry are also included in this class object.
		MHDLevelDataState state(pd,inputs.BoxSize*Point::Ones(),dx, dy, dz, inputs.gamma);
		(state.m_U).setToZero();  

		// This is used to find number of boxes in each processor.
		int count=0;
		for(DataIterator dit=state.m_U.begin(); *dit!=dit.end(); ++dit)
		{
			count++;
		}
		std::cout << "proc_id: " << pid << ";      num boxes: " << count << std::endl;

		

		for(DataIterator dit=(state.m_Jacobian_ave).begin(); *dit!=dit.end(); ++dit) {
			if (inputs.grid_type_global != 2){
				MHD_Mapping::Jacobian_Ave_calc((state.m_Jacobian_ave)[*dit],dx,dy,dz,state.m_U[*dit].box().grow(1));
			}
			MHD_Mapping::N_ave_f_calc_func((state.m_N_ave_f)[*dit],dx, dy, dz);
			// MHD_Output_Writer::WriteBoxData_Nsd_nocoord(state.m_N_ave_f[*dit], dx, dy, dz, "N_ave_sd_f_numeric");
		}
		if (inputs.grid_type_global != 2){
			(state.m_Jacobian_ave).exchange();
		}
		#if DIM == 3
		for(DataIterator dit=(state.m_detAA_avg).begin(); *dit!=dit.end(); ++dit) {
			MHD_Mapping::Spherical_map_calc_func((state.m_Jacobian_ave)[*dit], (state.m_A_1_avg)[*dit], (state.m_A_2_avg)[*dit], (state.m_A_3_avg)[*dit], (state.m_detAA_avg)[*dit], (state.m_detAA_inv_avg)[*dit], (state.m_r2rdot_avg)[*dit], (state.m_detA_avg)[*dit], (state.m_r2detA_1_avg)[*dit], (state.m_r2detAA_1_avg)[*dit], (state.m_r2detAn_1_avg)[*dit], (state.m_rrdotdetA_2_avg)[*dit], (state.m_rrdotdetAA_2_avg)[*dit], (state.m_rrdotd3ncn_2_avg)[*dit], (state.m_rrdotdetA_3_avg)[*dit], (state.m_rrdotdetAA_3_avg)[*dit], (state.m_rrdotncd2n_3_avg)[*dit],dx, dy, dz);
		}	
		#endif
		if (inputs.grid_type_global == 2 && (inputs.initialize_in_spherical_coords == 1)){
			for(DataIterator dit=state.m_U.begin(); *dit!=dit.end(); ++dit){
				MHD_Initialize::initializeState_Spherical((state.m_U)[*dit], (state.m_detAA_avg)[*dit], (state.m_detAA_inv_avg)[*dit], (state.m_r2rdot_avg)[*dit], (state.m_detA_avg)[*dit], state.m_dx, state.m_dy, state.m_dz,state.m_gamma);
			}
		} else {
			for(DataIterator dit=state.m_U.begin(); *dit!=dit.end(); ++dit){
				MHD_Initialize::initializeState((state.m_U)[*dit] ,state.m_dx, state.m_dy, state.m_dz, state.m_gamma);
			}
		}
		LevelBoxData<double,DIM+NUMCOMPS> OUT[3];
		// LevelBoxData<double,NUMCOMPS> OUT[3];
		if (inputs.grid_type_global == 2) {
			#if DIM == 1
			OUT[lev].define(DisjointBoxLayout(pd,Point(inputs.domainSizex)),Point::Zeros());
			#endif
			#if DIM == 2
			OUT[lev].define(DisjointBoxLayout(pd,Point(inputs.domainSizex, inputs.domainSizey)),{{0,1}});
			#endif
			#if DIM == 3
			OUT[lev].define(DisjointBoxLayout(pd,Point(inputs.domainSizex, inputs.domainSizey, inputs.domainSizez)), {{0,0,1}});
			#endif
		} else {
			#if DIM == 1
			OUT[lev].define(DisjointBoxLayout(pd,inputs.domainSizex*Point::Ones()),Point::Zeros());
			#endif
			#if DIM == 2
			OUT[lev].define(DisjointBoxLayout(pd,Point(inputs.domainSizex, inputs.domainSizey)),Point::Zeros());
			#endif
			#if DIM == 3
			OUT[lev].define(DisjointBoxLayout(pd,Point(inputs.domainSizex, inputs.domainSizey,  inputs.domainSizez)), Point::Zeros());
			#endif
		}
		LevelBoxData<double,NUMCOMPS> new_state(state.m_dbl,Point::Ones(NGHOST));
		LevelBoxData<double,DIM> phys_coords(state.m_dbl,Point::Ones(NGHOST));
		LevelBoxData<double,NUMCOMPS+DIM> out_data(state.m_dbl,Point::Ones(NGHOST));
		double dt_new = 0.;
		double time = 0.;
		if(pid==0) cout << "starting time loop, maxStep = "<< inputs.maxStep << endl;
		for (int k = 0; (k <= inputs.maxStep) && (time < inputs.tstop); k++)
		{		
			if (k!=0){
				if (inputs.convTestType == 0){
					double dt_temp = 1.0e10;
					for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
						if (inputs.grid_type_global == 2){
							MHD_Mapping::JU_to_W_Sph_ave_calc_func(new_state[*dit], state.m_U[*dit], (state.m_detAA_inv_avg)[*dit], (state.m_r2rdot_avg)[*dit], (state.m_detA_avg)[*dit], inputs.gamma);
						} else {
							MHD_Mapping::JU_to_W_bar_calc(new_state[*dit],state.m_U[*dit],(state.m_detAA_inv_avg)[*dit], (state.m_r2rdot_avg)[*dit], (state.m_detA_avg)[*dit],dx,dy,dz,inputs.gamma);
						}
						MHD_CFL::Min_dt_calc_func(dt_new, new_state[*dit], dx, dy, dz, inputs.gamma);
						if (dt_new < dt_temp) dt_temp = dt_new;
					}
					double mintime;
					#ifdef PR_MPI
						MPI_Reduce(&dt_temp, &mintime, 1, MPI_DOUBLE, MPI_MIN, 0,MPI_COMM_WORLD);
						MPI_Bcast(&mintime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					#endif
					dt = inputs.CFL*mintime;
					if ((inputs.tstop - time) < dt) dt = inputs.tstop - time;
				}
				// Below objects need to be created inside the time loop. Otherwise the init in dx keeps on eating memory
				// This is used to take rk4 step
				RK4<MHDLevelDataState,MHDLevelDataRK4Op,MHDLevelDataDX> rk4;
				// This will be used to take Euler step (Primarily used in convergence tests and debugging)
				EulerStep<MHDLevelDataState,MHDLevelDataRK4Op,MHDLevelDataDX> eulerstep;
				// Both Powell divergence cleaning and viscosity implementation need Euler steps at each state update
				EulerStep<MHDLevelDataState, MHDLevelDatadivBOp, MHDLevelDataDX> divBstep;
				EulerStep<MHDLevelDataState, MHDLevelDataViscosityOp, MHDLevelDataDX> viscositystep;
				auto start = chrono::steady_clock::now();
				if (inputs.convTestType == 1 || inputs.timeIntegratorType == 1) {
					eulerstep.advance(time,dt,state);
				} else {
					if (inputs.timeIntegratorType == 4){
						rk4.advance(time,dt,state);
					}
				}
				if (takeviscositystep) {
					// Take step for artificial viscosity
					viscositystep.advance(time,dt,state);
				}
				if (inputs.takedivBstep == 1) {
					// Take step for divB term
					divBstep.advance(time,dt,state);
				}
				auto end = chrono::steady_clock::now();
				time += dt;
				if(pid==0) cout <<"nstep = " << k << " time = " << time << " time step = " << dt << " Time taken: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms"  << endl;
			}
			
			if (inputs.convTestType == 0)
			{
				if(((inputs.outputInterval > 0) && ((k)%inputs.outputInterval == 0)) || time == inputs.tstop || ((inputs.outputInterval > 0) && (k == 0)))
				{
					for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {	
						if (inputs.grid_type_global == 2){
							MHD_Mapping::JU_to_W_Sph_ave_calc_func(new_state[*dit], state.m_U[*dit], (state.m_detAA_inv_avg)[*dit], (state.m_r2rdot_avg)[*dit], (state.m_detA_avg)[*dit], inputs.gamma);
						} else {
							//W_bar itself is not 4th order W. But it is calculated from 4th order accurate U for output.
							//JU_to_W_calc is not suitable here as m_U doesn't have ghost cells, and deconvolve doesn't work at boundaries.
						    MHD_Mapping::JU_to_W_bar_calc(new_state[*dit],state.m_U[*dit],(state.m_detAA_inv_avg)[*dit], (state.m_r2rdot_avg)[*dit], (state.m_detA_avg)[*dit],dx,dy,dz,inputs.gamma);
						}
						MHD_Mapping::phys_coords_calc(phys_coords[*dit],state.m_U[*dit].box(),dx,dy,dz);
						MHD_Mapping::out_data_calc(out_data[*dit],phys_coords[*dit],new_state[*dit]);
						// MHD_Mapping::out_data_calc(out_data[*dit],phys_coords[*dit],state.m_U[*dit]);
					}
					//Solution on a single patch
					(out_data).copyTo(OUT[lev]);
					// (state.m_U).copyTo(OUT[lev]);
					OUT[lev].exchange();
					std::string filename_data="Output_"+std::to_string(k);
					MHD_Output_Writer::WriteSinglePatchLevelData(OUT[lev], dx,dy,dz, filename_data);
					// MHD_Output_Writer::WriteSinglePatchLevelData_nocoord(OUT[lev], dx,dy,dz, filename_data);
					if(pid==0) cout << "Written data file after step "<< k << endl;		
				}
			}
		}

		if (inputs.convTestType != 0) {
			//Solution on a single patch
			#if DIM == 1
			U[lev].define(DisjointBoxLayout(pd,inputs.domainSizex*Point::Ones()),Point::Zeros());
			#endif
			#if DIM == 2
			U[lev].define(DisjointBoxLayout(pd,Point(inputs.domainSizex, inputs.domainSizey)),Point::Zeros());
			#endif
			#if DIM == 3
			U[lev].define(DisjointBoxLayout(pd,Point(inputs.domainSizex, inputs.domainSizey, inputs.domainSizez)), Point::Zeros());
			#endif
			(state.m_U).copyTo(U[lev]);
			inputs.domainSizex *= 2;
			inputs.domainSizey *= 2;
			inputs.domainSizez *= 2;
			inputs.BoxSize *= 2; //For debugging: if you want to keep the number of boxes the same
			if (inputs.convTestType == 2){
				inputs.maxStep *= 2;
			}
		}
	}
	//Until we have a coarsening operation for LevelBoxData, we perform the error calculations on a single patch.
	if(pid==0 && (inputs.convTestType != 0))
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
				std::string filename="Comp_"+std::to_string(varr)+"_err_"+std::to_string(ilev);
				//NOTE: this assumes that the domain length is 1.0, which is assumed throughout this code. May cause errors if this changes.
				double dx=1./(err.box().size(0));
				if (inputs.saveConvTestData){
					WriteBoxData(filename.c_str(),err,dx);
				}
				std::cout << "Lev: " << ilev << " , " << ErrMax[ilev] << std::endl;
			}
			double rate = log(abs(ErrMax[0]/ErrMax[1]))/log(2.0);
			std::cout << "order of accuracy = " << rate << std::endl;
		}
	}
	PR_TIMER_REPORT();
#ifdef PR_MPI
	MPI_Finalize();
#endif
}
