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


// For Chrono Timer (Talwinder)
#include <chrono>
#include <iostream>
//////////////////////////////

using namespace std;
using namespace Proto;

double time_globalll;
int grid_type_global;


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
void WriteSinglePatchLevelData(LevelBoxData<double,NUMCOMPS>& data,
                               const double& dx,
                               const string& filename)
{
    if(procID()==0)
    {
        DataIterator dit=data.begin();

#if DIM == 1
        const char* varnames[4]; 		
		varnames[0] = "density";
		varnames[1] = "Mom_x";
		varnames[2] = "e";
		varnames[3] = "Bx";
#endif
		
#if DIM == 2	
        const char* varnames[6]; 	
		varnames[0] = "density";
		varnames[1] = "Mom_x";
		varnames[2] = "Mom_y";
		varnames[3] = "e";
		varnames[4] = "Bx";
		varnames[5] = "By";
#endif	

#if DIM == 3		
        const char* varnames[8]; 
		varnames[0] = "density";
		varnames[1] = "Mom_x";
		varnames[2] = "Mom_y";
		varnames[3] = "Mom_z";
		varnames[4] = "e";
		varnames[5] = "Bx";
		varnames[6] = "By";
		varnames[7] = "Bz";
#endif		

		double origin[DIM];
		for (int ii = 0; ii < DIM; ii++)
		{
			origin[ii] = 0.0;
		}
        WriteBoxData(filename.c_str(),data[*dit],varnames,origin,dx);
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
	grid_type_global = 0;  // 0: 2D-Rectangular;  1: 2D-Wavy;  2: 2D-Polar 
    double tstop = 0.5, CFL  = 0.125, domsize = 1.0, gamma = 1.6666666666666666666667;
    int size1D = 128, maxStep = 500, outputInterval = 10;
    //double gamma = params.gamma;
	//grid_type_global = params.grid_type;
	//maxStep = params.nstepmax;
	//tstop = params.tmax;
	//CFL  = params.cfl;
    //outputInterval = params.outinterv;
    //int domainSize=params.nx;
	int domainSize=size1D;
    int sizeDomain=domainSize;
	
	
	
    LevelBoxData<double,NUMCOMPS> U[3];
	
	if (convTest){
		maxLev = 3;
	} else {
		maxLev = 1;
	}
	
    for (int lev=0; lev<maxLev; lev++)
    {
        Box domain(Point::Zeros(),Point::Ones()*(domainSize -1));
        array<bool,DIM> per;
        for(int idir = 0; idir < DIM; idir++) per[idir]=true;
		if (grid_type_global == 2){
			per[0]=true;
			per[1]=true;
		}
        ProblemDomain pd(domain,per);
        double dt, dx;
        dx = domsize/domainSize;
        //dx = params.domsize/domainSize;
        if (!convTest) dt = CFL*(1.0/domainSize);
		if (convTest) dt = CFL*(1.0/1024.);
        if(pid==0) std::cout << "domainSize: " << domainSize << std::endl;
        if(pid==0) std::cout << "dt: " << dt << std::endl;

        RK4<MHDLevelDataState,MHDLevelDataRK4Op,MHDLevelDataDX> rk4;
        //EulerStep<MHDLevelDataState,MHDLevelDataRK4Op,MHDLevelDataDX> rk4;
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
        for (int k = 0;(k < maxStep) && (time < tstop);k++)
        {
			
			if (!convTest)
            {
              //if(pid==0) cout <<"nstep = " << k << " time = " << time << " time step = " << dt << endl;
              if((outputInterval > 0) && (k == 0))
                {
				  //Solution on a single patch
                  U[lev].define(DisjointBoxLayout(pd,domainSize*Point::Ones()),Point::Zeros());
				  (state.m_U).copyTo(U[lev]);
				  
				  std::string filename="Output_"+std::to_string(k);
				  WriteSinglePatchLevelData(U[lev], dx, filename);  
				  
				  if(pid==0) cout << "Written .vtk file after step "<< k << endl;
                }
            }
			
			
			auto start = chrono::steady_clock::now();
			rk4.advance(time,dt,state);
			auto end = chrono::steady_clock::now();
			//cout << "Elapsed time in milliseconds in main rk4: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;
			// Take step for artificial viscosity
			viscositystep.advance(time,dt,state);
			// Take step for divB term
			divBstep.advance(time,dt,state);
            time += dt;
			time_globalll = time;
			if(pid==0) cout <<"nstep = " << k+1 << " time = " << time << " time step = " << dt << endl;
			if (!convTest)
            {
              //if(pid==0) cout <<"nstep = " << k << " time = " << time << " time step = " << dt << endl;
              if((outputInterval > 0) && ((k+1)%outputInterval == 0))
                {
				  //Solution on a single patch
                  U[lev].define(DisjointBoxLayout(pd,domainSize*Point::Ones()),Point::Zeros());
				  (state.m_U).copyTo(U[lev]);
			  
				  std::string filename="Output_"+std::to_string(k+1);
				  WriteSinglePatchLevelData(U[lev], dx, filename);
		  
				  if(pid==0) cout << "Written .vtk file after step "<< k+1 << endl;
                }
            }
        }

       
        
        if (convTest){
			//Solution on a single patch
			U[lev].define(DisjointBoxLayout(pd,domainSize*Point::Ones()),Point::Zeros());
			(state.m_U).copyTo(U[lev]);
			//std::string filename="rho_"+std::to_string(lev);
			//WriteSinglePatchLevelData(U[lev], 0, dx, filename);
			domainSize *= 2;
			sizeDomain *= 2; //For debugging: if you want to keep the number of boxes the same
			//maxStep *= 2;
		}
		
		
    }

    //Until we have a coarsening operation for LevelBoxData, we perform the error calculations on a single patch.
    if(pid==0 && convTest)
    {
		for (int varr = 0; varr < 6; varr++){
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
