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

// For Chrono Timer (Talwinder)
#include <chrono>
#include <iostream>
//////////////////////////////

using namespace std;
using namespace Proto;

void
parseCommandLine(double& a_tmax, int& a_nx, int& a_maxstep, int& a_outputinterval, int argc, char* argv[])
{
  a_tmax= 1.0;
  a_maxstep = 1;
  a_outputinterval = -1;
  a_nx = 128;
  for(int iarg = 0; iarg < argc-1; iarg++)
  {
    if(strcmp(argv[iarg],"-n") == 0)
    {
      a_nx = atoi(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg], "-m") == 0)
    {
      a_maxstep = atoi(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg], "-o") == 0)
    {
      a_outputinterval = atoi(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg],"-t") == 0)
    {
      a_tmax = atof(argv[iarg+1]);
    }
  }
}



void InitializeMHDLevelDataState(MHDLevelDataState& state)
{
    (state.m_U).setToZero();
    for(DataIterator dit=state.m_U.begin(); *dit!=dit.end(); ++dit)
        MHDOp::initializeState((state.m_U)[*dit],state.m_dx,state.m_gamma);
}

/**
 * @brief Write out component comp of data to filename.vtk. This routine is a no-op for any process
 * that is not process 0.
 * @param data LevelBoxData that is defined on a single box assigned to process 0
 * @param comp integer in the range [0,NUMCOMPS-1]
 */
void WriteSinglePatchLevelData(LevelBoxData<double,NUMCOMPS>& data,
                               const int comp,
                               const double& dx,
                               const string& filename)
{
    if(procID()==0)
    {
        DataIterator dit=data.begin();
        BoxData<double> data_slice = slice(data[*dit],comp);
		
		const char* varnames[0];
#if DIM == 1		
		if (comp == 0) varnames[0] = "density";
		if (comp == 1) varnames[0] = "Vx";
		if (comp == 2) varnames[0] = "p";
		if (comp == 3) varnames[0] = "Bx";
#endif
		
#if DIM == 2		
		if (comp == 0) varnames[0] = "density";
		if (comp == 1) varnames[0] = "Vx";
		if (comp == 2) varnames[0] = "Vy";
		if (comp == 3) varnames[0] = "p";
		if (comp == 4) varnames[0] = "Bx";
		if (comp == 5) varnames[0] = "By";
#endif	

#if DIM == 3		
		if (comp == 0) varnames[0] = "density";
		if (comp == 1) varnames[0] = "Vx";
		if (comp == 2) varnames[0] = "Vy";
		if (comp == 3) varnames[0] = "Vz";
		if (comp == 4) varnames[0] = "p";
		if (comp == 5) varnames[0] = "Bx";
		if (comp == 6) varnames[0] = "By";
		if (comp == 7) varnames[0] = "Bz";
#endif		
		double origin[DIM];
		for (int ii = 0; ii < DIM; ii++)
		{
			origin[ii] = 0.0;
		}
        WriteBoxData(filename.c_str(),data_slice,varnames,origin,dx);
    }
}

int main(int argc, char* argv[])
{
	//cerr << "here"<<endl;
#ifdef PR_MPI
    MPI_Init(&argc,&argv);
#endif

    int pid = procID();
    bool convTest = false;
	int maxLev;
	
    if(pid == 0) {
        cout << "MHD simulation of Orszag Tang vortex.  Periodic bcs." << endl;
        cout << "usage:  " << argv[0] << " -n nx  -t tmax -m maxstep -o output_interval" << endl;
    }

    double tstop;
    int size1D, maxStep, outputInterval;
    parseCommandLine(tstop, size1D, maxStep, outputInterval, argc, argv);
    double gamma = 1.6666666666666666666667;
    //double gamma = 1.4;
	

    int domainSize=size1D;
    //int sizeDomain=64;
    int sizeDomain=size1D;
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
        ProblemDomain pd(domain,per);
        double dt, dx;
        dx = 1.0/domainSize;
        if (!convTest) dt = (.25/domainSize);
		if (convTest) dt = (.25/1024.);
        if(pid==0) std::cout << "domainSize: " << domainSize << std::endl;
        if(pid==0) std::cout << "dt: " << dt << std::endl;

        //RK4<MHDLevelDataState,MHDLevelDataRK4Op,MHDLevelDataDX> rk4;
        EulerStep<MHDLevelDataState,MHDLevelDataRK4Op,MHDLevelDataDX> rk4;
		EulerStep<MHDLevelDataState, MHDLevelDataEulerOp, MHDLevelDataDX> eulerstep;
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
#if DIM == 1				  
				  std::string filename="rho_"+std::to_string(k);
				  WriteSinglePatchLevelData(U[lev], 0, dx, filename);
				  filename="Vx_"+std::to_string(k);
				  WriteSinglePatchLevelData(U[lev], 1, dx, filename);
				  filename="p_"+std::to_string(k);
				  WriteSinglePatchLevelData(U[lev], 2, dx, filename);
				  filename="Bx_"+std::to_string(k);
				  WriteSinglePatchLevelData(U[lev], 3, dx, filename);
#endif
				  
#if DIM == 2				  
				  std::string filename="rho_"+std::to_string(k);
				  WriteSinglePatchLevelData(U[lev], 0, dx, filename);
				  filename="Vx_"+std::to_string(k);
				  WriteSinglePatchLevelData(U[lev], 1, dx, filename);
				  filename="Vy_"+std::to_string(k);
				  WriteSinglePatchLevelData(U[lev], 2, dx, filename);
				  filename="p_"+std::to_string(k);
				  WriteSinglePatchLevelData(U[lev], 3, dx, filename);
				  filename="Bx_"+std::to_string(k);
				  WriteSinglePatchLevelData(U[lev], 4, dx, filename);
				  filename="By_"+std::to_string(k);
				  WriteSinglePatchLevelData(U[lev], 5, dx, filename);
#endif

#if DIM == 3				  
				  std::string filename="rho_"+std::to_string(k);
				  WriteSinglePatchLevelData(U[lev], 0, dx, filename);
				  filename="Vx_"+std::to_string(k);
				  WriteSinglePatchLevelData(U[lev], 1, dx, filename);
				  filename="Vy_"+std::to_string(k);
				  WriteSinglePatchLevelData(U[lev], 2, dx, filename);
				  filename="Vz_"+std::to_string(k);
				  WriteSinglePatchLevelData(U[lev], 3, dx, filename);
				  filename="p_"+std::to_string(k);
				  WriteSinglePatchLevelData(U[lev], 4, dx, filename);
				  filename="Bx_"+std::to_string(k);
				  WriteSinglePatchLevelData(U[lev], 5, dx, filename);
				  filename="By_"+std::to_string(k);
				  WriteSinglePatchLevelData(U[lev], 6, dx, filename);
				  filename="Bz_"+std::to_string(k);
				  WriteSinglePatchLevelData(U[lev], 7, dx, filename);
#endif					  
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
			eulerstep.advance(time,dt,state);
            time += dt;
			if(pid==0) cout <<"nstep = " << k+1 << " time = " << time << " time step = " << dt << endl;
			if (!convTest)
            {
              //if(pid==0) cout <<"nstep = " << k << " time = " << time << " time step = " << dt << endl;
              if((outputInterval > 0) && ((k+1)%outputInterval == 0))
                {
				  //Solution on a single patch
                  U[lev].define(DisjointBoxLayout(pd,domainSize*Point::Ones()),Point::Zeros());
				  (state.m_U).copyTo(U[lev]);
				  
#if DIM == 1				  
				  std::string filename="rho_"+std::to_string(k+1);
				  WriteSinglePatchLevelData(U[lev], 0, dx, filename);
				  filename="Vx_"+std::to_string(k+1);
				  WriteSinglePatchLevelData(U[lev], 1, dx, filename);
				  filename="p_"+std::to_string(k+1);
				  WriteSinglePatchLevelData(U[lev], 2, dx, filename);
				  filename="Bx_"+std::to_string(k+1);
				  WriteSinglePatchLevelData(U[lev], 3, dx, filename);
#endif
				  
#if DIM == 2				  
				  std::string filename="rho_"+std::to_string(k+1);
				  WriteSinglePatchLevelData(U[lev], 0, dx, filename);
				  filename="Vx_"+std::to_string(k+1);
				  WriteSinglePatchLevelData(U[lev], 1, dx, filename);
				  filename="Vy_"+std::to_string(k+1);
				  WriteSinglePatchLevelData(U[lev], 2, dx, filename);
				  filename="p_"+std::to_string(k+1);
				  WriteSinglePatchLevelData(U[lev], 3, dx, filename);
				  filename="Bx_"+std::to_string(k+1);
				  WriteSinglePatchLevelData(U[lev], 4, dx, filename);
				  filename="By_"+std::to_string(k+1);
				  WriteSinglePatchLevelData(U[lev], 5, dx, filename);
#endif

#if DIM == 3				  
				  std::string filename="rho_"+std::to_string(k+1);
				  WriteSinglePatchLevelData(U[lev], 0, dx, filename);
				  filename="Vx_"+std::to_string(k+1);
				  WriteSinglePatchLevelData(U[lev], 1, dx, filename);
				  filename="Vy_"+std::to_string(k+1);
				  WriteSinglePatchLevelData(U[lev], 2, dx, filename);
				  filename="Vz_"+std::to_string(k+1);
				  WriteSinglePatchLevelData(U[lev], 3, dx, filename);
				  filename="p_"+std::to_string(k+1);
				  WriteSinglePatchLevelData(U[lev], 4, dx, filename);
				  filename="Bx_"+std::to_string(k+1);
				  WriteSinglePatchLevelData(U[lev], 5, dx, filename);
				  filename="By_"+std::to_string(k+1);
				  WriteSinglePatchLevelData(U[lev], 6, dx, filename);
				  filename="Bz_"+std::to_string(k+1);
				  WriteSinglePatchLevelData(U[lev], 7, dx, filename);
#endif				  
				  
				  
				  
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
