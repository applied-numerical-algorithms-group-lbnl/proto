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
#include "EulerLevelDataRK4.H"
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"

#define PI 3.141592653589793

using namespace std;
using namespace Proto;

typedef Var<double,DIM> V;
typedef Var<double,NUMCOMPS> State;

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

PROTO_KERNEL_START 
unsigned int InitializeStateF(State& a_U,
                             const V& a_x)
{
    double gamma = 1.4;
    double rho0 = gamma;
    double p0 = 1.;
    double umag = 0.;
    double rho = rho0;
    rho += .01*rho0*sin(2*2*PI*a_x(0));
    double p = p0*pow(rho/rho0,gamma);
    a_U(0) = rho;
    double c0 = sqrt(gamma*p0/rho0);
    double c = sqrt(gamma*p/rho);
    umag = 2*(c-c0)/(gamma-1.);
    double ke = 0.;
    // FIX
    for (int dir = 1; dir <= 1; dir++)
    {
        ke += umag*umag;
        a_U(dir) = rho*umag;
    }
    ke *= .5;
    a_U(NUMCOMPS-1) = p/(gamma-1.0) + rho*ke;
    return 0;
}
PROTO_KERNEL_END(InitializeStateF, InitializeState)

//=================================================================================================
PROTO_KERNEL_START
void iotaFuncF(Point           & a_p,
               V               & a_X,
               double            a_h)
{
  for (int ii = 0; ii < DIM; ii++)
  {
    a_X(ii) = a_p[ii]*a_h + 0.5*a_h;
  }
}
PROTO_KERNEL_END(iotaFuncF,iotaFunc)


void InitializeEulerLevelDataState(EulerLevelDataState& state)
{
    for(DataIterator dit=state.m_U.begin(); *dit!=dit.end(); ++dit) {
        BoxData<double,NUMCOMPS>& UState = state.m_U[*dit];
        Box dbx0=UState.box();
        Box dbx = dbx0.grow(NGHOST);
        Box dbx1 = dbx.grow(1);
        BoxData<double,NUMCOMPS> UBig(dbx1);
        BoxData<double,DIM> x(dbx1);
        forallInPlace_p(iotaFunc, dbx1, x, EulerOp::s_dx);
        forallInPlace(InitializeState,dbx1,UBig,x);
        Stencil<double> Lap2nd = Stencil<double>::Laplacian();
        UState |= Lap2nd(UBig,dbx,1.0/24.0); 
        UState += UBig;
    }
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
        BoxData<double> rho=slice(data[*dit],comp);
        WriteBoxData(filename.c_str(),rho,dx);
    }
}

int main(int argc, char* argv[])
{
#ifdef PR_MPI
    MPI_Init(&argc,&argv);
#endif

    int pid=procID();

    if(pid==0) {
        cout << "Navier Stokes simulation of shear flow with sinusoidal perturbation.  Periodic bcs." << endl;
        cout << "usage:  " << argv[0] << " -n nx  -t tmax -m maxstep -o output_interval" << endl;
    }

    double tstop;
    int size1D, maxStep, outputInterval;
    parseCommandLine(tstop, size1D, maxStep, outputInterval, argc, argv);
    EulerOp::s_gamma = 1.4;

    int domainSize=size1D;
    int sizeDomain=64;
    //int sizeDomain=size1D;
    LevelBoxData<double,NUMCOMPS> U[3];
    for (int lev=0; lev<3; lev++)
    {
        Box domain(Point::Zeros(),Point::Ones()*(domainSize -1));
        array<bool,DIM> per;
        for(int idir = 0; idir < DIM; idir++) per[idir]=true;
        ProblemDomain pd(domain,per);

        double dx = 1.0/domainSize;
        EulerOp::s_dx=dx;
        double dt = .25/domainSize;
        std::cout << "domainSize: " << domainSize << std::endl;
        std::cout << "dt: " << dt << std::endl;

        RK4<EulerLevelDataState,EulerLevelDataRK4Op,EulerLevelDataDX> rk4;
        EulerLevelDataState state(pd,sizeDomain*Point::Ones());
        InitializeEulerLevelDataState(state);

        int count=0;
        for(DataIterator dit=state.m_U.begin(); *dit!=dit.end(); ++dit)
        {
            count++;
        }
        std::cout << "proc_id, num boxes " << pid << ", " << count << std::endl;
/*
  double max_init_val=0.0;
  int count=0;
  for(DataIterator dit=state.m_U.begin(); *dit!=dit.end(); ++dit) {
  max_init_val=std::max(state.m_U[*dit].absMax(),max_init_val);
  count++;
  }
  std::cout << "Max init val: " << max_init_val << ", count: " << count << std::endl;
*/

        double time = 0.;
        if(pid==0)
            cout << "starting time loop, maxStep = "<< maxStep << endl;
        //maxStep=0;
        for (int k = 0;(k < maxStep) && (time < tstop);k++)
        {
            rk4.advance(time,dt,state);
            time += dt;
        }

        //Solution on a single patch
        U[lev].define(DisjointBoxLayout(pd,domainSize*Point::Ones()),Point::Zeros());
        (state.m_U).copyTo(U[lev]);
        std::string filename="rho_"+std::to_string(lev);
        WriteSinglePatchLevelData(U[lev], 0, dx, filename);

        /*
          double max_val=0.0;
          for(DataIterator dit=state.m_U.begin(); *dit!=dit.end(); ++dit) {
          max_val=std::max(slice(state.m_U[*dit],0).absMax(),max_val);
          }
          std::cout << "Max val: " << max_val << std::endl;
        */

        domainSize *= 2;
        //sizeDomain *= 2;
        maxStep *= 2;
    }

    //Until we have a coarsening operation for LevelBoxData, we perform the error calculations on a single patch.
    if(pid==0)
    {
        double rhoErrMax[2];
        for(int ilev=0; ilev<2; ilev++)
        {
            DataIterator dit_lev=U[ilev].begin();
            DataIterator dit_levp1=U[ilev+1].begin();
            BoxData<double,1> err=slice(U[ilev][*dit_lev],0);
            err-=Stencil<double>::AvgDown(2)(slice(U[ilev+1][*dit_levp1],0));
            rhoErrMax[ilev]=err.absMax();
            std::string filename="rho_err"+std::to_string(ilev);
//NOTE: this assumes that the domain length is 1.0, which is assumed throughout this code. May cause errors if this changes.
            double dx=1./(err.box().size(0));
            WriteBoxData(filename.c_str(),err,dx);
            std::cout << "Lev: " << ilev << " , " << rhoErrMax[ilev] << std::endl;
        }
        double rate = log(abs(rhoErrMax[0]/rhoErrMax[1]))/log(2.0);
        std::cout << "order of accuracy = " << rate << std::endl;
    }

#ifdef PR_MPI
    MPI_Finalize();
#endif
}
