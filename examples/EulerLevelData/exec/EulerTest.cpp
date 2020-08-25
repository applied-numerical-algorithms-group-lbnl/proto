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
  cout << "Navier Stokes simulation of shear flow with sinusoidal perturbation.  Periodic bcs." << endl;
  cout << "usage:  " << argv[0] << " -n nx  -t tmax -m maxstep -o output_interval" << endl;
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

int main(int argc, char* argv[])
{
#ifdef PR_MPI
    MPI_Init(&argc,&argv);
#endif

    double tstop;
    int size1D, maxStep, outputInterval;
    parseCommandLine(tstop, size1D, maxStep, outputInterval, argc, argv);
    EulerOp::s_gamma = 1.4;

    int domainSize=size1D;
    int sizeDomain=64;
    Box domain(Point::Zeros(),Point::Ones()*(domainSize -1));
    array<bool,DIM> per;
    for(int idir = 0; idir < DIM; idir++) per[idir]=true;
    double dx = 1.0/domainSize;
    EulerOp::s_dx=dx;
    ProblemDomain pd(domain,per);

    RK4<EulerLevelDataState,EulerLevelDataRK4Op,EulerLevelDataDX> rk4;
    EulerLevelDataState state(pd,sizeDomain*Point::Ones());
    InitializeEulerLevelDataState(state);
    double max_init_val=0.0;
    int count=0;
    for(DataIterator dit=state.m_U.begin(); *dit!=dit.end(); ++dit) {
        max_init_val=std::max(state.m_U[*dit].absMax(),max_init_val);
        count++;
    }
    std::cout << "Max init val: " << max_init_val << ", count: " << count << std::endl;
    double dt = .25/domainSize;
    double time = 0.;
    cout << "starting time loop, maxStep = "<< maxStep << endl;
    for (int k = 0;(k < maxStep) && (time < tstop);k++)
    {
        rk4.advance(time,dt,state);
        time += dt;
    }
    double max_val=0.0;
    for(DataIterator dit=state.m_U.begin(); *dit!=dit.end(); ++dit) {
        max_val=std::max(slice(state.m_U[*dit],0).absMax(),max_val);
    }
    std::cout << "Max val: " << max_val << std::endl;

#ifdef PR_MPI
    MPI_Finalize();
#endif
}
