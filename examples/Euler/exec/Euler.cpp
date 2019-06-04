#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>

#include "Proto.H"
#include "EulerRK4.H"
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"

#define PI 3.141592653589793
#define NUMCOMPS DIM+2
#define NUMCELLS 32

using namespace std;
using namespace Proto;

typedef Var<double,DIM> V;
typedef Var<double,NUMCOMPS> State;

void WriteData(int                              a_iter,
               BoxData<double,NUMCOMPS>&    a_state,
               double                           a_dx)
{
    char basename[1024];
    sprintf(basename,"euler.%06d",a_iter);

    const char* varnames[DIM+2];
    varnames[0] = "density";

    varnames[1] = "momentum-x";
    #if DIM > 1
    varnames[2] = "momentum-y";
    #endif
    #if DIM > 2
    varnames[3] = "momentum-z";
    #endif

    varnames[DIM+1] = "energy";
    double origin[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        origin[ii] = 0.0;
    }
    WriteBoxData(basename,a_state,varnames,origin,a_dx);
};  

PROTO_KERNEL_START 
unsigned int InitializeStateF(State& a_U,
                             const V& a_x)
{
    double gamma = 1.4;
    double rho0 = gamma;
    double p0 = 1.;
    double umag = .5/sqrt(1.*(DIM));
    double rho = rho0;
    for (int dir = 0; dir < DIM; dir++)
    {
        rho += .1*rho0*sin(2*2*PI*a_x(0));
    }
    double p = p0 + (rho - rho0)*gamma*p0/rho0;
    a_U(0) = rho;
    double ke = 0.;
    for (int dir = 1; dir <= DIM; dir++)
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
/***/
void
parseCommandLine(double& a_tmax, int& a_nx, int& a_maxstep, int& a_outputinterval, int argc, char* argv[])
{
  cout << "Navier Stokes simulation of shear flow with sinusoidal perturbation.  Periodic bcs." << endl;
  cout << "usage:  " << argv[0] << " -n nx  -t tmax -m maxstep -o output_interval" << endl;
  a_tmax= 1.0;
  a_maxstep = 10;
  a_outputinterval = -1;
  a_nx = NUMCELLS;
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

#define offset2(i,j,M) ((j)+(i)*(M))
#define offset3(i,j,k,M,N) ((k)+((j)+(i)*(M))*(N))
#define offset4(i,j,k,l,M,N,P) ((l)+((k)+((j)+(i)*(M))*(N))*(P))

template <unsigned D>
void printBox(const BoxData<double,D>& x) {
    unsigned size = x.size();
    int xmin = (x.box().low())[0];
    int ymin = (x.box().low())[1];
    int xmax = (x.box().high())[0];
    int ymax = (x.box().high())[1];
    int xsize = (xmax - xmin) + 1;
    int ysize = (ymax - ymin) + 1;
    const double *data = x.data();
    for (int x = xmin; x <= xmax; x++) {
        int i = x - xmin;
        for (int y = ymin; y <= ymax; y++) {
            int j = y - ymin;
            for (int d = 0; d < D; d++) {
                double xd = data[offset3(i, j, d, ysize, D)];
                fprintf(stderr, "(x,y,i,j,d,v) = (%d,%d,%d,%d,%d,%g)\n", x, y, i, j, d, xd);
            }
        }
    }
}

/***/
int main(int argc, char* argv[])
{
  //have to do this to get a time table
  PR_TIMER_SETFILE("proto.time.table");
  {
    PR_TIME("main");
    double tstop;
    int size1D, maxStep, outputInterval;
    parseCommandLine(tstop, size1D, maxStep, outputInterval, argc, argv);

    int nGhost = NGHOST;
    EulerOp::s_gamma = 1.4;
    EulerRK4Op::s_count = 0;
    Point lo = Point::Zeros();
    Point hi = Point::Ones(size1D - 1);
    Box dbx0(lo,hi);
    EulerOp::s_dx = 1./size1D;
    EulerState state(dbx0);
    RK4<EulerState,EulerRK4Op,EulerDX> rk4;
    Box dbx = dbx0.grow(nGhost);
    Box dbx1 = dbx.grow(1);
    BoxData<double,NUMCOMPS> UBig(dbx1);
    BoxData<double,DIM> x(dbx1);
    forallInPlace_p(iotaFunc, dbx1, x, EulerOp::s_dx);

    BoxData<double,NUMCOMPS>& U = state.m_U;
    //iota(x,EulerOp::s_dx);
    double dt = .25/size1D;
    Stencil<double> Lap2nd = Stencil<double>::Laplacian();
    cout << "before initializestate"<< endl;
    forallInPlace(InitializeState,dbx1,UBig,x);
    cout << "after initializestate"<< endl;

    U |= Lap2nd(UBig,dbx,1.0/24.0); 
    U += UBig;
    //printBox<NUMCOMPS>(U);

    double time = 0.;
    string resStr = "_"+std::to_string(size1D);
    string fileRoot = "outfile";
    cout << "starting time loop"<< endl;
    for (int k = 0;(k < maxStep) && (time < tstop);k++)
    {
      rk4.advance(time,dt,state);
      time += dt;
      dt = std::min(1.1*dt,.8/size1D/state.m_velSave);
      state.m_velSave = 0.; 
      cout <<"nstep = " << k << " time = " << time << " time step = " << dt << endl;
      if((outputInterval > 0) && (k%outputInterval == 0))
      {
        WriteData(k,state.m_U,EulerOp::s_dx);
      }
    }
  }    
  PR_TIMER_REPORT();

}
