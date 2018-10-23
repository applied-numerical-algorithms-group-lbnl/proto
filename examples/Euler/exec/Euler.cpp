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
#include "WriteBoxData.H"
#include "Proto_Timer.H"

#define PI 3.141592653589793
#define NUMCOMPS DIM+2

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


int main(int argc, char* argv[])
{
    PR_TIME("main");

    int size1D, nGhost,  maxStep;
    double tstop;

    tstop  = 0.1;
    size1D  = 64;
    maxStep = 2;

    nGhost = NGHOST;
    EulerOp::s_gamma = 1.4;
    EulerRK4Op::s_count = 0;
//    for (int level = 0; level < maxLev;level++)
//    {
        Point lo = Point::Zeros();
        Point hi = Point::Ones(size1D - 1);
        Bx dbx0(lo,hi);
        EulerOp::s_dx = 1./size1D;
        EulerState state(dbx0);
        RK4<EulerState,EulerRK4Op,EulerDX> rk4;
        Bx dbx = dbx0.grow(nGhost);
        Bx dbx1 = dbx.grow(1);
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
        double time = 0.;
        string resStr = "_"+std::to_string(size1D);
        string fileRoot = "outfile";
        cout << "starting time loop"<< endl;
        for (int k = 0;(k < maxStep) && (time < tstop);k++)
        {
            rk4.advance(time,dt,state);
            time += dt;
            dt = min(1.1*dt,.8/size1D/state.m_velSave);
            state.m_velSave = 0.; 
            cout <<"nstep = " << k << " time = " << time << " time step = " << dt << endl;
//            WriteData(k,state.m_U,EulerOp::s_dx);
        }
        string fileString = fileRoot + resStr;
        size1D*=2;
//    }
    
    TraceTimer::report();

}
