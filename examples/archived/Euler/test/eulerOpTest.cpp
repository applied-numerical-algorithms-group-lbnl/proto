#include <cstdio>
#include <cstring>
#include <cmath>

#include <vector>

#include <iostream>
#include <fstream>
#include <sstream>

#include "Proto.H"
#include "EulerOp.H"

#define PI 3.141592653589793

using namespace std;
using namespace Proto;

typedef Var<double,DIM> V;
typedef Var<double,NUMCOMPS> State;

unsigned int InitializeState(State& a_U,
                             const V& a_x)
{
    double gamma = 1.4;
    double rho0 = gamma;
    double p0 = 1.;
    double umag = 1.0;
    double rho = rho0;
    for (int dir = 0; dir < DIM; dir++)
    {
       rho += .1*rho0*sin(4*2*PI*a_x(0));
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
    a_U(DIM+1) = p*gamma/(gamma-1.0) + rho*ke;
    return 0;
}   
  
int main(int argc, char* argv[])
{
    #ifdef CH_MPI
    MPI_Init(&argc, &argv);
    #endif
  
    //CH_TIMERS("main");
    int size1D, nGhost; 
    size1D = atoi(argv[1]);
    int size1D_0 = size1D;
    nGhost = 4;
    BoxData<double,NUMCOMPS> mDivF[3];
    EulerOp eop;
    int numLevel = 3;
    for (int level = 0; level < numLevel; level++)
    {
        Box b = Box::Cube(size1D);
        EulerOp::s_dx = 1.0/size1D;
        Box B0 = b.grow(nGhost); 
        mDivF[level].define(b);
        Box B1 = B0.grow(1);
        BoxData<double,NUMCOMPS> UBig(B1),U(B0);
      
        BoxData<double,DIM> x = BoxData<double,DIM>::iota(B1,EulerOp::s_dx);
        forall(InitializeState,B1,UBig,x);
        U |= Stencil<double>::Laplacian()(UBig,B0); 
        U *= (1.0/24.0);
        U += UBig;
        //U.copyTo(mDivF[level]);
        double maxvel = eop(mDivF[level],U,b);
        size1D*=2;
    }
    
    Stencil<double> avgDown = Stencil<double>::AvgDown(2);
    Box printBox = Box(Point::Zeros(),Point::Basis(0,15));
    double error[numLevel-1];
    for (int ii = 1; ii < numLevel; ii++)
    {
        auto coarse = alias(mDivF[ii-1]);
        auto fine = alias(mDivF[ii]);
        BoxData<double,NUMCOMPS> fineToCoarse(coarse.box());
        fineToCoarse |= avgDown(fine,coarse.box());
        fineToCoarse -= coarse;
        error[ii-1] = fineToCoarse.absMax();
    }
    cout << "Coarse Error: " << error[0] << endl;
    cout << "Fine Error: " << error[1] << endl;
    cout << "Convergence Rate: " << log(error[0]/error[1])/log(2.0) << endl;
    #ifdef CH_MPI
    //CH_TIMER_REPORT();
    MPI_Finalize();
    #endif
    return 0;
}
