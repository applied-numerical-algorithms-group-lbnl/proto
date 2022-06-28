#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <memory>
#include <stdio.h>
#include <fstream>
#include "Proto.H"
#include <iostream>
#include <cstring>
#include <memory>
#include "LevelMultigrid.H"
#include "InputParser.H"

using namespace std;
using namespace Proto;

PROTO_KERNEL_START void f_initT(const Point& a_pt, Var<double> a_rho,double a_h)
{
    a_rho(0) = 1.;
    for (int idir = 0; idir < DIM; idir++)
    {
        a_rho(0) = a_rho(0)*sin(M_PI*2*(a_pt[idir]*a_h + .5*a_h + .125));
    }
}
PROTO_KERNEL_END(f_initT, f_init);

int main(int argc, char* argv[])
{
#ifdef PR_MPI
    MPI_Init (&argc, &argv);
#endif
    // DEFAULT PARAMETERS
    int domainSize = 256;
    int boxSize = 64;
    int numLevels = 8;
    int maxIter = 20;
    double tolerance = 1e-10;

    // PARSE COMMAND LINE
    InputArgs args;
    args.add("domainSize",    domainSize);
    args.add("boxSize",       boxSize);
    args.add("numLevels",     numLevels);
    args.add("maxIter",       maxIter);
    args.add("tolerance",     tolerance);
    args.parse(argc, argv);
    args.print();

    // INITIALIZE TIMERS
    PR_TIMER_SETFILE(to_string(domainSize) + ".DIM" + to_string(DIM) + ".LevelMultigrid.time.table");
    PR_TIMERS("main");

    // INITIALIZE DOMAIN
    auto domain = Box::Cube(domainSize);
    array<bool,DIM> per;
    per.fill(true);
    double dx = 1.0/domainSize;
    ProblemDomain pd(domain,per);
    DisjointBoxLayout layout(pd,Point::Ones(boxSize));

    // INITIALIZE DATA
    LevelBoxData<double > rho(layout,Point::Zeros());
    LevelBoxData<double > phi(layout,Point::Ones());
    rho.setToZero();
    phi.setToZero();
    for (auto dit : layout)
    {
        BoxData<double>& rho_i = rho[dit];
        forallInPlace_p(f_init, rho_i, dx);
    }
    //phi.exchange();
    
    // SOLVE
    LevelMultigrid mg(layout, dx, numLevels);
    double resmax0 = mg.resnorm(phi, rho);
    Proto::pout() << "initial residual = " << resmax0 << endl;
    for (int iter = 0; iter < maxIter; iter++)
    {
        PR_TIMERS("main::solver");
        mg.vCycle(phi,rho);
#ifdef PR_HDF5
        HDF5Handler h5;
        h5.writeLevel(phi, "MG_Phi_I%i.hdf5", iter);
#endif
        double resmax = mg.resnorm(phi, rho);
        Proto::pout() << "\titer: " << iter << " | resnorm (max-norm): " << resmax << endl;
        if (resmax < tolerance*resmax0)
        {
            Proto::pout() << "Converged in " << iter+1 << " iterations." << std::endl;
            break;
        }
    }
    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
}

