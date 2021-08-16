#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <memory>
#include <stdio.h>
#include <fstream>
#include "Proto.H"
#include "Proto_SPMD.H"
#include <iostream>
#include <cstring>
#include <memory>
#include "LevelMultigrid.H"

using namespace std;
using namespace Proto;
//inline 

PROTO_KERNEL_START void rhsPointT(const Point& a_pt, Var<double> a_rho,double a_h)
{
    a_rho(0) = 1.;
    for (int idir = 0; idir < DIM; idir++)
    {
        a_rho(0) = a_rho(0)*sin(M_PI*2*(a_pt[idir]*a_h + .5*a_h + .125));
    }
}
PROTO_KERNEL_END(rhsPointT, rhsPoint);

//Compute the max of the residual across all processes.
//The max is then broadcast to all the processes.
double computeMaxResidualAcrossProcs(LevelMultigrid& mg,
        LevelBoxData<double>& phi,
        LevelBoxData<double>& rho)
{
    std::cout << "Computing residual... ";
    double ret_val;//=0;
    double resnorm = mg.resnorm(phi,rho);
#ifdef PR_MPI
    MPI_Reduce(&resnorm,&ret_val,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Bcast(&ret_val, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
    ret_val=resnorm;
#endif
    barrier();
    std::cout << "done." << std::endl;
    return ret_val;
}


int main(int argc, char* argv[])
{
#ifdef PR_MPI
    MPI_Init (&argc, &argv);
#endif
    int logDomainSize;
    int numLevels;
    int maxiter;
    double tol;
    int myproc = procID();
    if (myproc == 0)
    {
        cout << "input log_2(domainSize), number of multigrid levels" << endl;
        cin >> logDomainSize >> numLevels; 
        cout << "input max number of iterations, convergence tolerance " << endl;
        cin >> maxiter >> tol;
    }
#ifdef PR_MPI
    MPI_Bcast(&logDomainSize, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numLevels, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxiter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    barrier();
    cout << "log_domain: " << logDomainSize
        << " | numLevels: " << numLevels
        << " | maxiter: " << maxiter
        << " | tol: " << tol
        << " | proc: " << myproc << endl;
    int domainSize = ipow(2,logDomainSize);
    PR_TIMER_SETFILE(to_string(domainSize) + ".forall.proto.time.table");
    PR_TIMERS("main");

    auto domain = Box::Cube(domainSize);

    array<bool,DIM> per;
    for(int idir = 0; idir < DIM; idir++) { per[idir]=true; }
    double dx = 1.0/domainSize;
    int scalarBoxSize = 64;
    int modulus = domainSize % scalarBoxSize;
    PROTO_ASSERT((modulus == 0), "Domain not nested: %i mod %i != 0", domainSize, scalarBoxSize);
    ProblemDomain pd(domain,per);
    DisjointBoxLayout dbl(pd,Point::Ones(scalarBoxSize));

    LevelBoxData<double > rho(dbl,Point::Zeros());
    LevelBoxData<double > phi(dbl,Point::Ones());

    rho.setToZero();
    phi.setToZero();
    double resmax;
    for (auto dit = phi.begin();*dit != dit.end();++dit)
    {
        BoxData<double>& rhoPatch = rho[*dit];
        forallInPlace_p(rhsPoint,rhoPatch,dx);
    }
    cout << "Building LevelMultigrid: " << Proto::procID() << endl;
    LevelMultigrid mg(dbl,dx,numLevels);
    cout << "Finished Building LevelMultigrid: address: " << &mg << " | process: " << Proto::procID() << endl;
    cout << "Computing Residual: " << Proto::procID() << endl;
    double resmax0=computeMaxResidualAcrossProcs(mg,phi,rho);
    cout << "Finished Computing Residual: " << Proto::procID() << endl;
    if (myproc==0) 
    {
        cout << "initial residual = " << resmax0 << endl;
    }
    for (int iter = 0; iter < maxiter; iter++)
    {
        PR_TIMERS("MG top level");
        mg.vCycle(phi,rho);
        // WriteData(phi,iter,dx,"phi");
        double resmax=computeMaxResidualAcrossProcs(mg,phi,rho);
        if (myproc==0) 
        {
            //cout << "iter = " << iter << ", resmax = " << resmax << endl;
        }
        if (resmax < tol*resmax0) break;
    }

    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
}

