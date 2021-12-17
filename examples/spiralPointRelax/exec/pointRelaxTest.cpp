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
#include "Proto_DisjointBoxLayout.H"
#include "Proto_LevelBoxData.H"

using namespace std;
using namespace Proto;
//inline 
PROTO_KERNEL_START void jacobiUpdateT(Var<double> a_phi, 
				    Var<double> a_Lphi, 
				    Var<double> a_rhs, 
				    double a_lambda)
{
	a_phi(0) = a_phi(0) + a_Lphi(0) - a_lambda*a_rhs(0);
};
PROTO_KERNEL_END(jacobiUpdateT, jacobiUpdate);

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
double computeMaxResidualAcrossProcs(
                                     LevelBoxData<double>& a_phi,
                                     LevelBoxData<double>& a_rho,
                                     const double& a_h)
{
  Reduction<double,Abs> rxn;
  PR_TIMERS("resnorm");
  double hsqi = 1./a_h/a_h;
  a_phi.exchange();
  double maxnorm = 0.;
  
  for (auto dit=a_phi.begin();*dit != dit.end();++dit)
    {
      BoxData<double>& phi = a_phi[*dit];
      BoxData<double>& rho = a_rho[*dit];
      BoxData<double> res(dit.box());
      
      res.setVal(0.);
      res -= rho;
      res += Stencil<double>::Laplacian()(phi,hsqi);
      //Reduction doesn't seem to work properly for me
      //maxnorm = max(m_rxn.fetch(),maxnorm);
      maxnorm = max(res.absMax(),maxnorm);
    }
  return maxnorm;
}


int main(int argc, char* argv[])
{
#ifdef PR_MPI
    MPI_Init (&argc, &argv);
#endif
    int logDomainSize = 8;
    int maxiter = 20;
    double tol = 1e-10;
    int myproc = procID();
    if (myproc == 0)
    {
        cout << "input log_2(domainSize)" << endl;
        cin >> logDomainSize; 
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
    double wgt = 1./(4.*DIM);
    double lambda = dx*dx/(4.0*DIM);
    int scalarBoxSize = 64;

    PROTO_ASSERT((domainSize % scalarBoxSize == 0), "Domain not nested: %i mod %i != 0", 
            domainSize, scalarBoxSize);

    ProblemDomain pd(domain,per);
    DisjointBoxLayout dbl(pd,Point::Ones(scalarBoxSize));

    LevelBoxData<double > rho(dbl,Point::Zeros());
    LevelBoxData<double > phi(dbl,Point::Ones());

    rho.setToZero();
    phi.setToZero();
    for (auto dit = phi.begin();*dit != dit.end();++dit)
    {
        BoxData<double>& rhoPatch = rho[*dit];
        forallInPlace_p(rhsPoint,rhoPatch,dx);
    }
    double resmax0=computeMaxResidualAcrossProcs(phi,rho,dx);
    
    if (myproc==0) 
    {
        cout << "initial residual = " << resmax0 << endl;
    }
    for (int iter = 0; iter < maxiter; iter++)
    {
        PR_TIMERS("PointRelax");

         phi.exchange();
         for (auto dit=phi.begin();*dit != dit.end();++dit)
           {
             BoxData<double>& phiPatch = phi[*dit];
             BoxData<double>& rhoPatch = rho[*dit];
             BoxData<double> temp = Stencil<double>::Laplacian()(phiPatch,wgt);
             forallInPlace(jacobiUpdate,phiPatch,temp,rhoPatch,lambda);
           }
         double resmax=computeMaxResidualAcrossProcs(phi,rho,dx);
         if (myproc==0) 
           {
             cout << "iter = " << iter << ", resmax = " << resmax << endl;
           }
         if (resmax < tol*resmax0) break;
    }
#ifdef PR_HDF5
        HDF5Handler h5;
        h5.writeLevel(phi, "MG_PHI.hdf5", 0);
#endif
    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
}

