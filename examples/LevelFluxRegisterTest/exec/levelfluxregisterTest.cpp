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
#include "Proto_LevelFluxRegister.H"

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
#ifdef PR_MPI
  if (myproc == 0)
#endif
    {
      cout << "input log_2(domainSize)"  
        //, number of multigrid levels"
           << endl;
      cin >> logDomainSize; //>> numLevels; 
      //cout << "input max number of iterations, convergence tolerance " << endl;
      //cin >> maxiter >> tol;
    }
#ifdef PR_MPI
  MPI_Bcast(&logDomainSize, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
  MPI_Bcast(&numLevels, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
  MPI_Bcast(&maxiter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
  barrier();
  //cout << "numLevels = " << numLevels << ", myproc = " << myproc << endl;
  int domainSize = ipow(2,logDomainSize);
  PR_TIMER_SETFILE(to_string(domainSize) + "LFR.time.table");
  PR_TIMERS("main");
  
  Box coarseDomain(Point::Zeros(),Point::Ones()*(domainSize -1));
  Point refRatio = 4*Point::Ones();
  Box fineDomain = coarseDomain.refine(refRatio);
  
  array<bool,DIM> per;
  for(int idir = 0; idir < DIM; idir++) per[idir]=true;
  
  double dxCoarse = 1.0/domainSize;
  double dxFine = dxCoarse/4;
  double dxCinv = 1.0/dxCoarse;
  double dxFinv = 1.0/dxFine;

  int scalarBoxSize = 64;
  int modulus = domainSize % scalarBoxSize;
  PROTO_ASSERT((modulus == 0),"Domain not nested");
  ProblemDomain pdCoarse(coarseDomain,per);
  ProblemDomain pdFine(fineDomain,per);
  Box fineBitMap =fineDomain.coarsen(scalarBoxSize*Point::Ones());
  vector<Point > refinedPatches;
  Point pt = (fineBitMap.high() + Point::Ones())/2;
  refinedPatches.push_back(pt);
  DisjointBoxLayout dblCoarse(pdCoarse,scalarBoxSize*Point::Ones());
  DisjointBoxLayout dblFine(pdFine,refinedPatches,scalarBoxSize*Point::Ones());
  LevelBoxData<double,1,MEMTYPE_DEFAULT> phiCoarse(dblCoarse,Point::Ones());
  LevelBoxData<double,1,MEMTYPE_DEFAULT> phiFine(dblFine,Point::Ones());
  LevelBoxData<double,1,MEMTYPE_DEFAULT> LOfPhiCoarse(dblCoarse,Point::Zeros());
  LevelBoxData<double,1,MEMTYPE_DEFAULT> LOfPhiFine(dblFine,Point::Zeros());
  phiFine.setToZero();
  phiCoarse.setToZero();
  LOfPhiFine.setToZero();
  LOfPhiCoarse.setToZero();

   LevelFluxRegister<double> lfr(dblCoarse,dblFine,refRatio);
   lfr.reset();

   array<Stencil<double>,DIM > fluxStencil;
   array<Stencil<double>,DIM > divStencil;

   for (int d = 0; d< DIM;d++)
     {
       fluxStencil[d] = 1.0*Shift(Point::Zeros()) + (-1.0)*Shift(Point::Basis(d,-1));
       divStencil[d] = 1.0*Shift(Point::Basis(d,1)) + (-1.0)*Shift(Point::Zeros());
     }
  for (auto dit = phiCoarse.begin();*dit != dit.end();++dit)
    {
      BoxData<double>& phiPatch = phiCoarse[*dit];
      phiPatch.setVal(1.0);
    }
for (auto dit = phiFine.begin();*dit != dit.end();++dit)
    {
      BoxData<double>& phiPatch = phiFine[*dit];
      phiPatch.setVal(1.0);
    }


 for (auto dit = phiCoarse.begin();*dit!=dit.end();++dit)
   {
     for (int d = 0; d < DIM; d++)
       {
         BoxData<double,1> flux = fluxStencil[d](phiCoarse[*dit],dxCinv);
         lfr.incrementCoarse(flux,*dit,-1.0,d);
         LOfPhiCoarse[*dit] += divStencil[d](flux,dxCinv);
       }     
   }
 for (auto dit = phiFine.begin();*dit!=dit.end();++dit)
   {
     for (int d = 0; d < DIM; d++)
       {         
         auto flux = fluxStencil[d](phiFine[*dit],dxFinv);
         lfr.incrementCoarse(flux,*dit,1.0,d);
         LOfPhiFine[*dit] += divStencil[d](flux,dxFinv);
       }
   }
 lfr.reflux(LOfPhiCoarse,dxCinv);
  PR_TIMER_REPORT();
#ifdef PR_MPI
  MPI_Finalize();
#endif
}

