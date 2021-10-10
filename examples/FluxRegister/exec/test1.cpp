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
#include "Proto.H"
#include "Proto_DBLInternals.H"
#include "Proto_parstream.H"
#include "Proto_LevelFluxRegister.H"

using namespace std;
using namespace Proto;
int main(int argc, char* argv[])
{
#ifdef PR_MPI
  MPI_Init (&argc, &argv);
#endif
  int logDomainSize;
  int refRatio1D = 2;
  int myproc = procID();
  logDomainSize = 6;
  
#ifdef PR_MPI
  // MPI_Bcast(&logDomainSize, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
  setPoutBaseName("testOut");
  //cout << poutFileName() << endl;
  // barrier();
#endif
      int domainSize = ipow(2,logDomainSize);
      PR_TIMER_SETFILE(to_string(domainSize) + "proto.time.table");
      PR_TIMERS("main");
      Box domain(Point::Zeros(),Point::Ones()*(domainSize -1));
      array<bool,DIM> per;
      for(int idir = 0; idir < DIM; idir++) per[idir]=true;
      double dx = 1.0/domainSize;
      int boxSize = 8;
      int modulus = domainSize % boxSize;
      PROTO_ASSERT((modulus == 0),"Domain not nested");
      Point refRatio = Point::Ones(refRatio1D);
      ProblemDomain pdCoarse(domain,per);
      ProblemDomain pdFine = pdCoarse.refine(refRatio);
      //cout << "bx = " << domain <<  ", coarsened box = " << 
      // domain.coarsen(Point::Ones()*scalarBoxSize) << endl;
      DisjointBoxLayout coarseDBL(pdCoarse,Point::Ones(boxSize));
      Box bxFine = Box(Point::Ones(domainSize/4),Point::Ones(3*domainSize/4 - 1)).refine(refRatio).coarsen(Point::Ones(boxSize));
      cout << bxFine<< endl;
      std::vector<Point> finePts;
      BoxIterator bxit(bxFine);
      for (bxit.begin();!bxit.done();++bxit)
        {
          finePts.push_back(*bxit);
        }
      DisjointBoxLayout fineDBL(pdFine,finePts,Point::Ones(boxSize));
      LevelBoxData<double > phiCoarse(coarseDBL,Point::Zeros());
      LevelBoxData<double > phiFine(fineDBL,Point::Ones());
      phiCoarse.setToZero();
      phiFine.setToZero();
      LevelFluxRegister<double> lfr(coarseDBL,fineDBL,refRatio);
      for (auto dit = phiCoarse.begin();*dit!=dit.end();++dit)
        {
          Box bxCoarse = dit.box();
          for (int dir = 0; dir < DIM; dir++)
            {
              Box bxFlux = bxCoarse.growHi(dir,1);
              BoxData<double> flux(bxFlux);
              flux.setVal(1.0);
              lfr.incrementCoarse(flux,*dit,dir,1.0);
            }
        }
      for (auto dit = phiFine.begin();*dit!=dit.end();++dit)
        {
          Box bxFine = dit.box();
          for (int dir = 0; dir < DIM; dir++)
            {
              Box bxFlux = bxFine.growHi(dir,1);
              BoxData<double> flux(bxFlux);
              flux.setVal(1.0);
              lfr.incrementFine(flux,*dit,dir,1.0);
            }
        }
      lfr.reflux(phiCoarse,1.0);
      PR_TIMER_REPORT();
#ifdef PR_MPI
      MPI_Finalize();
#endif
}
