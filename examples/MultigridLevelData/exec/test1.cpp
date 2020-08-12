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

using namespace std;
using namespace Proto;
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
          //cout << "input log_2(domainSize), number of multigrid levels" << endl;
          //cin >> logDomainSize >> numLevels; 
          //cout << "input max number of iterations, convergence tolerance " << endl;
          //cin >> maxiter >> tol;
          logDomainSize = 4;
          numLevels = 3;
          maxiter = 1;
          tol = 1.e-10;
        }
#ifdef PR_MPI
      MPI_Bcast(&logDomainSize, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
      MPI_Bcast(&numLevels, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
      MPI_Bcast(&maxiter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
      MPI_Bcast(&tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      setPoutBaseName("testOut");
      //cout << poutFileName() << endl;
#endif
      barrier();
      //cout << "numLevels = " << numLevels << ", myproc = " << myproc << endl;
      int domainSize = 3*4*5; //ipow(2,logDomainSize);
      PR_TIMER_SETFILE(to_string(domainSize) + "proto.time.table");
      PR_TIMERS("main");
      Box domain(Point::Zeros(),Point::Ones()*(domainSize -1));
      array<bool,DIM> per;
      for(int idir = 0; idir < DIM; idir++) per[idir]=true;
      double dx = 1.0/domainSize;
      int phiBoxSize = 3*5;
      int rhoBoxSize = 2*5;
      int modulus = domainSize % phiBoxSize;
      PROTO_ASSERT((modulus == 0),"Domain not nested");
      ProblemDomain pd(domain,per);
      //cout << "bx = " << domain <<  ", coarsened box = " << 
      // domain.coarsen(Point::Ones()*scalarBoxSize) << endl;
      DisjointBoxLayout dblRho(pd,rhoBoxSize*Point::Ones());
      DisjointBoxLayout dbl(pd,phiBoxSize*Point::Ones());
#if 0
      {
        PR_TIMERS("dbltest");
        DataIterator dit(dbl);
        for (dit = dit.begin();*dit!=dit.end();++dit)
          {
            Point pt = dit.point();
            int index = (*dit).intIndex();
            pout() << "test:: DataIterator Point = " << pt << ", index = " << index << endl;
            Box bx = dit.box();
            bx = bx.grow(Point::Ones());
            pout() << "test:: Neighbor Box = " << bx << endl;
            NeighborIterator nit(dbl,bx);
            for (nit.begin();*nit != nit.end(); ++nit)
              {
                Box srcbx = nit.srcBox();
                Box destbx = nit.destBox();
                Point shift=nit.shift();
                pout() << "test::       src, size = " << srcbx << " , " 
                       << srcbx.size() << ", dest, size = "<< destbx  << " , " 
                       << destbx.size() << ", shift = " << shift << ", point = "
                       << nit.point() <<  endl;
              }
          }
      }
#endif
      LevelBoxData<double > rho(dblRho,Point::Zeros());
      LevelBoxData<double > phi(dbl,Point::Ones());
      rho.setToZero();
      phi.setToZero();
      phi.copyTo(rho);
      rho.copyTo(phi);
      //cout << "before exchange" << endl;
      //phi.exchange();
      PR_TIMER_REPORT();
      MPI_Finalize();
}
