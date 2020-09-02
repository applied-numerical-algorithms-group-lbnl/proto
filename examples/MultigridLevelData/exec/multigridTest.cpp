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
  int logDomainSize;
  int numLevels;
  int maxiter;
  double tol;
  cout << "input log_2(domainSize), number of multigrid levels" << endl;
  cin >> logDomainSize >> numLevels;
  cout << "input max number of iterations, convergence tolerance " << endl;
  cin >> maxiter >> tol;
  int domainSize = ipow(2,logDomainSize);
  PR_TIMER_SETFILE(to_string(domainSize) + "proto.time.table");
  PR_TIMERS("main");
  
  Box domain(Point::Zeros(),Point::Ones()*(domainSize -1));
  array<bool,DIM> per;
  for(int idir = 0; idir < DIM; idir++) per[idir]=true;
  double dx = 1.0/domainSize;
  DisjointBoxLayout bl(domain,min(MAXBOXSIZE,domainSize),per);

  LevelData<BoxData<double> > rho(bl,Point::Zeros());
  LevelData<BoxData<double> > phi(bl,Point::Ones());
  rho.setToZero();
  phi.setToZero();
  for (int i = 0; i < bl.size();i++)
    {
      BoxData<double>& rhoPatch = rho[i];
      forallInPlace_p(rhsPoint,rhoPatch,dx);
    }
  LevelMultigrid mg(domain,dx,numLevels);
  {
    double resnorm0 = mg.resnorm(phi,rho);
    cout << "initial residual = " << resnorm0 << endl;
    for (int iter = 0; iter < maxiter; iter++)
      {
        PR_TIMERS("MG top level");
        mg.vCycle(phi,rho);
        // WriteData(phi,iter,dx,"phi");
        double resnorm = mg.resnorm(phi,rho);
        cout << "iter = " << iter << ", resnorm = " << resnorm << endl;
        if (resnorm < tol*resnorm0) break;
      }
  }
  PR_TIMER_REPORT();
}
