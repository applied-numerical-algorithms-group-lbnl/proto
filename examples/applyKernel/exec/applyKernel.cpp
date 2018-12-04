#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>

#include "SGMultigrid.H"
#include "Proto_DebugHooks.H"
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"
using std::cout;
using std::endl;
using namespace Proto;

/**/
void
parseCommandLine(int & a_nx, int & a_numapplies, int argc, char* argv[])
{
  cout << "kernel timings of various laplacians" << endl;
  cout << "usage:  " << argv[0] << " -n nx -m num_iterations" << endl;
  for(int iarg = 0; iarg < argc-1; iarg++)
  {
    if(strcmp(argv[iarg],"-n") == 0)
    {
      a_nx = atoi(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg], "-m") == 0)
    {
      a_numapplies = atoi(argv[iarg+1]);
    }
  }
}

/**/

void
applyStuff(int  a_nx, int a_numapplies)
{
  PR_TIME("applyStuff");

  Point lo = Point::Zeros();
  Point hi = Point::Ones(a_nx - 1);
  Box domain(lo, hi);

  Stencil<double> loOrderLap = Stencil<double>::Laplacian();
#if DIM==2
  Stencil<double> hiOrderLap = Stencil<double>::Laplacian_9();
#else 
  Stencil<double> hiOrderLap = Stencil<double>::Laplacian_27();
#endif
  Point ghostPt = hiOrderLap.ghost();
  Box   ghostBx = domain.grow(ghostPt);
  
  BoxData<double> phi(ghostBx);
  BoxData<double> lap(domain);
  
  //remember this is just for timings
  phi.setVal(0.);
  lap.setVal(0.);
  double dx = 1.0/(std::max(double(a_nx), 1.));
  cout << "apply standard laplacian " << a_numapplies << " times" << endl;
  for(int iapp = 0; iapp < a_numapplies; iapp++)
  {
    PR_TIME("STD  laplacian");
    loOrderLap.apply(phi, lap, domain, true, 1.0/(dx*dx));
  }
  cout << "apply dense laplacian " << a_numapplies << " times" << endl;
  for(int iapp = 0; iapp < a_numapplies; iapp++)
  {
    PR_TIME("DENSE laplacian");
    hiOrderLap.apply(phi, lap, domain, true, 1.0/(dx*dx));
  }

}
/**/
int main(int argc, char* argv[])
{
  //have to do this to get a time table
  PR_TIMER_SETFILE("proto.time.table");
  int nx, niter;
  parseCommandLine(nx, niter, argc, argv);
  applyStuff(nx, niter);


  PR_TIMER_REPORT();

}  
