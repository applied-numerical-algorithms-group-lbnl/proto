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
#include "Proto_DebugHooks.H"
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"
using std::cout;
using std::endl;
using namespace Proto;
constexpr unsigned int NUMCOMPS=DIM+2;
typedef Var<double,NUMCOMPS> State;


/**/
void
parseCommandLine(int & a_nx, int & a_numapplies, int argc, char* argv[])
{
  //defaults
  a_nx = 64;
  a_numapplies = 100;
  cout << "kernel timings of various laplacians" << endl;
  cout << "usage:  " << argv[0] << " -n nx[default:64] -m num_iterations[default:100]" << endl;
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

PROTO_KERNEL_START
void upwindStateF(State& a_out,
                  const State& a_low,
                  const State& a_high,
                  int   a_dir,
                  double a_gamma)
{
  const double& rhol = a_low(0);
  const double& rhor = a_high(0);
  const double& ul = a_low(a_dir+1);
  const double& ur = a_high(a_dir+1);
  const double& pl = a_low(NUMCOMPS-1);
  const double& pr = a_high(NUMCOMPS-1);
  double gamma = a_gamma;
  double rhobar = (rhol + rhor)*.5;
  double pbar = (pl + pr)*.5;
  double ubar = (ul + ur)*.5;
  double cbar = sqrt(gamma*pbar/rhobar);
  double pstar = (pl + pr)*.5 + rhobar*cbar*(ul - ur)*.5;
  double ustar = (ul + ur)*.5 + (pl - pr)/(2*rhobar*cbar);
  int sign;
  if (ustar > 0) 
  {
    sign = -1;
    for (int icomp = 0;icomp < NUMCOMPS;icomp++)
    {
      a_out(icomp) = a_low(icomp);
    }
  }
  else
  {
    sign = 1;
    for (int icomp = 0;icomp < NUMCOMPS;icomp++)
    {
      a_out(icomp) = a_high(icomp);
    }
  }
  if (cbar + sign*ubar > 0)
  {
    a_out(0) += (pstar - a_out(NUMCOMPS-1))/(cbar*cbar);
    a_out(a_dir+1) = ustar;
    a_out(NUMCOMPS-1) = pstar;
  }
}
PROTO_KERNEL_END(upwindStateF, upwindState)


PROTO_KERNEL_START
void doNothingF(State& a_out,
                  const State& a_low,
                  const State& a_high,
                  int   a_dir,
                  double a_gamma)
{
}
PROTO_KERNEL_END(doNothingF, doNothing)

/**/
inline void sync()
{
  #ifdef PROTO_CUDA
    {
      PR_TIME("device sync");
      cudaDeviceSynchronize();
    }
#endif
}

template <class T> void
doSomeForAlls(int  a_nx, int a_numapplies,
              BoxData<T, NUMCOMPS>& out,
              BoxData<T, NUMCOMPS>& low,
              BoxData<T, NUMCOMPS>& hig,
              Box domain)
{

  PR_TIME("doSomeForAlls");

  //remember this is just for timings
  out.setVal(1.);
  low.setVal(1.);
  hig.setVal(1.);
  T gamma = 1.4;
  int idir = 0;
  cout << "do riemann problem " << a_numapplies << " times" << endl;
  {
    PR_TIME("riemann problem");
    for(int iapp = 0; iapp < a_numapplies; iapp++)
    {
      PR_TIME("actual forall");
      forallInPlace(upwindState, domain, out, low, hig, idir, gamma);
    }
    sync();
  }
  cout << "do nothing " << a_numapplies << " times" << endl;
  {
    for(int iapp = 0; iapp < a_numapplies; iapp++)
    {
      PR_TIME("empty forall");
      forallInPlace(doNothing, domain, out, low, hig, idir, gamma);
    
    }
  }
}
/**/
int main(int argc, char* argv[])
{
  //have to do this to get a time table
  PR_TIMER_SETFILE("proto.time.table");
  int nx, niter;
  parseCommandLine(nx, niter, argc, argv);

  Point lo = Point::Zeros();
  Point hi = Point::Ones(nx - 1);
  Box domain(lo, hi);
  {
    PR_TIME("forall test");

    BoxData<double, NUMCOMPS> outd, lowd, higd;
    {
      PR_TIME("dataholder definition");
      outd.define(domain);
      lowd.define(domain);
      higd.define(domain);
    }
    doSomeForAlls<double>(nx, niter, outd, lowd, higd, domain);
  }



  PR_TIMER_REPORT();

}  
