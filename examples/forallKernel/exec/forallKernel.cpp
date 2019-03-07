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
#include "Proto_DisjointBoxLayout.H"
#include "Proto_LevelData.H"
using std::cout;
using std::endl;
using namespace Proto;
constexpr unsigned int NUMCOMPS=DIM+2;
typedef Var<double,NUMCOMPS> State;
#define NMULT 10
/**/
void
parseCommandLine(int & a_nx, int & a_numapplies, int & a_maxgrid, int& a_numstreams, int argc, char* argv[])
{
  //defaults
  a_nx = 128;
  a_numapplies = 100;
  a_maxgrid = 32;
  a_numstreams = 8;
  cout << "kernel timings of riemann and empty forall" << endl;
  cout << "usage:  " << argv[0] << "-s numstreams -m a_maxgrid -n nx[default:128] -a num_iterations[default:10]" << endl;
  for(int iarg = 0; iarg < argc-1; iarg++)
  {
    if(strcmp(argv[iarg],"-n") == 0)
    {
      a_nx = atoi(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg], "-a") == 0)
    {
      a_numapplies = atoi(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg], "-m") == 0)
    {
      a_maxgrid = atoi(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg], "-s") == 0)
    {
      a_numstreams = atoi(argv[iarg+1]);
    }
  }

  cout << "nx          = " << a_nx         << endl;
  cout << "numapplies  = " << a_numapplies << endl;
  cout << "maxgrid     = " << a_maxgrid    << endl;
  cout << "numstreams  = " << a_numstreams << endl;
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
  //took this one out for a bunch of multiplies so
  //I can have flops I can count
//  double cbar = sqrt(gamma*pbar/rhobar);
  double cbar = gamma*pbar/rhobar;
  //NMULT
  for(int iter = 0; iter < NMULT; iter++)
  {
    cbar *= gamma;
  }

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
//2
  double cond =  (cbar + sign*ubar > 0);
 //took out conditional to be sure flops happen
 //   if(cond > 0)
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
doSomeForAlls(int  a_nx, int a_numapplies, int a_maxgrid, int a_numstream,
              LevelData< BoxData<T, NUMCOMPS> > & a_out,
              LevelData< BoxData<T, NUMCOMPS> > & a_low,
              LevelData< BoxData<T, NUMCOMPS> > & a_hig,
              const Box                         & a_domain,
              const DisjointBoxLayout           & a_dbl)
{

  PR_TIME("doSomeForAlls");

  using namespace Proto;
  //remember this is just for timings
  for(int idx = 0; idx < a_dbl.size(); idx++)
  {
    PR_TIME("setVal");
    a_out[idx].setVal(1.);
    a_low[idx].setVal(1.);
    a_hig[idx].setVal(1.);
  }
  double gamma = 1.4;
  int idir = 0;

#ifdef PROTO_CUDA
  vector<cudaStream_t> streams(a_numstream);
  for(unsigned int ibox = 0; ibox < a_numstream; ibox++)
  {
    cudaStreamCreate(&streams[ibox]);
  }
#endif
  

  {
    cout << "doing riemann problems " << endl;
    for(unsigned int ibox = 0; ibox < a_dbl.size(); ibox++)
    {
      int istream = ibox%a_numstream;
      for(unsigned int iapp = 0; iapp < a_numapplies; iapp++)
      {
        PR_TIME("riemann_on_level_multiStream");
        Box appBox       = a_dbl[ibox];

        unsigned long long int count = (28 + NMULT)*appBox.size();
        PR_FLOPS(count);
#ifdef PROTO_CUDA
        cudaForallStream(streams[istream], upwindState, appBox, a_out[ibox], a_low[ibox], a_hig[ibox], idir, gamma, a_nmult);
#else
        forallInPlace(upwindState, appBox, a_out[ibox], a_low[ibox], a_hig[ibox], idir, gamma);
#endif

      }
    }
    sync();
  }

  {
    cout << "doing empty foralls " << endl;
    for(unsigned int ibox = 0; ibox < a_dbl.size(); ibox++)
    {
      int istream = ibox%a_numstream;
      for(unsigned int iapp = 0; iapp < a_numapplies; iapp++)
      {
        PR_TIME("do_nothing_on_level_multiStream");
        Box appBox       = a_dbl[ibox];

#ifdef PROTO_CUDA
        cudaForallStream(streams[istream], doNothing  , appBox, a_out[ibox], a_low[ibox], a_hig[ibox], idir, gamma, a_nmult);
#else
        forallInPlace(doNothing, appBox, a_out[ibox], a_low[ibox], a_hig[ibox], idir, gamma);
#endif
      }
    }
    sync();
  }


#ifdef PROTO_CUDA
  for(unsigned int ibox = 0; ibox < a_numstream; ibox++)
  {
    cudaStreamDestroy(streams[ibox]);
  }
#endif
}

/**/
int main(int argc, char* argv[])
{
  //have to do this to get a time table
  PR_TIMER_SETFILE("proto.time.table");
  int nx, niter, numstreams, maxgrid;

  parseCommandLine(nx, niter, maxgrid, numstreams, argc, argv);

  Point lo = Point::Zeros();
  Point hi = Point::Ones(nx - 1);
  Box domain(lo, hi);
  {
    PR_TIME("forall test");

    LevelData<BoxData<double, NUMCOMPS> > outd, lowd, higd;
    array<bool, DIM> periodic;
    for(int idir = 0 ; idir < DIM; idir++) periodic[idir] = true;
    DisjointBoxLayout dbl(domain, maxgrid, periodic);
    {
      PR_TIME("dataholder definition");
      outd.define(dbl, Point::Zeros());
      lowd.define(dbl, Point::Zeros());
      higd.define(dbl, Point::Zeros());
    }
    doSomeForAlls<double>(nx, niter, maxgrid, numstreams, outd, lowd, higd, domain, dbl);
  }



  PR_TIMER_REPORT();

}  
