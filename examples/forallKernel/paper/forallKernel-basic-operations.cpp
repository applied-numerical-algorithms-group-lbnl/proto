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
#include "Proto_LevelBoxData.H"
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
  a_nx = 512;
  a_numapplies = 1;
  a_maxgrid = 512;
  a_numstreams = 1;
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
void doAddF(State& a_out,
                  const State& a_in)
{
#pragma unroll
    for (int icomp = 0;icomp < NUMCOMPS;icomp++)
    {
	    auto in = a_in(icomp);
	    auto &out = a_out(icomp);
	    out += in;
     // a_out(icomp) += a_in(icomp);
    }
}
PROTO_KERNEL_END(doAddF, doAdd)


PROTO_KERNEL_START
void doSubF(State& a_out,
                  const State& a_in)
{
      unsigned int size = 512*512*512;
      double* out = &a_out(0);
      const double* in = &a_in(0);
#pragma unroll
      for (int icomp = 0;icomp < 5;icomp++)
      {
        *out -= *in;
        out += size; in += size;
      }
}
PROTO_KERNEL_END(doSubF, doSub)

PROTO_KERNEL_START
void doMulF(State& a_out,
                  const State& a_in)
{
#pragma unroll
    for (int icomp = 0;icomp < NUMCOMPS;icomp++)
    {
      a_out(icomp) *= a_in(icomp);
    }
}
PROTO_KERNEL_END(doMulF, doMul)

PROTO_KERNEL_START
void doDivF(State& a_out,
                  const State& a_in)
{
#pragma unroll
    for (int icomp = 0;icomp < NUMCOMPS;icomp++)
    {
      a_out(icomp) /= a_in(icomp);
    }
}
PROTO_KERNEL_END(doDivF, doDiv)
/**/
inline void sync()
{
  #ifdef PROTO_CUDA
    {
      PR_TIME("device sync");
      protoDeviceSynchronize(DEVICE);
    }
#endif
}

template <class T> void
doSomeForAlls(int  a_nx, int a_numapplies, int a_maxgrid, int a_numstream,
              LevelBoxData<T, NUMCOMPS> & a_out,
              LevelBoxData<T, NUMCOMPS> & a_low,
              LevelBoxData<T, NUMCOMPS> & a_hig,
              const Box                         & a_domain,
              const DisjointBoxLayout           & a_dbl)
{

  PR_TIME("doSomeForAlls");

  using namespace Proto;
  //remember this is just for timings
  for(DataIterator dit = a_out.begin(); (*dit)!=dit.end(); ++dit)
  {
    PR_TIME("setVal");
    a_out[*dit].setVal(1.);
    a_low[*dit].setVal(1.);
  }
  
  {
    cout << "do basics operations - forAll " << endl;
    PR_TIME("do basics operations - forAll" );
    //for(DataIterator dit(a_out.getDBL()); *dit!=dit.end(); ++dit)
    for(DataIterator dit = a_out.begin(); (*dit)!=dit.end(); ++dit)
    {
      for(unsigned int iapp = 0; iapp < a_numapplies; iapp++)
      {
        Box appBox       = a_out[*dit].box();
	
	{
          PR_TIME("do add - forAll" );
          unsigned int count = 0; //appBox.size();
          PR_FLOPS(count);
	  protoForall(doAdd, appBox, a_out[*dit], a_low[*dit]);
	}
	{
          PR_TIME("do sub - forAll" );
          unsigned int count = 0;// appBox.size();
          PR_FLOPS(count);
	  protoForall(doSub, appBox, a_out[*dit], a_low[*dit]);
	}
	{
          PR_TIME("do mul - forAll" );
          unsigned int count = 0; //appBox.size();
          PR_FLOPS(count);
	  protoForall(doMul, appBox, a_out[*dit], a_low[*dit]);
	}
	{
          PR_TIME("do div - forAll" );
          unsigned int count = 0; //appBox.size();
          PR_FLOPS(count);
	  protoForall(doDiv, appBox, a_out[*dit], a_low[*dit]);
	}
      }
    }
  }
  {
    PR_TIME("do basics operations - no forAll" );
    for(DataIterator dit = a_out.begin(); (*dit)!=dit.end(); ++dit)
    {
      for(unsigned int iapp = 0; iapp < a_numapplies; iapp++)
      {
        Box appBox       = a_out[*dit].box();
	
	{
          PR_TIME("do add - no forAll" );
          unsigned int count = 0;//appBox.size();
          PR_FLOPS(count);
	  a_out[*dit] += a_low[*dit];
	}
	{
          PR_TIME("do sub - no forAll" );
          unsigned int count = 0;//appBox.size();
          PR_FLOPS(count);
	  a_out[*dit] -= a_low[*dit];
	}
	{
          PR_TIME("do mul - no forAll" );
          unsigned int count = 0;//appBox.size();
          PR_FLOPS(count);
	  a_out[*dit] *= a_low[*dit];
	}
	{
          PR_TIME("do div - no forAll" );
          unsigned int count = 0;//appBox.size();
          PR_FLOPS(count);
//	  a_out[*dit] \= a_low[*dit];
	}
      }
    }

  }
  sync();
}

/**/
int main(int argc, char* argv[])
{
#ifdef PR_MPI
    MPI_Init(&argc,&argv);
#endif
  //have to do this to get a time table
  PR_TIMER_SETFILE("proto.time.table");
  int nx, niter, numstreams, maxgrid;

  parseCommandLine(nx, niter, maxgrid, numstreams, argc, argv);

  Point lo = Point::Zeros();
  Point hi = Point::Ones(nx - 1);
  Box domain(lo, hi);
  {
    PR_TIME("forall test");

    LevelBoxData<double, NUMCOMPS> outd, lowd, higd;
    array<bool, DIM> periodic;
    for(int idir = 0 ; idir < DIM; idir++) periodic[idir] = true;
    ProblemDomain probDom(domain,periodic);
    DisjointBoxLayout dbl(probDom, maxgrid*Point::Ones()); //This will create a disjoint layout with maxgrid size boxes
    {
      PR_TIME("dataholder definition");
      outd.define(dbl, Point::Zeros());
      lowd.define(dbl, Point::Zeros());
      higd.define(dbl, Point::Zeros());
    }
    doSomeForAlls<double>(nx, niter, maxgrid, numstreams, outd, lowd, higd, domain, dbl);
  }



  PR_TIMER_REPORT();

#ifdef PR_MPI
  MPI_Finalize();
#endif

}  
