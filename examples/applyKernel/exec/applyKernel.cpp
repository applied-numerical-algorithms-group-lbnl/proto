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
#include "implem/Proto_LevelData.H"
#include "implem/Proto_DisjointBoxLayout.H"
#include "Proto_DebugHooks.H"
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"
using std::cout;
using std::endl;
using namespace Proto;
constexpr unsigned int NUMCOMPS=DIM+2;

/**/
void
parseCommandLine(int & a_nx, int & a_numapplies, int & a_maxgrid, int& a_numstreams, int argc, char* argv[])
{
  //defaults
  a_nx = 128;
  a_numapplies = 10;
  a_maxgrid = 32;
  a_numstreams = 8;
  cout << "kernel timings of stencil applies" << endl;
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
#ifdef PROTO_CUDA
__global__ void empty(){ ;}
#endif

inline void emptyKernel(int a_nx)
{
#ifdef PROTO_CUDA
  empty<<<a_nx*a_nx, a_nx>>>();
#endif
}
inline void sync()
{
#ifdef PROTO_CUDA
  {
    PR_TIME("device sync");
    protoDeviceSynchronize();
  }
#endif
}
/**/

template <class T> void
applyLaplacians(int  a_nx, int a_numapplies, int a_numstream, 
                LevelData< BoxData<T> > & a_phi, 
                LevelData< BoxData<T> > & a_lap, 
                const DisjointBoxLayout& a_dbl)
{

  PR_TIME("applyLaplacians");

#ifdef PROTO_CUDA
  vector<protoStream_t> streams(a_numstream);
  for(unsigned int ibox = 0; ibox < a_numstream; ibox++)
  {
    protoStreamCreate(&streams[ibox]);
  }
#endif

  Stencil<T> emptySten;

  Stencil<T> loOrderLap = Stencil<T>::Laplacian();
#if DIM==2
  Stencil<T> hiOrderLap = Stencil<T>::Laplacian_9();
#else 
  Stencil<T> hiOrderLap = Stencil<T>::Laplacian_27();
#endif

  int maxnx = std::max(a_nx, 1);
  T dx = 1.0/maxnx;
  cout << "laplacian tests" << endl;
  for(unsigned int ibox = 0; ibox < a_dbl.size(); ibox++)
  {
    int istream = ibox%a_numstream;
    BoxData<T>& phi  = a_phi[ibox];
    BoxData<T>& lap  = a_lap[ibox];
    Box domain = a_dbl[ibox];
    //remember this is just for timings
    phi.setVal(0.);
    lap.setVal(0.);
    {
      PR_TIME("STD  laplacian with sync");
#ifdef PROTO_CUDA 
      for(int iapp = 0; iapp < a_numapplies; iapp++)
      {
        PR_TIME("actual apply");
        unsigned long long int flops;
        loOrderLap.cudaApplyStream(phi, lap, domain, true, 1.0/(dx*dx), streams[istream], flops);
      }
      sync();
#else
      for(int iapp = 0; iapp < a_numapplies; iapp++)
      {
        PR_TIME("actual apply");
        loOrderLap.apply(phi, lap, domain, true, 1.0/(dx*dx));
      }
#endif
    }

    {
      PR_TIME("DENSE  laplacian with sync");
#ifdef PROTO_CUDA 
      for(int iapp = 0; iapp < a_numapplies; iapp++)
      {
        PR_TIME("actual apply");
        unsigned long long int flops;
        hiOrderLap.cudaApplyStream(phi, lap, domain,  true, 1.0/(dx*dx), streams[istream], flops);
        PR_FLOPS(flops);
      }
      sync();
#else
      for(int iapp = 0; iapp < a_numapplies; iapp++)
      {
        PR_TIME("actual apply");
        hiOrderLap.apply(phi, lap, domain, true, 1.0/(dx*dx));
      }
#endif
    
    }

    {
      PR_TIME("empty stencil");
#ifdef PROTO_CUDA 
      for(int iapp = 0; iapp < a_numapplies; iapp++)
      {
        PR_TIME("actual apply");
        unsigned long long int flops;
        emptySten.cudaApplyStream(phi, lap, domain,  true, 1.0/(dx*dx), streams[istream], flops);
        PR_FLOPS(flops);
      }
      sync();
#else
      for(int iapp = 0; iapp < a_numapplies; iapp++)
      {
        PR_TIME("actual apply");
        //emptyKernel(a_nx);
        emptySten.apply(phi, lap, domain, true, 1.0/(dx*dx));
      }
      sync();
#endif

    } 
    {
      PR_TIME("empty kernel"); 
      for(int iapp = 0; iapp < a_numapplies; iapp++)
      { 
        PR_TIME("actual apply");
        emptyKernel(a_nx); 
      }
      sync();
    }

  }
#ifdef PROTO_CUDA 
  for(unsigned int ibox = 0; ibox < a_numstream; ibox++)
  {
    protoStreamDestroy(streams[ibox]);
  }
#endif
}
template <class T> void
applyEulerish(int  a_nx, int a_numapplies, int a_numstream,
              LevelData< BoxData<T,NUMCOMPS> > & a_U, 
              LevelData< BoxData<T,NUMCOMPS> > & a_W, 
              const DisjointBoxLayout& a_dbl)
{
  PR_TIME("applyEulerish");

  
#ifdef PROTO_CUDA
  vector<protoStream_t> streams(a_numstream);
  for(unsigned int ibox = 0; ibox < a_numstream; ibox++)
  {
    protoStreamCreate(&streams[ibox]);
  }
#endif


  Stencil<T> m_laplacian;
  Stencil<T> m_deconvolve;
  Stencil<T> m_laplacian_f[DIM];
  Stencil<T> m_deconvolve_f[DIM];
  Stencil<T> m_interp_H[DIM];
  Stencil<T> m_interp_L[DIM];
  Stencil<T> m_divergence[DIM];

  
  {
    PR_TIME("stencil definition");
    m_laplacian = Stencil<T>::Laplacian();
    m_deconvolve = ((T)(-1.0/24.0))*m_laplacian + ((T)1.0)*Shift(Point::Zeros());
    for (int dir = 0; dir < DIM; dir++)
    {
      m_laplacian_f[dir] = Stencil<T>::LaplacianFace(dir);
      m_deconvolve_f[dir] = ((T)(-1.0/24.0))*m_laplacian_f[dir] + ((T)1.0)*Shift(Point::Zeros());
      m_interp_H[dir]   = Stencil<T>::CellToEdgeH(dir);
      m_interp_L[dir]   = Stencil<T>::CellToEdgeL(dir);
      m_divergence[dir] = Stencil<T>::FluxDivergence(dir);
    }
  }

  cout << "Euler proxy stencils"<<endl;
  for(unsigned int ibox = 0; ibox < a_dbl.size(); ibox++)
  {
    BoxData<T, NUMCOMPS>& U = a_U[ibox];
    BoxData<T, NUMCOMPS>& W = a_W[ibox];
    Box domain = a_dbl[ibox];
    Box ghostBx = U.box();
    int istream = ibox%a_numstream;

    Box facedom[DIM];
    for(int idir = 0; idir < DIM; idir++)
    {
      facedom[idir] = domain.extrude(idir,1,true);
    }
    BoxData<T, NUMCOMPS> W_f[DIM]; 
    for(int idir = 0; idir < DIM; idir++)
    {
      Box faceBx = ghostBx.extrude(idir, 1, true);
      W_f[idir].define(faceBx);
    }
    {
      PR_TIME("setValZero");
      U.setVal(0.);
      W.setVal(0.);
      for(int idir = 0; idir < DIM; idir++)
      {
        W_f[idir].setVal(0.);
      }
    }

    for(int iapp = 0; iapp < a_numapplies; iapp++)
    {
      {
        PR_TIME("deconvolve");
#ifdef PROTO_CUDA
        unsigned long long int flops;
        m_deconvolve.cudaApplyStream(U, W, domain, true, 1.0, streams[istream], flops);
        PR_FLOPS(flops);
        sync();
#else
        m_deconvolve.apply(U, W, domain, true, 1.0);
#endif
      }
  
      {
        PR_TIME("laplacian");
#ifdef PROTO_CUDA
        unsigned long long int flops;
        m_laplacian.cudaApplyStream(U, W, domain, true, 1.0, streams[istream], flops);
        PR_FLOPS(flops);
        sync();
#else
        m_laplacian.apply(U, W, domain, true, 1.0);
#endif
      }

      {
        PR_TIME("interpLandH");
#ifdef PROTO_CUDA
        for(int idir = 0; idir < DIM; idir++)
        {
          PR_TIME("actual apply");
          unsigned long long int flops;
          m_interp_L[idir].cudaApplyStream(U, W_f[idir], facedom[idir], true, 1.0, streams[istream], flops);
          PR_FLOPS(flops);
          m_interp_H[idir].cudaApplyStream(U, W_f[idir], facedom[idir], true, 1.0, streams[istream], flops);
          PR_FLOPS(flops);
        }
        sync();
#else
        for(int idir = 0; idir < DIM; idir++)
        {
          PR_TIME("actual apply");
          m_interp_L[idir].apply(U, W_f[idir], facedom[idir], true, 1.0);
          m_interp_H[idir].apply(U, W_f[idir], facedom[idir], true, 1.0);
        }
#endif
      }
      {
        PR_TIME("deconvolve_f");
#ifdef PROTO_CUDA
        for(int idir = 0; idir < DIM; idir++)
        {
          PR_TIME("actual apply");
          unsigned long long int flops;
          m_deconvolve_f[idir].cudaApplyStream(U, W_f[idir], facedom[idir], true, 1.0, streams[istream], flops);
        }
        sync();
#else
        for(int idir = 0; idir < DIM; idir++)
        {
          PR_TIME("actual apply");
          m_deconvolve_f[idir].apply(U, W_f[idir], facedom[idir], true, 1.0);
        }
#endif
      }

      {
        PR_TIME("divergence");
#ifdef PROTO_CUDA
        for(int idir = 0; idir < DIM; idir++)
        {
          PR_TIME("actual apply");
          unsigned long long int flops;
          m_divergence[idir].cudaApplyStream(W_f[idir], U, domain,  true, 1.0, streams[istream], flops);
        }
        sync();
#else
        for(int idir = 0; idir < DIM; idir++)
        {
          PR_TIME("actual apply");
          m_divergence[idir].apply(W_f[idir], U, domain, true, 1.0);
        }
#endif

      }
    }
  }
  cout << "done with Euler proxy stencils"<<endl;

#ifdef PROTO_CUDA
  for(unsigned int ibox = 0; ibox < a_numstream; ibox++)
  {
    protoStreamDestroy(streams[ibox]);
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
  std::array<bool, DIM> periodic;
  for(int idir = 0; idir < DIM; idir++) periodic[idir]=true;

  DisjointBoxLayout dbl(domain, maxgrid, periodic);
  {
    PR_TIME("laplacian test");

    LevelData< BoxData<float,  1> > phif,lapf;
    LevelData< BoxData<double, 1> > phid,lapd;
    {
      PR_TIME("dataholder definition");
      phif.define(dbl, Point::Ones());
      phid.define(dbl, Point::Ones());
      lapf.define(dbl, Point::Zeros());
      lapd.define(dbl, Point::Zeros());
    }
    {
      PR_TIME("SINGLE_precision_laplacian");
      applyLaplacians<float >(nx, niter, numstreams, phif, lapf, dbl);
    }

    {
      PR_TIME("DOUBLE_precision_laplacian");
      applyLaplacians<double>(nx, niter, numstreams, phid, lapd, dbl);
    }
  }


  {
    PR_TIME("Euler_stencil_test");

    //float data
    LevelData< BoxData<float,NUMCOMPS> > Uf_c;
    LevelData< BoxData<float,NUMCOMPS> > Wf_c;


    LevelData< BoxData<double, NUMCOMPS> > Ud_c;
    LevelData< BoxData<double, NUMCOMPS> > Wd_c;
    {
      PR_TIME("dataholder definition");
      Uf_c.define(dbl, Point::Ones(4));
      Ud_c.define(dbl, Point::Ones(4));
      Wf_c.define(dbl, Point::Ones(4));
      Wd_c.define(dbl, Point::Ones(4));
    }

    {
      PR_TIME("SINGLE_precision_euler");
      applyEulerish<float>(nx, niter, numstreams, Uf_c, Wf_c,  dbl);
    }
    {
      PR_TIME("DOUBLE_precision_euler");
      applyEulerish<double>(nx, niter, numstreams, Ud_c, Wd_c, dbl);

    }
  }

  PR_TIMER_REPORT();

}  
