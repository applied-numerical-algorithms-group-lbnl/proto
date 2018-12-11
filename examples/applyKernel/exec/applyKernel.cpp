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
constexpr unsigned int NUMCOMPS=DIM+2;

/**/
void
parseCommandLine(int & a_nx, int & a_numapplies, int argc, char* argv[])
{
  //defaults
  a_nx = 8;
  a_numapplies = 1;
  cout << "kernel timings of various laplacians" << endl;
  cout << "usage:  " << argv[0] << " -n nx[default:8] -m num_iterations[default:1]" << endl;
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
      cudaDeviceSynchronize();
    }
#endif
}
/**/

template <class T> void
applyLaplacians(int  a_nx, int a_numapplies, BoxData<T>& phi, BoxData<T>& lap, Box domain, Box ghostBx)
{

  PR_TIME("applyLaplacians");

  Stencil<T> emptySten;

  Stencil<T> loOrderLap = Stencil<T>::Laplacian();
#if DIM==2
  Stencil<T> hiOrderLap = Stencil<T>::Laplacian_9();
#else 
  Stencil<T> hiOrderLap = Stencil<T>::Laplacian_27();
#endif

  
  
  //remember this is just for timings
  phi.setVal(0.);
  lap.setVal(0.);
  int maxnx = std::max(a_nx, 1);
  T dx = 1.0/maxnx;
  cout << "apply standard laplacian " << a_numapplies << " times" << endl;
  {
    PR_TIME("STD  laplacian with sync");
    for(int iapp = 0; iapp < a_numapplies; iapp++)
    {
      PR_TIME("actual apply");
      loOrderLap.apply(phi, lap, domain, true, 1.0/(dx*dx));
    }
    sync();
  }
  cout << "apply dense laplacian " << a_numapplies << " times" << endl;
  {
    PR_TIME("DENSE  laplacian with sync");
    for(int iapp = 0; iapp < a_numapplies; iapp++)
    {
      PR_TIME("actual apply");
      hiOrderLap.apply(phi, lap, domain, true, 1.0/(dx*dx));
    }
    sync();
    
  }

  cout<<" empty stencil launches"<<endl;
  {
    PR_TIME("empty stencil");
    for(int iapp = 0; iapp < a_numapplies; iapp++)
      {
        PR_TIME("actual apply");
        //emptyKernel(a_nx);
        emptySten.apply(phi, lap, domain, true, 1.0/(dx*dx));
      }
    sync();
    
  } 
  cout << "actual empty kernel launches"<<endl;
  {
   PR_TIME("empty kernel"); 
    for(int iapp = 0; iapp < a_numapplies; iapp++)
      { 
        emptyKernel(a_nx); 
      }
    sync();
  }

}
template <class T> void
applyEulerish(int  a_nx, int a_numapplies, BoxData<T,NUMCOMPS>& U, BoxData<T,NUMCOMPS>& W, BoxData<T,NUMCOMPS> W_f[DIM], Box domain, Box ghostBx)
{
  PR_TIME("applyEulerish");

  
  cout << "Euler proxy stencils"<<endl;
  Box facedom[DIM];
  for(int idir = 0; idir < DIM; idir++)
  {
    facedom[idir] = domain.extrude(idir,1,true);
  }

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
      m_deconvolve.apply(U, W, domain, true, 1.0);
    }
  
    {
      PR_TIME("laplacian");
      m_laplacian.apply(U, W, domain, true, 1.0);
    }
    for(int idir = 0; idir < DIM; idir++)
    {
      PR_TIME("interpLandH");
      m_interp_L[idir].apply(U, W_f[idir], facedom[idir], true, 1.0);
      m_interp_H[idir].apply(U, W_f[idir], facedom[idir], true, 1.0);
    }

    for(int idir = 0; idir < DIM; idir++)
    {
      PR_TIME("deconvolve_f");
      m_deconvolve_f[idir].apply(U, W_f[idir], facedom[idir], true, 1.0);
    }
    for(int idir = 0; idir < DIM; idir++)
    {
      PR_TIME("divergence");
      m_divergence[idir].apply(W_f[idir], U, domain, true, 1.0);
    }
  }
  cout << "done with Euler proxy stencils"<<endl;

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
    PR_TIME("laplacian test");

    Point ghostPt = Point::Ones();
    Box   ghostBx = domain.grow(ghostPt);
    BoxData<float,  1> phif,lapf;
    BoxData<double, 1> phid,lapd;
    {
      PR_TIME("dataholder definition");
      phif.define(ghostBx);
      lapf.define(domain);
      phid.define(ghostBx);
      lapd.define(domain);
    }
    {
      PR_TIME("SINGLE_precision_laplacian");
      applyLaplacians<float>(nx, niter, phif, lapf, domain, ghostBx);
    }

    {
      PR_TIME("DOUBLE_precision_laplacian");
      applyLaplacians<double>(nx, niter, phid, lapd, domain, ghostBx);
    }
  }


  {
    PR_TIME("Euler_stencil_test");

    Point ghostPt = Point::Ones(4);
    Box   ghostBx = domain.grow(ghostPt);

    //float data
    BoxData<float,NUMCOMPS> Uf_c;
    BoxData<float,NUMCOMPS> Wf_c;
    BoxData<float,  NUMCOMPS> Uf_f[DIM];


    BoxData<double, NUMCOMPS> Ud_c;
    BoxData<double, NUMCOMPS> Wd_c;
    BoxData<double, NUMCOMPS> Ud_f[DIM];
    {
      PR_TIME("dataholder definition");
      Uf_c.define(ghostBx);
      Ud_c.define(ghostBx);
      Wf_c.define(ghostBx);
      Wd_c.define(ghostBx);
      for(int idir = 0; idir < DIM; idir++)
      {
        Box faceBx = ghostBx.extrude(idir, 1, true);
        Uf_f[idir].define(faceBx);
        Ud_f[idir].define(faceBx);
      }
    }

    {
      PR_TIME("SINGLE_precision_euler");

      applyEulerish<float>(nx, niter, Uf_c, Wf_c, Uf_f, domain, ghostBx);
    }
    {
      PR_TIME("DOUBLE_precision_euler");
      applyEulerish<double>(nx, niter, Ud_c, Wd_c, Ud_f, domain, ghostBx);

    }
  }

  PR_TIMER_REPORT();

}  
