/*
  Copyright Marcin Krotkiewski, University of Oslo, 2012


  updated by Brian Van Straalen 2018
   updated for the loss of cutil.h
   updated to use 6 streams of execution
   updated to use std::chrono

   updated by dtg 2018
   added timers and forced double precision
   moved bulk of computation out of main for purely aesthetic reasons
   added some minimal documentation of command line options 
   changed from stencil calc to forall calc
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <time.h>
#include <errno.h>
#include <chrono>
#include <algorithm>
#include <iostream>

//#include <cutil_inline.h>
#include <cuda_runtime_api.h>
#include <vector_types.h>
#include <vector_functions.h>
#include "../../include/Proto_Timer.H"

__global__ void forall_riemann_proxy(double *rout, double *rhig, double *rlow,
                                     double *uout, double *uhig, double *ulow,
                                     double *pout, double *phig, double *plow,
                                     uint pitch, uint pitchy,
                                     uint kstart, uint kend)
{
  const  int ix = blockIdx.x*blockDim.x + threadIdx.x;	
  const  int iy = blockIdx.y*blockDim.y + threadIdx.y;
  double gamma = 1.4;

  for(unsigned int kk=kstart; kk<kend; kk++)
  {

    double& rhoo = rout[ix + iy*pitch + kk*pitch*pitchy];
    double& rhol = rlow[ix + iy*pitch + kk*pitch*pitchy];
    double& rhor = rhig[ix + iy*pitch + kk*pitch*pitchy];
    double& uo   = uout[ix + iy*pitch + kk*pitch*pitchy];
    double& ul =   ulow[ix + iy*pitch + kk*pitch*pitchy];
    double& ur =   uhig[ix + iy*pitch + kk*pitch*pitchy];
    double& po   = pout[ix + iy*pitch + kk*pitch*pitchy];
    double& pl =   plow[ix + iy*pitch + kk*pitch*pitchy];
    double& pr =   phig[ix + iy*pitch + kk*pitch*pitchy];


    //2
    double rhobar = (rhol + rhor)*.5;
    //2
    double pbar = (pl + pr)*.5;
    //2
    double ubar = (ul + ur)*.5;
    //took this one out for a bunch of multiplies so
    //I can have flops I can count
//  double cbar = sqrt(gamma*pbar/rhobar);
    //2
    double cbar = gamma*pbar/rhobar;
    //NMULT
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;
    cbar *= pbar;

    //7
    double pstar = (pl + pr)*.5 + rhobar*cbar*(ul - ur)*.5;
    //7
    double ustar = (ubar + ur)*.5 + (pl - pr)/(2*rhobar*cbar);

    rhoo  = rhol;
    uo    =   ul;
    po    =   pl;
    //4
    rhoo += (pstar - po)/(cbar*cbar);
    uo    = ustar;
    po = pstar;

  }
}


#define cutilSafeCall(err)     __cudaSafeCall   (err, __FILE__, __LINE__)
#define cutilCheckError(err)   __cutilCheckError(err, __FILE__, __LINE__)
inline void __cudaSafeCall(cudaError err,
                           const char *file, const int line){
  if(cudaSuccess != err) {
    printf("%s(%i) : cutilSafeCall() Runtime API error : %s.\n",
           file, line, cudaGetErrorString(err) );
    exit(-1);
  }
}
inline void __cutilCheckError( bool err, const char *file, const int line )
{
    if( true != err) {
        fprintf(stderr, "CUTIL CUDA error in file <%s>, line %i.\n",
                file, line);
        exit(-1);
    }
}




dim3 get_grid(dim3 block, int nx, int ny, int nz, int thrdim_x, int thrdim_y)
{
  int modx = nx%thrdim_x;
  int mody = ny%thrdim_y;
  
  dim3 grid(nx/block.x, ny/block.y, 1);
  if(modx){
    grid.x++;
  }
  if(mody){
    grid.y++;
  }
  return grid;
}

void GetCmdLineArgumenti(int argc, const char** argv, const char* name, int* rtn)
{
  size_t len = strlen(name);
  for(int i=1; i<argc; i+=2)
    {
      if(strcmp(argv[i]+1,name) ==0)
        {
         *rtn = atoi(argv[i+1]);
         std::cout<<name<<"="<<" "<<*rtn<<std::endl;
          break;
        }
    }
}

int runTest(int argc, char*argv[])
{
  
  int device = 0;
  int nx = 32;
  int ny, nz;
  int iters = 100;


  cudaExtent gridExtent;
  cudaPitchedPtr p_T1r, p_T2r, p_T3r;
  cudaPitchedPtr p_T1u, p_T2u, p_T3u;
  cudaPitchedPtr p_T1p, p_T2p, p_T3p;
  double *d_T1r, *d_T2r, *d_T3r;
  double *d_T1u, *d_T2u, *d_T3u;
  double *d_T1p, *d_T2p, *d_T3p;



  int  thrdim_x = 32, thrdim_y = 6;
  int texsize;

  /* -------------------- */
  /* command-line parameters */
  /* -------------------- */
  GetCmdLineArgumenti(argc, (const char**)argv, "nx", &nx);
  ny = nx;
  nz = nx;
  printf("(nx, ny, nz)= (%d, %d, %d)\n", nx, ny, nz);
  GetCmdLineArgumenti(argc, (const char**)argv, "ny", &ny);
  GetCmdLineArgumenti(argc, (const char**)argv, "nz", &nz);
  GetCmdLineArgumenti(argc, (const char**)argv, "device", &device);
  GetCmdLineArgumenti(argc, (const char**)argv, "iters", &iters);

  /* choose device */
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, device);
  cudaSetDevice(device);
  if(strstr(deviceProp.name, "1060")){
    texsize = 22;
  } else {
    texsize = 24;
  }
  GetCmdLineArgumenti(argc, (const char**)argv, "texsize", &texsize);

  int pitch  = nx;
  int pitchy = ny;

  /* special case - much worse performance with default pitch on tesla1060 */
  if(nx==352) pitch = nx + 64;

  GetCmdLineArgumenti(argc, (const char**)argv, "pitch", &pitch);
  GetCmdLineArgumenti(argc, (const char**)argv, "pitchy", &pitchy);
  printf("pitch = %d \n", pitch);
  printf("pitchy = %d\n", pitchy);
  /* -------------------- */
  /* Initialization */
  /* -------------------- */

  /* allocate alligned 3D data on the GPU */
  gridExtent = make_cudaExtent(pitch*sizeof(double), pitchy, nz);

  cutilSafeCall(cudaMalloc3D(&p_T1r, gridExtent));
  cutilSafeCall(cudaMalloc3D(&p_T2r, gridExtent));
  cutilSafeCall(cudaMalloc3D(&p_T3r, gridExtent));


  cutilSafeCall(cudaMalloc3D(&p_T1u, gridExtent));
  cutilSafeCall(cudaMalloc3D(&p_T2u, gridExtent));
  cutilSafeCall(cudaMalloc3D(&p_T3u, gridExtent));

  cutilSafeCall(cudaMalloc3D(&p_T1p, gridExtent));
  cutilSafeCall(cudaMalloc3D(&p_T2p, gridExtent));
  cutilSafeCall(cudaMalloc3D(&p_T3p, gridExtent));

  d_T1r  = (double*)p_T1r.ptr;
  d_T2r  = (double*)p_T2r.ptr;
  d_T3r  = (double*)p_T3r.ptr;


  d_T1u  = (double*)p_T1u.ptr;
  d_T2u  = (double*)p_T2u.ptr;
  d_T3u  = (double*)p_T3u.ptr;

  d_T1p  = (double*)p_T1p.ptr;
  d_T2p  = (double*)p_T2p.ptr;
  d_T3p  = (double*)p_T3p.ptr;

  pitch = p_T1r.pitch/sizeof(double);
  printf("pitch %li, xsize %li, ysize %li\n", p_T1r.pitch/sizeof(double), p_T1r.xsize/sizeof(double), p_T1r.ysize);

  cutilSafeCall(cudaMemset(d_T1r, 1, pitch*pitchy*nz*sizeof(double)));
  cutilSafeCall(cudaMemset(d_T2r, 1, pitch*pitchy*nz*sizeof(double)));
  cutilSafeCall(cudaMemset(d_T3r, 1, pitch*pitchy*nz*sizeof(double)));
  cutilSafeCall(cudaMemset(d_T1u, 1, pitch*pitchy*nz*sizeof(double)));
  cutilSafeCall(cudaMemset(d_T2u, 1, pitch*pitchy*nz*sizeof(double)));
  cutilSafeCall(cudaMemset(d_T3u, 1, pitch*pitchy*nz*sizeof(double)));
  cutilSafeCall(cudaMemset(d_T1p, 1, pitch*pitchy*nz*sizeof(double)));
  cutilSafeCall(cudaMemset(d_T2p, 1, pitch*pitchy*nz*sizeof(double)));
  cutilSafeCall(cudaMemset(d_T3p, 1, pitch*pitchy*nz*sizeof(double)));

  /* -------------------- */
  /* performance tests    */
  /* -------------------- */
  
  cudaStream_t streams[6];
  cudaStreamCreate(streams);
  cudaStreamCreate(streams+1);
  cudaStreamCreate(streams+2);
  cudaStreamCreate(streams+3);
  cudaStreamCreate(streams+4);
  cudaStreamCreate(streams+5);
  
  {
    PR_TIME("riemann_calls");

    for(int it=0; it<iters; it++)
    {
      dim3 block(thrdim_x, thrdim_y, 1);
      dim3 grid = get_grid(block, nx, ny, nz, thrdim_x, thrdim_y);
    
      int kstep  = std::min((1<<texsize)/(pitch*pitchy), nz);
      //printf("kstep %d\n", kstep);
    
      int kstart = 0;
      int kstop = nz;
      
      
      forall_riemann_proxy<<<grid, block, 2*(block.x)*(block.y)*sizeof(double),streams[it%6]>>>
        (d_T1r, d_T2r, d_T3r, d_T1u, d_T2u, d_T3u, d_T1p, d_T2p, d_T3p, 
         pitch, pitchy, kstart, kstop);
      
    }

    /* finalize */
    cudaDeviceSynchronize();
    unsigned long long int count = 48*nx*ny*nz;
    PR_FLOPS(count);
    
  }
  return 0;
}
int main(int argc, char* argv[]) 
{
  PR_TIMER_SETFILE("proto.time.table");
  printf("usage: baroqueforallexe -nx nx -ny ny -nz nz -routine routine -device device -iters iters -pitch pitch x -pitchy pitchy -texsize texsize");

  runTest(argc, argv);

  
  PR_TIMER_REPORT();
  return 0;
}
