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
                                     double *pout, double *phig, double *plow)
{
  const  int ix = blockIdx.x*blockDim.x + threadIdx.x;	
  double gamma = 1.4;
  double& rhoo = rout[ix];
  double& rhol = rlow[ix];
  double& rhor = rhig[ix];
  double& uo   = uout[ix];
  double& ul =   ulow[ix];
  double& ur =   uhig[ix];
  double& po   = pout[ix];
  double& pl =   plow[ix];
  double& pr =   phig[ix];


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
  
  int nx = 256;
  int ny, nz;
  int iters = 512;


  double *d_T1r, *d_T2r, *d_T3r;
  double *d_T1u, *d_T2u, *d_T3u;
  double *d_T1p, *d_T2p, *d_T3p;



  /* -------------------- */
  /* command-line parameters */
  /* -------------------- */
  GetCmdLineArgumenti(argc, (const char**)argv, "nx", &nx);
  ny = nx;
  nz = nx;
  printf("(nx, ny, nz)= (%d, %d, %d)\n", nx, ny, nz);
  GetCmdLineArgumenti(argc, (const char**)argv, "ny", &ny);
  GetCmdLineArgumenti(argc, (const char**)argv, "nz", &nz);
  GetCmdLineArgumenti(argc, (const char**)argv, "iters", &iters);
  printf("num_iters = %d\n", iters);
  int numstreams = 8;
  GetCmdLineArgumenti(argc, (const char**)argv, "ns", &numstreams);
  printf("numstreams = %d\n", numstreams);

  /* -------------------- */
  /* Initialization */
  /* -------------------- */
  size_t memsize = nx*ny*nz*sizeof(double);
  cutilSafeCall(cudaMalloc(&d_T1r, memsize));
  cutilSafeCall(cudaMalloc(&d_T2r, memsize));
  cutilSafeCall(cudaMalloc(&d_T3r, memsize));
  cutilSafeCall(cudaMalloc(&d_T1u, memsize));
  cutilSafeCall(cudaMalloc(&d_T2u, memsize));
  cutilSafeCall(cudaMalloc(&d_T3u, memsize));
  cutilSafeCall(cudaMalloc(&d_T1p, memsize));
  cutilSafeCall(cudaMalloc(&d_T2p, memsize));
  cutilSafeCall(cudaMalloc(&d_T3p, memsize));

  cutilSafeCall(cudaMemset(d_T1r, 1, memsize));
  cutilSafeCall(cudaMemset(d_T2r, 1, memsize));
  cutilSafeCall(cudaMemset(d_T3r, 1, memsize));
  cutilSafeCall(cudaMemset(d_T1u, 1, memsize));
  cutilSafeCall(cudaMemset(d_T2u, 1, memsize));
  cutilSafeCall(cudaMemset(d_T3u, 1, memsize));
  cutilSafeCall(cudaMemset(d_T1p, 1, memsize));
  cutilSafeCall(cudaMemset(d_T2p, 1, memsize));
  cutilSafeCall(cudaMemset(d_T3p, 1, memsize));

  /* -------------------- */
  /* performance tests    */
  /* -------------------- */
  
  std::vector<cudaStream_t> streams(numstreams);
  for(int istr = 0; istr < numstreams; istr++)
  {
    cudaStreamCreate(&streams[istr]);
  }

  
  {
    PR_TIME("riemann_calls");

    for(int it=0; it<iters; it++)
    {
      //printf("kstep %d\n", kstep);
    
      
      int stride = nx;
      int blocks = ny*nz;
      size_t smem = 0;
      forall_riemann_proxy<<<blocks, stride, smem, streams[it%numstreams]>>>
        (d_T1r, d_T2r, d_T3r, d_T1u, d_T2u, d_T3u, d_T1p, d_T2p, d_T3p);
      
    }

    /* finalize */
    cudaDeviceSynchronize();
    unsigned long long int count = 48*nx*ny*nz;
    PR_FLOPS(count);
    
  }
  for(int istr = 0; istr < numstreams; istr++)
  {
    cudaStreamDestroy(streams[istr]);
  }
  cutilSafeCall(cudaFree(d_T1r));
  cutilSafeCall(cudaFree(d_T2r));
  cutilSafeCall(cudaFree(d_T3r));
  cutilSafeCall(cudaFree(d_T1u));
  cutilSafeCall(cudaFree(d_T2u));
  cutilSafeCall(cudaFree(d_T3u));
  cutilSafeCall(cudaFree(d_T1p));
  cutilSafeCall(cudaFree(d_T2p));
                cutilSafeCall(cudaFree(d_T3p));
  return 0;
}
int main(int argc, char* argv[]) 
{
  PR_TIMER_SETFILE("proto.time.table");
  printf("usage: baroqueforallexe -nx nx -ny ny -nz nz -iters iters  -ns num_streams\n");

  runTest(argc, argv);
  
  PR_TIMER_REPORT();
  return 0;
}
