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

 using namespace std::chrono;

void ctoc(high_resolution_clock::time_point timer, uint iters, float unit_mem, int nrw_center, int nro_halo, int thrdim_x, int thrdim_y, int nx, int ny, int nz,
          unsigned long long int numflops, size_t numbytesread, size_t numbyteswritten)
{  
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  duration<double> time_span = duration_cast<duration<double>>(t2-timer);
  double fps = ((double)iters) / (time_span.count());
  double halo_overhead = (double)(2*thrdim_x + 2*(thrdim_y+2)*32/sizeof(double))/(double)(thrdim_x*thrdim_y);
  double effmembwd = (nrw_center+nro_halo)*unit_mem / (time_span.count()/ (double)iters);
  double hwmembwd  = (nrw_center+nro_halo+halo_overhead)*unit_mem / (time_span.count()/ (double)iters);
  double ptsthrough = (double)nx*ny*nz*iters/(double)(time_span.count());
  double mflops = ((double) numflops)/((double)time_span.count());
  double bpersecread = ((double) numbytesread)/((double)time_span.count());
  double bpersecwritten = ((double) numbyteswritten)/((double)time_span.count());
  printf("numbytesread = %d, bytespersec_read %e, bytespersec_written %e\n", numbytesread, bpersecread, bpersecwritten);

  printf("(%d, %d, %d): (TX, TY) = (%d, %d), mflop %e, fps %e, time %e, pts/s %e, effmembwd GB/s %3.1e, overhead %3.3e hwmembwd (GB/s) %3.3e)\n", 
         nx, ny, nz, thrdim_x, thrdim_y,mflops, fps, time_span.count(), ptsthrough, effmembwd/1e9, halo_overhead, hwmembwd/1e9);
}

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
inline void __cudaSafeCall(protoError err,
                           const char *file, const int line){
  if(protoSuccess != err) {
    printf("%s(%i) : cutilSafeCall() Runtime API error : %s.\n",
           file, line, protoGetErrorString(err) );
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
  int iters = 1;


  protoExtent gridExtent;
  protoPitchedPtr p_T1r, p_T2r, p_T3r;
  protoPitchedPtr p_T1u, p_T2u, p_T3u;
  protoPitchedPtr p_T1p, p_T2p, p_T3p;
  double *d_T1r, *d_T2r, *d_T3r;
  double *d_T1u, *d_T2u, *d_T3u;
  double *d_T1p, *d_T2p, *d_T3p;



  int  thrdim_x = 32, thrdim_y = 6;

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
  int pitch  = nx;
  int pitchy = ny;

  /* special case - much worse performance with default pitch on tesla1060 */
  if(nx==352) pitch = nx + 64;

  GetCmdLineArgumenti(argc, (const char**)argv, "pitch", &pitch);
  printf("pitch = %d \n", pitch);
  printf("pitchy = %d\n", pitchy);
  pitchy = pitch;
  GetCmdLineArgumenti(argc, (const char**)argv, "pitchy", &pitchy);
  /* -------------------- */
  /* Initialization */
  /* -------------------- */

  /* allocate alligned 3D data on the GPU */
  gridExtent = make_protoExtent(pitch*sizeof(double), pitchy, nz);

  cutilSafeCall(protoMalloc3D(&p_T1r, gridExtent));
  cutilSafeCall(protoMalloc3D(&p_T2r, gridExtent));
  cutilSafeCall(protoMalloc3D(&p_T3r, gridExtent));


  cutilSafeCall(protoMalloc3D(&p_T1u, gridExtent));
  cutilSafeCall(protoMalloc3D(&p_T2u, gridExtent));
  cutilSafeCall(protoMalloc3D(&p_T3u, gridExtent));

  cutilSafeCall(protoMalloc3D(&p_T1p, gridExtent));
  cutilSafeCall(protoMalloc3D(&p_T2p, gridExtent));
  cutilSafeCall(protoMalloc3D(&p_T3p, gridExtent));

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


  cutilSafeCall(protoMemset(d_T1r, 1, pitch*pitchy*nz*sizeof(double)));
  cutilSafeCall(protoMemset(d_T2r, 1, pitch*pitchy*nz*sizeof(double)));
  cutilSafeCall(protoMemset(d_T3r, 1, pitch*pitchy*nz*sizeof(double)));
  cutilSafeCall(protoMemset(d_T1u, 1, pitch*pitchy*nz*sizeof(double)));
  cutilSafeCall(protoMemset(d_T2u, 1, pitch*pitchy*nz*sizeof(double)));
  cutilSafeCall(protoMemset(d_T3u, 1, pitch*pitchy*nz*sizeof(double)));
  cutilSafeCall(protoMemset(d_T1p, 1, pitch*pitchy*nz*sizeof(double)));
  cutilSafeCall(protoMemset(d_T2p, 1, pitch*pitchy*nz*sizeof(double)));
  cutilSafeCall(protoMemset(d_T3p, 1, pitch*pitchy*nz*sizeof(double)));

  /* -------------------- */
  /* performance tests    */
  /* -------------------- */
  
  std::vector<protoStream_t> streams(numstreams);
  for(int istr = 0; istr < numstreams; istr++)
  {
    protoStreamCreate(&streams[istr]);
  }

  
  {
    PR_TIME("riemann_calls");

  high_resolution_clock::time_point timer = high_resolution_clock::now(); 
    for(int it=0; it<iters; it++)
    {
      dim3 block(thrdim_x, thrdim_y, 1);
      dim3 grid = get_grid(block, nx, ny, nz, thrdim_x, thrdim_y);
    
      //printf("kstep %d\n", kstep);
    
      int kstart = 0;
      int kstop = nz;
      
      
      forall_riemann_proxy<<<grid, block, 2*(block.x)*(block.y)*sizeof(double),streams[it%numstreams]>>>
        (d_T1r, d_T2r, d_T3r, d_T1u, d_T2u, d_T3u, d_T1p, d_T2p, d_T3p, 
         pitch, pitchy, kstart, kstop);
      
    }

    /* finalize */
    protoDeviceSynchronize();
    unsigned long long int count = 48*nx*ny*nz;
    PR_FLOPS(count);
    
    //8 is for double
    unsigned long long int numbytesread = 6*iters*nx*ny*nz*8;
    size_t numbyteswritten = iters*nx*ny*nz*sizeof(double);
    std::cout << "numbytesread = " << numbytesread << std::endl;
    ctoc(timer, iters, nx*ny*nz*sizeof(double), 1, 1, thrdim_x, thrdim_y, nx, ny, nz, count, numbytesread, numbyteswritten);   
  }
  for(int istr = 0; istr < numstreams; istr++)
  {
    protoStreamDestroy(streams[istr]);
  }
//  cutilSafeCall(protoFree(p_T1r.ptr));
//  cutilSafeCall(protoFree(p_T2r.ptr));
//  cutilSafeCall(protoFree(p_T3r.ptr));
//  cutilSafeCall(protoFree(p_T1u.ptr));
//  cutilSafeCall(protoFree(p_T2u.ptr));
//  cutilSafeCall(protoFree(p_T3u.ptr));
//  cutilSafeCall(protoFree(p_T1p.ptr));
//  cutilSafeCall(protoFree(p_T2p.ptr));
//  cutilSafeCall(protoFree(p_T3p.ptr));
  return 0;
}
int main(int argc, char* argv[]) 
{
  PR_TIMER_SETFILE("proto.time.table");
  printf("usage: baroqueforallexe -nx nx -ny ny -nz nz -iters iters -pitch pitch x -pitchy pitchy -ns num_streams\n");

  runTest(argc, argv);
  
  PR_TIMER_REPORT();
  return 0;
}
