
#include <cstdlib>
#include <cstdio>
#include <functional>
#include <iostream>
#include <vector>
#include "../../include/Proto_Timer.H"
#include "../../include/Proto.H"

/* stencil apply  header material ============================ */
using std::cout;
using std::endl;
using std::vector;

constexpr unsigned int NUMCOMPS=3;
constexpr unsigned int NUMMULTS=10;

/**/
void
parseCommandLine(unsigned int & a_nx,     unsigned int & a_napplies,
                 unsigned int & a_nboxes, unsigned int & a_nstreams, 
                 int argc, char* argv[])
{
  //defaults
  a_nx        = 2048;
  a_napplies  = 100;
  a_nboxes    = 64;
  a_nstreams  = 8;
  cout << "simple forall breadboard implementation--[] indicates default values" << endl;
  cout << "usage:  " << argv[0] << " -n nx[2048] -a num_applies[100] -b num_boxes[64] -s num_streams[8]" << endl;
  for(int iarg = 0; iarg < argc-1; iarg++)
  {
    if(strcmp(argv[iarg],"-n") == 0)
    {
      a_nx = atoi(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg], "-a") == 0)
    {
      a_napplies = atoi(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg], "-b") == 0)
    {
      a_nboxes = atoi(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg], "-s") == 0)
    {
      a_nstreams = atoi(argv[iarg+1]);
    }
    cout << "nx          = " << a_nx       << endl;
    cout << "napplies    = " << a_napplies << endl;
    cout << "nstreams    = " << a_nstreams << endl;
    cout << "nboxes      = " << a_nboxes   << endl;
  }
}
/**/
__device__
void riemannProxy(double&        a_rout,
                  double&        a_uout,
                  double&        a_pout,
                  const double&  a_rlow,
                  const double&  a_ulow,
                  const double&  a_plow,
                  const double&  a_rhig,
                  const double&  a_uhig,
                  const double&  a_phig)
{
  constexpr double gamma = 1.4;

  const double& rhol = a_rlow;
  const double& rhor = a_rhig;
  const double& ul   = a_ulow;
  const double& ur   = a_uhig;
  const double& pl   = a_plow;
  const double& pr   = a_phig;

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
  for(int iter = 0; iter < NUMMULTS; iter++)
  {
    cbar *= pbar;
  }
  //7
  double pstar = (pl + pr)*.5 + rhobar*cbar*(ul - ur)*.5;
  //7
  double ustar = (ul + ur)*.5 + (pl - pr)/(2*rhobar*cbar);

  //took out conditionals to make flop count more honest

  a_rout = a_rlow;
  a_uout = a_ulow;
  a_pout = a_plow;

  //5
  a_rout += (pstar - a_pout)/(cbar*cbar);
  a_uout = ustar + ubar;
  a_pout = pstar;

  //I get 27 + NMULT
}
/**/
__global__
void proxyIndexer(double*       a_out, 
                  double*       a_low,
                  double*       a_hig, 
                  unsigned int  a_nx)
{
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  if(idx < a_nx)
  {
    //variable order = rho, u, p
    unsigned int rind = idx;
    unsigned int uind = idx +   a_nx;
    unsigned int pind = idx + 2*a_nx;
    riemannProxy(a_out[rind],a_out[uind],a_out[pind],
                 a_low[rind],a_low[uind],a_low[pind],
                 a_hig[rind],a_hig[uind],a_hig[pind]);
  }
}
/**/
void forallProxy(double*       a_outptr, 
                 double*       a_lowptr,
                 double*       a_higptr, 
                 unsigned int  a_nx, 
                 cudaStream_t& a_stream)
{
  unsigned int stride = 64;
  unsigned int blocks = a_nx/64;
  size_t smem = 0;
  proxyIndexer<<<blocks, stride, smem, a_stream>>>(a_outptr, a_lowptr, a_higptr, a_nx);
}

inline void sync()
{
  {
    PR_TIME("device sync");
    cudaDeviceSynchronize();
  }
}
/**/
void runTest(unsigned int a_nx    , unsigned int a_napplies, 
             unsigned int a_nboxes, unsigned int a_nstreams)
{
  PR_TIME("test_function");
  //This is the data for *all* the boxes with one allocation per dataholder.
  //Yes, I could have put all this into one horrible pointer.
  //I chose to make this at least somewhat readable.
  double* low;
  double* hig;
  double* out;

  cout << "allocating memory" << endl;
  {
    PR_TIME("cudamalloc");
    size_t datsize = a_nboxes*NUMCOMPS*a_nx*sizeof(double);
    cudaMalloc(&low, datsize);
    cudaMalloc(&hig, datsize);
    cudaMalloc(&out, datsize);
  }

  cout << "setting values" << endl;
  {
    PR_TIME("thrust::fill (setval)");
    double val = 1.0;
    unsigned int len = a_nboxes*a_nx*NUMCOMPS;
    thrust::device_ptr<double> outptr = thrust::device_pointer_cast(out);
    thrust::device_ptr<double> lowptr = thrust::device_pointer_cast(low);
    thrust::device_ptr<double> higptr = thrust::device_pointer_cast(hig);
    thrust::fill(thrust::device, outptr, outptr+len, val);
    thrust::fill(thrust::device, lowptr, lowptr+len, val);
    thrust::fill(thrust::device, higptr, higptr+len, val);
  }

  
  cout << "making streams" << endl;
  vector<cudaStream_t> streams(a_nstreams);
  {
    PR_TIME("stream_create");
    for(unsigned int ibox = 0; ibox < a_nstreams; ibox++)
    {
      cudaStreamCreate(&streams[ibox]);
    }
  }
  
  cout << "doing forall proxy" << endl;
  for(unsigned int iter = 0; iter  < a_napplies; iter++)
  {
    PR_TIME("forallproxy");

    for(unsigned int ibox = 0; ibox < a_nstreams; ibox++)
    {
      int istream = ibox%a_nboxes;
      double* outptr = out + ibox*a_nx;
      double* lowptr = low + ibox*a_nx;
      double* higptr = hig + ibox*a_nx;
      forallProxy(outptr, lowptr, higptr, a_nx, streams[istream]);
      unsigned long long int nflop = (27+NUMMULTS)*a_nx;
      PR_FLOPS(nflop);
    }
    sync();
  }

  cout << "freeing memory" << endl;
  {
    PR_TIME("cudafree");
    cudaFree(low);
    cudaFree(hig);
    cudaFree(out);
  }

  cout << "destroying streams" << endl;
  {
    PR_TIME("cudaStreamDestroy");
    for(unsigned int ibox = 0; ibox < a_nstreams; ibox++)
    {
      cudaStreamDestroy(streams[ibox]);
    }
  }
}
/**/
int main(int argc, char* argv[]) 
{
  PR_TIMER_SETFILE("proto.time.table");
  unsigned int nx, napplies, nboxes, nstreams;
  parseCommandLine(nx,  napplies, nboxes, nstreams,argc, argv);

  runTest(nx, napplies, nboxes, nstreams);
  

  PR_TIMER_REPORT();
  
  return 0;
}
