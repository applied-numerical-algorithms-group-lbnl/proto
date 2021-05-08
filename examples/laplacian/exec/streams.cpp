
#include "Proto.H"
#include "implem/Proto_LevelData.H"

#ifdef PROTO_BRICK
#include "brick.h"
#include "bricksetup.h"
#include "multiarray.h"

namespace {

#if DIM == 3
#define SCRIPT "laplacian3d.py"
#else
#define SCRIPT "laplacian2d.py"
#endif

  void optimized_laplacian(Proto_brick &src,
                           Proto_brick &dest,
                           int brick_idx,
                           bool initToZero,
                           vector<double> &coefs) {
    brick(SCRIPT, BVEC, (BDIM), (BFOLD), brick_idx);
  }
}
#endif

namespace Proto
{
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
  
  void testFunc(int argc, char* argv[])
  {
      int nx = 256;
    int ny = 256;
    int nz = 256;
    int maxbox = 64;
    int niters = 10;

    /* -------------------- */
    /* command-line parameters */
    /* -------------------- */
    GetCmdLineArgumenti(argc, (const char**)argv, "nx", &nx);
    ny = nx;
    nz = nx;
    maxbox = 64;
    GetCmdLineArgumenti(argc, (const char**)argv, "ny", &ny);
    GetCmdLineArgumenti(argc, (const char**)argv, "nz", &nz);
    GetCmdLineArgumenti(argc, (const char**)argv, "maxbox", &maxbox);
    GetCmdLineArgumenti(argc, (const char**)argv, "niters", &niters);
    int nstream = 8;
#ifdef PROTO_CUDA
    GetCmdLineArgumenti(argc, (const char**)argv, "nstream", &nstream);
    Proto::DisjointBoxLayout::setNumStreams(nstream);
#endif
  
#if DIM==3
    static Stencil<double> sten = Stencil<double>::Laplacian_27();
#else
    static Stencil<double> sten = Stencil<double>::Laplacian();
#endif
#if defined(PROTO_BRICK) && !defined(PROTO_BRICK_JIT)
    sten.host_optimized = (void*)&optimized_laplacian;
#endif
    printf("nx = %d, ny = %d, nz= %d\n", nx, ny, nx);
    printf("maxbox = %d, niters = %d, nstream = %d\n", maxbox, niters, nstream);
    Box domain(Point::Zeros(), Point::Ones(nx-1));
    std::array<bool, DIM> periodic;
    for(int idir = 0; idir < DIM; idir++) periodic[idir]=true;
    DisjointBoxLayout   dbl(domain, maxbox, periodic);

    // Extended Dbl should use the same type of layout. However, extended /= non-extended
    LevelData<BoxData<double, 2>> phild(dbl, Point::Unit());
    LevelData<BoxData<double, 2>> lphld(dbl, Point::Zero());

    for(unsigned int i=0; i<dbl.size(); i++)
      {
        auto& phi = phild[i];
        phi.setVal(1.5);
      }
    {
      PR_TIME("apply_laplacian_current");
      for(unsigned int iter = 0; iter < niters; iter++)
        {
          for(unsigned int i=0; i<dbl.size(); i++)
            {
              
              auto& phi = phild[i];
              auto& lph = lphld[i];
              sten.apply(phi, lph, dbl[i], true);
              sten.apply(phi, lph, dbl[i], true);
            } 
        }
#ifdef PROTO_CUDA    
      protoDeviceSynchronize();
      protoError err = protoGetLastError();
      if (err != protoSuccess)
        {
          fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                  __FILE__, __LINE__, protoGetErrorString(err));
        }
#endif    
    }
  
    {
#ifdef PROTO_CUDA
    PR_TIME("apply_laplacian_update");
    for(unsigned int iter = 0; iter < niters; iter++)
      {
        for(unsigned int i=0; i<dbl.size(); i++)
          {
            
            auto& phi = phild[i];
            auto& lph = lphld[i];
            sten.cudaApply2(phi, lph, dbl[i], true, 1.0);
            sten.cudaApply2(phi, lph, dbl[i], true, 1.0);
          }
        
      }
#endif

#ifdef PROTO_CUDA    
    protoDeviceSynchronize();
    protoError err = protoGetLastError();
    if (err != protoSuccess)
      {
        fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                __FILE__, __LINE__, protoGetErrorString(err));
      }
#endif    
    }
  }

}
int main(int argc, char* argv[])
{
  
  //have to do this to get a time table
  PR_TIMER_SETFILE("proto.time.table");
  Proto::testFunc(argc, argv);
  PR_TIMER_REPORT();
  return 0;
}
  
      
