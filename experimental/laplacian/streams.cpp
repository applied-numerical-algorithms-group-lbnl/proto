
#include "Proto.H"
#include "Proto_LevelData.H"
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
  
    printf("nx = %d, ny = %d, nz= %d\n", nx, ny, nx);
    printf("maxbox = %d, niters = %d, nstream = %d\n", maxbox, niters, nstream);
    Box domain(Point::Zeros(), Point::Ones(nx-1));
    std::array<bool, DIM> periodic;
    for(int idir = 0; idir < DIM; idir++) periodic[idir]=true;
    DisjointBoxLayout   dbl(domain, maxbox, periodic);

    LevelData<BoxData<double, 1>> phild(dbl, Point::Unit());
    LevelData<BoxData<double, 1>> lphld(dbl, Point::Zero());

    for(unsigned int i=0; i<dbl.size(); i++)
    {
      BoxData<double>& phi = phild[i];
      BoxData<double>& lph = lphld[i];
      phi.setVal(0.);
      lph.setVal(0.);
    }
#if DIM==3
    Stencil<double> sten = Stencil<double>::Laplacian_27();
#else
    Stencil<double> sten = Stencil<double>::Laplacian();
#endif

    for(unsigned int iter = 0; iter < niters; iter++)
    {
      PR_TIME("apply_laplacian");
      for(unsigned int i=0; i<dbl.size(); i++)
      {
        BoxData<double>& phi = phild[i];
        BoxData<double>& lph = lphld[i];
        sten.apply(phi, lph, dbl[i], true);
      }
#ifdef PROTO_CUDA    
      cudaDeviceSynchronize();
      cudaError err = cudaGetLastError();
      if (err != cudaSuccess)
      {
        fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                __FILE__, __LINE__, cudaGetErrorString(err));
      }
#endif    
    }

    printf("out of loop --- writing report\n");
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
  
      
