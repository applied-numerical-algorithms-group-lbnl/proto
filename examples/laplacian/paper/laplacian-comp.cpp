
#include "Proto.H"
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
    int maxbox = 256;
    int niters = 10;

    /* -------------------- */
    /* command-line parameters */
    /* -------------------- */
    GetCmdLineArgumenti(argc, (const char**)argv, "nx", &nx);
    int ny = nx;
    int nz = nx;
    GetCmdLineArgumenti(argc, (const char**)argv, "ny", &ny);
    GetCmdLineArgumenti(argc, (const char**)argv, "nz", &nz);
    GetCmdLineArgumenti(argc, (const char**)argv, "maxbox", &maxbox);
    GetCmdLineArgumenti(argc, (const char**)argv, "niters", &niters);
    int nstream = 1;
#ifdef PROTO_CUDA
    GetCmdLineArgumenti(argc, (const char**)argv, "nstream", &nstream);
#endif
  
#if DIM==3
    static Stencil<double> sten = Stencil<double>::Laplacian_27();
#else
    static Stencil<double> sten = Stencil<double>::Laplacian();
#endif
    printf("nx = %d, ny = %d, nz= %d\n", nx, ny, nx);
    printf("maxbox = %d, niters = %d, nstream = %d\n", maxbox, niters, nstream);
    Box domain(Point::Zeros(), Point::Ones(nx-1));
    std::array<bool, DIM> periodic;
    for(int idir = 0; idir < DIM; idir++) periodic[idir]=true;
    ProblemDomain probDom(domain,periodic);
    DisjointBoxLayout dbl(probDom, maxbox*Point::Ones()); //This will create a disjoint layout with maxgrid size boxes

    constexpr unsigned int comp = 2;
    LevelBoxData<double, comp> phild, lphld;
#ifdef PROTO_CUDA
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
#endif
    {
      PR_TIME("dataholder definition");
      phild.define(dbl, Point::Unit());
      lphld.define(dbl, Point::Zero());
    }

    for(DataIterator dit = phild.begin(); (*dit)!=dit.end(); ++dit)
    {
      auto& phi = phild[*dit];
      phi.setVal(1.5);
    }
    {
      PR_TIME("apply_laplacian_current");
      for(unsigned int iter = 0; iter < niters; iter++)
        {
          for(DataIterator dit = phild.begin(); (*dit)!=dit.end(); ++dit)
            {
              
              auto& phi = phild[*dit];
              auto& lph = lphld[*dit];
#ifdef PROTO_CUDA
	      cudaEventRecord(start);
#endif
	      sten.apply(phi, lph, dbl[*dit], true);
#ifdef PROTO_CUDA
	      cudaEventRecord(stop);
	      cudaEventSynchronize(stop);
	      float milliseconds = 0;
              cudaEventElapsedTime(&milliseconds, start, stop);

	      double s = milliseconds * 1E-3;
	      double Gflops = (comp) * (nx*nx*nx) * (27*2+1) * 1E-9;
	      std::cout << " laplacian:  " << Gflops/s << " GFLOPS"<<std::endl;
#endif
            } 
        }
      protoDeviceSynchronizeGPU();
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
  
      
