
#include "Proto.H"
#include "Proto_LevelData.H"
#include "EulerOp.H"

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
  
int main(int argc, char* argv[])
{
  
  //have to do this to get a time table
  PR_TIMER_SETFILE("proto.time.table");

  int nx = 32;
  int ny = 32;
  int nz = 32;
  int maxbox = 32;
  int niters = 1;

  /* -------------------- */
  /* command-line parameters */
  /* -------------------- */
  GetCmdLineArgumenti(argc, (const char**)argv, "nx", &nx);
  ny = nx;
  nz = nx;
  GetCmdLineArgumenti(argc, (const char**)argv, "ny", &ny);
  GetCmdLineArgumenti(argc, (const char**)argv, "nz", &nz);
  GetCmdLineArgumenti(argc, (const char**)argv, "maxbox", &maxbox);
  GetCmdLineArgumenti(argc, (const char**)argv, "niters", &niters);
#ifdef PROTO_CUDA
  int nstream = 1;
  GetCmdLineArgumenti(argc, (const char**)argv, "nstream", &nstream);
  DisjointBoxLayout::setNumStreams(nstream);
#endif
  
  printf("nx = %d, ny = %d, nz= %d\n", nx, ny, nx);
  printf("maxbox = %d, niters = %d, nstream = %d\n", maxbox, niters, nstream);
  Box domain(Point::Zeros(), Point::Ones(nx-1));
  std::array<bool, DIM> periodic;
  for(int idir = 0; idir < DIM; idir++) periodic[idir]=true;
  DisjointBoxLayout   dbl(domain, maxbox, periodic);

  LevelData<BoxData<double, NUMCOMPS>> U(dbl, NGHOST*Point::Unit());
  LevelData<BoxData<double, NUMCOMPS>> RHS(dbl, Point::Zero());

  for(unsigned int iter = 0; iter < niters; iter++)
  {
    PR_TIME("full_euler_iteration");
    for(unsigned int i=0; i<dbl.size(); i++)
    {
      auto u = U[i];
      auto rhs = RHS[i];
      Box rbox = dbl[i];
      printf("before step i = %d\n", iter);
      double wave = EulerOp::step(rhs, u, rbox);
      printf("after step i = %d\n", iter);
    }
#ifdef PROTO_CUDA    
    printf("before cudadevice sync\n");
    cudaDeviceSynchronize();
    cudaError err = cudaGetLastError();
    if (err != cudaSuccess)
    {
      fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
              __FILE__, __LINE__, cudaGetErrorString(err));
    }
    printf("after cudadevice sync\n");
#endif    
  }

  PR_TIMER_REPORT();
  return 0;
}
  
      
