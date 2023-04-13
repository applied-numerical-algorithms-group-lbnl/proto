
#include "Proto.H"
#include "implem/Proto_LevelData.H"
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

  int nx = 128;
  int ny = 128;
  int nz = 128;
  int maxbox = 32;
  int niters = 10;

  /* -------------------- */
  /* command-line parameters */
  /* -------------------- */
  GetCmdLineArgumenti(argc, (const char**)argv, "nx", &nx);
  ny = nx;
  nz = nx;
  maxbox = 32;
  GetCmdLineArgumenti(argc, (const char**)argv, "ny", &ny);
  GetCmdLineArgumenti(argc, (const char**)argv, "nz", &nz);
  GetCmdLineArgumenti(argc, (const char**)argv, "maxbox", &maxbox);
  GetCmdLineArgumenti(argc, (const char**)argv, "niters", &niters);
  int nstream = 8;
#ifdef PROTO_ACCEL
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
  
  Reduction<double,Abs> Rxn;

  for(unsigned int i=0; i<dbl.size(); i++)
  {
    auto& u = U[i];
    auto& rhs = RHS[i];
    u.setVal(1.);
    rhs.setVal(0.);
  }
  //brace here just for timers.
  Reduction<double> rxn;  //Abs, init to zero
  double maxwave = 1.0;
  {
    PR_TIME("full_euler_iteration");
    for(unsigned int iter = 0; iter < niters; iter++)
    {
      U.exchange();
      if(iter != 0)
        {
          maxwave = rxn.fetch();
          rxn.reset();
        }
      else{
        rxn.reset();
      }
      for(unsigned int i=0; i<dbl.size(); i++)
      {
       auto& u = U[i];
       auto& rhs = RHS[i];
       Box rbox = dbl[i];

       EulerOp::step( u, rhs, rbox, rxn, false, false);
      }

    }
#ifdef PROTO_ACCEL    
      protoDeviceSynchronize();
      protoError err = protoGetLastError();
      if (err != protoSuccess)
      {
        fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                __FILE__, __LINE__, protoGetErrorString(err));
      }
#endif
  }    

  printf("out of loop --- writing report\n");
  PR_TIMER_REPORT();
  return 0;
}
  
      
