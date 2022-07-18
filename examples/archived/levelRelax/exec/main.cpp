
#include "Proto.H"
#include "Proto_LevelBoxData.H"
#include "Proto_DisjointBoxLayout.H"
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

PROTO_KERNEL_START void fusedRelaxT(Var<double> u, Var<double> temp, Var<double> rhs, double lambda)
{
  u(0) += lambda*(rhs(0)-temp(0));
}
PROTO_KERNEL_END(fusedRelaxT, fusedRelax);

#define PFLOPS 3
  
PROTO_KERNEL_START void initRHST(const Point& p, Var<double> rhs)
{
  rhs(0) = 0.0;
  if(p==Point(8,8,8)) rhs(0) = 1.0;
  if(p==Point(12,12,8)) rhs(0) = -1.0;
}

PROTO_KERNEL_END(initRHST, initRHS);

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
#ifdef PROTO_CUDA
  //GetCmdLineArgumenti(argc, (const char**)argv, "nstream", &nstream);
 // DisjointBoxLayout::setNumStreams(nstream);
#endif
  
  printf("nx = %d, ny = %d, nz= %d\n", nx, ny, nz);
  printf("maxbox = %d, niters = %d, nstream = %d\n", maxbox, niters, nstream);
  Box domain(Point::Zeros(), Point::Ones(nx-1));
  double dx = 1.0/(nx-1);
  std::array<bool, DIM> periodic;
  Point boxSize(maxbox, maxbox, maxbox);
  for(int idir = 0; idir < DIM; idir++) periodic[idir]=false;
  DisjointBoxLayout   dbl(ProblemDomain(domain, periodic),boxSize );

  LevelBoxData<double, 1> U(dbl, NGHOST*Point::Unit());
  LevelBoxData<double, 1> RHS(dbl, Point::Zero());      

  const double alpha = 0.1;
  const double beta = 2.5;


  for(auto dit=U.begin();*dit != dit.end();++dit)
  {
    BoxData<double>& u = U[*dit];
    BoxData<double>& rhs = RHS[*dit];
    u.setVal(0.);
    forallInPlace_p(initRHS, rhs);
  }
 
  const double diag = alpha + (beta*(-2.*DIM)/(dx*dx));
  const double lambda = 1/diag;
  Stencil<double> alphaSten = (-alpha)*Shift(Point::Zeros());
  Stencil<double> betaSten =  (-beta/(dx*dx))*(Stencil<double>::Laplacian());
  Stencil<double> negoperator = (-alpha)*Shift(Point::Zeros()) + (-beta/(dx*dx))*(Stencil<double>::Laplacian());
  {
    PR_TIME("Expanded version");
    for(unsigned int iter = 0; iter < niters; iter++)
    {
     // U.exchange();
 
      for(auto dit=U.begin();*dit != dit.end();++dit)
        {
          auto& u = U[*dit];
          auto& rhs = RHS[*dit];
          BoxData<double> temp = alphaSten(u);
          temp += betaSten(u);
          temp += rhs;
          temp*= lambda;
          u+= temp;
        }
    }
#ifdef PROTO_CUDA    
      {  PR_TIME("deviceSynch");
      protoDeviceSynchronize(MEMTYPE_DEFAULT);
      protoError err = protoGetLastError();
      if (err != protoSuccess)
      {
        fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                __FILE__, __LINE__, protoGetErrorString(err));
      }}
#endif
   }
   {
    PR_TIME("Folded version");
    for(unsigned int iter = 0; iter < niters; iter++)
      {
      //  U.exchange();
        for(auto dit=U.begin();*dit != dit.end();++dit)
          {
            auto& u = U[*dit];
            auto& rhs = RHS[*dit];
            BoxData<double> temp = negoperator(u);
          // forallInPlaceOp(PFLOPS, "algebra parts", fusedRelax, u, temp, rhs, lambda);
           forallInPlaceOp(PFLOPS, "algebra parts", [] PROTO_LAMBDA (Var<double> u, Var<double> temp, Var<double> rhs, double lambda)
                            {    u(0) += lambda*(rhs(0)-temp(0));}, u, temp, rhs, lambda);	
            
          }
      }
#ifdef PROTO_CUDA 
      { PR_TIME("deviceSynch");   
      protoDeviceSynchronize(MEMTYPE_DEFAULT);
      protoError err = protoGetLastError();
      if (err != protoSuccess)
      {
        fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                __FILE__, __LINE__, protoGetErrorString(err));
      }}
#endif
  }    

  printf("out of loop --- writing report\n");
  PR_TIMER_REPORT();
  return 0;
}
  
      
