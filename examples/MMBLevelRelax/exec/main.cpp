#include "Proto.H"
#include "Proto_LevelBoxData.H"
#include "Proto_DisjointBoxLayout.H"
#include "TestMapping.H"

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
  // Setting up simple example of Poisson solver where the
  // input and output data is allocated by the user.
  
  int nx = 128;
  /* -------------------- */
  /* command-line parameters */
  /* -------------------- */
  GetCmdLineArgumenti(argc, (const char**)argv, "nx", &nx);
  std::cout << "command line input:" << std::endl;
  std::cout << argv[0] << " -nx " << nx  << std::endl;

  
  Box bx(Point::Zeros(),(nx-1)*Point::Ones());

  // Compute mapping evaluated at corners, rows of NT at faces, Jacobian.
  double length = 1.0;
  double h = length/nx;
  int nGhost = 4;
  
  PatchMap mapping(h);
  BoxData<double,DIM> X = mapping.map(bx,nGhost);
  std::array<BoxData<double,DIM>,DIM> NT;
  
  for (int dir = 0; dir < DIM;dir++)
    {
      NT[dir] = Operator::_cofactor(X,dir);
    }
  BoxData<double> J = Operator::_jacobian(X,NT);

  // compute divergence of a flux.
  
  BoxData<double> divNonNorm(bx);
  divNonNorm.setToZero();
  for (int dir = 0; dir < DIM; dir++)
    {
      BoxData<double,DIM> FAvDir = fAv(bx,X,dir);
      BoxData<double> fluxdir =
        Operator::_faceMatrixProductATB(NT[dir],FAvDir,NT[dir],FAvDir,dir);
      divNonNorm += Stencil<double>::FluxDivergence(dir)(fluxdir);
    }
  cout << "Non-norm Box = " << divNonNorm.box() << endl;
  auto divF = Operator::_cellTensorQuotient(J,divNonNorm,J,divNonNorm);
  auto divFError = divFExact(divF.box(),X);
  divFError -= divF;
  auto errnorm = divFError.absMax();
  cout << "max error = " << errnorm << endl;
  cout << "Box = " << divF.box() << endl;
}
