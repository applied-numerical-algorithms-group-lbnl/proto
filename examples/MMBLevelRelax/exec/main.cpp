#include "Proto.H"
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
  HDF5Handler h5;
  
  Box bx(Point::Zeros(),(nx-1)*Point::Ones());
  // Compute mapping evaluated at corners, rows of NT at faces, Jacobian.
  double length = 1.0;
  double h = length/nx;
  int nGhost = 4;
  Point waveNumber(1,0);
  cout << waveNumber << endl;
  PatchMap mapping(h);
  BoxData<double,DIM> X = mapping.map(bx,nGhost);
  string strX = "X";
 h5.writePatch(1.0/nx,X,strX);
  std::array<BoxData<double,DIM>,DIM> NT;
  
  for (int dir = 0; dir < DIM;dir++)
    {
      NT[dir] = Operator::_cofactor(X,dir);
      string strloc = "NT" + to_string(dir);
      h5.writePatch(1.0/nx,NT[dir],strloc);
    }
  BoxData<double> J = Operator::_jacobian(X,NT);
  string strJ = "J";
  h5.writePatch(1.0/nx,J,strJ);
  // compute divergence of a flux.
  
  BoxData<double> divNonNorm(bx);
  divNonNorm.setToZero();
  for (int dir = 0; dir < DIM; dir++)
    {
      BoxData<double,DIM> FAvDir = fAv(X,waveNumber,dir);
      string strAvdir = "fAvDir" + to_string(dir);
       h5.writePatch(1.0/nx,FAvDir,strAvdir);
      BoxData<double> fluxdir =
        Operator::_faceMatrixProductATB(NT[dir],FAvDir,NT[dir],FAvDir,dir);
      string strdir = "fluxdir" + to_string(dir);
       h5.writePatch(1.0/nx,fluxdir,strdir);
      divNonNorm += Stencil<double>::FluxDivergence(dir)(fluxdir);
      string str="divFNonNormDir"+to_string(dir);
      h5.writePatch(1.0/nx,divNonNorm,str);
    }
  string str="divFNonNorm";
  h5.writePatch(1.0/nx,divNonNorm,str);
  cout << "Non-norm Box = " << divNonNorm.box() << endl;
  auto divF = Operator::_cellTensorQuotient(divNonNorm,J,divNonNorm,J);
  //auto divF = Operator::_cellQuotient(divNonNorm,J,divNonNorm,J);
  str="divF";
  h5.writePatch(1.0/nx,divF,str);
  auto divFError = divFExact(divF.box(),X,waveNumber);
  str="divFExact";
  h5.writePatch(1.0/nx,divFError,str);
  divFError -= divF;
  str="divFError";
  h5.writePatch(1.0/nx,divFError,str);
  auto errnorm = divFError.absMax();
  cout << "max error = " << errnorm << endl;
  cout << "Box = " << divF.box() << endl;
}
