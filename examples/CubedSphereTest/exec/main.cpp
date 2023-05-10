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
  PROTO_ASSERT(DIM==3,"cubed-sphere works only for DIM=3");
  int nx = 8;
    GetCmdLineArgumenti(argc, (const char**)argv, "nx", &nx);
    std::cout << "command line input:" << std::endl;
    std::cout << argv[0] << " -nx " << nx  << std::endl << endl;
    PR_TIMER_SETFILE(to_string(nx) + ".DIM" + to_string(DIM) + ".cubedSphere.fastSpread.time.table");
    PR_TIMERS("main");
    int nGhost = 3;  
    int maxRef = 3;
    int testCase = 0;
    Point waveNumber(1,1,1);
    Array<BoxData<double>,3> JRefs;
    Array<double,3> divNTNorm;
    HDF5Handler h5;
    for (int ref=0;ref < maxRef; ref++)
      {
        cout << endl << "nx : " << nx << endl;
        Array<BoxData<double>,DIM> dfdxi;
        BoxData<double,DIM> X;
        BoxData<double> divNonNorm;
        BoxData<double> divF;
        Array<BoxData<double,DIM>,DIM> NT;
        Box bx(Point::Zeros(),Point::Ones(nx-1));
        bx = bx.grow(nGhost);
        BoxData<double>& J = JRefs[ref];
        Box bxnode = bx.extrude(Point::Ones());
        BoxData<double> radius(bxnode);
        double h = 1.0/nx;
        forallInPlace_p([]PROTO_LAMBDA(Point a_pt,                            
                                       Var<double,1>& a_rad,
                                       double a_h)
                        {
                          double xi0 = 1.0 + a_pt[0]*a_h;
                          a_rad(0) = xi0*xi0/2;
                        },bxnode,radius,h);
        
        // Compute almost analytic X,NT,J.
        Operator::cubedSphereGeometry(X,NT,J,radius,bx,h);
        
        for (int dir = 0; dir < DIM; dir++)         
        {
          PR_TIMERS("Flux calculation");
          // cout << "dir = " << dir << endl;
          BoxData<double,1,MEMTYPE_DEFAULT,DIM> FAvDir;
          switch (testCase)
            {
            case 0:
              {
                //testing just the divergence calculation.
                FAvDir = fAv(X,waveNumber,dir);               
                break;
              }
            case 1:
              // Testing Scalar Laplacian.
              {
                BoxData<double> phiAv;
                {
                  PR_TIMERS("initialize phi");
                  auto phiE = phiExact(X,waveNumber);
                  phiAv = ((1.0)*Shift::Zeros())(phiE,bx.grow(nGhost));
                }
                // Average of J on face.
                BoxData<double> JFace = Stencil<double>::CellToFace(dir)(J);
                
                // Face-centered cofactor matrix N.
                auto NTMatrix = Operator::cofactorMatrix(NT,dir);
                
                // FAvDir is the gradient of phi with respect to x variables.
                {
                  PR_TIMERS("gradient for Laplacian");
                  FAvDir = 
                    Operator::_faceGradxPhi(phiAv,phiAv,
                                            NTMatrix,NTMatrix,
                                            JFace,JFace,
                                            dir);
                  FAvDir *= 1./h;
                }
                break;
              }
            default:
              {
                cout << "testCase = "<< testCase << " is not a valid test case"<< endl;
                break;
              }
            }
          BoxData<double> fluxdir =
              Operator::_faceMatrixProductAB(FAvDir,NT[dir],FAvDir,NT[dir],dir);
          dfdxi[dir] = Stencil<double>::FluxDivergence(dir)(fluxdir,1.0/h);       
        }
        // diagnostics :
        
        if (testCase < 2)
          {
            PR_TIMERS("Divergence");
            
            divNonNorm = dfdxi[0] + dfdxi[1] + dfdxi[2];
            //cout << "norm dfdx0: " << dfdxi[0].absMax() << endl;
            //cout << "norm dfdx1: " << dfdxi[1].absMax() << endl;
            //cout << "norm dfdx2: " << dfdxi[2].absMax() << endl;
            divF = Operator::_cellQuotient(divNonNorm,J,divNonNorm,J);            
            h5.writePatch(1./nx,divF,"divF"+to_string(nx));
            auto divfexact = divFExact(divF.box(),X,waveNumber);
            h5.writePatch(1./nx,divfexact,"divFExact"+to_string(nx));
            divF -= divfexact;
            h5.writePatch(1./nx,divF,"divFError"+to_string(nx));            
            auto erroldnorm = divF.absMax();
            cout << "max error = " << erroldnorm << endl;    
            cout << "divF Box = " << divF.box() << endl;
          }

        BoxData<double,DIM> divNT(bx);
        
        divNT.setToZero();
        for (int dir = 0; dir < DIM; dir++)
          {
            divNT += Stencil<double>::FluxDivergence(dir)(NT[dir],1.0/h);
            //cout << NT[dir].absMax() << endl;
          }
        divNTNorm[ref] = divNT.absMax() / JRefs[ref].absMax();
        cout << "norm of div(NT): " << divNTNorm[ref] << endl;
        nx*=2;
      }
        
    if (maxRef >= 3)
      {
        /*cout << endl;
        double normalize = JRefs[maxRef-1].absMax();
        BoxData<double> JAv(Box(Point::Zeros(),Point::Ones(nx/4-1)));
        JRefs[0].copyTo(JAv);
        JAv += Stencil<double>::AvgDown(2)(JRefs[1],-1.0);
        double jerrnorm0 = JAv.absMax()/normalize;
        cout <<"Jerr0 norm: "<<  jerrnorm0 <<endl;
        h5.writePatch(4./nx,JAv,"Jerr0"+to_string(nx));
        JAv.define(Box(Point::Zeros(),Point::Ones(nx/2-1)));       
        JAv += Stencil<double>::AvgDown(2)(JRefs[2],-1.0);
        JAv += JRefs[1];
        double jerrnorm1 = JAv.absMax()/normalize;
        cout <<"Jerr1 norm: "<<  jerrnorm1 <<endl;
        h5.writePatch(2./nx,JAv,"Jerr1"+to_string(nx));
        cout << "error exponent: "
             << log(jerrnorm0/jerrnorm1)/log(2.0) << endl;*/
      }
    
    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
}

