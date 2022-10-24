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
    GetCmdLineArgumenti(argc, (const char**)argv, "nx", &nx);
    std::cout << "command line input:" << std::endl;
    std::cout << argv[0] << " -nx " << nx  << std::endl << endl;
    PR_TIMER_SETFILE(to_string(nx) + ".DIM" + to_string(DIM) + ".MMBOperator.time.table");
    PR_TIMERS("main");
    int nGhost = 5;  
    int numLevels = 1;
    array<array<double,DIM > , DIM> arr;
    arr[0][0] = 1.0;
    arr[0][1] = 1.0; 
    arr[1][0] = 0.0;
    arr[1][1] = 1.0;
#if DIM==2
    //array<double,DIM> coef = {0.025,0.025};
    array<double,DIM> coef = {0.0,0.0};
    Point waveNumber(1,1);
#endif  

#if DIM==3
    arr[0][2] = 0.0;
    arr[1][2] = 0.0;
    arr[2][0] = 0.0;
    arr[2][1] = 0.0;
    arr[2][2] = 1.0;
    array<double,DIM> coef = {0.025,0.025,.025};
    Point waveNumber(1,1,1);
    //array<double,DIM> coef = {0.0,0.0,0.0};
#endif

    cout << "waveNumber = " << waveNumber << endl << endl;
    double length = 1.0;
    /* -------------------- */
    /* command-line parameters */
    /* -------------------- */
    int testCase = 1;
    HDF5Handler h5;
    for (int refiter = 0;refiter < numLevels;refiter++)
    {
        Box bx(Point::Zeros(),(nx-1)*Point::Ones());
        // Compute mapping evaluated at corners, rows of NT at faces, Jacobian.
        double h = length/nx;            
        PatchMap mapping(arr,coef,h);
        BoxData<double,DIM> X = mapping.map(bx,nGhost);
        std::array<BoxData<double,DIM>,DIM> NT;

        for (int dir = 0; dir < DIM;dir++)
        {
            PR_TIMERS("NT");
            NT[dir] = Operator::cofactor(X,dir);
        }
        BoxData<double> J;
        {
            PR_TIMERS("Jacobian");
            J = Operator::jacobian(X,NT);
        }
        BoxData<double> divNonNorm(bx);
        divNonNorm.setToZero();
        for (int dir = 0; dir < DIM; dir++)
        {
          BoxData<double,DIM> FAvDir;
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
                PR_TIMERS("gradient for Laplacian");
                BoxData<double> phiAv;
                {
                  PR_TIMERS("initialize phi");
                  phiAv = phiExact(X,waveNumber);
                }
                // Average of J on face.
                BoxData<double> JFace = Stencil<double>::CellToFace(dir)(J);
                
                // Face-centered cofactor matrix N.
                auto NTMatrix = Operator::_cofactorMatrix(NT,dir);
                cout << NTMatrix.box() << " NTMatrix Box " << endl;
                
                // Gradient of phi with respect to xi variables
                auto gradxiphi = Operator::_faceGradient(phiAv,dir);
                cout << gradxiphi.box() << " gradxiphi Box " << endl;
                 
                // FAVDir is the gradient of phi with respect to x variables.
                auto temp = Operator::_faceMatrixProductABT
                  (NTMatrix,gradxiphi,NTMatrix,gradxiphi,dir);
                cout << temp.box() << " temp Box " << endl;
                cout << JFace.box() << " JFace Box " << endl;
                FAvDir =
                   Operator::_faceTensorQuotient(temp,JFace,temp,JFace,dir);
                cout << FAvDir.box() << " FAvDir Box " << endl;
                break;
              }
            default:
              {
                cout << "testCase = "<< testCase << " is not a valid test case"<< endl;
                abort;
              }
            }
          {
            PR_TIMERS("Divergence");
            BoxData<double> fluxdir =
              Operator::_faceMatrixProductATB(NT[dir],FAvDir,NT[dir],FAvDir,dir);
            divNonNorm += Stencil<double>::FluxDivergence(dir)(fluxdir);
          }
        }
        auto divFOld = Operator::_cellQuotient(divNonNorm,J,divNonNorm,J);
        h5.writePatch(1./nx,divFOld,"divF"+to_string(nx));
        auto divfexact = divFExact(divFOld.box().grow(1),X,waveNumber);
        h5.writePatch(1./nx,divfexact,"divFExact"+to_string(nx));
        divFOld -= divfexact;
        h5.writePatch(1./nx,divFOld,"divFError"+to_string(nx));
        auto erroldnorm = divFOld.absMax();
        cout << "max error = " << erroldnorm << endl;    
        cout << "divF Box = " << divFOld.box() << endl;
        nx*=2;
    }
    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
}

