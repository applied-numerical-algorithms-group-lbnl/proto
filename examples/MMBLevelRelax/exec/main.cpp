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

    int nx = 64;
    GetCmdLineArgumenti(argc, (const char**)argv, "nx", &nx);
    std::cout << "command line input:" << std::endl;
    std::cout << argv[0] << " -nx " << nx  << std::endl << endl;
    PR_TIMER_SETFILE(to_string(nx) + ".DIM" + to_string(DIM) + ".MMBOperator.time.table");
    PR_TIMERS("main");
    int nGhost = 4;  
    int numLevels = 3;
    Array<Array<double,DIM > , DIM> arr;
    arr[0][0] = 1.0;
    arr[0][1] = 0.0; 
    arr[1][0] = 0.0;
    arr[1][1] = 1.0;
#if DIM==2
    Array<double,DIM> coef = {0.025,0.025};
    //Array<double,DIM> coef = {0.0,0.0};
    //Array<double,DIM> coef = {0.0,0.025};
    Point waveNumber(1,1);
#endif  

#if DIM==3
    arr[0][2] = 0.0;
    arr[1][2] = 0.0;
    arr[2][0] = 0.0;
    arr[2][1] = 0.0;
    arr[2][2] = 1.0;
    Array<double,DIM> coef = {0.025,0.025,.025};
    Point waveNumber(1,1,1);
    
    //Array<double,DIM> coef = {0.0,0.0,0.0};
#endif

    cout << "waveNumber = " << waveNumber << endl;
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
        BoxData<double,DIM> X = mapping.map(bx,nGhost+2);
        Array<BoxData<double,DIM>,DIM> NT;
        
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
        Array<BoxData<double,1>,DIM> dfdxi;
        
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
                auto NTMatrix = Operator::_cofactorMatrix(NT,dir);
                
                // FAVDir is the gradient of phi with respect to x variables.
                {
                  PR_TIMERS("gradient for Laplacian");
                FAvDir = 
                  Operator::_faceGradxPhi(phiAv,phiAv,NTMatrix,NTMatrix,JFace,JFace,dir);
                }
                break;
              }
            case 2:
              {
                // Test MHD flux calculation.
                PROTO_ASSERT(DIM==3,"MHD spherical flux works only for DIM=3");
                Box bxFace = bx.grow(nGhost).extrude(dir);
                BoxData<double,8> primvars4(bxFace);
                BoxData<double,8> primvars2(bxFace);
                BoxData<double,DIM,MEMTYPE_DEFAULT,DIM> A4(bxFace);
                BoxData<double,DIM,MEMTYPE_DEFAULT,DIM> A2(bxFace);
                BoxData<double,DIM> drAdjA4(bxFace);
                BoxData<double,DIM> drAdjA2(bxFace);
                BoxData<double,1> drDetA4(bxFace);
                BoxData<double,1> drDetA2(bxFace);
                double gamma = 5./3.;                 
                auto flux = Operator::MHDSphericalFlux<double,8,8,MEMTYPE_DEFAULT>
                  (primvars4,primvars2,A4,A2,drDetA4,drDetA2,drAdjA4, drAdjA2,gamma,dir);
                cout << "dir = " << dir << ", input Box = " << bxFace << ", return Box = " << flux.box() << endl;
                break;
              }
            default:
              {
                cout << "testCase = "<< testCase << " is not a valid test case"<< endl;
                abort;
              }
            }
          if (testCase != 2)
            {
              PR_TIMERS("Divergence");
              BoxData<double> fluxdir =
                Operator::_faceMatrixProductAB(FAvDir,NT[dir],FAvDir,NT[dir],dir);
              dfdxi[dir] = Stencil<double>::FluxDivergence(dir)(fluxdir);
            }
        }
        if (testCase != 2)
          { 
            if (DIM == 2)
              {
                divNonNorm = dfdxi[0] + dfdxi[1];
              }
            else
              {
                divNonNorm = dfdxi[0] + dfdxi[1] + dfdxi[2];
              }
            
            auto divFOld = Operator::_cellQuotient(divNonNorm,J,divNonNorm,J);
            h5.writePatch(1./nx,divFOld,"divF"+to_string(nx));
            auto divfexact = divFExact(divFOld.box(),X,waveNumber);
            h5.writePatch(1./nx,divfexact,"divFExact"+to_string(nx));
            divFOld -= divfexact;
            h5.writePatch(1./nx,divFOld,"divFError"+to_string(nx));            
            auto erroldnorm = divFOld.absMax();
            cout << "max error = " << erroldnorm << endl;    
            cout << "divF Box = " << divFOld.box() << endl;
          }
          nx*=2;     
        }
    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
}

