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
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    // Setting up simple example of Poisson solver where the
    // input and output data is allocated by the user.

    int nx = 128;
    GetCmdLineArgumenti(argc, (const char**)argv, "nx", &nx);
    std::cout << "command line input:" << std::endl;
    std::cout << argv[0] << " -nx " << nx  << std::endl << endl;
    PR_TIMER_SETFILE(to_string(nx) + ".DIM" + to_string(DIM) + ".MMBOperator.time.table");
    PR_TIMERS("main");
    int nGhost = 4;  
    int numLevels = 1;
    array<array<double,DIM > , DIM> arr;
    arr[0][0] = 1.0;
    arr[0][1] = 1.0; 
    arr[1][0] = 0.0;
    arr[1][1] = 1.0;
#if DIM==2
    array<double,DIM> coef = {0.025,0.025};
    //array<double,DIM> coef = {0.0,0.0};
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



    HDF5Handler h5;
    for (int refiter = 0;refiter < numLevels;refiter++)
    {
        Box bx(Point::Zeros(),(nx-1)*Point::Ones());
        // Compute mapping evaluated at corners, rows of NT at faces, Jacobian.

        double h = length/nx;            
        PatchMap mapping(arr,coef,h);
        BoxData<double,DIM> X = mapping.map(bx,nGhost);
        //h5.writePatch(h,X,
        //            "X"+to_string(nx));
        std::array<BoxData<double,DIM>,DIM> NT;

        for (int dir = 0; dir < DIM;dir++)
        {
            PR_TIMERS("NT");
            NT[dir] = Operator::cofactor(X,dir);
            //h5.writePatch(h,NT[dir],
            //           "NT" + to_string(nx) + "_" + to_string(dir));
            cout << "NT Box for " << dir << " direction = " << NT[dir].box() << endl;
        }
        BoxData<double> J;
        {
            PR_TIMERS("Jacobian");
            J = Operator::jacobian(X,NT);
        }
        // h5.writePatch(1.0/nx,J,"J" + to_string(nx));
        // compute divergence of a flux.

        BoxData<double> divNonNorm(bx);
        divNonNorm.setToZero();
        for (int dir = 0; dir < DIM; dir++)
        {
            PR_TIMERS("Divergence");
            BoxData<double,DIM> FAvDir = fAv(X,waveNumber,dir);
            //h5.writePatch(h,FAvDir,
            //              "fAvDir" + to_string(dir)+"_"+ to_string(nx));
            BoxData<double> fluxdir =
                Operator::_faceMatrixProductATB(NT[dir],FAvDir,NT[dir],FAvDir,dir);
            //h5.writePatch(h,fluxdir,
            //              "fluxdir" + to_string(dir)+"_" +to_string(nx));

            divNonNorm += Stencil<double>::FluxDivergence(dir)(fluxdir);
            //h5.writePatch(h,divNonNorm,
            //              "divFNonNormDir"+to_string(dir)+"_"+to_string(nx));

        }
        auto divFOld = Operator::_cellQuotient(divNonNorm,J,divNonNorm,J);
        //h5.writePatch(h,divFOld,
        //               "divFOld"+to_string(nx));
        /*
           h5.writePatch(h,divNonNorm,
           "divFNonNorm"+to_string(nx));
           cout << "Non-norm Box = " << divNonNorm.box() << endl; */

        //J.setVal(h*h);

        /*auto divF = forall<double>
          ([] PROTO_LAMBDA(Var<double>& a_ret,
          Var<double>& a_ql,
          Var<double>& a_qr)
          {a_ret(0) = a_ql(0)/a_qr(0);},divNonNorm,J);   */                     
        //BoxData<double> divF(divNonNorm.box());
        //divNonNorm.copyTo(divF);
        //divF *= (1./h/h);
        /*    
              BoxData<double> dJdx =  Stencil<double>::Derivative(1,0,2)(J);                   BoxData<double> dJdy =  Stencil<double>::Derivative(1,1,2)(J);
              BoxData<double> dDFdx =  Stencil<double>::Derivative(1,0,2)(divNonNorm);
              BoxData<double> dDFdy =  Stencil<double>::Derivative(1,1,2)(divNonNorm);
              auto divF = forall<double>
              ([] PROTO_LAMBDA(Var<double>& a_ret,
              Var<double>& a_J,
              Var<double>& a_divF,
              Var<double>& a_dJdx,
              Var<double>& a_dJdy,
              Var<double>& a_dDFdx,
              Var<double>& a_dDFdy)
              {
              a_ret(0) = a_divF(0)/a_J(0)*
              (1.0 + (a_dJdx(0)*a_dJdx(0) + a_dJdy(0)*a_dJdy(0))
         *(1./(12.0*a_J(0)*a_J(0)))) + 
         (a_dDFdx(0)*a_dJdx(0) + a_dDFdy(0)*a_dJdy(0))*
         (-1./(12.0*a_J(0)*a_J(0)));
         },J,divNonNorm,dJdx,dJdy,dDFdx,dDFdy);

         h5.writePatch(h,divF,
         "divF"+to_string(nx)); */

        //BoxData<double> id(divNonNorm.box());
        //id.setVal(1.0);
        //auto Jinv = Operator::_cellQuotient(id,J,id,J);
        auto divfexact = divFExact(divFOld.box().grow(1),X,waveNumber);
        //divF -= divfexact;
        divFOld -= divfexact;
        /* BoxData<double> tmpJ(divfexact.box());
           J.copyTo(tmpJ);
           auto divFError = Operator::_cellProduct(divfexact,tmpJ,divfexact,tmpJ);
           divFError -= divNonNorm;
        //divFError *= 1./h/h;
        auto divFError2 = Operator::_cellQuotient(divFError,J,divFError,J);
        h5.writePatch(h,divFError2,
        "divFError2_"+to_string(nx)); */
        //h5.writePatch(h,divFOld,
        //              "divFError"+to_string(nx));
        //divFError -= divF;
        //h5.writePatch(h,divFError,
        //            "divFError"+to_string(nx));

        //auto errnorm = divF.absMax();
        //cout << "max error = " << errnorm << endl;
        auto erroldnorm = divFOld.absMax();
        cout << "max error = " << erroldnorm << endl;
        //auto errnorm2 = divFError2.absMax();     
        cout << "divF Box = " << divFOld.box() << endl;
        nx*=2;
    }
    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
}

