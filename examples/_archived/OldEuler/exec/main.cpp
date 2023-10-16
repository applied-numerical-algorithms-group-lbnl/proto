#include "Proto.H"
#include "RK4.H"
#include "InputParser.H"
#include "EulerLevelDataRK4.H"
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"

using namespace std;
using namespace Proto;
/*
void
parseCommandLine(double& a_tmax, int& a_nx, int& a_maxstep, int& a_outputinterval, int argc, char* argv[])
{
  a_tmax= 1.0;
  a_maxstep = 1;
  a_outputinterval = -1;
  a_nx = 128;
  for(int iarg = 0; iarg < argc-1; iarg++)
  {
    if(strcmp(argv[iarg],"-n") == 0)
    {
      a_nx = atoi(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg], "-m") == 0)
    {
      a_maxstep = atoi(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg], "-o") == 0)
    {
      a_outputinterval = atoi(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg],"-t") == 0)
    {
      a_tmax = atof(argv[iarg+1]);
    }
  }
}
*/



void InitializeEulerLevelDataState(EulerLevelDataState& state)
{
    (state.m_U).setToZero();
    for(auto dit : state.m_U) {
        EulerOp::initializeState((state.m_U)[dit],state.m_dx,state.m_gamma);
    }
}

/**
 * @brief Write out component comp of data to filename.vtk. This routine is a no-op for any process
 * that is not process 0.
 * @param data LevelBoxData that is defined on a single box assigned to process 0
 * @param comp integer in the range [0,NUMCOMPS-1]
 */
void WriteSinglePatchLevelData(LevelBoxData<double,NUMCOMPS>& data,
                               const int comp,
                               const double& dx,
                               const string& filename)
{
    if(procID()==0)
    {
        auto dit=data.begin();
        BoxData<double> rho=slice(data[*dit],comp);
        WriteBoxData(filename.c_str(),rho,dx);
    }
}
void WriteSinglePatchLevelData(LevelBoxData<double,NUMCOMPS>& data,
                               const int comp,
                               const double& dx,
                               const string& fileroot,
                               const string& filename)
{
    if(procID()==0)
    {
        auto dit=data.begin();
        BoxData<double> rho=slice(data[*dit],comp);
        std::array<double,DIM> corner;
        corner.fill(0.);
        const char* filerootArg[1];
        filerootArg[0] = fileroot.c_str();
        WriteBoxData(filename.c_str(),rho,filerootArg,&corner[0],dx);
    }
}
int main(int argc, char* argv[])
{
#ifdef PR_MPI
    MPI_Init(&argc,&argv);
#endif
    /*    
    int pid=procID();

    if(pid==0) {
        cout << "Navier Stokes simulation of shear flow with sinusoidal perturbation.  Periodic bcs." << endl;
        cout << "usage:  " << argv[0] << " -n nx  -t tmax -m maxstep -o output_interval" << endl;
    }
    */

    // DEFAULT PARAMETERS
    int domainSize = 64;
    int boxSize = 32;
    double maxTime = 1.0;
    int maxStep = 10;
    int outputInterval = 1;
    double gamma = 1.4;

    // PARSE COMMAND LINE
    InputArgs args;
    args.add("domainSize",     domainSize);
    args.add("boxSize",        boxSize);
    args.add("maxTime",        maxTime);
    args.add("maxStep",        maxStep);
    args.add("outputInterval", outputInterval);
    args.parse(argc, argv);
    args.print();
    
    PR_TIMER_SETFILE(to_string(domainSize) + "OldEuler.time.table");
    PR_TIMERS("main");
    LevelBoxData<double,NUMCOMPS> U[3];
    for (int lev=0; lev<1 ; lev++)
    {
        Box domainBox = Box::Cube(domainSize);
        array<bool,DIM> per;
        per.fill(true);
        ProblemDomain pd(domainBox,per);

        double dx = 1.0/domainSize;
        double dt = .25/domainSize;
        if(procID() == 0)
        {
            std::cout << "domainSize: " << domainSize << std::endl;
            std::cout << "dt: " << dt << std::endl;
        }

        RK4<EulerLevelDataState,EulerLevelDataRK4Op,EulerLevelDataDX> rk4;
        EulerLevelDataState state(pd, Point::Ones(domainSize), dx, gamma);
        InitializeEulerLevelDataState(state);

        double time = 0.0;
        if(procID() == 0)
        {
            cout << "starting time loop, maxStep = "<< maxStep << endl;
        }
#ifdef PR_HDF5
        HDF5Handler h5;
        std::vector<std::string> varnames(NUMCOMPS);
        varnames[0] = "rho";
        for (int ii = 1; ii <= DIM; ii++) { varnames[ii] = ("rho_v" + std::to_string(ii-1)); }
        varnames[NUMCOMPS-1] = "rho_E";
        h5.writeLevel(varnames, dx, state.m_U, "U_D%i_I0", domainSize);
#endif
        for (int k = 0; (k < maxStep) && (time < maxTime); k++)
        {
            rk4.advance(time, dt, state);
            time += dt;
            if ((k+1) % outputInterval == 0)
            {
#ifdef PR_HDF5
                h5.writeLevel(varnames, dx, state.m_U, "U_D%i_I%i", domainSize, k+1);
#else
                LevelBoxData<double,NUMCOMPS>
                    UOut(DisjointBoxLayout(pd,domainSize*Point::Ones()),Point::Zeros());
                (state.m_U).copyTo(UOut);
                std::string fileroot="rho_"+std::to_string(domainSize);
                std::string filename=fileroot+"_"+to_string(k+1);
                WriteSinglePatchLevelData(UOut, 0, dx, fileroot,filename);
#endif
            }
        }

        //Solution on a single patch
        //U[lev].define(DisjointBoxLayout(pd,domainSize*Point::Ones()),Point::Zeros());
        //(state.m_U).copyTo(U[lev]);
        //std::string filename="rho_"+std::to_string(lev);
        //WriteSinglePatchLevelData(U[lev], 0, dx, filename);

        domainSize *= 2;
        //sizeDomain *= 2; //For debugging: if you want to keep the number of boxes the same
        maxStep *= 2;
        dt *= .5;
    }

    //Until we have a coarsening operation for LevelBoxData, we perform the error calculations on a single patch.
#if 0
    if(pid==0)
    {
        Reduction<double> rxn;
        double rhoErrMax[2];
        for(int ilev=0; ilev<2; ilev++)
        {
            rxn.reset();
            DataIterator dit_lev=U[ilev].begin();
            DataIterator dit_levp1=U[ilev+1].begin();
            BoxData<double,1> err=slice(U[ilev][*dit_lev],0);
            err-=Stencil<double>::AvgDown(2)(slice(U[ilev+1][*dit_levp1],0));
            err.absMax(rxn);
            rhoErrMax[ilev]=rxn.fetch();
            //rhoErrMax[ilev]=err.absMax();
            std::string filename="rho_err"+std::to_string(ilev);
//NOTE: this assumes that the domain length is 1.0, which is assumed throughout this code. May cause errors if this changes.
            double dx=1./(err.box().size(0));
            WriteBoxData(filename.c_str(),err,dx);
            std::cout << "Lev: " << ilev << " , " << rhoErrMax[ilev] << std::endl;
        }
        double rate = log(abs(rhoErrMax[0]/rhoErrMax[1]))/log(2.0);
        std::cout << "order of accuracy = " << rate << std::endl;
    }
#endif
    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
}
