#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "ProtoAMR.H"
#include "AdvectionOp.H"
#include "AMRSubcycleExplicit.H"
#include "AdvectionTest.H"

using namespace std;
using namespace Proto;
namespace Proto
{
  void writeAMRGridPlot(double a_dx,AMRGrid a_grids,string a_label)
  {
    AMRData<short,1> data(a_grids,Point::Zeros());
    for (int ii = 0; ii < a_grids.numLevels();ii++)
      {
        auto & dataLev = data[ii];
        for (auto dit=dataLev.begin();*dit != dit.end();++dit)
          {
            dataLev[*dit].setVal(ii);
          }
      }
    HDF5Handler h5;
    h5.writeAMRData(a_dx,data,a_label);
  }
}
void
parseCommandLine(double& a_tmax, int& a_nx, int& a_maxstep, int& a_outputinterval,
                 int& a_refratio, int argc, char* argv[])
{
  a_tmax= 1.0;
  a_maxstep = 12;
  a_outputinterval = 1;
  a_nx = 32;
  a_refratio = 4;
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
    else if(strcmp(argv[iarg],"-r") == 0)
    {
      a_refratio = atof(argv[iarg+1]);
    }
  }
  cout << "nx = " << a_nx << endl;
  cout << "maxstep = " << a_maxstep << endl;
  cout << "output interval = " << a_outputinterval << endl;
  cout << "max time = " << a_tmax << endl;
  cout << "refinement ratio = " << a_refratio << endl;
}
int main(int argc, char* argv[])
{
#ifdef PR_MPI
    MPI_Init(&argc,&argv);
#endif
    

    int pid=procID();
   
    if (pid==0) {
      cout << "Two-level test for RK4 Explicit using EulerOp" << endl;
      cout << "usage:  " << argv[0] << " -n nx  -t tmax -m maxstep -o output_interval -r refratio" << endl;
       }
    double tstop;
    int size1D, maxStep, outputInterval, refratio;
    parseCommandLine(tstop, size1D, maxStep, outputInterval, refratio, argc, argv);
    Point spaceRefRatio = Point::Ones(refratio);
    int boxsize1D = min(size1D,32);
    int bitmapsize1D = size1D/boxsize1D;
    Box bigbox0 = Box(Point::Zeros(),Point::Ones(size1D-1));
    PR_assert(boxsize1D*bitmapsize1D == size1D);
    
    Point boxsize = boxsize1D*Point::Ones();
    Point coarseSize = size1D*Point::Ones();
    //int sizeDomain=size1D;
    PR_TIMER_SETFILE(to_string(size1D) + ".DIM="+ to_string(DIM) + ".numProc=" +to_string(numProc())+ ".AMRSubcycleAdvection.time.table");
    PR_TIMERS("main");
    
    Box domain(Point::Zeros(),Point::Ones()*(size1D -1));
    array<bool,DIM> per;
    for(int idir = 0; idir < DIM; idir++) per[idir]=true;
    ProblemDomain pd(domain,per);

    double dx = 1.0/size1D;
    double dt = .5/size1D;
    cout << "dx = " << dx << endl;
    
    if (pid==0)
      {
        std::cout << "Coarsest domain: " << domain << std::endl;
        std::cout << "Coarsest dt: " << dt << std::endl;
      }
    DisjointBoxLayout dblCoarse(pd,boxsize);
    ProblemDomain pdfine = pd.refine(spaceRefRatio);

    int numLevels = 3;
    int buffersize = 2;
   
    // Initialize grids. All refinement ratios equal, isotropic.
    // Start by initializing the coarsest dbl, then fill the other levels
    // using regrid.
    vector<Point> refRatios(numLevels-1,Point::Ones()*refratio);
    AMRGrid initGrids(dblCoarse,refRatios,numLevels);
    Advection op;    
    double dxLevel = dx;

     for (int lev = 0; lev < numLevels-1;lev++)
      {
        // Compute the exact solution at this level.
        LevelBoxData<double,NUMCOMPS,MEMTYPE_DEFAULT>
                             initData(initGrids[lev],Point::Zeros());
        advectionExact<double>(initData,dxLevel,0);
        // Generate tags.
        LevelTagData levtags(initGrids[lev],Point::Ones(buffersize));
        levtags.setToZero();
        op.generateTags(levtags,initData,dxLevel,0);
        AMRGrid::buffer(levtags,buffersize);
        cout << "end buffering in initialization in proc " << pid << endl;
        initGrids.regrid(levtags,lev,boxsize);
        dxLevel = dxLevel/refratio;
      }
    
    // Make another pass through the hierarchy to insure proper nesting.
    for (int lev = 2;lev < numLevels;lev++)
      {
        cout << "enforce nesting" << endl;
        initGrids.enforceNesting2(lev);
      }

    // Define pointer to solution data.
    
    Point ghostsize = Advection::ghostSize();
    auto amrdataPtr = shared_ptr<AMRData<double,NUMCOMPS,MEMTYPE_DEFAULT> >
      (new AMRData<double,NUMCOMPS,MEMTYPE_DEFAULT>(initGrids,ghostsize));
    auto globalStep = shared_ptr<vector<int> >(new vector<int >);
    // Initialize solution data.
    dxLevel = dx;
    for (int lev = 0; lev < numLevels; lev++)
      {
        advectionExact((*amrdataPtr)[lev],dxLevel,0.0);
        dxLevel /= refratio;
        globalStep -> push_back(0);
      }
    
    AMRSubcycleExplicit<Advection,double,NUMCOMPS,MEMTYPE_DEFAULT>
      amrAdvection(amrdataPtr,globalStep,dx,refratio,0);
    // average down to obtain consistent conservation sums at all levels.
    amrdataPtr->averageDown();
    
    // Compute initial conservation sum.
    double sum0 = (*amrdataPtr)[0].sum();
    
#if 1
    if(pid==0)
      {
        cout <<"Initial level 0 conservation sum = " << sum0  << endl;
        cout << "starting time loop, maxStep = "<< maxStep << endl;
      }
#endif
    double time = 0;
    {
      HDF5Handler h5;
      h5.writeAMRData(dx, *amrdataPtr,"U_N0");
    }
    for (int k = 0;((k < maxStep) && (time < tstop));k++)
      {
        {
          PR_TIMERS("main advance");
          LevelFluxRegister<double,NUMCOMPS,MEMTYPE_DEFAULT> lfrdummy;
          amrAdvection.advance(lfrdummy,dt,0,false);
          time += dt;
        }     
        if ((k+1)%outputInterval == 0)
          {
#ifdef PR_HDF5
            if(pid==0) cout << "writing data before time step = " << k+1 << " , time = " << time << endl;
            HDF5Handler h5;
            h5.writeAMRData(dx, *amrdataPtr,"U_N%i", k+1);           
#endif
          }
      }
    /*
    double dxFine = dx/spaceRefRatio[0];
    advectionError(error[0],(*amrdataPtr)[0],dx,1.,time);
    advectionError(error[1],(*amrdataPtr)[1],dxFine,1.,time);
#ifdef PR_HDF5
    HDF5Handler h5;
    h5.writeAMRData(dx, error,"errorFinal");
#endif    
    double errmax = error[0].absMax(); 
    */
    sum0 = (*amrdataPtr)[0].sum();
    if (pid==0)
      {
        //cout << "max error = " << errmax << endl;
        cout << "Final level 0 conservation sum = " << sum0 << endl;
      }
    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
}
