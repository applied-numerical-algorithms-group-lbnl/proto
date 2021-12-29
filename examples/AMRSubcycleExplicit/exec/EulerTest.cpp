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

#include "Proto.H"
#include "EulerOpWrapper.H"
#include "AMRSubcycleExplicit.H"
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"
#include "AdvectionTest.H"

using namespace std;
using namespace Proto;
  
 void
parseCommandLine(double& a_tmax, int& a_nx, int& a_maxstep, int& a_outputinterval, int argc, char* argv[])
{
  a_tmax= 1.0;
  a_maxstep = 1;
  a_outputinterval = 100000;
  a_nx = 64;
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
int main(int argc, char* argv[])
{
#ifdef PR_MPI
    MPI_Init(&argc,&argv);
#endif
    
    int pid=procID();

    if(pid==0) {
        cout << "Two-level test for RK4 Explicit using EulerOp" << endl;
        cout << "usage:  " << argv[0] << " -n nx  -t tmax -m maxstep -o output_interval" << endl;
    }
    double tstop;
    int size1D, maxStep, outputInterval;
    parseCommandLine(tstop, size1D, maxStep, outputInterval, argc, argv);
    double gamma = 1.4;

    int boxsize1D = size1D;
    int bitmapsize1D = size1D/boxsize1D;
    PR_assert(boxsize1D*bitmapsize1D == size1D);
    
    Point boxsize = boxsize1D*Point::Ones();
    Point coarseSize = size1D*Point::Ones();
    //int sizeDomain=size1D;
    PR_TIMER_SETFILE(to_string(size1D) + ".DIM = "+ to_string(DIM) + ".AMRSubcycleEuler.time.table");
    PR_TIMERS("main");
    
    Box domain(Point::Zeros(),Point::Ones()*(size1D -1));
    array<bool,DIM> per;
    for(int idir = 0; idir < DIM; idir++) per[idir]=true;
    ProblemDomain pd(domain,per);

    double dx = 1.0/size1D;
    double dt = .25/size1D;
    if(pid==0)
      {
        std::cout << "Coarsest domain size: " << size1D << std::endl;
        std::cout << "Coarsest dt: " << dt << std::endl;
      }
    DisjointBoxLayout dblCoarse(pd,boxsize);

    ProblemDomain pdfine = pd.refine(PR_AMR_REFRATIO*Point::Ones());
    DisjointBoxLayout dblFine(pdfine,boxsize);
#if 0
    int bitmapsize1DFine = size1D*PR_AMR_REFRATIO/boxsize1D;
    vector<Point> finePatches = {(bitmapsize1DFine/2-1)*Point::Ones()};
    DisjointBoxLayout dblFine(pdfine,finePatches,boxsize);
#endif
    int numLevels = 1;
    
    // amrDataPtr points to data holders for AMR calculation.
    
    vector<DisjointBoxLayout> dbls = {dblCoarse};
    vector<double> dxlevel ={dx,dx/PR_AMR_REFRATIO};
    AMRGrid amrgrid(dbls,numLevels);
    Point ghostsize = Advection::ghostSize();
    auto amrdataPtr = shared_ptr<AMRData<double,DIM+2,MEMTYPE_DEFAULT> >
      (new AMRData<double,DIM+2,MEMTYPE_DEFAULT>(amrgrid,ghostsize));

    // Error for debugging.
    AMRData<double,DIM+2,MEMTYPE_DEFAULT> error(amrgrid,Point::Zeros());
    
    AMRSubcycleExplicit<Advection,double,DIM+2,MEMTYPE_DEFAULT> amreuler;
    amreuler.define(amrdataPtr,dx,PR_AMR_REFRATIO,0);

    double time = 0.;
    double dxLevel = dx;
    for (int level = 0;level < numLevels;level++)
      {
        advectionExact<double>((*amrdataPtr)[level],dxLevel,time);
        dxLevel /= PR_AMR_REFRATIO;
      }
    Stencil<double> interp= Stencil<double>::CellToFace(0,Side::Lo,5);
    interp.invert(0);
    cout << "Stencil Coefs, offsets" << endl;
    double coefTot = 0.;
    for (int isten; isten < interp.size();isten++)
      {
        cout << interp.coefs()[isten] <<  " , " << interp.offsets()[isten] << endl;
        coefTot += interp.coefs()[isten];
      }
    cout << "sum of coefs = " << coefTot << endl;
    if(pid==0)
      cout << "starting time loop, maxStep = "<< maxStep << endl;
    for (int k = 0;(k < maxStep) && (time < tstop);k++)
      {
        // writeDivergenceError<double>((*amrdataPtr)[0],dx,time,-1);
        {
          PR_TIMERS("main advance");
          LevelFluxRegister<double,DIM+2,MEMTYPE_DEFAULT> lfrdummy;
          amreuler.advance(lfrdummy,dt,0,false);
          time += dt;
        }
        //writeDivergenceError<double>((*amrdataPtr)[0],dx,time,4);
        
        if (k%outputInterval == 0)
          {
#ifdef PR_HDF5
            if(pid==0) cout << "writing data for time step = " << k << endl;
            HDF5Handler h5;
            std::array<double, DIM> dx_vect;
            advectionError(error[0],(*amrdataPtr)[0],dx,time,time);
            h5.writeAMRData(dx, *amrdataPtr,"U_N%i", k);
            h5.writeAMRData(dx, error,"error_N%i", k);
#endif
          }
      }
    
    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
}
