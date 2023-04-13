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
#include "EulerRK4.H"
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"

#define NUMCOMPS DIM+2

using namespace std;
using namespace Proto;

void WriteData(int                              a_iter,
               BoxData<double,NUMCOMPS>&    a_state,
               double                           a_dx)
{
    char basename[1024];
    sprintf(basename,"euler.%06d",a_iter);

    const char* varnames[DIM+2];
    varnames[0] = "density";

    varnames[1] = "momentum-x";
    #if DIM > 1
    varnames[2] = "momentum-y";
    #endif
    #if DIM > 2
    varnames[3] = "momentum-z";
    #endif

    varnames[DIM+1] = "energy";
    double origin[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        origin[ii] = 0.0;
    }
    WriteBoxData(basename,a_state,varnames,origin,a_dx);
};  

void  viewDataNC(BoxData<double, NUMCOMPS>* a_dataPtr)
{
  if(a_dataPtr != NULL)
  {
    WriteData<NUMCOMPS>(*a_dataPtr, -1, 1.0, string("var"), string("debugdat"));
    int sysret = system("visit -o debugdat.vtk");
    std::cout << "system call returned " << sysret << std::endl;;
  }
}

void printDataNC(BoxData<double, NUMCOMPS>* a_dataPtr, int icomp)
{
  std::cout    << setprecision(9)
               << setiosflags(ios::showpoint)
               << setiosflags(ios::scientific);
  if(a_dataPtr != NULL)
  {
    Point ivlo = a_dataPtr->box().low();
    Point ivhi = a_dataPtr->box().high();
    for(int j = ivlo[1]; j <= ivhi[1]; j++)
    {
      for(int i = ivlo[0]; i <= ivhi[0]; i++)
      {
        Point pt(i,j);
        std::cout << (*a_dataPtr)(pt, icomp) << "  ";
      }
      std::cout << std::endl;
    }
  }
}
/***/
void
parseCommandLine(double& a_tmax, int& a_nx, int& a_maxstep, int& a_outputinterval, int argc, char* argv[])
{
  cout << "Navier Stokes simulation of shear flow with sinusoidal perturbation.  Periodic bcs." << endl;
  cout << "usage:  " << argv[0] << " -n nx  -t tmax -m maxstep -o output_interval" << endl;
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
/***/
int main(int argc, char* argv[])
{
  //have to do this to get a time table
  PR_TIMER_SETFILE("proto.time.table");
  viewDataNC(NULL);
  printDataNC(NULL, 0);
  
  PR_TIME("main");
  double tstop;
  int size1D, maxStep, outputInterval;
  parseCommandLine(tstop, size1D, maxStep, outputInterval, argc, argv);
  
  // IMPORTANT: if convTest = false, does a single grid run with imput Grid size size1D. 
  // if convTest = true, performs a Richardson error run at fixed time step  
  // for grid sizes size1D, 2*size1D, 4 size1D, for maxStep = input, 2*input, 4*input.
  // in the latter case, only the final solutions, errors in density are output to .vtk fines, and
  // the max norm of the errors themselves order of accuracy output to cout.

  bool convTest = true;
  
  double gamma = 1.4;
  EulerRK4Op::s_count = 0;
  BoxData<double,NUMCOMPS> err[2],U[3];
  int maxLev = 1;
  if (convTest)
    {
      int sizeLev = size1D;
      tstop = 1.e+10;
      maxLev = 3;
      for (int lev = 0; lev < 3 ; lev++)
        {
          Point lo = Point::Zeros();
          Point hi = Point::Ones(sizeLev - 1);
          Box dbx0(lo,hi);
          U[lev].define(dbx0);
          if (lev < 2) err[lev].define(dbx0);
          sizeLev *=2;
        }    
    }
  int sizeLev = size1D;
  for (int lev = 0; lev < maxLev;lev++)
    {
      Point lo = Point::Zeros();
      Point hi = Point::Ones(sizeLev - 1);
      Box dbx0(lo,hi);
      double dx = 1./sizeLev;
      double dt = .25/sizeLev;
      EulerState state(dbx0,dx,gamma);
      RK4<EulerState,EulerRK4Op,EulerDX> rk4;
      EulerOp::initializeState(state.m_U,dx,gamma);

      double time = 0.;
      string resStr = "_"+std::to_string(size1D);
      string fileRoot = "outfile";
      cout << "starting time loop, maxStep = "<< maxStep << endl;
      Reduction<double>& rxn = state.m_Rxn;
      for (int k = 0;(k < maxStep) && (time < tstop);k++)
        {
          rk4.advance(time,dt,state);
          time += dt;
          if (!convTest)
            {
              dt = min(1.1*dt,.8/size1D/rxn.fetch());
              rxn.reset();
              cout <<"nstep = " << k << " time = " << time << " time step = " << dt << endl;
              if((outputInterval > 0) && (k%outputInterval == 0))
                {
                  WriteData(k,state.m_U,dx);
                }
            }
        }
      if (convTest)
        {
          (state.m_U).copyTo(U[lev]);
          sizeLev *= 2;
          maxStep *= 2;
        }
    }
  if (convTest)
    {
      Reduction<double> rxn;
      double rhoErrMax[2];
      for (int lev = 0; lev < maxLev-1; lev++)
        {
          rxn.reset();
          err[lev] |= Stencil<double>::AvgDown(2)(U[lev+1]);
          err[lev] -= U[lev];
          int size = U[lev].box().size(0);
          auto errRho = slice(err[lev],0);
          auto rhoLev = slice(U[lev],0);
	  errRho.absMax(rxn);
	  rhoErrMax[lev] = rxn.fetch();
          //rhoErrMax[lev] = errRho.absMax();
          string errStr = "error_"+std::to_string(size);
          string rhoStr = "rho_"+std::to_string(size);
          WriteBoxData(errStr.c_str(),errRho,1./size);
          WriteBoxData(rhoStr.c_str(),rhoLev,1./size);
          std::cout << "size = " << size << " , max error = " << rhoErrMax[lev] << std::endl;
        } 
      double rate = log(abs(rhoErrMax[0]/rhoErrMax[1]))/log(2.0);
      std::cout << "order of accuracy = " << rate<< std::endl;
    }    
  PR_TIMER_REPORT();
}
