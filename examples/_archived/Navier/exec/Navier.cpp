#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>

#include "Proto.H"
#include "Proto_Timer.H"
#include "Proto_DebugHooks.H"
#include "GodunovAdvectionOp.H"
#include "BCG_Integrator.H"
#include "Proto_WriteBoxData.H"
#include "CommonTemplates.H"

#define PI 3.141592653589793

using namespace std;
using namespace Proto;

class RunParams
{
public:
  RunParams()
  {
    nstepmax = 10;
    nx       = 128;
    outinterv= -10;
    tmax     = 2.0e10;
    domsize  = 1.0;
//    vortrad = 0.375*domsize;
    vortrad = 0.125*domsize;
    vortloc = 0.5*domsize;
    cfl      = 0.1;
    viscosity= 0.0002;
    presiter = 1;
    resetDx();
  }

  int nstepmax;
  int presiter;
  int nx;
  int outinterv;
  double viscosity;
  double tmax;
  double domsize;
  double vortrad;
  double vortloc;
  double dx;
  double cfl;
  void resetDx()
  {
    dx = domsize/nx;
    vortrad = 0.125*domsize;
    vortloc = 0.5*domsize;
  }

  void coarsen(int a_refrat)
  {
    dx *= a_refrat;
    nx /= a_refrat;;
    
  }
  void print() const
  {

    cout << "parameters: "                           << endl;
    cout << "nx                  =  "   << nx        << endl;
    cout << "output interval     =  "   << outinterv << endl;
    cout << "nstepmax            =  "   << nstepmax  << endl;
    cout << "tmax                =  "   << tmax      << endl;
    cout << "domain size         =  "   << domsize   << endl;
    cout << "viscosity           =  "   << viscosity << endl;
    cout << "pressure iterations =  "   << presiter  << endl;
    cout << "CFL number          =  "   << cfl       << endl;
   }
};                  
////////////
void
parseCommandLine(RunParams& a_params, int argc, char* argv[])
{
  cout << "Navier Stokes simulation weird litttle vortex.  Periodic bcs." << endl;
  cout << "usage:  " << argv[0] << " -n nx -i pressure_iterations -c cfl -t tmax -m maxstep  -v viscosity -d domain_size -o output_interval" << endl;
  for(int iarg = 0; iarg < argc-1; iarg++)
  {
    if(strcmp(argv[iarg],"-n") == 0)
    {
      a_params.nx = atoi(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg], "-i") == 0)
    {
      a_params.presiter = atoi(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg], "-m") == 0)
    {
      a_params.nstepmax = atoi(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg], "-o") == 0)
    {
      a_params.outinterv = atoi(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg],"-t") == 0)
    {
      a_params.tmax = atof(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg],"-v") == 0)
    {
      a_params.viscosity = atof(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg],"-c") == 0)
    {
      a_params.cfl = atof(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg],"-d") == 0)
    {
      a_params.domsize = atof(argv[iarg+1]);
    }
  }
  a_params.resetDx();
  a_params.print();
}


CUDA_DECORATION
double getVel(double x[DIM], int idir)
{
  double sinx = sin(2.*PI*x[0]);
  double siny = sin(2.*PI*x[1]);
  double retval = 0;
  if(idir == 0)
  {
    retval = siny;
  }
  else
  {
    retval = sinx;
  }
  return retval;
}


CUDA_DECORATION
double getGrad(double x[DIM], int idir)
{
  double sinx = 2.*PI*sin(2.*PI*x[0]);
  double siny = 2.*PI*sin(2.*PI*x[1]);
  double retval = 0;
  if(idir == 0)
  {
    retval = sinx;
  }
  else
  {
    retval = siny;
  }
  return retval;
}

PROTO_KERNEL_START void IncrementVelWithGradF(Point&            a_p,
                                              Pvector         & a_U,
                                              RunParams         a_params)
{
  
  double x[DIM];
  for(int jdir = 0; jdir <DIM; jdir++)
  {
    x[jdir] = (a_p[jdir] + 0.5)*a_params.dx;
  }

  double sinx = sin(2.*PI*x[0]);
  double siny = sin(2.*PI*x[1]);
  a_U(0) += getGrad(x, 0);
  a_U(1) += getGrad(x, 1);
}
PROTO_KERNEL_END(IncrementVelWithGradF, IncrementVelWithGrad)


PROTO_KERNEL_START void InitializeVelSinesF(Point&            a_p,
                                            Pvector         & a_U,
                                            RunParams         a_params)
{
  
  double x[DIM];
  for(int jdir = 0; jdir <DIM; jdir++)
  {
    x[jdir] = (a_p[jdir] + 0.5)*a_params.dx;
    a_U(jdir) = 0.;
  }

  double sinx = sin(2.*PI*x[0]);
  double siny = sin(2.*PI*x[1]);
  a_U(0) = getVel(x, 0);
  a_U(1) = getVel(x, 1);
}
PROTO_KERNEL_END(InitializeVelSinesF, InitializeVelSines)


PROTO_KERNEL_START void ExactGradSinesF(Point&            a_p,
                                        Pvector         & a_U,
                                        RunParams         a_params)
{
  
  double x[DIM];
  for(int jdir = 0; jdir <DIM; jdir++)
  {
    x[jdir] = (a_p[jdir] + 0.5)*a_params.dx;
    a_U(jdir) = 0.;
  }

  double sinx = sin(2.*PI*x[0]);
  double siny = sin(2.*PI*x[1]);
  a_U(0) = getGrad(x, 0);
  a_U(1) = getGrad(x, 1);
}
PROTO_KERNEL_END(ExactGradSinesF, ExactGradSines)


PROTO_KERNEL_START void ExactUdelUSinesF(Point&            a_p,
                                         Pvector         & a_udelu,
                                         RunParams         a_params)
{
  
  double x[DIM];
  for(int jdir = 0; jdir <DIM; jdir++)
  {
    x[jdir] = (a_p[jdir] + 0.5)*a_params.dx;
    a_udelu(jdir) = 0.;
  }

  double sinx = sin(2.*PI*x[0]);
  double siny = sin(2.*PI*x[1]);
  double cosx = cos(2.*PI*x[0]);
  double cosy = cos(2.*PI*x[1]);

  a_udelu(0) = 2.*PI*sinx*cosy;
  a_udelu(1) = 2.*PI*siny*cosx;
}
PROTO_KERNEL_END(ExactUdelUSinesF, ExactUdelUSines)

///
void getDtAndMaxWave(double                     & a_dt,
                     double                     & a_maxwave,
                     const BoxData<double, DIM> & a_velocity,
                     RunParams                    a_params)
{
  a_maxwave = 0;
  for(int idir = 0; idir < DIM; idir++)
  {
    //get the maximum wave speed
    a_maxwave = std::max(a_maxwave, a_velocity.absMax(idir));
  }
  //just in case it comes out zero.
  a_maxwave = std::max(a_maxwave, 1.0e-1);
  a_dt = a_params.cfl*a_params.dx/a_maxwave;
}
///
void
compareErrorMaxNorm(const BoxData<double, DIM>& a_errrFine, 
                    const BoxData<double, DIM>& a_errrCoar)
{
  cout << " comparing errr using max norm: " << endl;
  for(int ivar = 0; ivar < DIM; ivar++)
  {
    double coarnorm = a_errrCoar.absMax(ivar);
    double finenorm = a_errrFine.absMax(ivar);
    double order = 0;
    if(std::abs(finenorm) > 1.0e-12)
    {
      order = log(std::abs(coarnorm/finenorm))/log(2.0);
    }
    cout << "var = " << ivar << ", coarnorm = " << coarnorm << ", finenorm = "  << finenorm << ", order = " << order << endl;
  }     
}
void
enforceBoundaryConditions(BoxData<double, DIM>& a_phi, 
                          int a_numghost)
{
  for(int idir = 0; idir < DIM; idir++)
  {
    protocommon::enforceSGBoundaryConditions<double, DIM>(a_phi, a_numghost, idir);
  }
};

///
void udeluConvergence(const RunParams& a_params)
{
  int nghost = 3;
  Point lo = Point::Zeros();
  Point hi = Point::Ones(a_params.nx - 1);
  Box domainFine(lo, hi);
  Box domainCoar = domainFine.coarsen(2);

  RunParams paramFine = a_params;
  RunParams paramCoar = a_params;
  paramCoar.coarsen(2);
  

  
  Box ghostFine = domainFine.grow(nghost);
  Box ghostCoar = domainCoar.grow(nghost);
  //define and initialize velocity 
  BoxData<double, DIM> veloFine(ghostFine);
  BoxData<double, DIM> veloCoar(ghostFine);

  forallInPlace_p(InitializeVelSines, ghostFine, veloFine, paramFine);
  forallInPlace_p(InitializeVelSines, ghostCoar, veloCoar, paramCoar);

  BoxData<double, DIM> calcFine(domainFine);
  BoxData<double, DIM> calcCoar(domainCoar);
  BoxData<double, DIM> exacFine(domainFine);
  BoxData<double, DIM> exacCoar(domainCoar);

  forallInPlace_p(ExactUdelUSines, domainFine, exacFine, paramFine);
  forallInPlace_p(ExactUdelUSines, domainCoar, exacCoar, paramCoar);


  double dt = 0;  //do not want a real dt here.

  BCG_Integrator opFine(domainFine, paramFine.dx, paramFine.viscosity, nghost);
  BCG_Integrator opCoar(domainCoar, paramCoar.dx, paramCoar.viscosity, nghost);

  opFine.getUDotDelU(calcFine, veloFine, dt);
  opCoar.getUDotDelU(calcCoar, veloCoar, dt);

  BoxData<double, DIM> errrFine(domainFine);
  BoxData<double, DIM> errrCoar(domainCoar);
  errrFine.setVal(0.);
  errrCoar.setVal(0.);
  errrFine  += calcFine;
  errrCoar  += calcCoar;
  errrFine  -= exacFine;
  errrCoar  -= exacCoar;

  compareErrorMaxNorm(errrFine, errrCoar);
  WriteData<DIM>(errrCoar, -1, paramCoar.dx, string("erudu"), string("errrCoar"));
  WriteData<DIM>(errrFine, -1, paramFine.dx, string("erudu"), string("errrFine"));
  WriteData<DIM>(exacCoar, -1, paramCoar.dx, string("udelu"), string("exacCoar"));
  WriteData<DIM>(exacFine, -1, paramFine.dx, string("udelu"), string("exacFine"));
  WriteData<DIM>(calcCoar, -1, paramCoar.dx, string("udelu"), string("calcCoar"));
  WriteData<DIM>(calcFine, -1, paramFine.dx, string("udelu"), string("calcFine"));


  cout << "finished udelu convergence test" << endl;
}


///
void projectionConvergence(const RunParams& a_params)
{
  int nghost = 3;
  Point lo = Point::Zeros();
  Point hi = Point::Ones(a_params.nx - 1);
  Box domainFine(lo, hi);
  Box domainCoar = domainFine.coarsen(2);

  RunParams paramFine = a_params;
  RunParams paramCoar = a_params;
  paramCoar.coarsen(2);
  

  
  Box ghostFine = domainFine.grow(nghost);
  Box ghostCoar = domainCoar.grow(nghost);
  //define and initialize velocity 
  BoxData<double, DIM> calcVeloFiGh(ghostFine);
  BoxData<double, DIM> calcVeloCoGh(ghostCoar);
  BoxData<double, DIM> exacVeloFine(domainFine);
  BoxData<double, DIM> exacVeloCoar(domainCoar);

  BoxData<double, DIM> calcGradFiGh(ghostFine);
  BoxData<double, DIM> calcGradCoGh(ghostCoar);
  BoxData<double, DIM> exacGradFine(domainFine);
  BoxData<double, DIM> exacGradCoar(domainCoar);


  //initialize velocity as exact + grad scalar
  forallInPlace_p(InitializeVelSines, ghostFine, calcVeloFiGh, paramFine);
  forallInPlace_p(InitializeVelSines, ghostCoar, calcVeloCoGh, paramCoar);


  
  //put exact and grad scalar into their holders
  forallInPlace_p(IncrementVelWithGrad, ghostFine, calcVeloFiGh, paramFine);
  forallInPlace_p(IncrementVelWithGrad, ghostCoar, calcVeloCoGh, paramCoar);

  forallInPlace_p(ExactGradSines, domainFine, exacGradFine, paramFine);
  forallInPlace_p(ExactGradSines, domainCoar, exacGradCoar, paramCoar);


  forallInPlace_p(InitializeVelSines, domainFine, exacVeloFine, paramFine);
  forallInPlace_p(InitializeVelSines, domainCoar, exacVeloCoar, paramCoar);

  BCG_Integrator opFine(domainFine, paramFine.dx, paramFine.viscosity, nghost);
  BCG_Integrator opCoar(domainCoar, paramCoar.dx, paramCoar.viscosity, nghost);

  opCoar.ccProject(calcVeloCoGh, calcGradCoGh);
  opFine.ccProject(calcVeloFiGh, calcGradFiGh);

  BoxData<double, DIM> calcVeloFine(domainFine);
  BoxData<double, DIM> calcVeloCoar(domainCoar);
  BoxData<double, DIM> calcGradFine(domainFine);
  BoxData<double, DIM> calcGradCoar(domainCoar);

  calcVeloFiGh.copyTo(calcVeloFine);
  calcVeloCoGh.copyTo(calcVeloCoar);
  calcGradFiGh.copyTo(calcGradFine);
  calcGradCoGh.copyTo(calcGradCoar);

  BoxData<double, DIM> errrFine(domainFine);
  BoxData<double, DIM> errrCoar(domainCoar);
  errrFine.setVal(0.);
  errrCoar.setVal(0.);
  errrFine  += calcVeloFine;
  errrCoar  += calcVeloCoar;
  errrFine  -= exacVeloFine;
  errrCoar  -= exacVeloCoar;

  cout << "comparing error in velocity" << endl;
  compareErrorMaxNorm(errrFine, errrCoar);
  WriteData<DIM>(errrCoar, -1, paramCoar.dx, string("ervel"), string("errrVeloCoar"));
  WriteData<DIM>(errrFine, -1, paramFine.dx, string("ervel"), string("errrVeloFine"));
  WriteData<DIM>(exacVeloCoar, -1, paramCoar.dx, string("veloc"), string("exacVeloCoar"));
  WriteData<DIM>(exacVeloFine, -1, paramFine.dx, string("veloc"), string("exacVeloFine"));
  WriteData<DIM>(calcVeloCoar, -1, paramCoar.dx, string("veloc"), string("calcVeloCoar"));
  WriteData<DIM>(calcVeloFine, -1, paramFine.dx, string("veloc"), string("calcVeloFine"));


  errrFine.setVal(0.);
  errrCoar.setVal(0.);
  errrFine  += calcGradFine;
  errrCoar  += calcGradCoar;
  errrFine  -= exacGradFine;
  errrCoar  -= exacGradCoar;

  cout << "comparing error in gradient of scalar" << endl;
  compareErrorMaxNorm(errrFine, errrCoar);
  WriteData<DIM>(errrCoar, -1, paramCoar.dx, string("ergrad"), string("errrGradCoar"));
  WriteData<DIM>(errrFine, -1, paramFine.dx, string("ergrad"), string("errrGradFine"));
  WriteData<DIM>(exacGradCoar, -1, paramCoar.dx, string("gradsc"), string("exacGradCoar"));
  WriteData<DIM>(exacGradFine, -1, paramFine.dx, string("gradsc"), string("exacGradFine"));
  WriteData<DIM>(calcGradCoar, -1, paramCoar.dx, string("gradsc"), string("calcGradCoar"));
  WriteData<DIM>(calcGradFine, -1, paramFine.dx, string("gradsc"), string("calcGradFine"));
  cout << "finished projection convergence test" << endl;
}
///
void navierRun(const RunParams& a_params)
{
  int nghost = 3;
  Point lo = Point::Zeros();
  Point hi = Point::Ones(a_params.nx - 1);
  Box domain(lo, hi);
  Box ghostBox = domain.grow(nghost);
  

  //define and initialize velocity 
  BoxData<double, DIM> gradpres(ghostBox);
  BoxData<double, DIM> velocity(ghostBox);
  BoxData<double, DIM> scratch(domain);
                              
  RunParams params = a_params;
  forallInPlace_p(InitializeVelSines, ghostBox, velocity, params);
#if DIM==2
  BoxData<double,   1> vorticity(domain);
#else
  BoxData<double, DIM> vorticity(domain);
#endif
  BoxData<double, DIM> gradpress(ghostBox);
  
  //run the simulation
  int   nstep = 0; 
  double time = 0;
  BCG_Integrator op(domain, a_params.dx, a_params.viscosity, nghost);

  cout << "initial projection" << endl;
  op.ccProject(velocity, gradpress);
  gradpress.setVal(0.);
  cout << "iterating for initial pressure" << endl;
  double maxwave = 0;
  double dt = 0;
  getDtAndMaxWave(dt, maxwave, velocity, a_params);

//  for(int iiter = 0; iiter < a_params.presiter; iiter++)
  if(1)
  {
    velocity.copyTo(scratch);
    op.advanceSolution(velocity, gradpress, dt);
    scratch.copyTo(velocity);
  }
  cout << "starting navier stokes run" << endl;
  if(a_params.outinterv >= 0)
  {
    velocity.copyTo(scratch);
    WriteData<DIM>(scratch , nstep, a_params.dx, string("vel"), string("velocity"));
    gradpress.copyTo(scratch);
    WriteData<DIM>(scratch, nstep, a_params.dx, string("gradp"), string("gradp"));
    op.computeVorticity(vorticity, velocity);
#if DIM==2
    WriteData<1>(vorticity, nstep, a_params.dx, string("vort"), string("vorticity"));
#else
    WriteData<DIM>(vorticity, nstep, a_params.dx, string("vort"), string("vorticity"));
#endif
  }

  while((nstep < a_params.nstepmax) && (time < a_params.tmax))
  {

    cout << "start of step = "     << nstep << ", old time = " << time << ", maxwave = " << maxwave << " dt = "  << dt  << endl; 

    enforceBoundaryConditions(velocity, nghost);

    op.advanceSolution(velocity, gradpress, dt);
    op.computeVorticity(vorticity, velocity);
    
    nstep++;
    time += dt;
    getDtAndMaxWave(dt, maxwave, velocity, a_params);

    if((a_params.outinterv >= 0) && ((nstep %a_params.outinterv == 0) || (nstep == a_params.nstepmax-1)))
    {
      velocity.copyTo(scratch);
      WriteData<DIM>(scratch , nstep, a_params.dx, string("vel"), string("velocity"));
      gradpress.copyTo(scratch);
      WriteData<DIM>(scratch, nstep, a_params.dx, string("gradp"), string("gradp"));
      op.computeVorticity(vorticity, velocity);
#if DIM==2
      WriteData<1>(vorticity, nstep, a_params.dx, string("vort"), string("vorticity"));
#else
      WriteData<DIM>(vorticity, nstep, a_params.dx, string("vort"), string("vorticity"));
#endif
    }
  }

  cout << "finished navier run" << endl;

}


////////////
int main(int a_argc, char* a_argv[])
{
  //have to do this to get a time table
  PR_TIMER_SETFILE("proto.time.table");

  RunParams params;
  parseCommandLine(params, a_argc, a_argv);
  navierRun(params);
  //udeluConvergence(params);
  //projectionConvergence(params);

  PR_TIMER_REPORT();
}

