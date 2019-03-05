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
#include "GodunovAdvectionOp.H"
#include "CommonTemplates.H"
#include "Proto_WriteBoxData.H"

#define PI 3.141592653589793

using namespace std;
using namespace Proto;

class RunParams
{
public:
  RunParams()
  {
    nstepmax = 0;
    nx       = 64;
    outinterv= -10;
    tmax     = 1.0;
    blobrad  = 0.1;
    blobscal = 0.1;
    vortrad  = 0.325;
    vortloc  = 0.5;
    domsize  = 1.0;
    cfl      = 0.1;
    for(int idir = 0; idir < DIM; idir++)
    {
      blobloc[idir] = 0.;
    }
    blobloc[0] = 0.1;
    resetDx();
  }

  int nstepmax;
  int nx;
  int outinterv;
  double tmax;
  double vortrad;
  double vortloc;
  double blobrad;
  double domsize;
  double dx;
  double cfl;
  double blobloc[DIM];
  double blobscal;
  void resetDx()
  {
    dx = domsize/nx;

    for(int idir = 0; idir < DIM; idir++)
    {
      blobloc[idir] = domsize/2.;
    }
    blobloc[1] += blobscal;
  }

  void print() const
  {

    cout << "parameters: " << endl;
    cout << "nx        =  "   << nx << endl;
    cout << "output interv =  "   << outinterv << endl;
    cout << "nstepmax  =  "   << nstepmax << endl;
    cout << "tmax      =  "   << tmax << endl;
    cout << "blobrad   =  "   << blobrad << endl;
    cout << "vortrad   =  "   << vortrad << endl;
    cout << "vortloc   =  "   << vortloc << endl;
    cout << "blobscal   =  "   << blobscal << endl;
    cout << "dom size  =  "   << domsize << endl;
    cout << "cfl number=  "   << cfl     << endl;
    for(int idir = 0; idir < DIM; idir++)
    {
      cout << "blobloc[" << idir << "] = " << blobloc[idir] << endl;
    }
  }
};                  
////////////
void
parseCommandLine(RunParams& a_params, int argc, char* argv[])
{
  cout << "Godunov advection around center of domain with a made up vortex.  Blob starts on x axis.   Periodic bcs." << endl;
  cout << "usage:  " << argv[0] << " -n nx -c cfl -t tmax -m maxstep -b blobrad -d domain_size -l blobloc -o output_interval" << endl;
  for(int iarg = 0; iarg < argc-1; iarg++)
  {
    if(strcmp(argv[iarg],"-n") == 0)
    {
      a_params.nx = atoi(argv[iarg+1]);
      a_params.resetDx();
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
    else if(strcmp(argv[iarg],"-l") == 0)
    {
      a_params.blobscal = atof(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg],"-c") == 0)
    {
      a_params.cfl = atof(argv[iarg+1]);
    }
    else if(strcmp(argv[iarg],"-d") == 0)
    {
      a_params.domsize = atof(argv[iarg+1]);
      a_params.vortloc = a_params.domsize/2.;
      a_params.vortrad = 0.375*a_params.domsize;
      a_params.resetDx();
    }
    else if(strcmp(argv[iarg],"-b") == 0)
    {
      a_params.blobrad = atof(argv[iarg+1]);
    }
  }
  a_params.resetDx();
  a_params.print();
}


PROTO_KERNEL_START unsigned int InitializePhiF(Point&             a_p,
                                               Scalar           & a_phi,
                                               RunParams          a_params)
{
  double x[DIM];
  for(int idir = 0; idir < DIM; idir++)
  {    
    x[idir] = (a_p[idir] + 0.5)*a_params.dx;
  }
  double rad0sq = a_params.blobrad*a_params.blobrad;
  const double* const x0 = a_params.blobloc;
  double radsq = 0;
  for(int idir = 0; idir < DIM; idir++)
  {  
    radsq += (x[idir]-x0[idir])*(x[idir]- x0[idir]);
  }
  if(radsq < rad0sq)
  {
    double cosval = cos(0.5*PI*(radsq/rad0sq));
    double phival =  cosval*cosval;

    a_phi(0) = phival;
  }
  else
  {
    a_phi(0) = 0.0;
  }
//  a_phi(0) = a_p[0];
  return 0;
}
PROTO_KERNEL_END(InitializePhiF, InitializePhi)

PROTO_KERNEL_START unsigned int InitializeVelF(Point&           a_p,
                                               Scalar        & a_U,
                                               RunParams       a_params,
                                               int             a_idir)
{
  
  double xrel[DIM];
  double radsq = 0;
  double rad0sq = a_params.vortrad*a_params.vortrad;
  {
    double  x[DIM];
    double x0[DIM];
    for(int idir = 0; idir < DIM; idir++)
    {    
      x[idir] = (a_p[idir] + 0.5)*a_params.dx;
      x0[idir] = a_params.vortloc;
    }
    x[a_idir] = a_p[a_idir]*a_params.dx;


    for(int idir = 0; idir < DIM; idir++)
    {    
      xrel[idir] = x[idir] - x0[idir];
    }
    for(int idir = 0; idir < DIM; idir++)
    {  
      radsq += xrel[idir]*xrel[idir];
    }
  }

  double uval = 0;
//  if(1)
  if(radsq < rad0sq)
  {
    double cosval = cos(0.5*PI*(radsq/rad0sq));
    double phival = cosval*cosval;
    if(a_idir == 0)
    {
      uval = -(xrel[1])*phival;
    }
    else if(a_idir == 1)
    {
      uval =  xrel[0]*phival;
    }
  }
  else
  {
    uval = 0.0;
  }
  a_U(0) = uval;
  return 0;
}
PROTO_KERNEL_END(InitializeVelF, InitializeVel)

//cheerfully stolen from the euler example
void
enforceBoundaryConditions(BoxData<double, 1>& a_phi, 
                          const Box         & a_domain,        
                          const Box         & a_ghostBox,
                          int a_numghost)
{
  for(int idir = 0; idir < DIM; idir++)
  {
    protocommon::enforceSGBoundaryConditions<double, 1>(a_phi, a_numghost, idir);
  }
};
void godunovRun(const RunParams& a_params)
{
  int nghost = 3;
  Point lo = Point::Zeros();
  Point hi = Point::Ones(a_params.nx - 1);
  Box domain(lo, hi);
  Box ghostBox = domain.grow(nghost);

  //define and initialize scalar phi
  BoxData<double, 1> phi(ghostBox);
  
  RunParams params = a_params;
  forallInPlace_p(InitializePhi, ghostBox, phi, params);


  //define and initialize face-centered velocities
  BoxData<double, DIM> velCell(ghostBox);
  BoxData<double, 1> velFace[DIM];
  BoxData<double, 1> src(ghostBox);
  src.setVal(0.);
  for(int idir = 0; idir < DIM; idir++)
  {
    Box facedom = domain.extrude(idir);
    Box ghostFace = facedom.grow(nghost);
    velFace[idir].define(ghostFace);
    BoxData<double,1> velcomp = slice(velCell, idir);

    forallInPlace_p(InitializeVel, ghostFace, velFace[idir], params, idir);
    forallInPlace_p(InitializeVel, ghostBox,  velcomp      , params, idir);

    double facemax = velFace[idir].max();
    double cellmax = velCell.max(idir);
    double facemin = velFace[idir].min();
    double cellmin = velCell.min(idir);
    cout << idir << " max vel face = " << facemax << ", max vel cell = " << cellmax << endl;
    cout << idir << " min vel face = " << facemin << ", min vel cell = " << cellmin << endl;
  }
  
  double pmax = phi.max();
  double pmin = phi.min();
  cout << "after initialize,  phi max = " << pmax   << ", phi min  = " << pmin << endl;
  //get the maximum wave speed
  double maxwave = 0;
  for(int idir = 0; idir < DIM; idir++)
  {
    maxwave = std::max(maxwave, velFace[idir].absMax());
  }

  //just in case it comes out zero.
  maxwave = std::max(maxwave, 1.23456789e-10);

  double dt = a_params.cfl*a_params.dx/maxwave;

  cout << "maxwave = " << maxwave << ", dt = " << dt << endl;
  //run the simulation
  int   nstep = 0; 
  double time = 0;
  GodunovAdvectionOp op(a_params.dx);


  cout << "starting godunov advection run" << endl;
  if(a_params.outinterv >= 0)
  {
    BoxData<double,1> phiPrint(domain);
    phi.copyTo(phiPrint);
    WriteData<1>(phiPrint, 0, a_params.dx, string("phi"), string("phi"));

    BoxData<double,DIM> velPrint(domain);
    velCell.copyTo(velPrint);
    WriteData<DIM>(velPrint, -1, a_params.dx, string("vel"), string("velCell"));
    for(int idir = 0; idir < DIM; idir++)
    {
      Box facedom = domain.extrude(idir);
      BoxData<double,1> velFacePrint(domain);
      velFace[idir].copyTo(velFacePrint);
      WriteData<1>(velFacePrint, idir, a_params.dx, string("vel"), string("velFace"));
    }
  }
  while((nstep < a_params.nstepmax) && (time < a_params.tmax))
  {
    double max = phi.max();
    double min = phi.min();
    cout << "start of step = "     << nstep << ", old time = " << time 
         << ", phi max = " << max   << ", phi min  = " << min << endl;

    BoxData<double, 1> divF(ghostBox);
    bool doingVel = false;

    enforceBoundaryConditions(phi, domain, ghostBox, nghost);
    op.divFluxNPH(divF, velFace, phi, src, velCell, domain, doingVel, dt);
    
    divF *= (-dt);

    double divfmax = divF.max();
    double divfmin = divF.min();
    cout << "divfmax = " << divfmax << ", divfmin = " << divfmin << endl;

    phi += divF;


    nstep++;
    time += dt;

    if((a_params.outinterv >= 0) && (nstep %a_params.outinterv == 0))
    {
      BoxData<double,1> phiPrint(domain);
      phi.copyTo(phiPrint);
      WriteData<1>(phiPrint, nstep, a_params.dx, string("phi"), string("phi"));
    }
  }

  cout << "finished godunov advection run" << endl;

}
////////////
int main(int a_argc, char* a_argv[])
{
  //have to do this to get a time table
  PR_TIMER_SETFILE("proto.time.table");

  RunParams params;
  parseCommandLine(params, a_argc, a_argv);
  godunovRun(params);

  PR_TIMER_REPORT();
}

