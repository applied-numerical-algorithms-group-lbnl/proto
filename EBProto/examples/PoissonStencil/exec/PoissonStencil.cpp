#include <cmath>
#include <cstdio>
#include <iostream>


#include "EBProto.H"

#include <iomanip>

#define MAX_ORDER 2

using std::cout;
using std::endl;
using std::shared_ptr;
//why does this not work?
//static const int numCompVol = Proto::s_max_indmom_sizes[DIM-1][MAX_ORDER];
#if DIM==2
#define numCompVol 6
#else
#define numCompVol 10
#endif
using std::setw;
using std::setprecision;
using std::setiosflags;
using std::ios;
namespace Proto
{
/////
  void GetCmdLineArgumenti(int argc, const char** argv, const char* name, int* rtn)
  {
    size_t len = strlen(name);
    for(int i=1; i<argc; i+=2)
    {
      if(strcmp(argv[i]+1,name) ==0)
      {
        *rtn = atoi(argv[i+1]);
        break;
       }
    }
  }
  void GetCmdLineArgumentd(int argc, const char** argv, const char* name, double* rtn)
  {
    size_t len = strlen(name);
    for(int i=1; i<argc; i+=2)
    {
      if(strcmp(argv[i]+1,name) ==0)
      {
        *rtn = atof(argv[i+1]);
        break;
       }
    }
  }
/////

/***************/
  int
  runTest(int a_argc, char* a_argv[])
  {
    int nx      = 32;
    int maxGrid = 32;
    double x0 = 0.5;
    double y0 = 0.5;
    double z0 = 0.5;
    double A = 1.0;
    double B = 1.0;
    double C = 1.0;
    double R = 0.25;
    int nIter       = 10;
    int nStream    = 8;
    GetCmdLineArgumenti(a_argc, (const char**)a_argv, "nx"     , &nx);
    GetCmdLineArgumenti(a_argc, (const char**)a_argv, "niter"  , &nIter);
    GetCmdLineArgumenti(a_argc, (const char**)a_argv, "nstream", &nStream);
    GetCmdLineArgumenti(a_argc, (const char**)a_argv, "maxGrid", &maxGrid);
    GetCmdLineArgumentd(a_argc, (const char**)a_argv, "x0"     , &x0);
    GetCmdLineArgumentd(a_argc, (const char**)a_argv, "y0"     , &y0);
    GetCmdLineArgumentd(a_argc, (const char**)a_argv, "z0"     , &z0);
    GetCmdLineArgumentd(a_argc, (const char**)a_argv, "A"      , &A);
    GetCmdLineArgumentd(a_argc, (const char**)a_argv, "B"      , &B);
    GetCmdLineArgumentd(a_argc, (const char**)a_argv, "C"      , &C);
    GetCmdLineArgumentd(a_argc, (const char**)a_argv, "R"      , &R);         

    cout << "nx      = " << nx       << endl;
    cout << "maxGrid = " << maxGrid  << endl;
    cout << "x0      = " << x0       << endl;
    cout << "y0      = " << y0       << endl;
    cout << "z0      = " << z0       << endl;
    cout << "A       = " << A        << endl;
    cout << "B       = " << B        << endl;
    cout << "C       = " << C        << endl;
    cout << "R       = " << R        << endl;

    cout << "nIter   = " << nIter    << endl;
    cout << "nstream = " << nStream  << endl;

    RealVect ABC, X0;
    ABC[0] = A;
    ABC[1] = B;
    X0[0] = x0;
    X0[1] = y0;
#if DIM==3
    ABC[2] = C;
    X0[2] = z0;
#endif
    Box domain(Point::Zeros(), Point::Ones(nx-1));
    double dx = 1.0/domain.size(0);
    std::array<bool, DIM> periodic;
    for(int idir = 0; idir < DIM; idir++) periodic[idir]=true;
    DisjointBoxLayout grids(domain, maxGrid, periodic);
    Point dataGhost = Point::Ones(1);
    Point geomGhost = Point::Ones(2);
    RealVect origin = RealVect::Zero();
    shared_ptr<BaseIF>                       impfunc(new SimpleEllipsoidIF(ABC, X0, R, false));
    shared_ptr<GeometryService<MAX_ORDER> >  geoserv(new GeometryService<MAX_ORDER>(impfunc, origin, dx, domain, grids, geomGhost, 0));

    EBDictionary<2, double, CELL, CELL> dictionary(geoserv, grids, dataGhost, dataGhost, dx, true);
    typedef EBStencil<2, double, CELL, CELL> ebstencil_t;
    string stenname("Second_Order_Poisson");
    string dombcname("Periodic");
    string  ebbcname("Neumann");

    dictionary.registerStencil(stenname, dombcname, ebbcname);
    LevelData<EBBoxData<CELL,  double, 1> > srcData(grids, dataGhost);
    LevelData<EBBoxData<CELL,  double, 1> > dstData(grids, dataGhost);
    for(int ibox = 0; ibox < grids.size(); ibox++)
    {
      double val = 0;
      EBBoxData<CELL, double, 1>& srcebbd = srcData[ibox];
      EBBoxData<CELL, double, 1>& dstebbd = dstData[ibox];
      srcebbd.setVal(val);
      dstebbd.setVal(val);

      srcData[ibox].setVal(0.);
      dstData[ibox].setVal(0.);
    }

    for(int iiter = 0; iiter < nIter; iiter++)
    {    
      for(int ibox = 0; ibox < grids.size(); ibox++)
      {
        shared_ptr<ebstencil_t> stencil = dictionary.getEBStencil(stenname, ebbcname, ibox);
        stencil->apply(dstData[ibox], srcData[ibox]);
      }
    }

  }
}

int main(int a_argc, char* a_argv[])
{
  int retval = Proto::runTest(a_argc, a_argv);
  return retval;

}
