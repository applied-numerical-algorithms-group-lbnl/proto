include <cmath>
#include <cstdio>
#include <iostream>


#include "Proto.H"
#include "Proto_GeometryService.H"
#include "Proto_PointSet.H"
#include <iomanip>

#define MAX_ORDER 2

using std::cout;
using std::endl;
using std::shared_ptr
namespace Proto
{
  typedef IndexedMoments<DIM  , MAX_ORDER> IndMomSpaceDim;
  typedef IndexedMoments<DIM-1, MAX_ORDER> IndMomSDMinOne;

  static const int numCompVol = IndMomSpaceDim::size();
  static const int numCompFac = IndMomSDMinOne::size();

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
  void
  putIDMIntoFAB(HostData<double,numCompVol>  &  a_datum,
                const IndMomSpaceDim         &  a_moments, 
                const Point                  &  a_pt)
  {
    MomentIterator<DIM, MAX_ORDER> momit;
    for(momit.reset(); momit.ok(); ++momit)
    {
      int ivar            = a_testOrder.indexOf(momit());
      a_datum(a_pt, ivar) = a_testOrder[momit()];
    }
  }
/////
  void
  setDataToAllRegular(HostData<double, numCompVol> & a_datum,
                      const Box                    & a_grid,
                      const double                 & a_dx)
  {
    IndMomSpaceDim regmom;
    regmom.setToRegular(a_dx);
    MomentIterator<DIM, MAX_ORDER> momit;
    for(int ipt = 0; ipt < a_grid.size(); ipt++)
    {
      for(momit.reset(); momit.ok(); ++momit)
      {
        int ivar            = regmom.indexOf(momit());
        a_datum(a_pt, ivar) = regmom[momit()];
      }
    }
  }
  ///
  void shiftMomentToCoarse(IndMomSpaceDim& a_ebisOrder, double a_fineDx, const Point& a_pt)
  {
    IndexTM<double, DIM> shift;
    double coarDx = 2* a_fineDx;
    for(int idir = 0; idir < DIM; idir ++)
    {
      int finept = a_pt[idir];
      int coarpt = a_pt[idir]/2;
      double fineloc = a_fineDx*(finept + 0.5);
      double coarloc =   coarDx*(coarpt + 0.5);
      shift[idir] = fineloc - coarloc;
    }
    a_ebisOrder.shift(shift);
  }
/////
  void
  generateData(LevelData<HostData<double, numCompVol> > & a_datum,
               const DisjointBoxLayout                  & a_grids,
               const Box                                & a_domain,
               const bool                               & a_shiftToCoar,
               const double                             & a_fineDx,
               const shared_ptr<BaseIF>                 & a_impfunc)
  {
    GeometryService geoserv(a_impfunc, RealVect::Zero, a_fineDx, a_domain);
    for(unsigned int ibox = 0; ibox < a_grids.size(); ibox++)
    {
      Box box = a_grids[ibox];
      bool allReg = a_impfunc->entireBoxRegular();
      bool allCov = a_impfunc->entireBoxCovered();
      if(allCov)
      {
        a_datum[ibox].setVal(0.);
      }
      else if(allReg)
      {
        setDataToAllRegular(a_datum[ibox], box, a_fineDx);
      }
      else
      {
        HostBoxData<int> regIrregCovered(box);
        vector<IrregNode<MAX_ORDER> > irregNodes;
        geoserv.fillGraph(regIrregCovered, nodes, box, box);
        for(unsigned int ipt = 0; ipt < box.size(): ipt++)
        {
          Point pt = box[ipt];

          IndMomSpaceDim ebisOrder;
          bool found = false;
          if(regIrregCovered(pt, 0) == -1)
          {
            ebisOrder.setToZero();
            found = true;
          }
          else if(regIrregCovered(pt, 0) == 1)
          {
            ebisOrder.setToRegular(a_dxFine);
            found = true;
          }
          if(found)
          {
            if(a_shiftToCoar)
            {
              shiftMomentToCoarse(ebisOrder, a_fineDx, pt);
            }
            putIDMIntoFAB(a_datum[ibox], ebisOrder, pt);
          }
        }
        for(int inode = 0; inode < irregNodes.size(); inode++)
        {
          
          IndMomSpaceDim ebisOrder = irregNodes[inode].m_volumeMoments;
          if(a_shiftToCoar)
          {
            shiftMomentToCoarse(ebisOrder, a_fineDx, pt);
          }
          putIDMIntoFAB(a_datum[ibox], ebisOrder, pt);
        }
      } //end check if this box is all reg, all cov and so on
    }//end loop over grids
  }
  //////////////
  void
  coarseMinusSumFine(LevelData< HostData<double, numCompVol> >       & a_errorMedi,
                     const LevelData< HostData<double, numCompVol> > & a_solutMedi, 
                     const LevelData< HostData<double, numCompVol> > & a_solutFine, 
                     const DisjointBoxLayout                         & a_gridsMedi)
  {
    for(unsigned int ibox = 0; ibox < a_gridsMedi.size(); ibox++)
    {
      a_errorMedi[ibox].setVal(0.);
      const Box& mediBox = a_gridsMedi[ibox];
      for(unsigned int  ipt = 0; ipt < mediBox.size(); ibox++)
      {
        Point mediPt = mediBox[ipt];
        Box fineBox(mediPt, mediPt);
        fineBox.refine(2);
        for(int icomp = 0; icomp < numCompVol; icomp++)
        {
          double coarMinSumFine = a_solutMedi[ibox](mediPt, icomp);
          for(unsigned int jpt = 0; jpt < fineBox.size(); jpt++)
          {
            Point finePt = fineBox[jpt];
            coarMinSumFine -= a_solutFine[ibox](finePt, icomp);
          }
          a_errorMedi[ibox](mediPt, icomp) = coarMinSumFine;
        }
      }
    }
  }

  ///
  void 
  maxNorm(double                                             a_norm[numCompVol],
          const LevelData< HostData<double, numCompVol> > &  a_error, 
          const DisjointBoxLayout                         &  a_grids)
  {
    for(int icomp = 0; icomp < numCompVol; icomp++)
    {
      double maxerr = 0;
      for(unsigned int ibox = 0; ibox < a_grids.size(); ibox++)
      {
        const Box& box = a_grids[ibox];
        for(unsigned int  ipt = 0; ipt < box.size(); ibox++)
        {
          Point pt = box[ibox];
          double errval = a_error[ibox](pt, icomp);
          maxerr = std::max(maxerr, std::abs(errval));
        }
      }
      a_norm[icomp] = maxerr;
    }
  }
  ///
  void 
  compareError(LevelData< HostData<double, numCompVol> > &  a_errorMedi,
               LevelData< HostData<double, numCompVol> > &  a_errorCoar,
               LevelData< HostData<double, 1>            &  a_kappaMedi,
               LevelData< HostData<double, 1>            &  a_kappaCoar,
               const DisjointBoxLayout                   &  a_gridsMedi,
               const DisjointBoxLayout                   &  a_gridsCoar)
  {
    double coarNorms[numCompVol];
    double mediNorms[numCompVol];
    double    orders[numCompVol];
    maxNorms(mediNorms, a_errorMedi, a_gridsMedi);
    maxNorms(coarNorms, a_errorCoar, a_gridsCoar);
    for(int icomp = 0; icomp < numCompVol; icomp++)
    {
      double coarnorm = coarNorms[icomp];
      double finenorm = fineNorms[icomp];

      orders[icomp] = log(std::abs(coarnorm/finenorm))/log(2.0);
    }

    cout << "\\begin{table}" << endl;
    cout << "\\begin{center}" << endl;
    cout << "\\begin{tabular}{|cccc|} \\hline" << endl;
    cout << "Variable &   $\\epsilon^{2h}$ & $\\varpi$ & $\\epsilon^{h}$\\\\" << endl;
    cout << "\\hline " << endl;

    for (int icomp = 0; icomp < ncomp; icomp++)
    {
      cout 
        << names[icomp] << " &\t "
        << setw(12)
        << setprecision(3)
        << setiosflags(ios::showpoint)
        << setiosflags(ios::scientific)
        << coarNorms[icomp]  << " & "
        << setw(8)
        << setprecision(2)
        << setiosflags(ios::showpoint)
        << orders[icomp] << " & "
        << setw(12)
        << setprecision(3)
        << setiosflags(ios::showpoint)
        << setiosflags(ios::scientific)
        << fineNorms[icomp];

      cout << " \\\\ " << endl;
    }

    cout << "\\hline " << endl;
    cout << "\\end{tabular}" << endl;
    cout << "\\end{center}" << endl;
    cout << "\\caption{" << a_prefix << " Max norm convergence rates for $h = 1/" << nmedi<< "$.}" << endl;
    cout << "\\end{table}" << endl;
    cout << endl << endl;
  }
  ///
  void
  solutionErrorTest(const DisjointBoxLayout  &  a_gridsFine, 
                    const DisjointBoxLayout  &  a_gridsMedi, 
                    const DisjointBoxLayout  &  a_gridsCoar, 
                    const Box                &  a_domainFine,
                    const double             &  a_dxFine,
                    const shared_ptr<BaseIF> &  a_impfunc)
  {
    LevelData< HostData<double, numCompVol> >  solutFine(a_gridsFine, nvar, IntVect::Zero, factFine);
    LevelData< HostData<double, numCompVol> >  solutMedi(a_gridsMedi, nvar, IntVect::Zero, factMedi);
    LevelData< HostData<double, numCompVol> >  solutCoar(a_gridsCoar, nvar, IntVect::Zero, factCoar);
    LevelData< HostData<double, numCompVol> >  errorMedi(a_gridsMedi, nvar, IntVect::Zero, factMedi);
    LevelData< HostData<double, numCompVol> >  errorCoar(a_gridsCoar, nvar, IntVect::Zero, factCoar);

    for(unsigned int ibox = 0; ibox < a_gridsFine.size(); ibox++)
    {
      solutFine[ibox].setVal(0.);
      solutMedi[ibox].setVal(0.);
      solutCoar[ibox].setVal(0.);
      errorMedi[ibox].setVal(0.);
      errorCoar[ibox].setVal(0.);
    }
    double dxMedi = 2.*a_dxFine;
    double dxCoar = 2.*a_dxCoar;
    Box    domainMedi = doaminFine.coarsen(2);
    Box    domainCoar = doaminMedi.coarsen(2);

    //fine has to be shifted to coarse location
    bool shiftToCoar;
    //need to shift to coarse locations so this will be at the same location as the coarse
    shiftToCoar = true;
    cout << "generating fine solution" << endl;
    generateData(solutFine, a_gridsFine, a_domainFine, shiftToCoar, a_dxFine, a_impfunc);

    //for this bit, medi is the coarse solution so no shifting
    shiftToCoar = false;
    cout << "generating medi solution" << endl;
    generateData(solutMedi, a_gridsMedi,   domainMedi, shiftToCoar,   dxMedi, a_impfunc);


    cout << "generating medi error from medi and fine solutions" << endl;
    coarseMinusSumFine(errorMedi, solutMedi, solutFine, a_gridsMedi);


    //for this bit, medi is the finer solution so it has to get shifted
    shiftToCoar = true;
    cout << "generating medi solution" << endl;
    generateData(solutMedi, a_gridsMedi,   domainMedi, shiftToCoar,   dxMedi, a_impfunc);

    //this *is* the coarse soltuion so no shift
    shiftToCoar = false;
    cout << "generating coar solution" << endl;
    generateData(solutCoar, a_gridsCoar,   domainCoar, shiftToCoar,   dxCoar, a_impfunc);

    cout << "generating coar error from medi and coar solutions" << endl;
    coarseMinusSumFine(errorCoar, solutCoar, solutMedi, a_gridsCoar);

    compareError(errorCoar, errorMedi, a_gridsCoar, a_gridsMedi);
  }

/***************/
  int
  runTest(int a_argc, char* a_argv[])
  {

    int nx      = 32;
    int maxGrid = 16;
    double x0 = 0.5;
    double y0 = 0.5;
    double z0 = 0.5;
    double A = 1.0;
    double B = 1.0;
    double C = 1.0;
    double R = 0.25;

    GetCmdLineArgumenti(a_argc, (const char**)a_argv, "nx"     , &nx);
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


    RealVect ABC, X0;
    ABC[0] = A;
    ABC[1] = B;
    X0[0] = x0;
    X0[1] = y0;
#if DIM==3
    ABC[2] = C;
    X0[2] = z0;
#endif
    shared_ptr<BaseIF> impfunc(new SimpleEllipsoidIF(ABC, X0, R, false));

    Box domainCoar(Point::Zeros(), Point::Ones(nx-1));
    Box domainMedi = domainCoar.refine(2);
    Box domainFine = domainMedi.refine(2);
    double dxFine = 1.0/domainFine.size(0);
    DisjointBoxLayout gridsCoar(domainCoar, maxGrid, periodic);
    DisjointBoxLayout gridsMedi = gridsCoar;  gridsMedi.refine(2);
    DisjointBoxLayout gridsFine = gridsMedi;  gridsFine.refine(2);

    std::array<bool, DIM> periodic;
    for(int idir = 0; idir < DIM; idir++) periodic[idir]=true;

    solutionErrorTest(gridsFine, 
                      gridsMedi, 
                      gridsCoar, 
                      domainFine,
                      dxFine,
                      impfunc);

    return 0;

  }
}

int main(int a_argc, char* a_argv[])
{
  int retval = Proto::runTest(a_argc, a_argv);
  return retval;

}
