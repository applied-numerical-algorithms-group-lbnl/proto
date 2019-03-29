include <cmath>
#include <cstdio>
#include <iostream>


#include "Proto.H"
#include "Proto_GeometryService.H"
#include "Proto_PointSet.H"

#define MAX_ORDER 2

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
        std::cout<<name<<"="<<" "<<*rtn<<std::endl;
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
                      const Box                 & a_grid,
                      const double              & a_dx)
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

/*****/
  void
  solutionErrorTest(const GridParameters &  a_paramFine, 
                    const GridParameters &  a_paramMedi, 
                    const GridParameters &  a_paramCoar,
                    const DisjointBoxLayout &  a_gridsFine, 
                    const DisjointBoxLayout &  a_gridsMedi, 
                    const DisjointBoxLayout &  a_gridsCoar, 
                    const EBISLayout        &  a_ebislFine, 
                    const EBISLayout        &  a_ebislMedi, 
                    const EBISLayout        &  a_ebislCoar, 
                    const EBISLayout        &  a_ebislCoFi, 
                    const EBISLayout        &  a_ebislCoMe, 
                    const TestType          &  a_type, 
                    const int               &  a_idir)
  {
    IndMomSpaceDim idmproxy;
    int nvar = idmproxy.size();

    string prefix;
    if(a_type == VOL_MOM)
    {
      prefix = string("volume_moment");
    }
    else if(a_type == EB_MOM)
    {
      prefix = string("eb_moment");
    }
    else if(a_type == EB_NORMAL_MOM)
    {
      prefix =  string("ebNormalMoment_") + convertInt(a_idir);
    }
    else
    {
      MayDay::Error("bogus type");
    }

    EBCellFactory       factFine(a_ebislFine);
    EBCellFactory       factMedi(a_ebislMedi);
    EBCellFactory       factCoar(a_ebislCoar);

    LevelData<EBCellFAB>  solutFine(a_gridsFine, nvar, IntVect::Zero, factFine);
    LevelData<EBCellFAB>  solutMedi(a_gridsMedi, nvar, IntVect::Zero, factMedi);
    LevelData<EBCellFAB>  solutCoar(a_gridsCoar, nvar, IntVect::Zero, factCoar);
    LevelData<EBCellFAB>  errorMedi(a_gridsMedi, nvar, IntVect::Zero, factMedi);
    LevelData<EBCellFAB>  errorCoar(a_gridsCoar, nvar, IntVect::Zero, factCoar);

    EBLevelDataOps::setToZero(solutFine);
    EBLevelDataOps::setToZero(solutMedi);
    EBLevelDataOps::setToZero(solutCoar);

    //fine has to be shifted to coarse location
    bool shiftToCoar;
    //need to shift to coarse locations so this will be at the same location as the coarse
    shiftToCoar = true;
    pout() << "generating fine solution" << endl;
    generateData(solutFine, a_gridsFine, a_ebislFine, a_paramFine, a_type, a_idir, shiftToCoar);

    //for this bit, medi is the coarse solution so no shifting
    shiftToCoar = false;
    pout() << "generating medi solution" << endl;
    generateData(solutMedi, a_gridsMedi, a_ebislMedi, a_paramMedi, a_type, a_idir, shiftToCoar);

    pout() << "generating medi error from medi and fine solutions" << endl;
    sumFineMinusCoarse(errorMedi,
                       solutMedi, a_gridsMedi, a_ebislCoFi, a_paramMedi,
                       solutFine, a_gridsFine, a_ebislFine, a_paramFine, a_type);

    //for this bit, medi is the finer solution so it has to get shifted
    shiftToCoar = true;
    pout() << "generating medi solution" << endl;
    generateData(solutMedi, a_gridsMedi, a_ebislMedi, a_paramMedi, a_type, a_idir, shiftToCoar);

    //this *is* the coarse soltuion so no shift
    shiftToCoar = false;
    pout() << "generating coar solution" << endl;
    generateData(solutCoar, a_gridsCoar, a_ebislCoar, a_paramCoar, a_type, a_idir, shiftToCoar);

    pout() << "generating coar error from medi and coar solutions" << endl;
    sumFineMinusCoarse(errorCoar,
                       solutCoar, a_gridsCoar, a_ebislCoMe, a_paramCoar,
                       solutMedi, a_gridsMedi, a_ebislMedi, a_paramMedi, a_type);


    Vector<Real> orders;
    Box domCoar = a_paramCoar.coarsestDomain.domainBox();
    if(a_type == VOL_MOM)
    {
  
      Vector<string> names(nvar);
      getMomentNames<CH_EBIS_ORDER>(names, string("m"));
      //the 1 is a verbosity flag.  leave it at one.  trust me.  
      EBArith::compareError(orders, errorMedi, errorCoar,
                            a_gridsMedi, a_gridsCoar, 
                            a_ebislMedi, a_ebislCoar,
                            domCoar, 1, NULL, names, prefix);
    }
    else if((a_type == EB_MOM) || (a_type == EB_NORMAL_MOM))
    {
      
      Vector<string> names(nvar);
      getMomentNames<CH_EBIS_ORDER>(names, string("b"));
      BaseIVFactory<Real> bivrFactMedi(a_ebislMedi);
      BaseIVFactory<Real> bivrFactCoar(a_ebislCoar);
      LevelData<BaseIVFAB<Real> > sparseErrorMedi(a_gridsMedi, nvar, IntVect::Zero, bivrFactMedi);
      LevelData<BaseIVFAB<Real> > sparseErrorCoar(a_gridsCoar, nvar, IntVect::Zero, bivrFactCoar);

      copyDenseToSparse(sparseErrorMedi, errorMedi);
      copyDenseToSparse(sparseErrorCoar, errorCoar);
      EBArith::compareIrregError(orders, sparseErrorMedi, sparseErrorCoar,
                                 a_gridsMedi, a_gridsCoar, 
                                 a_ebislMedi, a_ebislCoar,
                                 domCoar, prefix, names);
    }
    else
    {
      MayDay::Error("bogus type");
    }
    /**/
    pout() << "Outputting moments to file" << endl;
    string solutFileFine =  prefix + string("_Fine.hdf5");
    string solutFileMedi =  prefix + string("_Medi.hdf5");
    string solutFileCoar =  prefix + string("_Coar.hdf5");
    writeEBLevelName(solutFine, solutFileFine);
    writeEBLevelName(solutMedi, solutFileMedi);
    writeEBLevelName(solutCoar, solutFileCoar);

    pout() << "Outputting error to file" << endl;
    string errorFileMedi =  prefix + string("_Error_Medi.hdf5");
    string errorFileCoar =  prefix + string("_Error_Coar.hdf5");
    writeEBLevelName(errorCoar, errorFileCoar);
    writeEBLevelName(errorMedi, errorFileMedi);
    /**/

  }

/***************/
  int
  runTest(int a_argc, char* a_argv[])
  {
#ifdef CH_MPI
    MPI_Init(&a_argc,&a_argv);
#endif
    {

      // Check for an input file
      char* inFile = NULL;

 

      if (a_argc > 1)
      {
        inFile = a_argv[1];
      }
      else
      {
        pout() << "Usage: <executable name> <inputfile>" << endl;
        pout() << "No input file specified" << endl;
        return -1;
      }
      ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

    
      GridParameters    paramFine, paramMedi, paramCoar;
      DisjointBoxLayout gridsFine, gridsMedi, gridsCoar;
      //read params from file
      getGridParameters(paramFine, true);
      paramMedi = paramFine;
      paramMedi.coarsen(2);
      paramCoar = paramMedi;
      paramCoar.coarsen(2);


      Vector<int> procs;
      Vector<Box> boxes; 
      domainSplit(paramFine.coarsestDomain, boxes,
                  paramFine.maxGridSize   , paramFine.blockFactor);
      LoadBalance(procs, boxes);
      gridsFine = DisjointBoxLayout(boxes, procs, paramFine.coarsestDomain);
      coarsen(gridsMedi, gridsFine, 2);
      coarsen(gridsCoar, gridsMedi, 2);

      pout() << "rct: defining FINE geometry" << endl;
      definePoissonGeometry(paramFine);
      pout() << "saving fine geometry into eblg" << endl;
      EBLevelGrid eblgFine(gridsFine, paramFine.coarsestDomain, 2, Chombo_EBIS::instance());
      pout() << "making CoFi info into eblg" << endl;
      EBLevelGrid eblgCoFi(gridsMedi, paramMedi.coarsestDomain, 2, Chombo_EBIS::instance());

      barrier();
      pout() << "clearing singleton" << endl;
      Chombo_EBIS::instance()->clear();
      pout() << "rct: defining MEDI geometry" << endl;
      definePoissonGeometry(paramMedi);
      pout() << "saving medi geometry into eblg" << endl;
      EBLevelGrid eblgMedi(gridsMedi, paramMedi.coarsestDomain, 2, Chombo_EBIS::instance());
      pout() << "making CoMe info into eblg" << endl;
      EBLevelGrid eblgCoMe(gridsCoar, paramCoar.coarsestDomain, 2, Chombo_EBIS::instance());

      barrier();
      pout() << "clearing singleton" << endl;
      Chombo_EBIS::instance()->clear();
      pout() << "rct: defining Coar geometry" << endl;
      definePoissonGeometry(paramCoar);
      pout() << "saving medi geometry into eblg" << endl;
      EBLevelGrid eblgCoar(gridsCoar, paramCoar.coarsestDomain, 2, Chombo_EBIS::instance());

      EBISLayout ebislFine = eblgFine.getEBISL();
      EBISLayout ebislMedi = eblgMedi.getEBISL();
      EBISLayout ebislCoar = eblgCoar.getEBISL();
      EBISLayout ebislCoFi = eblgCoFi.getEBISL();
      EBISLayout ebislCoMe = eblgCoMe.getEBISL();
      //all thse calls to geometry because we are testing the 
      //accuracy of geometry generation
      //the CoFi and CoMe stuff is because you cannot call refine on
      //coar and medi stuff because as far as it is concerned, it is the 
      //finest level.  They also might have slightly different graphs so this finesses that
      //problem as well
      barrier();

      pout() << "test of volume moments " << endl;
      solutionErrorTest(paramFine,paramMedi,paramCoar,
                        gridsFine,gridsMedi,gridsCoar,
                        ebislFine,ebislMedi,ebislCoar,
                        ebislCoFi, ebislCoMe,
                        VOL_MOM, 0);

      pout() << "test eb area moments" << endl;
      solutionErrorTest(paramFine,paramMedi,paramCoar,
                        gridsFine,gridsMedi,gridsCoar,
                        ebislFine,ebislMedi,ebislCoar,
                        ebislCoFi, ebislCoMe,
                        EB_MOM , 0);

      pout() << "test eb normal moments" << endl;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        solutionErrorTest(paramFine,paramMedi,paramCoar,
                          gridsFine,gridsMedi,gridsCoar,
                          ebislFine,ebislMedi,ebislCoar,
                          ebislCoFi, ebislCoMe,
                          EB_NORMAL_MOM, idir);
      }

      pout() << "clearing singleton " << endl;
      Chombo_EBIS::instance()->clear();
    }
#ifdef CH_MPI
    CH_TIMER_REPORT();
    MPI_Finalize();
#endif

  }

}
int main(int a_argc, char* a_argv[])
{
  int retval = Proto::runTest(a_argc, a_argv);
  return retval;

}
