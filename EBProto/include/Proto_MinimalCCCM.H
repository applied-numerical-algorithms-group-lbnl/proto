
#ifndef _Proto_MINIMALCCCM_H_
#define _Proto_MINIMALCCCM_H_

#include <map>
using std::map;
using std::vector;


#include "Proto_IndexTM.H"
#include "Proto_LSquares.H"
#include "Proto_IFData.H"
#include "Proto_CutCellMoments.H"


namespace Proto
{
  template <int dim> class LSProblem;

  template <int dim> class MinimalCCCM
  {
  public:
    typedef IndexTM<int,dim>  IvDim;
    typedef IndexTM<double,dim> RvDim;

    typedef map<IvDim,double >  PthMoment;

    typedef map<IndexTM<int,dim-1>,double > PthMomentLesserDimension;

    typedef map<IndexTM<int,1>,double >     OneDMoments;
    typedef map<int,IvDim>                LocPthMoment;
    typedef map<IvDim,int > PthMomentLoc;

    typedef map<Iv2,CutCellMoments<dim-1> > BdCutCellMoments;

    // Constructors
    MinimalCCCM();
    MinimalCCCM(const MinimalCCCM<dim>& a_MinimalCCCM);

    // This is used to build bd MinimalCCCM
    MinimalCCCM(const IFData<dim>& a_info);

    // Destructor
    ~MinimalCCCM();

    void computeMoments(const int              & a_orderPmax,
                        const int              & a_degreePmax);

    void computeMomentsRecursively(const int   & a_orderPmax,
                                   const int   & a_degreePmax);

    vector<double> computeRhs(LSProblem<dim> & a_lsp,
                              const int      & a_order);


    // Operators
    void operator=(const MinimalCCCM<dim>& a_MinimalCCCM);

    double factorial(const IvDim & a_multiIndex) const;

    // All the moments (and other data)
    CutCellMoments<dim> m_cutCellMoments;

    // Flag whether the boundary moments have been computed already
    bool m_boundaryMomentsComputed;
  };

// One dimensional MinimalCCCM
  template <> class MinimalCCCM<1>
  {
  public:
    typedef map<IndexTM<int,1>,double> OneDMoments;

    // Constructors
    MinimalCCCM();
    MinimalCCCM(const MinimalCCCM<1> & a_MinimalCCCM);
    MinimalCCCM(const IFData<1>& a_info);

    // Destructor
    ~MinimalCCCM();

    void computeMoments(const int              & a_orderPmax,
                        const int              & a_degreePmax);


    void simpleComputeMoments(const double & a_loPt,
                              const double & a_hiPt,
                              const int    & a_degreePmax);

    void computeMomentsUsingBinomial(const double & a_loPt,
                                     const double & a_hiPt,
                                     const int    & a_loSign,
                                     const int    & a_hiSign,
                                     const int    & a_degreePmax);

    // Operators
    void operator=(const MinimalCCCM<1>& a_MinimalCCCM);

    // Member data
    CutCellMoments<1> m_cutCellMoments;
  };


}

#include "Proto_MinimalCCCMImplem.H"

#endif
