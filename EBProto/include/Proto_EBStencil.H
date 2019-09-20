
#ifndef _Proto_EBStencil_H_
#define _Proto_EBStencil_H_


#include "Proto_RealVect.H"
#include "Proto_Box.H"
#include "Proto_Point.H"
#include "Proto_LSquares.H"
#include "Proto_BaseIF.H"
#include "Proto_IrregNode.H"
#include "Proto_CutCellMoments.H"
#include "Proto_HostBoxData.H"
#include "Proto_HostIrregData.H"
#include "Proto_PointSet.H"
#include "Proto_LevelData.H"
#include "Proto_EBGraph.H"
#include "Proto_EBStencilArchive.H"


using std::shared_ptr;
using std::vector;

namespace Proto
{
  ///
  /**
  */
  template <int order, typename data_t, CENTERING srcCenter, CENTERING dstCenter>
  class EBStencil
  {

  public:
    typedef EBStencilArchive<srcCenter, dstCenter, order, data_t> archive_t;

    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    typedef HostIrregData<CELL    ,  IndMomDIM , 1>  VoluData;
    typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

    ///
    /**
       regular stencil gets applied in the interior of the domain.
       The aggstencil overwrites it near the EB.
       The domainBC stencils fill ghost cells outside the domain.  
       The vector length needs to be the same as the one-sided span of the stencil.
       
       If you are using periodic boundary conditions, those are enforced via LevelData::exchange
    */
    EBStencil(string              a_stencilName,
              string              a_domainBCName,
              string              a_ebbcName,
              const Box         & a_validBox,
              const Box         & a_domain,
              const Point       & a_srcGhost,
              const Point       & a_dstGhost,
              const EBGraph     & a_graph,
              const VoluData    & a_voludata,
              const EBFaData    & a_ebfadata,
              const XFacData    & a_xfacdata,
              const YFacData    & a_yfacdata,
              const ZFacData    & a_zfacdata, 
              double              a_dx,
              bool periodic_domain_bcs_only);


    ///
    inline
    void apply(EBBoxData<       dstCenter, data_t, 1>  &   a_dst,
               const EBBoxData< srcCenter, data_t, 1>  &   a_src,
               bool   a_initToZero = false,
               double a_scale      = 1)  const;
               
    
    shared_ptr<EBBoxData<dstCenter, data_t, 1> > getDiagonalWeights() const
    {
      return m_diagonalWeights;
    }
  protected:
    void
    computeDiagonalWeights(const vector< EBIndex<dstCenter> >              & a_dstVoFs,
                           const vector< LocalStencil<srcCenter, data_t> > & a_stencil,
                           const EBGraph                                   & a_graph,
                           const Box                                       & a_validBox);

    ///these 
    Stencil<double> m_regStencilInterior;
    //placeholder.  not sure what to do with domain bcs.
    //vector because there can be multiple ghost  cells
    vector<Stencil<double> >  m_regStencilBCS[2*DIM];
    vector<Box>               m_BCApplyBoxes[2*DIM];
    bool                      m_periodic_domain_bcs_only;
    Box                       m_regApplyBox;
    shared_ptr< AggStencil<srcCenter, dstCenter, data_t> >  m_irregStencil;
    shared_ptr<  EBBoxData<dstCenter, data_t, 1> >          m_diagonalWeights;


    private:
    EBStencil();
    EBStencil(     const EBStencil<order, data_t, srcCenter, dstCenter>& a_inp);
    void operator=(const EBStencil<order, data_t, srcCenter, dstCenter>& a_inp);


  };
}
#include "Proto_EBStencilImplem.H"
#endif
