
#ifndef _Proto_EBStencilARchive_H_
#define _Proto_EBStencilARchive_H_


#include "EBProto.H"

namespace Proto
{
  ///
  /**
  */
  template <CENTERING srcCenter, CENTERING dstCenter>
  class EBStencilArchive
  {
  public:
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    typedef IrregData<CELL,  IndMomDIM     , 1>  VoluData;
    typedef IrregData<CELL,  IndMomDIM     , 1>  EBFaData;
    typedef IrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef IrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef IrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

    EBArchive() {}
    ~EBArchive() {}
    
    /// 
    static void
    getStencil(vector<EBIndex<dstCenter> >     & a_dstVoFs,
               vector<LocalStencil<srcCenter> >& a_stencil,
               Stencil<double>                 & a_regStencilInterior,
               vector<Stencil<double> >          a_regStencilBCS[2*DIM],
               Box                             & a_regApplyBox,
               string              a_stencilName,
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
               bool                a_periodiconly);
  };
}
#include "Proto_EBStencilArchiveImplem.H"
#endif
