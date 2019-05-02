
#ifndef _Proto_EBStencilARchive_H_
#define _Proto_EBStencilARchive_H_


#include "EBProto.H"

namespace Proto
{
  ///
  /**
  */
  template <CENTERING srcCenter, CENTERING dstCenter, int order, typename data_t>
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

    EBStencilArchive() {}
    ~EBStencilArchive() {}
    
    /// 
    inline static void
    getStencil(vector<EBIndex<dstCenter> >              & a_dstVoFs,
               vector<LocalStencil<srcCenter, data_t> > & a_stencil,
               Stencil<data_t>                          & a_regStencilInterior,
               vector<Stencil<data_t> >                   a_regStencilBCS[2*DIM],
               Box                                      & a_regApplyBox,
               string                                     a_stencilName,
               string                                     a_domainBCName,
               string                                     a_ebbcName,
               const Box                                & a_validBox,
               const Box                                & a_domain,
               const Point                              & a_srcGhost,
               const Point                              & a_dstGhost,
               const EBGraph                            & a_graph,
               const VoluData                           & a_voludata,
               const EBFaData                           & a_ebfadata,
               const XFacData                           & a_xfacdata,
               const YFacData                           & a_yfacdata,
               const ZFacData                           & a_zfacdata,
               data_t                                     a_dx,
               bool                                       a_periodiconly);

  };


  ///
  template<CENTERING facecent, typename data_t>
  class Base2ndOrderFlux
  {
  public:
    Base2ndOrderFlux()
    {
    }

    virtual ~Base2ndOrderFlux()
    {
    }

    void getFaceInterpolationStencil(LocalStencil<facecent, data_t>              & a_facesten,
                                     const EBIndex<facecent>                     & a_face,
                                     const EBGraph                               & a_graph,
                                     const IrregData<facecent, IndMomSDMinOne, 1>& a_facedat);


    void getCentroidStencil(LocalStencil<CELL, data_t>                  & a_facesten,
                            const EBIndex<facecent>                     & a_face,
                            const EBGraph                               & a_graph,
                            const IrregData<facecent, IndMomSDMinOne, 1>& a_facedat,
                            data_t                                        a_dx);

    virtual void 
    getFaceCenteredFluxStencil(LocalStencil<CELL, data_t>                  & a_facesten,
                               const EBIndex<facecent>                     & a_face,
                               const EBGraph                               & a_graph,
                               const IrregData<facecent, IndMomSDMinOne, 1>& a_facedat,
                               data_t                                        a_dx) = 0;
  };
  
  ///second order Poisson flux
  template<CENTERING facecent, typename data_t>
  class Poisson2ndOrder: public Base2ndOrderFlux<facecent, data_t>
  {
  public:
    inline void 
    getFaceCenteredFluxStencil(LocalStencil<CELL, data_t>                  & a_facesten,
                               const EBIndex<facecent>                     & a_face,
                               const EBGraph                               & a_graph,
                               const IrregData<facecent, IndMomSDMinOne, 1>& a_facedat,
                               data_t                                        a_dx);
 
    virtual ~Poisson2ndOrder()
    {
    }

    Poisson2ndOrder()
    {
    }
  };
  ///some basic infrastructure to support the second order stuff
  template<typename data_t>
  class SecondOrderStencil
  {
  public:
    inline static void
    get2ndOrderDivFStencil(LocalStencil< CELL,  data_t>          & a_vofStencil,
                           const EBIndex<CELL>                   & a_vof,
                           const EBGraph                         & a_graph,
                           const VoluData                        & a_voludata,
                           const EBFaData                        & a_ebfadata,
                           const XFacData                        & a_xfacdata,
                           const YFacData                        & a_yfacdata,
                           const ZFacData                        & a_zfacdata,
                           const Base2ndOrderFlux<XFACE, data_t> & a_xfaceflux,
                           const Base2ndOrderFlux<YFACE, data_t> & a_yfaceflux,
                           const Base2ndOrderFlux<ZFACE, data_t> & a_zfaceflux,
                           string                                  a_ebbcName,
                           data_t                                  a_dx);




    virtual ~FluxBase2ndOrder()
    {
    }

    FluxBase2ndOrder()
    {
    }
  };



}
#include "Proto_EBStencilArchiveImplem.H"
#endif
