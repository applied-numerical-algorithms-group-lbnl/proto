
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
  template<CENTERING cent,  int order>
  class Normalizor
  {
  public:
    typedef IndexTM<int   , SpaceDim>    IvSpaceDim;
    typedef IndexTM<int   , SpaceDim-1>  IvSDMinOne;
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef IrregData<CELL,  IndMomDIM     , 1>  VoluData;
    typedef IrregData<CELL,  IndMomDIM     , 1>  EBFaData;
    typedef IrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef IrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef IrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

    Normalizor()
    {
    }

    virtual ~Normalizor()
    {
    }

    ///
    inline static void 
    normalizedFaceMoments(double                                      & a_areaFrac,
                          RealVect                                  & a_centroid,
                          const IrregData<cent,     IndMomSDMinOne, 1>& a_facedat,
                          const EBIndex<cent>                         & a_face,
                          const EBGraph                               & a_graph,
                          const double                                & a_dx   );

    

    ///only makes sense for cent==cell but I had to put this somewhere
    inline static void 
    normalizedVoluMoments(double                                      & a_volFrac,
                          RealVect                                  & a_centroid,
                          const IrregData<cent, IndMomSpaceDim, 1>    & a_celldat,
                          const EBIndex<cent>                         & a_vof,
                          const EBGraph                               & a_graph,
                          const double                                & a_dx   )


    ///only makes sense for cent==boundary but I had to put this somewhere
    inline static void 
    normalizedBndryMoments(RealVect                                        & a_centroid,
                           const IrregData<cent, IndMomSpaceDim, 1>        & a_celldat,
                           const EBIndex<cent>                             & a_vof,
                           const EBGraph                                   & a_graph,
                           const double                                    & a_dx   );

    ///only makes sense for cent==cell but I had to put this somewhere
    inline static void 
    getNormalAndBoundaryArea(double                        & a_bndryArea,
                             RealVect                      & a_normal,
                             const VoluData                & a_voludata,
                             const XFacData                & a_xfacdata,
                             const YFacData                & a_yfacdata,
                             const ZFacData                & a_zfacdata,
                             const EBIndex<CELL>           & a_vof,
                             const EBGraph                 & a_graph,
                             const double                  & a_dx   );
  };


  template<CENTERING facecent, int order> //order is so we can send in geometry info
  class Base2ndOrderFlux
  {
  public:
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;

    Base2ndOrderFlux()
    {
    }

    virtual ~Base2ndOrderFlux()
    {
    }

    inline static
    void getFaceInterpolationStencil(LocalStencil<facecent, double>              & a_facesten,
                                     const EBIndex<facecent>                     & a_face,
                                     const EBGraph                               & a_graph,
                                     const IrregData<facecent, IndMomSDMinOne, 1>& a_facedat);


    inline static 
    void getCentroidStencil(LocalStencil<CELL, double>                  & a_facesten,
                            const EBIndex<facecent>                     & a_face,
                            const EBGraph                               & a_graph,
                            const IrregData<facecent, IndMomSDMinOne, 1>& a_facedat,
                            double                                        a_dx);

    virtual void 
    getFaceCenteredFluxStencil(LocalStencil<CELL, double>                  & a_facesten,
                               const EBIndex<facecent>                     & a_face,
                               const EBGraph                               & a_graph,
                               const IrregData<facecent, IndMomSDMinOne, 1>& a_facedat,
                               double                                        a_dx) = 0;



  };
  
  ///second order Poisson flux
  template<CENTERING facecent, int order>
  class Poisson2ndOrder: public Base2ndOrderFlux<facecent>
  {
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    typedef IrregData<CELL,  IndMomDIM     , 1>  VoluData;
    typedef IrregData<CELL,  IndMomDIM     , 1>  EBFaData;
    typedef IrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef IrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef IrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;
  public:
    inline void 
    getFaceCenteredFluxStencil(LocalStencil<CELL, double>                  & a_facesten,
                               const EBIndex<facecent>                     & a_face,
                               const EBGraph                               & a_graph,
                               const IrregData<facecent, IndMomSDMinOne, 1>& a_facedat,
                               double                                        a_dx);
 
    virtual ~Poisson2ndOrder()
    {
    }

    Poisson2ndOrder()
    {
    }
  };
  ///some basic infrastructure to support the second order stuff
  template<int order>  //this is so we can send in the geometry informaion sanely
  class SecondOrderStencil
  {
  public:
    inline static void
    get2ndOrderDivFStencil(LocalStencil< CELL,  double>          & a_vofStencil,
                           const EBIndex<CELL>                   & a_vof,
                           const EBGraph                         & a_graph,
                           const VoluData                        & a_voludata,
                           const EBFaData                        & a_ebfadata,
                           const XFacData                        & a_xfacdata,
                           const YFacData                        & a_yfacdata,
                           const ZFacData                        & a_zfacdata,
                           const Base2ndOrderFlux<XFACE, double> & a_xfaceflux,
                           const Base2ndOrderFlux<YFACE, double> & a_yfaceflux,
                           const Base2ndOrderFlux<ZFACE, double> & a_zfaceflux,
                           string                                  a_ebbcName,
                           double                                  a_dx);




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
