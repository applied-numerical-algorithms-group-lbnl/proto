
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
    EBStencilArchive() {}
    ~EBStencilArchive() {}
    
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    typedef HostIrregData<CELL,  IndMomDIM     , 1>  VoluData;
    typedef HostIrregData<CELL,  IndMomDIM     , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

    ///   
    inline static void
    getStencil(vector<EBIndex<dstCenter> >               & a_dstVoFs,
               vector<LocalStencil<srcCenter, data_t> >  & a_stencil,
               Stencil<double>                           & a_regStencilInterior,
               vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
               Box                                       & a_regApplyBox,
               string                                      a_stencilName,
               string                                      a_domainBCName,
               string                                      a_ebbcName,
               const Box                                 & a_validBox,
               const Box                                 & a_domain,
               const Point                               & a_srcGhost,
               const Point                               & a_dstGhost,
               const EBGraph                             & a_graph,
               const VoluData                            & a_voludata,
               const EBFaData                            & a_ebfadata,
               const XFacData                            & a_xfacdata,
               const YFacData                            & a_yfacdata,
               const ZFacData                            & a_zfacdata,
               data_t                                      a_dx,
               bool                                        a_periodiconly);


 
  };


  ///
  /**
   */
  template <int order>
  class EBStencilArchive<CELL, CELL, order, double>
  {
  public:
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    typedef HostIrregData<CELL,  IndMomDIM     , 1>  VoluData;
    typedef HostIrregData<CELL,  IndMomDIM     , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

    EBStencilArchive() {}
    ~EBStencilArchive() {}
    
    ///   
    inline static void
    getStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,
               vector<LocalStencil<CELL, double> >       & a_stencil,
               Stencil<double>                           & a_regStencilInterior,
               vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
               Box                                       & a_regApplyBox,
               string                                      a_stencilName,
               string                                      a_domainBCName,
               string                                      a_ebbcName,
               const Box                                 & a_validBox,
               const Box                                 & a_domain,
               const Point                               & a_srcGhost,
               const Point                               & a_dstGhost,
               const EBGraph                             & a_graph,
               const VoluData                            & a_voludata,
               const EBFaData                            & a_ebfadata,
               const XFacData                            & a_xfacdata,
               const YFacData                            & a_yfacdata,
               const ZFacData                            & a_zfacdata,
               double                                      a_dx,
               bool                                        a_periodiconly);

    inline static void
    getSecondOrderPoissonStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                                 vector<LocalStencil<CELL, double> >       & a_stencil,                    
                                 Stencil<double>                           & a_regStencilInterior,         
                                 vector<Stencil<double> >                    a_regStencilBCS[2*DIM],       
                                 Box                                       & a_regApplyBox,                
                                 string                                      a_stencilName,                
                                 string                                      a_domainBCName,               
                                 string                                      a_ebbcName,                   
                                 const Box                                 & a_validBox,                   
                                 const Box                                 & a_domain,                     
                                 const Point                               & a_srcGhost,                   
                                 const Point                               & a_dstGhost,                   
                                 const EBGraph                             & a_graph,                      
                                 const VoluData                            & a_voludata,                   
                                 const EBFaData                            & a_ebfadata,                   
                                 const XFacData                            & a_xfacdata,                   
                                 const YFacData                            & a_yfacdata,                   
                                 const ZFacData                            & a_zfacdata,                   
                                 double                                      a_dx,                         
                                 bool                                        a_periodiconly);               
  };


  ///
  namespace Normalizor
  {
    typedef IndexTM<int  , DIM>                  IvDIM;
    typedef IndexTM<int  , DIM-1>                IvSDMinOne;

    ///
    template <CENTERING cent, typename data_t, int order>
    static void
    normalizedFaceMoments(data_t                                                      & a_areaFrac,
                          RealVect                                                    & a_centroid,
                          const HostIrregData<cent, IndexedMoments<DIM-1, order>, 1 > & a_facedat,
                          const EBIndex<cent>                                         & a_face,
                          const EBGraph                                               & a_graph,
                          data_t                                                        a_dx   )
    {


      EBIndex<CELL> voflo = a_face.getVoF(Side::Lo);
      EBIndex<CELL> vofhi = a_face.getVoF(Side::Hi);
      if(a_graph.isRegular(voflo.m_pt) || a_graph.isRegular(vofhi.m_pt))
      {
        a_areaFrac = 1;
        a_centroid  = RealVect::Zero();
      }
      else if(a_graph.isCovered(voflo.m_pt) || a_graph.isCovered(vofhi.m_pt))
      {
        a_areaFrac = 0;
        a_centroid  = RealVect::Zero();
      }
      else
      {
        double fullfacearea = 1.;
        for(int idir = 0; idir < DIM-1; idir++) fullfacearea *= a_dx;
        const IndexedMoments<DIM-1, order> & facemom = a_facedat(a_face, 0); 
        double area = facemom[IvSDMinOne::Zero];
        a_areaFrac = area/fullfacearea;
        a_centroid = RealVect::Zero();
        if(area > 0)
        {
          int iindex = 0;
          for (int idir = 0; idir < DIM; ++idir)
          {
            if(idir != a_face.direction())
            {
              IvSDMinOne mono = BASISV_TM<int,DIM-1>(iindex);
              a_centroid[idir] = facemom[mono];
              //normalize  by area
              a_centroid[idir] /= area;
              //normalize by dx
              a_centroid[idir] /= a_dx;
              iindex++;
            }
          }
        }
      }
    }

    

    template<int order, typename data_t>
    inline static void
    normalizedVoluMoments(data_t                                                      & a_volFrac,
                          RealVect                                                    & a_centroid,
                          const HostIrregData<CELL, IndexedMoments<DIM  , order> , 1> & a_celldat,
                          const EBIndex<CELL>                                         & a_vof,
                          const EBGraph                                               & a_graph,
                          data_t                                                        a_dx   )
    {

      if(a_graph.isRegular(a_vof.m_pt))
      {
        a_volFrac = 1;
        a_centroid  = RealVect::Zero();
      }
      else if(a_graph.isCovered(a_vof.m_pt))
      {
        a_volFrac = 0;
        a_centroid  = RealVect::Zero();
      }
      else
      {
        double fullcellvol = 1.;
        for(int idir = 0; idir < DIM; idir++) fullcellvol *= a_dx;
        const IndexedMoments<DIM  , order> & moments = a_celldat(a_vof, 0); 
        double volume = moments[IvDIM::Zero];
        a_volFrac = volume/fullcellvol;
        a_centroid = RealVect::Zero();
        if(volume > 0)
        {
          for (int idir = 0; idir < DIM; ++idir)
          {
            IvDIM mono = BASISV_TM<int,DIM>(idir);
            a_centroid[idir] = moments[mono];
            //normalize  by volume
            a_centroid[idir] /= volume;
            //normalize by dx
            a_centroid[idir] /= a_dx;
          }
        }
      }
    }


    ///
    template <typename data_t, int order>
    inline void
    normalizedBndryMoments(RealVect                                                       & a_centroid,
                           const HostIrregData<BOUNDARY, IndexedMoments<DIM , order>,1 >  & a_celldat,
                           const EBIndex<BOUNDARY>                                        & a_vof,
                           const EBGraph                                                  & a_graph,
                           data_t                                                           a_dx   )
    {

      if(a_graph.isRegular(a_vof.m_pt))
      {
        a_centroid  = RealVect::Zero();
      }
      else if(a_graph.isCovered(a_vof.m_pt))
      {
        a_centroid  = RealVect::Zero();
      }
      else
      {
        double fullcellvol = 1.;
        for(int idir = 0; idir < DIM; idir++) fullcellvol *= a_dx;
        double fullarea = fullcellvol/a_dx;
        const IndexedMoments<DIM,order> & moments = a_celldat(a_vof, 0); 
        double area = moments[IvDIM::Zero];

        a_centroid = RealVect::Zero();
        if(area > 0)
        {
          for (int idir = 0; idir < DIM; ++idir)
          {
            IvDIM mono = BASISV_TM<int,DIM>(idir);
            a_centroid[idir] = moments[mono];
            //normalize  by area
            a_centroid[idir] /= area;
            //normalize by dx
            a_centroid[idir] /= a_dx;
          }
        }
      }
    }


//    ///only makes sense for cent==cell but I had to put this somewhere
//    inline static void 
//    getNormalAndBoundaryArea(double                        & a_bndryArea,
//                             RealVect                      & a_normal,
//                             const VoluData                & a_voludata,
//                             const XFacData                & a_xfacdata,
//                             const YFacData                & a_yfacdata,
//                             const ZFacData                & a_zfacdata,
//                             const EBIndex<CELL>           & a_vof,
//                             const EBGraph                 & a_graph,
//                             data_t                          a_dx   );
  }


  template<CENTERING facecent, int order> //order is so we can send in geometry info
  class Base2ndOrderFlux
  {
  public:
    typedef IndexedMoments<DIM  , order>         IndMomDIM;
    typedef IndexedMoments<DIM-1, order>         IndMomSDMinOne;
    typedef HostIrregData<CELL,  IndMomDIM     , 1>  VoluData;
    typedef HostIrregData<CELL,  IndMomDIM     , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

    Base2ndOrderFlux()
    {
    }

    virtual ~Base2ndOrderFlux()
    {
    }

    inline void
    getFaceInterpStencil(LocalStencil<facecent, double>                   &  a_sten,
                         const EBIndex<facecent>                          &  a_face,
                         const EBGraph                                    &  a_graph,
                         const HostIrregData<facecent, IndMomSDMinOne, 1> &  a_facdata,
                         double                                              a_dx) const;

    inline void
    getCentroidStencil(LocalStencil<CELL, double>                         &  a_facesten,
                       const EBIndex<facecent>                            &  a_face,
                       const EBGraph                                      &  a_graph,
                       const HostIrregData<facecent, IndMomSDMinOne, 1>   &  a_facedat,
                       double                                                a_dx) const;
 
    virtual void 
    getFaceCenteredFluxStencil(LocalStencil<CELL, double>                      & a_facesten,
                               const EBIndex<facecent>                         & a_face,
                               const EBGraph                                   & a_graph,
                               const HostIrregData<facecent, IndMomSDMinOne, 1>& a_facedat,
                               double                                            a_dx) const = 0;



  };
  
  ///second order Poisson flux
  template<CENTERING facecent, int order>
  class Poisson2ndOrder: public Base2ndOrderFlux<facecent, order>
  {
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    typedef HostIrregData<CELL,  IndMomDIM     , 1>  VoluData;
    typedef HostIrregData<CELL,  IndMomDIM     , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;
  public:
    inline void 
    getFaceCenteredFluxStencil(LocalStencil<CELL, double>                  & a_facesten,
                               const EBIndex<facecent>                     & a_face,
                               const EBGraph                               & a_graph,
                               const HostIrregData<facecent, IndMomSDMinOne, 1>& a_facedat,
                               double                                        a_dx) const;
 
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
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    typedef HostIrregData<CELL,  IndMomDIM     , 1>  VoluData;
    typedef HostIrregData<CELL,  IndMomDIM     , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;
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
                           const Base2ndOrderFlux<XFACE, order>  & a_xfaceflux,
                           const Base2ndOrderFlux<YFACE, order>  & a_yfaceflux,
                           const Base2ndOrderFlux<ZFACE, order>  & a_zfaceflux,
                           string                                  a_ebbcName,
                           double                                  a_dx);




    virtual ~SecondOrderStencil()
    {
    }

    SecondOrderStencil()
    {
    }
  };



}
#include "Proto_EBStencilArchiveImplem.H"
#endif
