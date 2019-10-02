
#ifndef _Proto_EBStencilArchive_H_
#define _Proto_EBStencilArchive_H_


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
    typedef HostIrregData<CELL,     IndMomDIM  , 1>  VoluData;
    typedef HostIrregData<BOUNDARY, IndMomDIM  , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

    ///   returns true if there are domain boundary conditions (not all ops have them)
    inline static bool
    getStencil(vector<EBIndex<dstCenter> >               & a_dstVoFs,
               vector<LocalStencil<srcCenter,  data_t> > & a_stencil,
               Stencil<data_t>                           & a_regStencilInterior,
               vector<Stencil<data_t> >                    a_regStencilBCS[2*DIM],
               vector<Box>                                 a_BCApplyBoxes[2*DIM],       
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
               bool                                        a_periodiconly)               ;


 
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
    typedef HostIrregData<CELL,      IndMomDIM , 1>  VoluData;
    typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

    EBStencilArchive() {}
    ~EBStencilArchive() {}
    
    ///for prolongation and multicolor relaxation
    inline static Point getColor(unsigned int icolor)
    {
      Point retval;
      // Color offsets are grouped into "red"=even number of nonzeros (first 2^(DIM-1)) 
      // and "black= odd number of nonzeros (the rest).
#if DIM==2
      if(icolor == 0)
      {
        retval = Point::Zeros();//(0,0)
      }
      else if(icolor == 1)
      {
        retval = Point::Ones();//(1,1)
      }
      else if(icolor == 2)
      {
        retval = Point::Zeros() + Point::Basis(1);//(0,1)
      }
      else if(icolor == 3)
      {
        retval = Point::Zeros() + Point::Basis(0);//(1,0)
      }
      else
      {
        PROTO_ASSERT(false, "color not found");
      }
#elif DIM==3
      if(icolor == 0)
      {
        retval = Point::Zeros();//(0,0,0)
      }
      else if(icolor == 1)
      {
        retval = Point::Zeros() + Point::Basis(0) + Point::Basis(1);//(1,1,0)
      }
      else if(icolor == 2)
      {
        retval = Point::Zeros() + Point::Basis(1) + Point::Basis(2);//(0,1,1)
      }
      else if(icolor == 3)
      {
        retval = Point::Zeros() + Point::Basis(0) + Point::Basis(2);//(1,0,1)
      }
      else if(icolor == 4)
      {
        retval = Point::Zeros() + Point::Basis(1);//(0,1,0)
      }
      else if(icolor == 5)
      {
        retval = Point::Zeros() + Point::Basis(0);//(1,0,0)
      }
      else if(icolor == 6)
      {
        retval = Point::Zeros() + Point::Basis(2);//(0,0,1)
      }
      else if(icolor == 7)
      {
        retval = Point::Ones();//(1,1,1)
      }
      else
      {
        PROTO_ASSERT(false, "color not found");
      }
#else
      compiler_error_this_code_is_only_written_for_dim_2_or_3();
#endif
      return retval;
    }

    ///   returns true if there are domain boundary conditions
    inline static bool
    getStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,
               vector<LocalStencil<CELL, double> >       & a_stencil,
               Stencil<double>                           & a_regStencilInterior,
               vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
               vector<Box>                                 a_BCApplyBoxes[2*DIM],       
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
               bool                                        a_periodiconly)               ;

    inline static bool
    getSecondOrderPoissonStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                                 vector<LocalStencil<CELL, double> >       & a_stencil,                    
                                 Stencil<double>                           & a_regStencilInterior,         
                                 vector<Stencil<double> >                    a_regStencilBC[2*DIM],       
                                 vector<Box>                                 a_BCApplyBoxes[2*DIM],       
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

    inline static bool
    getRestrictionStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                          vector<LocalStencil<CELL, double> >       & a_stencil,                    
                          Stencil<double>                           & a_regStencilInterior,         
                          vector<Stencil<double> >                    a_regStencilBC[2*DIM],       
                          vector<Box>                                 a_BCApplyBoxes[2*DIM],       
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
                          bool                                        a_periodiconly,
                          bool                                        a_volumeWeighed)               ;

    ////
    inline static bool
    getProlongationStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                           vector<LocalStencil<CELL, double> >       & a_stencil,                    
                           Stencil<double>                           & a_regStencilInterior,         
                           vector<Stencil<double> >                    a_regStencilBC[2*DIM],       
                           vector<Box>                                 a_BCApplyBoxes[2*DIM],       
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
                           bool                                        a_periodiconly,
                           unsigned int                                a_icolor) ;              
  };


  ///
  template <int order>
  class Normalizor
  {
  public:
    typedef IndexTM<int  , DIM>                  IvDIM;
    typedef IndexTM<int  , DIM-1>                IvSDMinOne;
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
    template <CENTERING cent>
    inline static void
    normalizedFaceMoments(double                                                      & a_areaFrac,
                          RealVect                                                    & a_centroid,
                          const HostIrregData<cent, IndexedMoments<DIM-1, order>, 1 > & a_facedat,
                          const EBIndex<cent>                                         & a_face,
                          const EBGraph                                               & a_graph,
                          double                                                        a_dx   )
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

    

    inline static void
    normalizedVoluMoments(double                                                      & a_volFrac,
                          RealVect                                                    & a_centroid,
                          const HostIrregData<CELL, IndexedMoments<DIM  , order> , 1> & a_celldat,
                          const EBIndex<CELL>                                         & a_vof,
                          const EBGraph                                               & a_graph,
                          double                                                        a_dx   )
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
    inline static void
    normalizedBndryMoments(RealVect                                                       & a_centroid,
                           const EBFaData                                                 & a_ebfadat,
                           const EBIndex<BOUNDARY>                                        & a_vof,
                           const EBGraph                                                  & a_graph,
                           double                                                           a_dx   )
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
        //double fullarea = fullcellvol/a_dx;
        const IndexedMoments<DIM,order> & moments = a_ebfadat(a_vof, 0); 
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


    //sum of area fractions on a particular side of a vof
    inline static double
    sumArea(const VoluData                & a_voludata,
            const EBFaData                & a_ebfadata,
            const XFacData                & a_xfacdata,
            const YFacData                & a_yfacdata,
            const ZFacData                & a_zfacdata,
            const EBGraph                 & a_graph,
            const EBIndex<CELL>           & a_vof,
            const int                     & a_idir,
            const Side::LoHiSide          & a_sd,
            const double                  & a_dx)
    {
      
      double retval = 0;
      if(a_graph.isRegular(a_vof.m_pt))
      {
        retval = 1;
      }
      else if(a_graph.isCovered(a_vof.m_pt))
      {
        retval = 0;
      }
      else
      {
        if(a_idir == 0)
        {
          vector< EBIndex<XFACE> >faces = a_graph.getXFaces(a_vof, a_sd);
          for(int iface = 0; iface < faces.size(); iface++)
          {
            double   areaFrac;
            RealVect centroid;
            normalizedFaceMoments<XFACE>(areaFrac,
                                         centroid,
                                         a_xfacdata,
                                         faces[iface],
                                         a_graph,
                                         a_dx   );
            retval += areaFrac;
          }
        }
        else if(a_idir == 1)
        {
          vector< EBIndex<YFACE> >faces = a_graph.getYFaces(a_vof, a_sd);

          for(int iface = 0; iface < faces.size(); iface++)
          {
            double   areaFrac;
            RealVect centroid;
            normalizedFaceMoments<YFACE>(areaFrac,
                                         centroid,
                                         a_yfacdata,
                                         faces[iface],
                                         a_graph,
                                         a_dx   );
            retval += areaFrac;
          }
        }
#if DIM==3
        else if(a_idir == 2)
        {
          vector< EBIndex<ZFACE> >faces = a_graph.getZFaces(a_vof, a_sd);
          for(int iface = 0; iface < faces.size(); iface++)
          {
            double areaFrac;
            RealVect centroid;
            normalizedFaceMoments<ZFACE>(areaFrac,
                                         centroid,
                                         a_zfacdata,
                                         faces[iface],
                                         a_graph,
                                         a_dx   );
            retval += areaFrac;
          }
        }
#endif
        else
        {
          PROTO_ASSERT(false, "bogus idir");
        }
      } //else (is irregular)

      return retval;
    }

    template <typename data_t>
    inline static void 
    getNormalAndBoundaryArea(double                        & a_bndryArea,
                             RealVect                      & a_normal,
                             const VoluData                & a_voludata,
                             const EBFaData                & a_ebfadata,
                             const XFacData                & a_xfacdata,
                             const YFacData                & a_yfacdata,
                             const ZFacData                & a_zfacdata,
                             const EBIndex<CELL>           & a_vof,
                             const EBGraph                 & a_graph,
                             data_t                          a_dx   )
    {
      a_normal = BASISREALV(0); //give it something valid in case there is no area
      a_bndryArea = 0;
      if(a_graph.isIrregular(a_vof.m_pt))
      {

        data_t irregArea = 0.0;
        data_t hiArea[DIM];
        data_t loArea[DIM];
        for(int idir = 0; idir < DIM; idir++)
        {
          hiArea[idir] =  sumArea(a_voludata,
                                  a_ebfadata,
                                  a_xfacdata,
                                  a_yfacdata,
                                  a_zfacdata,
                                  a_graph,
                                  a_vof, idir,
                                  Side::Hi,
                                  a_dx);
          loArea[idir] =  sumArea(a_voludata,
                                  a_ebfadata,
                                  a_xfacdata,
                                  a_yfacdata,
                                  a_zfacdata,
                                  a_graph,
                                  a_vof, idir,
                                  Side::Lo,
                                  a_dx);
          irregArea += (hiArea[idir]-loArea[idir])*(hiArea[idir]-loArea[idir]);
        }
        a_bndryArea = sqrt(irregArea);
        if(a_bndryArea > 0)
        {
          for(int idir = 0; idir < DIM; idir++)
          {
            a_normal[idir] = (hiArea[idir]-loArea[idir])/a_bndryArea;
          }
        }
      }
    }
  };


  template<int order> //order is so we can send in geometry info
  class NeumannEBFlux
  {
  public:
    typedef IndexedMoments<DIM  , order>         IndMomDIM;
    typedef IndexedMoments<DIM-1, order>         IndMomSDMinOne;
    typedef HostIrregData<CELL    ,  IndMomDIM , 1>  VoluData;
    typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

    NeumannEBFlux()
    {
    }

    virtual ~NeumannEBFlux()
    {
    }

    inline static void
    getEBFluxStencil(LocalStencil<CELL, double>    &  a_facesten,
                     double                        &  a_weight, //for inhomogeneous bcs
                     const EBIndex<BOUNDARY>       &  a_face,
                     const EBGraph                 & a_graph,
                     const VoluData                & a_voludata,
                     const EBFaData                & a_ebfadata,
                     const XFacData                & a_xfacdata,
                     const YFacData                & a_yfacdata,
                     const ZFacData                & a_zfacdata,
                     const double                  & a_dx) 
    {

      //neumann has no stencil
    }
 
  };


  template<int order> //order is so we can send in geometry info
  class DirichletEBFlux
  {
  public:
    typedef IndexedMoments<DIM  , order>         IndMomDIM;
    typedef IndexedMoments<DIM-1, order>         IndMomSDMinOne;
    typedef HostIrregData<CELL    ,  IndMomDIM , 1>  VoluData;
    typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;

    DirichletEBFlux()
    {
    }

    virtual ~DirichletEBFlux()
    {
    }

    inline static void
    getEBFluxStencil(LocalStencil<CELL, double>    &  a_facesten,
                     double                        &  a_weight, //for inhomogeneous bcs
                     const EBIndex<BOUNDARY>       &  a_face,
                     const EBGraph                 & a_graph,
                     const VoluData                & a_voludata,
                     const EBFaData                & a_ebfadata,
                     const XFacData                & a_xfacdata,
                     const YFacData                & a_yfacdata,
                     const ZFacData                & a_zfacdata,
                     const double                  & a_dx); 

    
    inline static void
    calculateWeightingMatrix(RealVect           x0,
                             vector<RealVect>&  xp,
                             vector<RealVect>&  weightMatrix,
                             bool&              detZero);

  };
  ///stuff in common for all second order flux calculations
  template<CENTERING cent, int order>
  class Base2ndOrderFlux
  {
  public:
    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    
    virtual void
    getFaceCenteredFluxStencil(LocalStencil<CELL, double>                   &  a_facesten,
                               const EBIndex<cent>                          &  a_face,
                               const EBGraph                                &  a_graph,
                               const HostIrregData<cent, IndMomSDMinOne, 1> &  a_facdata,
                               double                                          a_dx) const = 0;

    void
    getFaceInterpStencil(LocalStencil<cent, double>                   &  a_sten,
                         const EBIndex<cent>                          &  a_face,
                         const EBGraph                                &  a_graph,
                         const HostIrregData<cent, IndMomSDMinOne, 1> &  a_facdata,
                         double                                          a_dx) const;

    void
    getCentroidStencil(LocalStencil<CELL, double>                   &  a_sten,
                       const EBIndex<cent>                          &  a_face,
                       const EBGraph                                &  a_graph,
                       const HostIrregData<cent, IndMomSDMinOne, 1> &  a_facdata,
                       double                                          a_dx) const;
    
  };
  
  ///second order Poisson flux
  template<CENTERING facecent, int order>
  class Poisson2ndOrder: public Base2ndOrderFlux<facecent, order>
  {

    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    typedef HostIrregData<CELL    ,  IndMomDIM , 1>  VoluData;
    typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
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
    typedef HostIrregData<CELL    ,  IndMomDIM , 1>  VoluData;
    typedef HostIrregData<BOUNDARY,  IndMomDIM , 1>  EBFaData;
    typedef HostIrregData<XFACE, IndMomSDMinOne, 1>  XFacData;
    typedef HostIrregData<YFACE, IndMomSDMinOne, 1>  YFacData;
    typedef HostIrregData<ZFACE, IndMomSDMinOne, 1>  ZFacData;
  public:

    ///
    inline static void
    get2ndOrderDivFStencil(LocalStencil< CELL,  double>             & a_vofStencil,
                           const EBIndex<CELL>                      & a_vof,
                           const EBGraph                            & a_graph,
                           const VoluData                           & a_voludata,
                           const EBFaData                           & a_ebfadata,
                           const XFacData                           & a_xfacdata,
                           const YFacData                           & a_yfacdata,
                           const ZFacData                           & a_zfacdata,
                           const Base2ndOrderFlux<XFACE   , order>  & a_xfaceflux,
                           const Base2ndOrderFlux<YFACE   , order>  & a_yfaceflux,
                           const Base2ndOrderFlux<ZFACE   , order>  & a_zfaceflux,
                           string                                  a_ebbcName,
                           double                                     a_dx);




    inline static void
    getDomainBCs(vector<Stencil<double> >                    a_regStencilBC[2*DIM],       
                 vector<Box>                                 a_BCApplyBoxes[2*DIM],       
                 string                                      a_domainBCName,               
                 const Box                                 & a_validBox,
                 const Box                                 & a_domain,
                 const Point                               & a_srcGhost);
    
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
  
