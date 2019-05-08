
#ifndef _Proto_EBStencilArchiveImplem_H_
#define _Proto_EBStencilArchiveImplem_H_

namespace Proto
{
  /*****************************/  
  template<CENTERING cent, int order>
  void
  Base2ndOrderFlux<cent, order>::
  getFaceInterpStencil(LocalStencil<cent, double> &  a_sten,
                       const EBIndex<cent>        &  a_face,
                       const EBGraph              &  a_graph,
                       const VoluData             &  a_voludata,
                       const EBFaData             &  a_ebfadata,
                       const XFacData             &  a_xfacdata,
                       const YFacData             &  a_yfacdata,
                       const ZFacData             &  a_zfacdata,
                       const double               &  a_dx)
  {
    a_sten.clear();
    double areaFrac;
    RealVect centroid;
    int  facedir = a_face.direction();
    int tandirs[DIM-1];
    int itan = 0;
    for(int idir = 0; idir < DIM; idir++)
    {
      if(idir != facedir)
      {
        tandirs[itan] = idir;
        itan++; 
      }
    }

    if(facedir == 0)
    {
      Normalizor<cent, order>::normalizedFaceMoments(a_areaFrac, centroid, a_xfacdata, a_graph, a_dx);
    }
    else if(facedir == 1)
    {
      Normalizor<cent, order>::normalizedFaceMoments(a_areaFrac, centroid, a_yfacdata, a_graph, a_dx);
    }
#if DIM==3
    else
    {
      Normalizor<cent, order>::normalizedFaceMoments(a_areaFrac, centroid, a_zfacdata, a_graph, a_dx);
    }
#endif
    else
    {
      PR_error("bogus direction");
    }
    Box domain = a_graph.getDomain();
    bool dropOrder = false;
    Side::LoHiSide sides[DIM-1];
    for(int itan = 0; itan < DIM-1; itan++)
    {
      if(centroid[tandirs[itan]] > 0)
      {
        sides[itan] = Side::Hi;
      }
      else
      {
        sides[itan] = Side::Lo;
      }
    }

#if DIM==2    
    EBIndex<cent> faces[2];
    faces[0] = a_face;
    bool hasFace = a_graph.getAdjacentFace<cent>(faces[1], a_face,tandirs[0], sides[0]);      
    if(!hasFace) dropOrder = true;
    if(dropOrder)
    {
      a_sten.add(a_face, 1.0);
    }
    else
    {
      double xbar = centroid[tandirs[0]];
      a_sten.add(faces[0], 1.0-xbar);
      a_sten.add(faces[1],     xbar);
    }
#else
    EBIndex<cent> faces[4];
    faces[0] = a_face;
    double xbar = centroid[tandirs[0]];
    double ybar = centroid[tandirs[1]];
    bool hasFacex = a_graph.getAdjacentFace<cent>(faces[1], a_face, tandirs[0], sides[0]);      
    bool hasFacey = a_graph.getAdjacentFace<cent>(faces[2], a_face, tandirs[1], sides[1]);      
    if((!hasFacex)  || (!hasFacey))
    {
      dropOrder = true;
    }
    else
    {
      EBIndex<cent> facexy, faceyx;
      bool hasFaceyx = a_graph.getAdjacentFace<cent>(faceyx, faces[1], tandirs[1], sides[1]);      
      bool hasFacexy = a_graph.getAdjacentFace<cent>(facexy, faces[2], tandirs[0], sides[0]);      
      if((!hasFaceyx) || (!hasFacexy) || (facexy != faceyx))
      {
        dropOrder = true;
        a_sten.add(a_face, 1.0);
      }
      else
      {
        faces[3] = facexy;
        double  f00coef = 1.0 - xbar - ybar + xbar*ybar;
        double  f10coef = xbar - xbar*ybar;
        double  f01coef = ybar - xbar*ybar;
        double  f11coef = xbar*ybar;
        a_sten.add(faces[0], f00coef);
        a_sten.add(faces[1], f10coef);
        a_sten.add(faces[2], f01coef);
        a_sten.add(faces[3], f11coef);
      }
    }
#endif
  }
/*******************************/
  template <CENTERING cent, typename data_t, int order>
  inline void
  Normalizor<cent, order>::
  normalizedFaceMoments(double                                      & a_areaFrac,
                        RealVect                                  & a_centroid,
                        const IrregData<facecent, IndMomSDMinOne, 1>& a_facedat,
                        const EBIndex<cent>                         & a_face,
                        const EBGraph                               & a_graph,
                        const double                                & a_dx   )
  {


    EBIndex<CELL> voflo = a_face.getVoF(Side::Lo);
    EBIndex<CELL> vofhi = a_face.getVoF(Side::Hi);
    if(a_graph.isRegular(voflo.m_pt) || a_graph.isRegular(vofhi.m_pt))
    {
      a_areaFrac = 1;
      a_centroid  = RealVect::Zero;
    }
    else if(a_graph.isCovered(voflo.m_pt) || a_graph.isCovered(vofhi.m_pt))
    {
      a_areaFrac = 0;
      a_centroid  = RealVect::Zero;
    }
    else
    {
      double fullfacearea = 1.;
      for(int idir = 0; idir < DIM-1; idir++) fullfacearea *= a_dx;
      IndMomSDMinOne facemom = a_facedat(a_face, 0); 
      double area = facemom[IvSDMinOne::Zero];
      a_areaFrac = area/fullfacearea;
      a_centroid = RealVect::Zero;
      if(area > 0)
      {
        int iindex = 0;
        for (int idir = 0; idir < SpaceDim; ++idir)
        {
          if(idir != a_face.direction())
          {
            IvSDMinOne mono = BASISV_TM<int,DIM-1>(iindex);
            a_faceCentroid[idir] = facemom[mono];
            //normalize  by area
            a_faceCentroid[idir] /= area;
            //normalize by dx
            a_faceCentroid[idir] /= a_dx;
            iindex++;
          }
        }
      }
    }
  }


/*******************************/
  template < typename int order>
  inline void
  Normalizor<CELL, order>::
  normalizedVoluMoments(double                                      & a_volFrac,
                        RealVect                                  & a_centroid,
                        const IrregData<CELL, IndMomSpaceDim, 1>    & a_celldat,
                        const EBIndex<CELL>                         & a_vof,
                        const EBGraph                               & a_graph,
                        const double                                & a_dx   )
  {
    if(a_graph.isRegular(a_vof.m_pt))
    {
      a_volFrac = 1;
      a_centroid  = RealVect::Zero;
    }
    else if(a_graph.isCovered(a_vof.m_pt))
    {
      a_volFrac = 0;
      a_centroid  = RealVect::Zero;
    }
    else
    {
      double fullcellvol = 1.;
      for(int idir = 0; idir < DIM; idir++) fullcellvol *= a_dx;
      IndMomSpaceDim moments = a_celldat(a_vof, 0); 
      double volume = moments[IvSpaceDim::Zero];
      a_volFrac = volume/fullcellvol;
      a_centroid = RealVect::Zero;
      if(volume > 0)
      {
        for (int idir = 0; idir < SpaceDim; ++idir)
        {
          IvSpaceDim mono = BASISV_TM<int,DIM>(idir);
          a_centroid[idir] = moments[mono];
          //normalize  by volume
          a_centroid[idir] /= volume;
          //normalize by dx
          a_centroid[idir] /= a_dx;
        }
      }
    }
  }

/*******************************/
  template < int order>
  inline void
  Normalizor<BOUNDARY,  order>::
  normalizedBndryMoments(RealVect                                        & a_centroid,
                         const IrregData<BOUNDARY, IndMomSpaceDim, 1>    & a_celldat,
                         const EBIndex<BOUNDARY>                         & a_vof,
                         const EBGraph                                   & a_graph,
                         const double                                    & a_dx   )
  {
    if(a_graph.isRegular(a_vof.m_pt))
    {
      a_areaFrac = 1;
      a_centroid  = RealVect::Zero;
    }
    else if(a_graph.isCovered(a_voflo.m_pt))
    {
      a_areaFrac = 0;
      a_centroid  = RealVect::Zero;
    }
    else
    {
      double fullcellvol = 1.;
      for(int idir = 0; idir < DIM; idir++) fullcellvol *= a_dx;
      double fullarea = fullcellvol/a_dx;
      IndMomSpaceDim moments = a_celldat(a_vof, 0); 
      double area = moments[IvSpaceDim::Zero];

      a_centroid = RealVect::Zero;
      if(area > 0)
      {
        for (int idir = 0; idir < SpaceDim; ++idir)
        {
          IvSpaceDim mono = BASISV_TM<int,DIM>(idir);
          a_centroid[idir] = moments[mono];
          //normalize  by area
          a_centroid[idir] /= area;
          //normalize by dx
          a_centroid[idir] /= a_dx;
        }
      }
    }
  }

/*******************************/
  template <  int order>
  inline void
  Normalizor<CELL, order>::
  getNormalAndBoundaryArea(double                        & a_bndryArea,
                           RealVect                      & a_normal,
                           const VoluData                & a_voludata,
                           const XFacData                & a_xfacdata,
                           const YFacData                & a_yfacdata,
                           const ZFacData                & a_zfacdata,
                           const EBIndex<CELL>           & a_vof,
                           const EBGraph                 & a_graph,
                           const double                  & a_dx   )
  {
    if(a_graph.isRegular(a_vof.m_pt))
    {
      a_bndryArea = 0;
      a_normal  = RealVect::Zero;
      a_normal[0] = 1;;
    }
    else if(a_graph.isCovered(a_voflo.m_pt))
    {
      a_bndryArea = 0;
      a_centroid  = RealVect::Zero;
      a_normal[0] = 1;
    }
    else
    {
      ///this is kind of a pain in the neck to do.
      //so I am leaving it for now
      PR_error("this bit as yet unfinished");
    }
  }
           
///
  template <CENTERING facecent, int order>
  Poisson2ndOrder<facecent>::
  void
  getFaceCenteredFluxStencil(LocalStencil<CELL, double>    &  a_facesten,
                             const EBIndex<facecent>       &  a_face,
                             const EBGraph                 &  a_graph,
                             const VoluData                &  a_voludata,
                             const EBFaData                &  a_ebfadata,
                             const XFacData                &  a_xfacdata,
                             const YFacData                &  a_yfacdata,
                             const ZFacData                &  a_zfacdata,
                             const double                  &  a_dx)
  {
    EBIndex<CELL> voflo = a_face.getVoF(Side::Lo);
    EBIndex<CELL> vofhi = a_face.getVoF(Side::Hi);
    a_facesten.clear();
    a_facesten.add(voflo, -1./a_dx);
    a_facesten.add(voflo, -1./a_dx);
  }

  template < int order>
  void
  SecondOrderStencil<order>::
  get2dOrderDivFStencil(LocalStencil<CELL, double>            & a_vofStencil,
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
                        double                                  a_dx)         
  {
    CH_TIME("NWOEBConductivityOp::getDivFStencil");
    a_vofStencil.clear();

    for (SideIterator sit; sit.ok(); ++sit)
    {
      int isign = sign(sit());
      vector< EBIndex<XFACE> > xfaces = a_graph.getXFaces(a_vof, sit());
      for (int iface = 0; iface < xfaces.size(); iface++)
      {
        LocalStencil<CELL, double> fluxStencil;
        a_xfaceflux.getCentroidStencil(fluxStencil, xfaces[iface],
                                       a_graph,    
                                       a_voludata, 
                                       a_ebfadata, 
                                       a_xfacdata, 
                                       a_dx);       

        Real areaFrac = ebisBox.areaFrac(faces[iface]);
        fluxStencil *= double(isign)*areaFrac/m_dx;
        a_vofStencil += fluxStencil;
      }
      vector< EBIndex<YFACE> > yfaces = a_graph.getYFaces(a_vof, sit());
      for (int iface = 0; iface < yfaces.size(); iface++)
      {
        LocalStencil<CELL, double> fluxStencil;
        a_yfaceflux.getCentroidStencil(fluxStencil, yfaces[iface],
                                       a_graph,    
                                       a_voludata, 
                                       a_ebfadata, 
                                       a_yfacdata, 
                                       a_dx);       

        Real areaFrac = ebisBox.areaFrac(faces[iface]);
        fluxStencil *= double(isign)*areaFrac/m_dx;
        a_vofStencil += fluxStencil;
      }
#if DIM==3
      vector< EBIndex<ZFACE> > zfaces = a_graph.getZFaces(a_vof, sit());
      for (int iface = 0; iface < zfaces.size(); iface++)
      {
        LocalStencil<CELL, double> fluxStencil;
        a_yfaceflux.getCentroidStencil(fluxStencil, zfaces[iface],
                                       a_graph,    
                                       a_voludata, 
                                       a_ebfadata, 
                                       a_zfacdata, 
                                       a_dx);       

        Real areaFrac = ebisBox.areaFrac(faces[iface]);
        fluxStencil *= double(isign)*areaFrac/m_dx;
        a_vofStencil += fluxStencil;
      }
#endif      
    }
    if(a_ebbName == string("Neumann"))
    {
      //nothing to add here because Neumann has no stencil
    }
    else
    {
      PR_error("unknown EBBC");
    }
  }

  ///this deep down in the bowels, everything is a template specialization
  template <int order, typename data_t>
  void
  EBStencilArchive<CELL, CELL>::
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
             string              a_ebbcName,
             data_t              a_dx,
             bool                a_periodiconly)
  {
    if(a_stencilName == string("Second_Order_Poisson"))
    {
      Poisson2ndOrder<XFACE, data_t> xflux;
      Poisson2ndOrder<YFACE, data_t> yflux;
      Poisson2ndOrder<ZFACE, data_t> zflux;

      get2ndOrderDivFStencil(a_dstVoFs, a_stencil, a_regStencilInterior, a_applyBox,
                             a_ebbcName,
                             a_validBox,
                             a_domain,
                             a_srcGhost,
                             a_dstGhost,
                             a_graph,
                             a_voludata,
                             a_ebfadata,
                             a_xfacdata, a_yfacdata, a_zfacdata,
                             xflux, yflux, zflux,
                             a_ebbcName,
                             a_dx);

    }
    else
    {
      PR_error("stencil_archive_cell_cell::stencil string not found");
    }
  }


}

#endif
