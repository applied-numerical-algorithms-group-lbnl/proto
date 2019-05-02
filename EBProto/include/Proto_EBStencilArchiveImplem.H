
#ifndef _Proto_EBStencilArchiveImplem_H_
#define _Proto_EBStencilArchiveImplem_H_

namespace Proto
{
  typedef IndexTM<int   , SpaceDim>    IvSpaceDim;
  typedef IndexTM<int   , SpaceDim-1>  IvSDMinOne;
  typedef IndexTM<double, SpaceDim>    RvSpaceDim;
  typedef IndexTM<double, SpaceDim-1>  RvSDMinOne;

/*******************************/
  template <CENTERING cent>
  inline void
  normalizedFaceMoments(double                                      & a_areaFrac,
                        RvSpaceDim                                  & a_centroid,
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
      a_centroid  = RvSpaceDim::Zero;
    }
    else if(a_graph.isCovered(voflo.m_pt) || a_graph.isCovered(vofhi.m_pt))
    {
      a_areaFrac = 0;
      a_centroid  = RvSpaceDim::Zero;
    }
    else
    {
      double fullfacearea = 1.;
      for(int idir = 0; idir < DIM-1; idir++) fullfacearea *= a_dx;
      IndMomSDMinOne facemom = a_facedat(a_face, 0); 
      double area = Frac = facemom[IvSDMinOne::Zero];
      a_areaFrac = facemom[IvSDMinOne]/fullfacearea;
      a_centroid = RvSpaceDim::Zero;
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
  inline void
  normalizedVoluMoments(double                                      & a_volFrac,
                        RvSpaceDim                                  & a_centroid,
                        const IrregData<CELL, IndMomSpaceDim, 1>    & a_celldat,
                        const EBIndex<cent>                         & a_vof,
                        const EBGraph                               & a_graph,
                        const double                                & a_dx   )
  {
    if(a_graph.isRegular(a_vof.m_pt))
    {
      a_volFrac = 1;
      a_centroid  = RvSpaceDim::Zero;
    }
    else if(a_graph.isCovered(voflo.m_pt) || a_graph.isCovered(vofhi.m_pt))
    {
      a_volFrac = 0;
      a_centroid  = RvSpaceDim::Zero;
    }
    else
    {
      double fullcellvol = 1.;
      for(int idir = 0; idir < DIM; idir++) fullcellvol *= a_dx;
      IndMomSpaceDim moments = a_celldat(a_vof, 0); 
      double area = Frac = facemom[IvSDMinOne::Zero];
      a_volFrac = facemom[IvSDMinOne]/fullcellvol;
      a_centroid = RvSpaceDim::Zero;
      if(fullcellvol > 0)
      {
        for (int idir = 0; idir < SpaceDim; ++idir)
        {
          IvSpaceDim mono = BASISV_TM<int,DIM>(idir);
          a_centroid[idir] = moments[mono];
          //normalize  by volume
          a_centroid[idir] /= fullcellvol;
          //normalize by dx
          a_centroid[idir] /= a_dx;
        }
      }
    }
  }

  HERE NEED TO PUT IN A BIT THAT GETS NORMALS, BOUNDARYAREAFRAC AND BOUNDARY CENTROIDS
           
///
  template <CENTERING facecent, typename data_t>
  Poisson2ndOrder<facecent, data_t>::
  void
  getFaceCenteredFluxStencil(LocalStencil<CELL, data_t>    &  a_facesten,
                             const EBIndex<facecent>       &  a_face,
                             const EBGraph                 &  a_graph,
                             const VoluData                &  a_voludata,
                             const EBFaData                &  a_ebfadata,
                             const XFacData                &  a_xfacdata,
                             const YFacData                &  a_yfacdata,
                             const ZFacData                &  a_zfacdata,
                             const data_t                  &  a_dx)
  {
    EBIndex<CELL> voflo = a_face.getVoF(Side::Lo);
    EBIndex<CELL> vofhi = a_face.getVoF(Side::Hi);
    a_facesten.clear();
    a_facesten.add(voflo, -1./a_dx);
    a_facesten.add(voflo, -1./a_dx);
  }

 template <typename data_t>
 void
 SecondOrderStencil<data_t>::
 get2dOrderDivFStencil(LocalStencil<CELL, data_t>            & a_vofStencil,
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
                       data_t                                  a_dx)         
  {
    CH_TIME("NWOEBConductivityOp::getDivFStencil");
    a_vofStencil.clear();

      for (SideIterator sit; sit.ok(); ++sit)
      {
        int isign = sign(sit());
        vector< EBIndex<XFACE> > xfaces = a_graph.getXFaces(a_vof, sit());
        for (int iface = 0; iface < xfaces.size(); iface++)
        {
          EBIndex<CELL> fluxStencil;
          a_xfaceflux.getCentroidStencil(fluxStencil, xfaces[iface],
                                         a_graph,    
                                         a_voludata, 
                                         a_ebfadata, 
                                         a_xfacdata, 
                                         a_dx);       

          Real areaFrac = ebisBox.areaFrac(faces[iface]);
          fluxStencil *= data_t(isign)*areaFrac/m_dx;
          a_vofStencil += fluxStencil;
        }
        vector< EBIndex<YFACE> > yfaces = a_graph.getYFaces(a_vof, sit());
        for (int iface = 0; iface < yfaces.size(); iface++)
        {
          EBIndex<CELL> fluxStencil;
          a_yfaceflux.getCentroidStencil(fluxStencil, yfaces[iface],
                                         a_graph,    
                                         a_voludata, 
                                         a_ebfadata, 
                                         a_yfacdata, 
                                         a_dx);       

          Real areaFrac = ebisBox.areaFrac(faces[iface]);
          fluxStencil *= data_t(isign)*areaFrac/m_dx;
          a_vofStencil += fluxStencil;
        }
        vector< EBIndex<ZFACE> > zfaces = a_graph.getZFaces(a_vof, sit());
      }
    }
  }

  ///this deep down in the bowels, everything is a template specialization
  template < >
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
      SecondOrderStencil<data_t>::
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
