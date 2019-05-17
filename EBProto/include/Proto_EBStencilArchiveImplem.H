
#ifndef _Proto_EBStencilArchiveImplem_H_
#define _Proto_EBStencilArchiveImplem_H_

namespace Proto
{
  /*****************************/  
  template<CENTERING cent, int order>
  void
  Base2ndOrderFlux<cent, order>::
  getFaceInterpStencil(LocalStencil<cent, double>                   &  a_sten,
                       const EBIndex<cent>                          &  a_face,
                       const EBGraph                                &  a_graph,
                       const HostIrregData<cent, IndMomSDMinOne, 1> &  a_facdata,
                       double                                          a_dx) const
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

    Normalizor::normalizedFaceMoments<cent, double, order>(areaFrac, centroid, a_facdata, a_face, a_graph, a_dx);

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
  template<CENTERING cent, int order>
  void
  Base2ndOrderFlux<cent, order>::
  getCentroidStencil(LocalStencil<CELL, double>                   &  a_sten,
                     const EBIndex<cent>                          &  a_face,
                     const EBGraph                                &  a_graph,
                     const HostIrregData<cent, IndMomSDMinOne, 1> &  a_facdata,
                     double                                          a_dx) const
  {
    LocalStencil<cent, double> facesten;
    getFaceInterpStencil(facesten, a_face,
                         a_graph,    
                         a_facdata, 
                         a_dx);   
    a_sten.clear();
    for(int isten = 0; isten < facesten.size(); isten++)
    {
      EBIndex<cent>& face   = facesten.m_entries[isten].m_vof;
      double       & weight = facesten.m_entries[isten].m_weight;
      LocalStencil<CELL, double> centSten;
      getFaceCenteredFluxStencil(centSten, face, a_graph, a_facdata, a_dx);
      centSten *= weight;
      a_sten   += centSten;
    }
  }

           
///
  template <CENTERING cent, int order>
  void
  Poisson2ndOrder<cent,order>::
  getFaceCenteredFluxStencil(LocalStencil<CELL, double>                   &  a_facesten,
                             const EBIndex<cent>                          &  a_face,
                             const EBGraph                                &  a_graph,
                             const HostIrregData<cent, IndMomSDMinOne, 1> &  a_facdata,
                             double                                          a_dx) const
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
                         double                                  a_dx)
  {
    PR_TIME("get2ndOrderDivFStencil");
    a_vofStencil.clear();

    for (SideIterator sit; sit.ok(); ++sit)
    {
      int isign = sign(sit());
      vector< EBIndex<XFACE> > xfaces = a_graph.getXFaces(a_vof, sit());
      for (int iface = 0; iface < xfaces.size(); iface++)
      {
        LocalStencil<CELL, double> fluxStencil;
        double   areaFrac;
        RealVect centroid;
        Normalizor::normalizedFaceMoments<XFACE, double, order>(areaFrac, centroid, a_xfacdata, xfaces[iface], a_graph, a_dx);

        a_xfaceflux.getCentroidStencil(fluxStencil, xfaces[iface],
                                       a_graph,    
                                       a_xfacdata, 
                                       a_dx);       

        fluxStencil *= double(isign)*areaFrac/a_dx;
        a_vofStencil += fluxStencil;
      }
      vector< EBIndex<YFACE> > yfaces = a_graph.getYFaces(a_vof, sit());
      for (int iface = 0; iface < yfaces.size(); iface++)
      {
        LocalStencil<CELL, double> fluxStencil;
        double   areaFrac;
        RealVect centroid;
        Normalizor::normalizedFaceMoments<YFACE, double, order>(areaFrac, centroid, a_yfacdata, yfaces[iface], a_graph, a_dx);
        a_yfaceflux.getCentroidStencil(fluxStencil, yfaces[iface],
                                       a_graph,    
                                       a_yfacdata, 
                                       a_dx);       

        fluxStencil *= double(isign)*areaFrac/a_dx;
        a_vofStencil += fluxStencil;
      }
#if DIM==3
      vector< EBIndex<ZFACE> > zfaces = a_graph.getZFaces(a_vof, sit());
      for (int iface = 0; iface < zfaces.size(); iface++)
      {
        LocalStencil<CELL, double> fluxStencil;
        double   areaFrac;
        RealVect centroid;
        Normalizor::normalizedFaceMoments<ZFACE, double, order>(areaFrac, centroid, a_zfacdata, zfaces[iface], a_graph, a_dx);

        a_yfaceflux.getCentroidStencil(fluxStencil, zfaces[iface],
                                       a_graph,    
                                       a_zfacdata, 
                                       a_dx);       


        fluxStencil *= double(isign)*areaFrac/a_dx;
        a_vofStencil += fluxStencil;
      }
#endif      
    }
    if(a_ebbcName == string("Neumann"))
    {
      //nothing to add here because Neumann has no stencil
    }
    else
    {
      PR_error("unknown EBBC");
    }
  }

  ///this deep down in the bowels, everything is a template specialization
  template <int order>
  void
  EBStencilArchive<CELL, CELL,  order, double>::
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
             bool                                        a_periodiconly)               
  {
    if(a_stencilName == string("Second_Order_Poisson"))
    {
      getSecondOrderPoissonStencil(a_dstVoFs,                 
                                   a_stencil,                 
                                   a_regStencilInterior,      
                                   a_regStencilBCS,    
                                   a_regApplyBox,             
                                   a_stencilName,             
                                   a_domainBCName,            
                                   a_ebbcName,                
                                   a_validBox,                
                                   a_domain,                  
                                   a_srcGhost,                
                                   a_dstGhost,                
                                   a_graph,                   
                                   a_voludata,                
                                   a_ebfadata,                
                                   a_xfacdata,                
                                   a_yfacdata,                
                                   a_zfacdata,                
                                   a_dx,                      
                                   a_periodiconly);            


    }
    else
    {
      PR_error("stencil_archive_cell_cell::stencil string not found");
    }
  }


  ///this deep down in the bowels, everything is a template specialization
  template <int order>
  void
  EBStencilArchive<CELL, CELL,  order, double>::
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
                               bool                                        a_periodiconly)               
  {

    //this bunch of stuff should probably be put in a separate function
    if(a_periodiconly)
    {
      a_regStencilInterior = Stencil<double>::Laplacian();
      a_regApplyBox = a_validBox;
      a_dstVoFs = a_graph.getIrregVoFs(a_validBox);
      Poisson2ndOrder<XFACE, order> xflux;
      Poisson2ndOrder<YFACE, order> yflux;
      Poisson2ndOrder<ZFACE, order> zflux;
      a_stencil.resize(a_dstVoFs.size());
      for(int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
      {
        SecondOrderStencil<order>::
          get2ndOrderDivFStencil(a_stencil[ivof],
                                 a_dstVoFs[ivof],
                                 a_graph,
                                 a_voludata,
                                 a_ebfadata,
                                 a_xfacdata, a_yfacdata, a_zfacdata,
                                 xflux, yflux, zflux,
                                 a_ebbcName,
                                 a_dx);
      }

    }
    else
    {
      PR_error("boundary conditions not programmed yet");
    }
  }

}

#endif
