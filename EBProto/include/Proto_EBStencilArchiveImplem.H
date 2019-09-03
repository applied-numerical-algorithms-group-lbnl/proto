
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

    Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_facdata, a_face, a_graph, a_dx);

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
      a_sten.add(a_face, 1.0);
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
        Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_xfacdata, xfaces[iface], a_graph, a_dx);

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
        Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_yfacdata, yfaces[iface], a_graph, a_dx);
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
        Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_zfacdata, zfaces[iface], a_graph, a_dx);

        a_zfaceflux.getCentroidStencil(fluxStencil, zfaces[iface],
                                       a_graph,    
                                       a_zfacdata, 
                                       a_dx);       


        fluxStencil *= double(isign)*areaFrac/a_dx;
        a_vofStencil += fluxStencil;
      }
#endif      
    }

    LocalStencil<CELL, double> ebflux;
    EBIndex<BOUNDARY> face;
    face.m_pt = a_vof.m_pt;
    face.m_vofIDMe = a_vof.m_vofIDMe;
    face.m_vofIDLo = 0;
    face.m_isBoundary = false;
    
    if(a_ebbcName == string("Neumann"))
    {
      NeumannEBFlux<order>::
        getEBFluxStencil(ebflux, face,
                         a_graph,
                         a_voludata,
                         a_ebfadata,
                         a_xfacdata,
                         a_yfacdata,
                         a_zfacdata,
                         a_dx);
    }
    else if(a_ebbcName == string("Dirichlet"))
    {
      DirichletEBFlux<order>::
        getEBFluxStencil(ebflux, face,
                         a_graph,
                         a_voludata,
                         a_ebfadata,
                         a_xfacdata,
                         a_yfacdata,
                         a_zfacdata,
                         a_dx);
    }
    else
    {
      PROTO_ASSERT(false, "could not find EBBC string");
    }
      
    a_vofStencil += ebflux;
  }

  ///this deep down in the bowels, everything is a template specialization
  template <int order>
  void
  EBStencilArchive<CELL, CELL,  order, double>::
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
             bool                                        a_periodiconly)               
  {
    if(a_stencilName == string("Second_Order_Poisson"))
    {
      getSecondOrderPoissonStencil(a_dstVoFs,                 
                                   a_stencil,                 
                                   a_regStencilInterior,      
                                   a_regStencilBCS,
                                   a_BCApplyBoxes,
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
                               bool                                        a_periodiconly)               
  {

    a_regStencilInterior = Stencil<double>::Laplacian();
    a_regApplyBox = a_validBox;
    a_dstVoFs = a_graph.getIrregVoFs(a_validBox);
    Poisson2ndOrder<XFACE, order> xflux;
    Poisson2ndOrder<YFACE, order> yflux;
    Poisson2ndOrder<ZFACE, order> zflux;
    Poisson2ndOrder<BOUNDARY, order> ebflux;
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
                               xflux, yflux, zflux, a_ebbcName,
                               a_dx);
    }

    if(!a_periodiconly)
    {
      SecondOrderStencil<order>::
        getDomainBCs(a_regStencilBC, a_BCApplyBoxes, a_domainBCName, a_validBox, a_graph.getDomain(), a_srcGhost);
    }
  }
  /***/
  template < int order>
  void
  SecondOrderStencil<order>::
  getDomainBCs(vector<Stencil<double> >                    a_regStencilBC[2*DIM],       
               vector<Box>                                 a_BCApplyBoxes[2*DIM],       
               string                                      a_domainBCName,               
               const Box                                 & a_validBox,
               const Box                                 & a_domain,
               const Point                               & a_srcGhost)
    
  {
    Point lobox = a_validBox.low();
    Point hibox = a_validBox.high();
    Point lodom = a_domain.low();
    Point hidom = a_domain.high();
    Box srcBox = a_validBox.grow(a_srcGhost);
    double  stenSign;
    if(a_domainBCName == string("Dirichlet"))
    {
      stenSign = -1;
    }
    else if(a_domainBCName == string("Neumann"))
    {
      stenSign =  1;
    }
    else
    {
      PROTO_ASSERT(false, "unknown bc name");
    }

    for(int idir = 0; idir < DIM; idir++)
    {
        

      int numghost = a_srcGhost[idir]; //need a bc stencil for each ghost cell
      if(lobox[idir] == lodom[idir])
      {
        int loindex = ebp_index(idir, Side::Lo);
        Point loapp = srcBox.low();
        Point hiapp = srcBox.high();
        loapp[idir] = lodom[idir];
        hiapp[idir] = lodom[idir];

        Box appBox0(loapp, hiapp);

        a_regStencilBC[loindex].resize(numghost);
        a_BCApplyBoxes[loindex].resize(numghost);
        for(int igh = 0; igh < numghost; igh++)
        {
          a_BCApplyBoxes[loindex][igh] = appBox0.shift(idir, igh);
          a_regStencilBC[loindex][igh] = (stenSign)*Shift(Point::Basis(idir, -igh-1));
        }
      }
      if(hibox[idir] == hidom[idir])
      {
        int hiindex = ebp_index(idir, Side::Hi);
        Point loapp = srcBox.low();
        Point hiapp = srcBox.high();
        loapp[idir] = hidom[idir];
        hiapp[idir] = hidom[idir];

        Box appBox0(loapp, hiapp);
        
        a_regStencilBC[hiindex].resize(numghost);
        a_BCApplyBoxes[hiindex].resize(numghost);
        for(int igh = 0; igh < numghost; igh++)
        {
          a_BCApplyBoxes[hiindex][igh] = appBox0.shift(idir, -igh);
          a_regStencilBC[hiindex][igh] = (stenSign)*Shift(Point::Basis(idir, igh+1));
        }
      }
    }
  }
}
#endif
