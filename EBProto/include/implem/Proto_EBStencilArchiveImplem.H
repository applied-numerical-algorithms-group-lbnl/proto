
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
//debug
//    dropOrder = true;
//end debug
    if(dropOrder)
    {
      a_sten.add(a_face, 1.0);
    }
    else
    {
      double xbar = std::abs(centroid[tandirs[0]]);
      a_sten.add(faces[0], 1.0-xbar);
      a_sten.add(faces[1],     xbar);
    }
#else
    EBIndex<cent> faces[4];
    faces[0] = a_face;
    double xbar = std::abs(centroid[tandirs[0]]);
    double ybar = std::abs(centroid[tandirs[1]]);
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
#if 1
    LocalStencil<cent, double> facesten;
    getFaceInterpStencil(facesten, a_face,
                         a_graph,    
                         a_facdata, 
                         a_dx);   
#else
    LocalStencil<cent, double> facesten;
    facesten.add(a_face, 1.0);
#endif 
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
    a_facesten.add(vofhi,  1./a_dx);
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
      vector< EBIndex<XFACE> > xfaces = a_graph.getXFaces(a_vof, sit());
      int isign = sign(sit());
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
        double factor = double(isign)*areaFrac/a_dx;
        fluxStencil *= factor;
        a_vofStencil += fluxStencil;
      }
    }

    for (SideIterator sit; sit.ok(); ++sit)
    {
      vector< EBIndex<YFACE> > yfaces = a_graph.getYFaces(a_vof, sit());
      int isign = sign(sit());
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
    }
#if DIM==3
    for (SideIterator sit; sit.ok(); ++sit)
    {
      vector< EBIndex<ZFACE> > zfaces = a_graph.getZFaces(a_vof, sit());
      int isign = sign(sit());
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
    }
#endif      

    LocalStencil<CELL, double> ebflux;
    EBIndex<BOUNDARY> face = a_vof.getCutFace();
    double weight; //for inhomogeneous bcs
    if(a_ebbcName == string("Neumann"))
    {
      NeumannEBFlux<order>::
        getEBFluxStencil(ebflux, weight,
                         face,
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
        getEBFluxStencil(ebflux, weight,
                         face,
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
  ///returns whether there are domain boundary conditions
  template <int order>
  bool 
  EBStencilArchive<CELL, CELL,  order, double>::
  getStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,
             vector<LocalStencil<CELL, double> >       & a_stencil,
             Stencil<double>                           & a_regStencil,
             vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
             vector<Box>                                 a_BCApplyBoxes[2*DIM],       
             Box                                       & a_regApplyBox,
             string                                      a_stencilName,
             string                                      a_domainBCName,
             string                                      a_ebbcName,
             const Box                                 & a_srcValid,
             const Box                                 & a_dstValid,
             const Box                                 & a_srcDomain,
             const Box                                 & a_dstDomain,
             const Point                               & a_srcGhost,
             const Point                               & a_dstGhost,
             const EBGraph                             & a_srcGraph,
             const EBGraph                             & a_dstGraph,
             const VoluData                            & a_voludata,
             const EBFaData                            & a_ebfadata,
             const XFacData                            & a_xfacdata,
             const YFacData                            & a_yfacdata,
             const ZFacData                            & a_zfacdata,
             double                                      a_dx,
             bool                                        a_periodiconly, 
             Point                                       a_dstGrown)               
  {
    if(a_stencilName == string("Second_Order_Poisson"))
    {
      //for this one srcdomain == dstdomain
      //             srcgraph  == dstgraph
      return
        getSecondOrderPoissonStencil(a_dstVoFs,                 
                                     a_stencil,                 
                                     a_regStencil,      
                                     a_regStencilBCS,
                                     a_BCApplyBoxes,
                                     a_regApplyBox,             
                                     a_stencilName,             
                                     a_domainBCName,            
                                     a_ebbcName,                
                                     a_srcValid,                
                                     a_srcDomain,                  
                                     a_srcGhost,                
                                     a_dstGhost,                
                                     a_srcGraph,                   
                                     a_voludata,                
                                     a_ebfadata,                
                                     a_xfacdata,                
                                     a_yfacdata,                
                                     a_zfacdata,                
                                     a_dx,                      
                                     a_periodiconly);            

    }
    else if( (a_stencilName.find("Volume_Weighted_Averaging_rad_")) != std::string::npos)
    {

      //figure out the radius
      auto pos = a_stencilName.find_first_of("0123456789");
      unsigned long radius = 999;
      if (std::string::npos != pos)
      {
        radius = std::stoul(a_stencilName.substr(pos));
      }
      //for this one srcdomain == dstdomain
      //             srcgraph  == dstgraph
      return
        getNormalizationStencil(a_dstVoFs,                 
                                a_stencil,                 
                                a_regStencil,      
                                a_regApplyBox,             
                                a_srcValid,                
                                a_srcDomain,                  
                                a_srcGhost,                
                                a_dstGhost,                
                                a_srcGraph,               
                                a_voludata,
                                radius,
                                a_dx);                      
    }
    else if( (a_stencilName.find("Restriction")) != std::string::npos)
    {
      bool kappaWeighted = true;
      if( (a_stencilName.find("Multigrid")) != std::string::npos)
      {
        kappaWeighted = false;
      }

      return
        getRestrictionStencil(a_dstVoFs,                 
                              a_stencil,                 
                              a_regStencil,      
                              a_regApplyBox,             
                              a_srcValid,                
                              a_dstValid,                
                              a_srcGhost,                
                              a_dstGhost,                
                              a_srcGraph,                   
                              a_dstGraph,                   
                              a_voludata,
                              a_dx,
                              kappaWeighted);            

    }
    else if( (a_stencilName.find("PWC_Prolongation_")) != std::string::npos)
    {

      //figure out the color
      auto pos = a_stencilName.find_first_of("0123456789");
      unsigned long color = 999;
      if (std::string::npos != pos)
      {
       color = std::stoul(a_stencilName.substr(pos));
      }

      return
        getProlongationStencil(a_dstVoFs,                 
                               a_stencil,                 
                               a_regStencil,      
                               a_regApplyBox,             
                               a_srcValid,
                               a_dstValid,
                               a_srcGhost,                
                               a_dstGhost,                
                               a_srcGraph,                   
                               a_dstGraph,
                               color);
    }
    else if( (a_stencilName.find("Volume_Weighted_Redistribution_rad_")) != std::string::npos)
    {
      //figure out the radius
      auto pos = a_stencilName.find_first_of("0123456789");
      unsigned long radius = 999;
      if (std::string::npos != pos)
      {
        radius = std::stoul(a_stencilName.substr(pos));
      }
      //for this one srcdomain == dstdomain
      //             srcgraph  == dstgraph
      return
        getRedistributionStencil(a_dstVoFs,                 
                                 a_stencil,                 
                                 a_regStencil,      
                                 a_regApplyBox,             
                                 a_srcValid,                
                                 a_srcDomain,                  
                                 a_srcGhost,                
                                 a_dstGhost,                
                                 a_srcGraph,               
                                 a_voludata,
                                 radius,
                                 a_dx);                      

    }
    else if( (a_stencilName.find("Slope_Low_")) != std::string::npos)
    {
      //figure out the radius
      auto pos = a_stencilName.find_first_of("0123456789");
      unsigned long dir = 999;
      if (std::string::npos != pos)
      {
        dir = std::stoul(a_stencilName.substr(pos));
      }
      unsigned   int idir = dir;

      //for this one srcdomain == dstdomain
      //             srcgraph  == dstgraph
      return
        getSlopeStencil(a_dstVoFs,                 
                        a_stencil,                 
                        a_regStencil,      
                        a_regStencilBCS,
                        a_BCApplyBoxes,
                        a_regApplyBox,  
                        a_domainBCName,
                        a_srcValid,                
                        a_srcDomain,                  
                        a_srcGhost,                
                        a_dstGhost,                
                        a_srcGraph,               
                        a_voludata,
                        a_dstGrown,
                        idir,
                        Side::Lo);       //because it  is low

    }
    else if( (a_stencilName.find("Slope_High_")) != std::string::npos)
    {
      //figure out the radius
      auto pos = a_stencilName.find_first_of("0123456789");
      unsigned long dir = 999;
      if (std::string::npos != pos)
      {
        dir = std::stoul(a_stencilName.substr(pos));
      }
      unsigned   int idir = dir;
      

      //for this one srcdomain == dstdomain
      //             srcgraph  == dstgraph
      return
        getSlopeStencil(a_dstVoFs,                 
                        a_stencil,                 
                        a_regStencil,      
                        a_regStencilBCS,
                        a_BCApplyBoxes,
                        a_regApplyBox,  
                        a_domainBCName,
                        a_srcValid,                
                        a_srcDomain,                  
                        a_srcGhost,                
                        a_dstGhost,                
                        a_srcGraph,               
                        a_voludata,
                        a_dstGrown,
                        idir,
                        Side::Hi);//because it  is high

    }
    else
    {
      PR_error("stencil_archive_cell_cell::stencil string not found");
    }
  }
  ///
  template <int order>
  bool
  EBStencilArchive<CELL, CELL,  order, double>::
  getSlopeStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                  vector<LocalStencil<CELL, double> >       & a_stencil,                    
                  Stencil<double>                           & a_regStencil,
                  vector<Stencil<double> >                    a_regStencilBC[2*DIM],       
                  vector<Box>                                 a_BCApplyBoxes[2*DIM],       
                  Box                                       & a_regApplyBox,                
                  string                                      a_domainBCName,
                  const Box                                 & a_validBox,                   
                  const Box                                 & a_domain,                   
                  const Point                               & a_srcGhost,                   
                  const Point                               & a_dstGhost,                   
                  const EBGraph                             & a_graph,                      
                  const VoluData                            & a_voludata,
                  const Point                               & a_dstGrown,
                  unsigned int  a_slopeDir, Side::LoHiSide a_slopeSign)
  {
    Box grownBox = a_validBox.grow(a_dstGrown) & a_domain;

    a_regApplyBox = a_validBox.grow(a_dstGrown); //need to apply outside domain 
    if(a_slopeSign == Side::Hi)
    {
      a_regStencil = (1.0)*Shift(Point::Basis(a_slopeDir)) + (-1.0)*Shift(Point::Zeroes());
    }
    else
    {
      a_regStencil = (1.0)*Shift(Point::Zeroes())          +  (-1.0)*Shift(-Point::Basis(a_slopeDir));  
    }
    
    //this stencil has a span so we have to evaluate the stencil on a grown set of the irregular cells so the regular and irregular data actually talk
    getIrregLocations(a_dstVoFs, a_regStencil, grownBox, a_graph);

    a_stencil.resize(a_dstVoFs.size());
    for(unsigned int idst = 0; idst < a_dstVoFs.size(); idst++)
    {
      a_stencil[idst].clear();
      const EBIndex<CELL>&     dstVoF = a_dstVoFs[idst];

      int isign = sign(a_slopeSign);

      vector< EBIndex<CELL> >otherVoFs = a_graph.getNeighborVoFs(dstVoF, a_slopeDir, a_slopeSign);
      if(otherVoFs.size() > 0)
      {
        double myweight    = -double(isign);
        double otherweight =  double(isign)/otherVoFs.size();
        a_stencil[idst].add(dstVoF, myweight);
        for(int jvof = 0; jvof < otherVoFs.size(); jvof++)
        {
          a_stencil[idst].add(otherVoFs[jvof], otherweight);
        }
      }
      else
      {
        vector< EBIndex<CELL> >flipVoFs = a_graph.getNeighborVoFs(dstVoF, a_slopeDir, flip(a_slopeSign));
        if(flipVoFs.size() > 0)
        {
          double myweight    =  double(isign);
          double flipweight =  -double(isign)/flipVoFs.size();
          a_stencil[idst].add(dstVoF, myweight);
          for(int jvof = 0; jvof < flipVoFs.size(); jvof++)
          {
            a_stencil[idst].add(flipVoFs[jvof], flipweight);
          }
        }
      }
    }    
    
    SecondOrderStencil<order>::
      getDomainBCs(a_regStencilBC, a_BCApplyBoxes, a_domainBCName, a_validBox, a_graph.getDomain(), a_srcGhost);

    //bcs are needed here
    return true;
  }
/******/
  template <int order>
  bool
  EBStencilArchive<CELL, CELL,  order, double>::
  getRestrictionStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                        vector<LocalStencil<CELL, double> >       & a_stencil,                    
                        Stencil<double>                           & a_regStencil,         
                        Box                                       & a_regApplyBox,                
                        const Box                                 & a_srcValid,                   
                        const Box                                 & a_dstValid,                   
                        const Point                               & a_srcGhost,                   
                        const Point                               & a_dstGhost,                   
                        const EBGraph                             & a_srcGraph,                      
                        const EBGraph                             & a_dstGraph,                      
                        const VoluData                            & a_voludata,
                        double                                      a_dx,
                        bool                                        a_volumeWeighted)               
  {
#if DIM==2
    unsigned int ncolor = 4;
#else
    unsigned int ncolor = 8;
#endif    
    //this stencil has no span that I can understand so getIrrregLocations not appropriate
    a_dstVoFs = a_dstGraph.getIrregVoFs(a_dstValid);
    a_regApplyBox = a_dstValid;

    a_regStencil = Stencil<double>();
    double dnumpts = double(ncolor);
    double rweight = 1.0/dnumpts;
    for(unsigned long icolor = 0; icolor < ncolor; icolor++)
    {
      Point colorpt= getColor(icolor);
      a_regStencil += rweight*Shift(colorpt);
    }    
    a_regStencil.srcRatio() = Point::Ones(2);

    a_stencil.resize(a_dstVoFs.size());
    for(int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
    {
      const EBIndex<CELL>&   coarVoF  = a_dstVoFs[ivof];
      vector<EBIndex<CELL> > fineVoFs = a_dstGraph.refine(coarVoF);
      
      if(a_volumeWeighted)
      {
        double sumvol = 0;
        for(int ifine = 0; ifine <fineVoFs.size(); ifine++)
        {
          Point pt = fineVoFs[ifine].m_pt;
          double kappa = 1;
          if(a_srcGraph.isCovered(pt))
          {
            kappa = 0;
          }
          else if(a_voludata.contains(fineVoFs[ifine]))
          {
            RealVect centroid;
            Normalizor<order>::normalizedVoluMoments(kappa, centroid, a_voludata, fineVoFs[ifine], a_srcGraph, a_dx);
          }
          sumvol += kappa;
          a_stencil[ivof].add(fineVoFs[ifine], kappa);
        }
        if(sumvol > 0)
        {
          a_stencil[ivof] *= (1.0/sumvol);
        }
      }
      else
      {
        //in multigrid, the rhs is already volume weighted so it is just
        //coarse = (1/ncolors)*(sum(fine))
        for(int ifine = 0; ifine <fineVoFs.size(); ifine++)
        {
          a_stencil[ivof].add(fineVoFs[ifine], rweight);
        }
      }

    }
    return false;//no domain boundary conditions needed here
  }

/******/
  template <int order>
  bool
  EBStencilArchive<CELL, CELL,  order, double>::
  getProlongationStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                         vector<LocalStencil<CELL, double> >       & a_stencil,                    
                         Stencil<double>                           & a_regStencil,         
                         Box                                       & a_regApplyBox,                
                         const Box                                 & a_srcValid,                   
                         const Box                                 & a_dstValid,                   
                         const Point                               & a_srcGhost,                   
                         const Point                               & a_dstGhost,                   
                         const EBGraph                             & a_srcGraph,                      
                         const EBGraph                             & a_dstGraph,                      
                         unsigned long                               a_icolor)               
  {
    a_regApplyBox = a_srcValid;
    a_regStencil  = Stencil<double>();
    Point colorpt = getColor(a_icolor);
    a_regStencil  =  (1.0)*Shift(Point::Zeros());
    a_regStencil.destRatio() = Point::Ones(2);
    a_regStencil.destShift() = colorpt;
    //since the regular sencil has no span, we can just get the vofs
    a_dstVoFs = a_dstGraph.getIrregVoFs(a_dstValid);

    
    a_stencil.resize(a_dstVoFs.size());
    for(int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
    {
      const EBIndex<CELL>&   fineVoF = a_dstVoFs[ivof];
      EBIndex<CELL>          coarVoF = a_dstGraph.coarsen(fineVoF);
      double weight = 1;
      bool onColor = true;
      Point vofpt = fineVoF.m_pt;
      for(int idir = 0; idir < DIM; idir++)
      {
        int mod = (vofpt[idir])%2;
        if(mod != colorpt[idir])
        {
          onColor = false;
          weight  = 0;
        }
      }
      
      a_stencil[ivof].add(coarVoF, weight);
    }

    return false;//no domain boundary conditions needed here
  }

  ///
  template <int order>
  bool
  EBStencilArchive<CELL, CELL,  order, double>::
  getRedistributionStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                           vector<LocalStencil<CELL, double> >       & a_stencil,                    
                           Stencil<double>                           & a_regStencil,
                           Box                                       & a_regApplyBox,                
                           const Box                                 & a_validBox,                   
                           const Box                                 & a_domain,                   
                           const Point                               & a_srcGhost,                   
                           const Point                               & a_dstGhost,                   
                           const EBGraph                             & a_graph,                      
                           const VoluData                            & a_voludata,
                           int  a_redistRad,
                           double a_dx)                   
  {
    Box grownBox = a_validBox.grow(a_redistRad) & a_domain;
    //there is no reg stencil.  This should be called with  incrementonly = true (initToZero = false)
    a_regStencil = (0.0)*Shift(Point::Zeroes());
    a_regApplyBox = a_validBox;
    
    //since the regular sencil has no span, we can just get the vofs
    vector<EBIndex<CELL> > srcVoFs = a_graph.getIrregVoFs(grownBox);
    for(unsigned int ivof = 0; ivof < srcVoFs.size(); ivof++)
    {
      const EBIndex<CELL>     srcVoF = srcVoFs[ivof];
      vector<EBIndex<CELL > > neighbors = a_graph.getAllVoFsWithinRadius(srcVoF.m_pt, a_redistRad);
      
      LocalStencil<CELL, double> redistSten;
      double totvol = 0;
      for(unsigned int inei = 0; inei < neighbors.size(); inei++)
      {
        double kappa;    RealVect centroid;
        Normalizor<order>::normalizedVoluMoments(kappa, centroid, a_voludata, neighbors[inei], a_graph, a_dx);
        totvol += kappa;
        redistSten.add(neighbors[inei], kappa);
      }
      if(totvol > 0)
      {
        redistSten *= (1./totvol);
      }
      //now incrementally to a_stencil the redistSten.  Only do so if the destination is within the valid box

      for(int isten = 0; isten < redistSten.size(); isten++)
      {
        EBIndex<CELL> stenpt = redistSten.m_entries[isten].m_vof;
        double        weight = redistSten.m_entries[isten].m_weight;
        if(a_validBox.contains(stenpt.m_pt))
        {
          bool found = false;
          int index  = -1;
          for(int idst = 0; idst < a_dstVoFs.size(); idst++)
          {
            if(a_dstVoFs[idst] == stenpt)
            {
              index = idst;
              found = true;
            }
          }
          if(found)
          {
            a_stencil[index].add(stenpt, weight);
          }
          else
          {
            a_dstVoFs.push_back(stenpt);
            LocalStencil<CELL, double> dststen;
            dststen.add(stenpt, weight);
            a_stencil.push_back(dststen);
          }
        }
      }
    }    
    //no bcs to apply here
    return false;
  }
  ///
  template <int order>
  bool
  EBStencilArchive<CELL, CELL,  order, double>::
  getNormalizationStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                          vector<LocalStencil<CELL, double> >       & a_stencil,                    
                          Stencil<double>                           & a_regStencil,         
                          Box                                       & a_regApplyBox,                
                          const Box                                 & a_validBox,                   
                          const Box                                 & a_domain,                     
                          const Point                               & a_srcGhost,                   
                          const Point                               & a_dstGhost,                   
                          const EBGraph                             & a_graph,                      
                          const VoluData                            & a_voludata,                   
                          unsigned long                             & a_radius,
                          double                                      a_dx)
  {
    //regular stencil is the identity matrix
    a_regStencil = (1.0)*Shift(Point::Zeroes());
    a_regApplyBox = a_validBox;
    //since the regular sencil has no span, we can just get the vofs
    a_dstVoFs = a_graph.getIrregVoFs(a_validBox);
    a_stencil.resize(a_dstVoFs.size());
    int ideb = 0;
    for(unsigned int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
    {
      vector<EBIndex<CELL> > neighbors;
      Point pt = a_dstVoFs[ivof].m_pt;
      Box ivbox(pt, pt);
      Box neighborhood = ivbox.grow(a_radius) & a_graph.getDomain();
      for(auto bit = neighborhood.begin(); bit != neighborhood.end(); ++bit)
      {
        vector<EBIndex<CELL> > cellVoFs = a_graph.getVoFs(*bit);
        neighbors.insert(neighbors.end(), cellVoFs.begin(), cellVoFs.end());
      }

      //because stuff is kappa weighted already, all stencil weights are 1/sum(kappa)
      double totvol = 0;
      for(unsigned int inei = 0; inei < neighbors.size(); inei++)
      {
        double kappa;    RealVect centroid;
        Normalizor<order>::normalizedVoluMoments(kappa, centroid, a_voludata, neighbors[inei], a_graph, a_dx);
        totvol += kappa;
      }
      double weight = 0;
      if(totvol > 0)
      {
        weight = (1./totvol);
      }
      for(unsigned int inei = 0; inei < neighbors.size(); inei++)
      {
        a_stencil[ivof].add(neighbors[inei], weight);
      }
      ideb++;
    }

    //no bcs to apply here
    return false;
  }              
/***/           
  template <int order>
  void
  EBStencilArchive<CELL, CELL,  order, double>::    
  getIrregLocations(vector<EBIndex<CELL> >        & a_dstVoFs,                    
                    const Stencil<double>         & a_regStencil,         
                    const Box                     & a_boundingBox,
                    const EBGraph                 & a_dstGraph)
  {
    a_dstVoFs.clear();
    PointSet irregset = a_dstGraph.getIrregCells(a_boundingBox);
    Box span = a_regStencil.span();
    Point lo = span.low();
    Point hi = span.high();
    Point grow;
    for(int idir = 0; idir < DIM; idir++)
    {
      grow[idir] = std::max(std::abs(lo[idir]), std::abs(hi[idir]));
      irregset.grow(idir, grow[idir]);
    }

    irregset &= a_dstGraph.getDomain();
    irregset &= a_boundingBox;
    for(PointSetIterator ivsit(irregset); ivsit.ok(); ++ivsit)
    {
      vector<EBIndex<CELL> > vofs = a_dstGraph.getVoFs(ivsit());
      a_dstVoFs.insert(a_dstVoFs.end(), vofs.begin(), vofs.end());
    }
  }
  ///this deep down in the bowels, everything is a template specialization
  template <int order>
  bool
  EBStencilArchive<CELL, CELL,  order, double>::
  getSecondOrderPoissonStencil(vector<EBIndex<CELL> >                    & a_dstVoFs,                    
                               vector<LocalStencil<CELL, double> >       & a_stencil,                    
                               Stencil<double>                           & a_regStencil,         
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

    a_regStencil = Stencil<double>::Laplacian();
    a_regStencil *= 1.0/(a_dx*a_dx);
    a_regApplyBox = a_validBox;
    //this stencil has a span so we have to evaluate the stencil on a grown set of the irregular cells so the regular and irregular data actually talk
    getIrregLocations(a_dstVoFs, a_regStencil, a_validBox, a_graph);

    Poisson2ndOrder<XFACE, order> xflux;
    Poisson2ndOrder<YFACE, order> yflux;
    Poisson2ndOrder<ZFACE, order> zflux;
    Poisson2ndOrder<BOUNDARY, order> ebflux;
    a_stencil.resize(a_dstVoFs.size());
    for(int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
    {
      Point pt = a_dstVoFs[ivof].m_pt;
      SecondOrderStencil<order>::
        get2ndOrderDivFStencil(a_stencil[ivof],
                               a_dstVoFs[ivof],
                               a_graph,
                               a_voludata,
                               a_ebfadata,
                               a_xfacdata, a_yfacdata, a_zfacdata,
                               xflux, yflux, zflux, a_ebbcName,
                               a_dx);
//begin debug
//      int ideb = 0;
//      if( (pt[0]==18) && (pt[1]==30))
//      {
//        ideb = 1;
//      }
//      if(ideb == 1)
//      {
//        a_stencil[ivof].print();
//      }
//end debug
    }

    if(!a_periodiconly)
    {
      SecondOrderStencil<order>::
        getDomainBCs(a_regStencilBC, a_BCApplyBoxes, a_domainBCName, a_validBox, a_graph.getDomain(), a_srcGhost);
    }
    return true;
  }
  /***/
  template < int order>
  void
  DirichletEBFlux<order>::
  calculateWeightingMatrix(RealVect           x0,
                           vector<RealVect>&  xp,
                           vector<RealVect>&  weightMatrix,
                           bool&              detZero)
  {
    int stenSize = xp.size();

    vector<RealVect> deltaX = xp;
    for (int isten = 0; isten < stenSize; isten++)
    {
      deltaX[isten] -= x0;
    }

    vector<RealVect>    aTransA(DIM, RealVect::Zero());
    vector<RealVect> invATransA(DIM, RealVect::Zero());
    for (int idir = 0; idir < DIM; idir++)
    {
      for (int jdir = 0; jdir < DIM; jdir++)
      {
        for (int isten = 0; isten < stenSize; isten++)
        {
          aTransA[idir][jdir] = aTransA[idir][jdir]
            + deltaX[isten][idir]*deltaX[isten][jdir];
        }
      }
    }

    double det;
#if DIM == 2
    det = aTransA[0][0] * aTransA[1][1] - aTransA[0][1] * aTransA[1][0];
    if (det < 1.e-15 && det > -1.e-15)
    {
      detZero = true;
    }
    else
    {
      invATransA[0][0] =  aTransA[1][1] / det;
      invATransA[0][1] = -aTransA[0][1] / det;
      invATransA[1][0] = -aTransA[1][0] / det;
      invATransA[1][1] =  aTransA[0][0] / det;
    }
#elif DIM == 3
    det = aTransA[0][0] * ( aTransA[1][1] * aTransA[2][2]
                            - aTransA[1][2] * aTransA[2][1])
      + aTransA[0][1] * ( aTransA[1][2] * aTransA[2][0]
                          - aTransA[1][0] * aTransA[2][2])
      + aTransA[0][2] * ( aTransA[1][0] * aTransA[2][1]
                          - aTransA[1][1] * aTransA[2][0]);

    if (det < 1.e-15 && det > -1.e-15)
    {
      detZero = true;
    }
    else
    {
      invATransA[0][0] = ( aTransA[1][1] * aTransA[2][2]
                           - aTransA[1][2] * aTransA[2][1]) / det;
      invATransA[0][1] = ( aTransA[1][2] * aTransA[2][0]
                           - aTransA[1][0] * aTransA[2][2]) / det;
      invATransA[0][2] = ( aTransA[1][0] * aTransA[2][1]
                           - aTransA[1][1] * aTransA[2][0]) / det;
      invATransA[1][0] = ( aTransA[2][1] * aTransA[0][2]
                           - aTransA[2][2] * aTransA[0][1]) / det;
      invATransA[1][1] = ( aTransA[2][2] * aTransA[0][0]
                           - aTransA[2][0] * aTransA[0][2]) / det;
      invATransA[1][2] = ( aTransA[2][0] * aTransA[0][1]
                           - aTransA[2][1] * aTransA[0][0]) / det;
      invATransA[2][0] = ( aTransA[0][1] * aTransA[1][2]
                           - aTransA[0][2] * aTransA[1][1]) / det;
      invATransA[2][1] = ( aTransA[0][2] * aTransA[1][0]
                           - aTransA[0][0] * aTransA[1][2]) / det;
      invATransA[2][2] = ( aTransA[0][0] * aTransA[1][1]
                           - aTransA[0][1] * aTransA[1][0]) / det;
    }
#else
    THIS_IS_AN_ERROR_MESSAGE__THIS_WILL_ONLY_COMPILE_WHEN_DIM_IS_2_OR_3;
#endif

    //if (!detZero)
    {
      weightMatrix = vector<RealVect>(stenSize,RealVect::Zero());
      for (int idir = 0; idir < DIM; idir++)
      {
        for (int isten = 0; isten < stenSize; isten++)
        {
          for (int jdir = 0; jdir < DIM; jdir++)
          {
            weightMatrix[isten][idir] += invATransA[idir][jdir] * deltaX[isten][jdir];
          }
        }
      }
    }
  }
  
  /***/
  template < int order>
  void
  DirichletEBFlux<order>::
  getEBFluxStencil(LocalStencil<CELL, double>    &  a_stencil,
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
    EBIndex<CELL> vof = a_face.getVoF();
    unsigned int radius = 1;
    vector< EBIndex<CELL> > vols=  a_graph.getAllVoFsWithinRadius(vof.m_pt, radius);
    RealVect centroid, normal;
    double bndryArea;
    Normalizor<order>::normalizedBndryMoments(centroid, a_ebfadata, a_face, a_graph, a_dx);
    Normalizor<order>::getNormalAndBoundaryArea(bndryArea, normal,   
                                                a_voludata,
                                                a_ebfadata,
                                                a_xfacdata,
                                                a_yfacdata,
                                                a_zfacdata, 
                                                vof, a_graph, a_dx);

#if DIM==2
    unsigned int minStenSize = 3;
#else 
    unsigned int minStenSize = 7;
#endif

    bool dropOrder= (vols.size() < minStenSize);
    
    if (!dropOrder)
    {
      RealVect x0;
      Point iv0 = vof.m_pt;
      for (int idir = 0; idir < DIM; idir++)
      {
        x0[idir] = a_dx * (0.5 + centroid[idir] + iv0[idir]);
      }

      vector<RealVect> xp(vols.size());
      for (int isten = 0; isten < vols.size(); isten++)
      {
        for (int idir = 0; idir < DIM; idir++)
        {
          xp[isten][idir] = a_dx * (0.5 + vols[isten].m_pt[idir]);
        }
      }

      unsigned int volsize = vols.size();
      vector<RealVect> invATransAdeltaX(volsize, RealVect::Zero());


      bool detZero = false;
      calculateWeightingMatrix(x0, xp, invATransAdeltaX, detZero);

      a_stencil.clear();
      a_weight = 0.0;

      for (int isten = 0; isten < vols.size(); isten++)
      {
        double dphidnWeight = 0.0;
        for (int idir = 0; idir < DIM; idir++)
        {
          dphidnWeight -= invATransAdeltaX[isten][idir] * normal[idir];
        }
              
        a_stencil.add(vols[isten],dphidnWeight);
        a_weight -= dphidnWeight;
      }
    }
    else
    {
      a_stencil.clear();
      a_weight = 0.0;
    }
    
  }
  /***/
  template < int order>
  void
  SecondOrderStencil<order>::
  getDomainBCs(vector<Stencil<double> >    a_regStencilBC[2*DIM],       
               vector<Box>                 a_BCApplyBoxes[2*DIM],       
               string                      a_domainBCName,               
               const Box                 & a_validBox,
               const Box                 & a_domain,
               const Point               & a_srcGhost)
    
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
        loapp[idir] = lodom[idir]-1;
        hiapp[idir] = lodom[idir]-1;

        Box appBox0(loapp, hiapp);

        a_regStencilBC[loindex].resize(numghost);
        a_BCApplyBoxes[loindex].resize(numghost);
        for(int igh = 0; igh < numghost; igh++)
        {
          a_BCApplyBoxes[loindex][igh] = appBox0.shift(idir, -igh);
          a_regStencilBC[loindex][igh] = (stenSign)*Shift(Point::Basis(idir, 2*igh+1));
        }
      }
      if(hibox[idir] == hidom[idir])
      {
        int hiindex = ebp_index(idir, Side::Hi);
        Point loapp = srcBox.low();
        Point hiapp = srcBox.high();
        loapp[idir] = hidom[idir]+1;
        hiapp[idir] = hidom[idir]+1;

        Box appBox0(loapp, hiapp);
        
        a_regStencilBC[hiindex].resize(numghost);
        a_BCApplyBoxes[hiindex].resize(numghost);
        for(int igh = 0; igh < numghost; igh++)
        {
          a_BCApplyBoxes[hiindex][igh] = appBox0.shift(idir,  igh);
          a_regStencilBC[hiindex][igh] = (stenSign)*Shift(Point::Basis(idir, -2*igh-1));
        }
      }
    }
  }
/////////////end cell to cell

  ///begin CELL to fcent
  ///this deep down in the bowels, everything is a template specialization
  ///returns whether there are domain boundary conditions
  template <CENTERING fcent, int order>
  bool 
  EBStencilArchive<CELL,fcent,  order, double>::
  getStencil(vector<EBIndex<fcent> >                   & a_dstFaces,
             vector<LocalStencil<CELL, double> >       & a_stencil,
             Stencil<double>                           & a_regStencil,
             vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
             vector<Box>                                 a_BCApplyBoxes[2*DIM],       
             Box                                       & a_regApplyBox,
             string                                      a_stencilName,
             string                                      a_domainBCName,
             string                                      a_ebbcName,
             const Box                                 & a_srcValid,
             const Box                                 & a_dstValid,
             const Box                                 & a_srcDomain,
             const Box                                 & a_dstDomain,
             const Point                               & a_srcGhost,
             const Point                               & a_dstGhost,
             const EBGraph                             & a_srcGraph,
             const EBGraph                             & a_dstGraph,
             const VoluData                            & a_voludata,
             const EBFaData                            & a_ebfadata,
             const XFacData                            & a_xfacdata,
             const YFacData                            & a_yfacdata,
             const ZFacData                            & a_zfacdata,
             double                                      a_dx,
             bool                                        a_periodiconly, 
             Point                                       a_dstGrown)               
    {
      if(a_stencilName == string("AverageCellToFace"))
      {
        return
          getCellToFaceAveStencil(a_dstFaces,
                                  a_stencil,
                                  a_regStencil,
                                  a_regStencilBCS,
                                  a_BCApplyBoxes,       
                                  a_regApplyBox,
                                  a_domainBCName,
                                  a_srcValid,
                                  a_srcDomain,
                                  a_srcGhost,
                                  a_dstGhost,
                                  a_dstGraph,
                                  a_dx,
                                  a_dstGrown);

      }
      else if(a_stencilName == string("MAC_Gradient"))
      {
        return
          getMACGradientStencil(a_dstFaces,
                                a_stencil,
                                a_regStencil,
                                a_regStencilBCS,
                                a_BCApplyBoxes,       
                                a_regApplyBox,
                                a_domainBCName,
                                a_srcValid,
                                a_srcDomain,
                                a_srcGhost,
                                a_dstGhost,
                                a_dstGraph,
                                a_dx,
                                a_dstGrown);

      }
      else if(a_stencilName == string("Cell_To_Face_Low"))
      {
        return
          getCellToFaceStencil(a_dstFaces,
                               a_stencil,
                               a_regStencil,
                               a_regApplyBox,
                               a_srcValid,
                               a_srcDomain,
                               a_srcGhost,
                               a_dstGhost,
                               a_dstGraph,
                               a_dx,
                               a_dstGrown,
                               Side::Lo);
                               
      }
      else if(a_stencilName == string("Cell_To_Face_High"))
      {
                               
        return
          getCellToFaceStencil(a_dstFaces,
                               a_stencil,
                               a_regStencil,
                               a_regApplyBox,
                               a_srcValid,
                               a_srcDomain,
                               a_srcGhost,
                               a_dstGhost,
                               a_dstGraph,
                               a_dx,
                               a_dstGrown,
                               Side::Hi);
      }
      else
      {
        PROTO_ASSERT(false, "string not  found");
      }

      return true;
    }

///
  template <CENTERING fcent, int order>
  bool 
  EBStencilArchive<CELL,fcent,  order, double>::
  getCellToFaceStencil(vector<EBIndex<fcent> >                   & a_dstFaces,
                       vector<LocalStencil<CELL, double> >       & a_stencil,
                       Stencil<double>                           & a_regStencil,
                       Box                                       & a_regApplyBox,
                       const Box                                 & a_validBox,
                       const Box                                 & a_srcDomain,
                       const Point                               & a_srcGhost,
                       const Point                               & a_dstGhost,
                       const EBGraph                             & a_graph,
                       double                                      a_dx,
                       Point                                       a_dstGrown,
                       Side::LoHiSide                              a_side)               
  {
    int facedir;
    if(fcent == XFACE)
    {
      facedir = 0;
    }
    else if(fcent == YFACE)
    {
      facedir = 1;
    }
    else if(fcent == ZFACE)
    {
      facedir = 2;
    }
    else
    {
      PROTO_ASSERT(false, "bogus face direction");
    }
    if(a_side == Side::Lo)
    {
      a_regStencil = (1.0)*Shift(Point::Basis(facedir, -1));
    }
    else if(a_side == Side::Hi)
    {
      a_regStencil = (1.0)*Shift(Point::Zeros());
    }
    else
    {
      PROTO_ASSERT(false, "bogus side");
    }

    Box faceBox = a_validBox.growHi(facedir, 1);
    Box grownFace = faceBox;
    for(int jdir = 0; jdir < DIM; jdir++)
    {
      if(jdir != facedir)
      {
        grownFace = grownFace.grow(jdir, a_dstGrown[jdir]);
      }
    }
    a_regApplyBox = grownFace;

    Box grown = a_validBox.grow(a_dstGrown);
    a_dstFaces = getIrregFaces<fcent>(grown, a_graph);
    a_stencil.resize(a_dstFaces.size());
    for(unsigned int iface = 0; iface < a_dstFaces.size(); iface++)
    {
      EBIndex<CELL> vof = a_dstFaces[iface].getVoF(a_side);
      a_stencil[iface].add(vof, 1.0);
    }

    return false;
  }
///
  template <CENTERING fcent, int order>
  bool 
  EBStencilArchive<CELL,fcent,  order, double>::
  getCellToFaceAveStencil(vector<EBIndex<fcent> >                   & a_dstFaces,
                          vector<LocalStencil<CELL, double> >       & a_stencil,
                          Stencil<double>                           & a_regStencil,
                          vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
                          vector<Box>                                 a_BCApplyBoxes[2*DIM],       
                          Box                                       & a_regApplyBox,
                          string                                      a_domainBCName,
                          const Box                                 & a_validBox,
                          const Box                                 & a_srcDomain,
                          const Point                               & a_srcGhost,
                          const Point                               & a_dstGhost,
                          const EBGraph                             & a_graph,
                          double                                      a_dx,
                          Point                                       a_dstGrown)
  {
    int facedir;
    if(fcent == XFACE)
    {
      facedir = 0;
    }
    else if(fcent == YFACE)
    {
      facedir = 1;
    }
    else if(fcent == ZFACE)
    {
      facedir = 2;
    }
    else
    {
      PROTO_ASSERT(false, "bogus face direction");
    }
    a_regStencil = (0.5)*Shift(Point::Zeros()) + (0.5)*Shift(Point::Basis(facedir, -1));
    Box faceBox = a_validBox.growHi(facedir, 1);
    Box grownFace = faceBox;
    for(int jdir = 0; jdir < DIM; jdir++)
    {
      if(jdir != facedir)
      {
        grownFace = grownFace.grow(jdir, a_dstGrown[jdir]);
      }
    }
    a_regApplyBox = grownFace;

    Box grown = a_validBox.grow(a_dstGrown);
    a_dstFaces = getIrregFaces<fcent>(grown, a_graph);
    a_stencil.resize(a_dstFaces.size());
    for(unsigned int iface = 0; iface < a_dstFaces.size(); iface++)
    {
      EBIndex<CELL> loVoF = a_dstFaces[iface].getVoF(Side::Lo);
      EBIndex<CELL> hiVoF = a_dstFaces[iface].getVoF(Side::Hi);
      a_stencil[iface].add(loVoF, 0.5);
      a_stencil[iface].add(hiVoF, 0.5);
    }

   
    SecondOrderStencil<order>::
      getDomainBCs(a_regStencilBCS, a_BCApplyBoxes, a_domainBCName, a_validBox, a_graph.getDomain(), a_srcGhost);
    return true;
  }

///
  template <CENTERING fcent, int order>
  bool 
  EBStencilArchive<CELL,fcent,  order, double>::
  getMACGradientStencil(vector<EBIndex<fcent> >                   & a_dstFaces,
                        vector<LocalStencil<CELL, double> >       & a_stencil,
                        Stencil<double>                           & a_regStencil,
                        vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
                        vector<Box>                                 a_BCApplyBoxes[2*DIM],       
                        Box                                       & a_regApplyBox,
                        string                                      a_domainBCName,
                        const Box                                 & a_validBox,
                        const Box                                 & a_srcDomain,
                        const Point                               & a_srcGhost,
                        const Point                               & a_dstGhost,
                        const EBGraph                             & a_graph,
                        double                                      a_dx,
                        Point                                       a_dstGrown)
  {
    int facedir;
    if(fcent == XFACE)
    {
      facedir = 0;
    }
    else if(fcent == YFACE)
    {
      facedir = 1;
    }
    else if(fcent == ZFACE)
    {
      facedir = 2;
    }
    else
    {
      PROTO_ASSERT(false, "bogus face direction");
    }
    a_regStencil = (1.0/a_dx)*Shift(Point::Zeros()) - (1.0/a_dx)*Shift(Point::Basis(facedir, -1));
    Box faceBox = a_validBox.growHi(facedir, 1);
    Box grownFace = faceBox;
    for(int jdir = 0; jdir < DIM; jdir++)
    {
      if(jdir != facedir)
      {
        grownFace = grownFace.grow(jdir, a_dstGrown[jdir]);
      }
    }
    a_regApplyBox = grownFace;

    Box grown = a_validBox.grow(a_dstGrown);
    a_dstFaces = getIrregFaces<fcent>(grown, a_graph);
    a_stencil.resize(a_dstFaces.size());
    for(unsigned int iface = 0; iface < a_dstFaces.size(); iface++)
    {
      EBIndex<CELL> loVoF = a_dstFaces[iface].getVoF(Side::Lo);
      EBIndex<CELL> hiVoF = a_dstFaces[iface].getVoF(Side::Hi);
      a_stencil[iface].add(loVoF, -1.0/a_dx);
      a_stencil[iface].add(hiVoF,  1.0/a_dx);
    }


    SecondOrderStencil<order>::
      getDomainBCs(a_regStencilBCS, a_BCApplyBoxes, a_domainBCName, a_validBox, a_graph.getDomain(), a_srcGhost);
    return true;
  }

////end cell to face
  ///begin  fcent to cell

  ///this deep down in the bowels, everything is a template specialization
  ///returns whether there are domain boundary conditions
  template <CENTERING fcent, int order>
  bool 
  EBStencilArchive<fcent,  CELL, order, double>::
  getStencil(vector<EBIndex<CELL> >                   & a_dstVoFs,
             vector<LocalStencil<fcent, double> >       & a_stencil,
             Stencil<double>                           & a_regStencil,
             vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
             vector<Box>                                 a_BCApplyBoxes[2*DIM],       
             Box                                       & a_regApplyBox,
             string                                      a_stencilName,
             string                                      a_domainBCName,
             string                                      a_ebbcName,
             const Box                                 & a_srcValid,
             const Box                                 & a_dstValid,
             const Box                                 & a_srcDomain,
             const Box                                 & a_dstDomain,
             const Point                               & a_srcGhost,
             const Point                               & a_dstGhost,
             const EBGraph                             & a_srcGraph,
             const EBGraph                             & a_dstGraph,
             const VoluData                            & a_voludata,
             const EBFaData                            & a_ebfadata,
             const XFacData                            & a_xfacdata,
             const YFacData                            & a_yfacdata,
             const ZFacData                            & a_zfacdata,
             double                                      a_dx,
             bool                                        a_periodiconly, 
             Point                                       a_dstGrown)               
    {

      if(a_stencilName == string("Divergence"))
      {
        return
          getDivergenceStencil(a_dstVoFs,
                               a_stencil,
                               a_regStencil,
                               a_regApplyBox,
                               a_srcValid,
                               a_srcDomain,
                               a_srcGhost,
                               a_dstGhost,
                               a_voludata,
                               a_ebfadata,
                               a_xfacdata,
                               a_yfacdata,
                               a_zfacdata,
                               a_dstGraph,
                               a_dx);
                               
      }
      else
      {
        PROTO_ASSERT(false, "string not  found");
      }

      return true;
    }

  template <CENTERING fcent, int order>
  void
  EBStencilArchive<fcent, CELL,  order, double>::    
  getIrregLocations(vector<EBIndex<CELL> >        & a_dstVoFs,                    
                    const Stencil<double>         & a_regStencil,         
                    const Box                     & a_boundingBox,
                    const EBGraph                 & a_dstGraph,
                    const int                     & a_facedir)
  {
    a_dstVoFs.clear();
    PointSet irregset = a_dstGraph.getIrregCells(a_boundingBox);
    Box span = a_regStencil.span();
    Point lo = span.low();
    Point hi = span.high();
    Point grow = Point::Zero();
    int idir = a_facedir;
    grow[idir] = std::max(std::abs(lo[idir]), std::abs(hi[idir]));
    irregset.grow(idir, grow[idir]);


    irregset &= a_dstGraph.getDomain();
    irregset &= a_boundingBox;
    for(PointSetIterator ivsit(irregset); ivsit.ok(); ++ivsit)
    {
      vector<EBIndex<CELL> > vofs = a_dstGraph.getVoFs(ivsit());
      a_dstVoFs.insert(a_dstVoFs.end(), vofs.begin(), vofs.end());
    }
  }

  template <CENTERING fcent, int order>
  bool 
  EBStencilArchive<fcent, CELL, order, double>::
  getDivergenceStencil(vector<EBIndex<CELL > >                   & a_dstVoFs,
                       vector<LocalStencil<fcent, double> >      & a_stencil,
                       Stencil<double>                           & a_regStencil,
                       Box                                       & a_regApplyBox,
                       const Box                                 & a_validBox,
                       const Box                                 & a_srcDomain,
                       const Point                               & a_srcGhost,
                       const Point                               & a_dstGhost,
                       const VoluData                            & a_voludata,
                       const EBFaData                            & a_ebfadata,
                       const XFacData                            & a_xfacdata,
                       const YFacData                            & a_yfacdata,
                       const ZFacData                            & a_zfacdata,
                       const EBGraph                             & a_graph,
                       double                                      a_dx)
  {
    int facedir;
    if(fcent == XFACE)
    {
      facedir = 0;
    }
    else if(fcent == YFACE)
    {
      facedir = 1;
    }
    else if(fcent == ZFACE)
    {
      facedir = 2;
    }
    else
    {
      PROTO_ASSERT(false, "bogus face direction");
    }
    
    a_regStencil = (1.0/a_dx)*Shift(Point::Basis(facedir, 1)) +  (-1.0/a_dx)*Shift(Point::Zeros());
    a_regApplyBox = a_validBox;

    //this stencil has a span so we have to evaluate the stencil on a grown set of the irregular cells so the regular and irregular data actually talk
    getIrregLocations(a_dstVoFs, a_regStencil, a_validBox, a_graph, facedir);

    a_stencil.resize(a_dstVoFs.size());
    for(unsigned int ivof = 0; ivof < a_dstVoFs.size(); ivof++)
    {
      EBIndex<CELL> vof = a_dstVoFs[ivof];
      RealVect centroid;
      double kappa;
      RealVect normal;
      double bndryArea;
      Normalizor<order>::getNormalAndBoundaryArea(bndryArea, normal,   
                                                a_voludata,
                                                a_ebfadata,
                                                a_xfacdata,
                                                a_yfacdata,
                                                a_zfacdata, 
                                                vof, a_graph, a_dx);
      Normalizor<order>::normalizedVoluMoments(kappa, centroid, a_voludata, vof, a_graph, a_dx);
      //begin debug
//      int ideb = 0;
//      if((vof.m_pt[0] == 55) && (vof.m_pt[1]==15))
//      {
//        ideb = 1;
//      }
      //end debug
      
      for(SideIterator sit; sit.ok(); ++sit)
      {
        vector<EBIndex<fcent> > faces = getSideFaces<fcent>(vof, sit(), a_graph);
        for(int iface = 0; iface < faces.size(); iface++)
        {
          const EBIndex<fcent>& face = faces[iface];
          double areaFrac;
          RealVect centroid;
          Normalizor<order>::normalizedFaceMoments(areaFrac, centroid, a_xfacdata, a_yfacdata, a_zfacdata, face, a_graph, a_dx, facedir);
          double weight = sign(sit())*areaFrac/a_dx;
//debug
//          weight /= kappa;
// end debug
          a_stencil[ivof].add(face, weight);
        }
      }
    }
    
    //no bcs here.
    return false;
  }
///end face to cell



  ///begin  fcent to fcent
  ///this deep down in the bowels, everything is a template specialization
  ///returns whether there are domain boundary conditions
  template <CENTERING fcent, int order>
  bool 
  EBStencilArchive<fcent,  fcent, order, double>::
  getStencil(vector<EBIndex<fcent> >                   & a_dstVoFs,
             vector<LocalStencil<fcent, double> >       & a_stencil,
             Stencil<double>                           & a_regStencil,
             vector<Stencil<double> >                    a_regStencilBCS[2*DIM],
             vector<Box>                                 a_BCApplyBoxes[2*DIM],       
             Box                                       & a_regApplyBox,
             string                                      a_stencilName,
             string                                      a_domainBCName,
             string                                      a_ebbcName,
             const Box                                 & a_srcValid,
             const Box                                 & a_dstValid,
             const Box                                 & a_srcDomain,
             const Box                                 & a_dstDomain,
             const Point                               & a_srcGhost,
             const Point                               & a_dstGhost,
             const EBGraph                             & a_srcGraph,
             const EBGraph                             & a_dstGraph,
             const VoluData                            & a_voludata,
             const EBFaData                            & a_ebfadata,
             const XFacData                            & a_xfacdata,
             const YFacData                            & a_yfacdata,
             const ZFacData                            & a_zfacdata,
             double                                      a_dx,
             bool                                        a_periodiconly, 
             Point                                       a_dstGrown)
    {

      if(a_stencilName == string("InterpolateToFaceCentroid"))
      {
        return
          getFaceCentroidInterpStencil(a_dstVoFs,
                                       a_stencil,
                                       a_regStencil,
                                       a_regApplyBox,
                                       a_srcValid,
                                       a_srcDomain,
                                       a_srcGhost,
                                       a_dstGhost,
                                       a_xfacdata,
                                       a_yfacdata,
                                       a_zfacdata,
                                       a_dstGraph,
                                       a_dx);
      }
      else
      {
        PROTO_ASSERT(false, "string not  found");
      }

      return true;
    }

  template <CENTERING fcent, int order>
  inline LocalStencil<fcent, double>
  getFaceCentroidInterpStencil(const EBIndex<fcent>                           & a_face, 
                               const HostIrregData<XFACE, IndexedMoments<DIM-1, order>, 1>  & a_xfacdata,
                               const HostIrregData<YFACE, IndexedMoments<DIM-1, order>, 1>  & a_yfacdata,
                               const HostIrregData<ZFACE, IndexedMoments<DIM-1, order>, 1>  & a_zfacdata,
                               const EBGraph                                  & a_graph,
                               const double                                   & a_dx)
  {
    PROTO_ASSERT(false, "not implemented");
  }


  template <int order>
  inline LocalStencil<XFACE, double>
  getFaceCentroidInterpStencil(const EBIndex<XFACE>                           & a_face, 
                               const HostIrregData<XFACE, IndexedMoments<DIM-1, order>, 1>  & a_xfacdata,
                               const HostIrregData<YFACE, IndexedMoments<DIM-1, order>, 1>  & a_yfacdata,
                               const HostIrregData<ZFACE, IndexedMoments<DIM-1, order>, 1>  & a_zfacdata,
                               const EBGraph                                  & a_graph,
                               const double                                   & a_dx)
  {
    LocalStencil<XFACE, double> retval;
    Base2ndOrderFlux<XFACE, order> xflux;
    xflux.getFaceInterpStencil(retval, a_face, a_graph, a_xfacdata, a_dx);
//begin debug
//#if DIM==2
//      Point pt(16, 18);
//      Point fp = a_face.m_pt;
//      int ideb = 0;
//      if(fp == pt)
//      {
//        ideb = 1;
//      }
//#endif      
//end debug

    return retval;
  }
  

  template <int order>
  inline LocalStencil<YFACE, double>
  getFaceCentroidInterpStencil(const EBIndex<YFACE>                           & a_face, 
                               const HostIrregData<XFACE, IndexedMoments<DIM-1, order>, 1>  & a_xfacdata,
                               const HostIrregData<YFACE, IndexedMoments<DIM-1, order>, 1>  & a_yfacdata,
                               const HostIrregData<ZFACE, IndexedMoments<DIM-1, order>, 1>  & a_zfacdata,
                               const EBGraph                                  & a_graph,
                               const double                                   & a_dx)
  {
    LocalStencil<YFACE, double> retval;
    Base2ndOrderFlux<YFACE, order> yflux;
    yflux.getFaceInterpStencil(retval, a_face, a_graph, a_yfacdata, a_dx);
    return retval;
  }
  

  template <int order>
  inline LocalStencil<ZFACE, double>
  getFaceCentroidInterpStencil(const EBIndex<ZFACE>                           & a_face, 
                               const HostIrregData<XFACE, IndexedMoments<DIM-1, order>, 1>  & a_xfacdata,
                               const HostIrregData<YFACE, IndexedMoments<DIM-1, order>, 1>  & a_yfacdata,
                               const HostIrregData<ZFACE, IndexedMoments<DIM-1, order>, 1>  & a_zfacdata,
                               const EBGraph                                  & a_graph,
                               const double                                   & a_dx)
  {
    LocalStencil<ZFACE, double> retval;
    Base2ndOrderFlux<ZFACE, order> zflux;
    zflux.getFaceInterpStencil(retval, a_face, a_graph, a_zfacdata, a_dx);
    return retval;
  }
  
  template <CENTERING fcent, int order>
  bool 
  EBStencilArchive<fcent, fcent, order, double>::
  getFaceCentroidInterpStencil(vector<EBIndex<fcent > >                  & a_dstFaces,
                               vector<LocalStencil<fcent, double> >      & a_stencil,
                               Stencil<double>                           & a_regStencil,
                               Box                                       & a_regApplyBox,
                               const Box                                 & a_validBox,
                               const Box                                 & a_srcDomain,
                               const Point                               & a_srcGhost,
                               const Point                               & a_dstGhost,
                               const XFacData                            & a_xfacdata,
                               const YFacData                            & a_yfacdata,
                               const ZFacData                            & a_zfacdata,
                               const EBGraph                             & a_graph,
                               const double                              & a_dx)
  {
    int facedir;
    if(fcent == XFACE)
    {
      facedir = 0;
    }
    else if(fcent == YFACE)
    {
      facedir = 1;
    }
    else if(fcent == ZFACE)
    {
      facedir = 2;
    }
    else
    {
      PROTO_ASSERT(false, "bogus face direction");
    }
    //no regular stencil span so we don't need to do anything fancy with dstFaces
    a_regStencil = (1.0)*Shift(Point::Zeros());
    a_regApplyBox = a_validBox.growHi(facedir, 1);

    a_dstFaces = getIrregFaces<fcent>(a_validBox, a_graph);

    a_stencil.resize(a_dstFaces.size());
    for(unsigned int iface = 0; iface < a_dstFaces.size(); iface++)
    {
      const EBIndex<fcent>& face = a_dstFaces[iface];
      a_stencil[iface] =
        Proto::getFaceCentroidInterpStencil(face, a_xfacdata, a_yfacdata, a_zfacdata, a_graph, a_dx);
    }
    //no bcs here.
    return false;
  }
///end face to cell

}
#endif
