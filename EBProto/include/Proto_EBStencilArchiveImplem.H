
#ifndef _Proto_EBStencilArchiveImplem_H_
#define _Proto_EBStencilArchiveImplem_H_

namespace Proto
{

  /*****************************/
  template<CENTERING cent, int order>
  void
  Base2ndOrderFlux<cent, order>
  getFaceInterpStencil2D(LocalStencil<cent, double>    &  a_sten,
                         const EBIndex<cent>           &  a_face,
                         const EBGraph                 &  a_graph,
                         const VoluData                &  a_voludata,
                         const EBFaData                &  a_ebfadata,
                         const XFacData                &  a_xfacdata,
                         const YFacData                &  a_yfacdata,
                         const ZFacData                &  a_zfacdata,
                         const double                  &  a_dx)
{
  CH_TIME("EBArith::getInterpStencil2D");
  CH_assert(SpaceDim == 2);
  a_sten.clear();
  RealVect faceCentroid = a_ebisBox.centroid(a_face);
  int faceDir = a_face.direction();
  int tanDir = 1 - faceDir;
  Side::LoHiSide tanFaceSide;
  if (faceCentroid[tanDir] > 0.)
    {
      tanFaceSide = Side::Hi;
    }
  else if (faceCentroid[tanDir] < 0.)
    {
      tanFaceSide = Side::Lo;
    }
  else
    {
      tanFaceSide = Side::Lo;
      // MayDay::Error("EBArith: getInterpStencil2D faceCentroid[tanDir] = 0");
    }

  FaceIndex otherFace;
  bool uniqueFace = EBArith::getAdjacentFace(otherFace, a_face,
                                             a_ebisBox, a_domainBox,
                                             tanDir, tanFaceSide);

  //drop order of interpolation to zero if the other face is not unique
  //or if this face or the otherface is on the coarse fine interface
  bool dropOrder = false;
  dropOrder = !uniqueFace;
  if (!dropOrder)
    {
      for (SideIterator sit; sit.ok(); ++sit)
        {
          if (a_coarseFineIVS.contains(   a_face.gridIndex(sit())) ||
             a_coarseFineIVS.contains(otherFace.gridIndex(sit())))
            {
              dropOrder = true;
            }
        }
    }

  if (dropOrder)
    {
      a_sten.add(a_face, 1.0);
    }
  else
    {
      CH_assert(a_face.isDefined());
      CH_assert(otherFace.isDefined());
      Real dist = Abs(faceCentroid[tanDir]);
      Real thisWeight = 1.0 - dist;
      Real otherWeight =  dist;
      a_sten.add(a_face, thisWeight);
      a_sten.add(otherFace, otherWeight);
    }
}
/*****************************/
void
EBArith::
getInterpStencil3D(FaceStencil&         a_sten,
                   const FaceIndex&     a_face,
                   const IntVectSet&    a_coarseFineIVS,
                   const EBISBox&       a_ebisBox,
                   const ProblemDomain& a_domainBox)
{
  CH_TIME("EBArith::getInterpStencil3D");
  CH_assert(SpaceDim == 3);
  a_sten.clear();

  int faceDir = a_face.direction();
  RealVect faceCentroid = a_ebisBox.centroid(a_face);
  //first check to see if this is a genuinely 3d face
  //centroid > 0 in all directions
  bool really3D = true;
  int izerod = 7;
  Real tol = PolyGeom::getTolerance();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (idir != faceDir)
        {
          if (Abs(faceCentroid[idir]) < tol)
            {
              really3D= false;
              izerod = idir;
            }
        }
    }
  if (!really3D)
    {
      //use 2D stencil because we have a centroid that is zero
      //in the given direction
      int tanDir = 3 - izerod - faceDir;
      Side::LoHiSide tanFaceSide;
      if (faceCentroid[tanDir] > 0)
        tanFaceSide = Side::Hi;
      else
        tanFaceSide = Side::Lo;

      bool dropOrder = false;
      FaceIndex otherFace;
      bool uniqueFace = EBArith::getAdjacentFace(otherFace, a_face,
                                                 a_ebisBox, a_domainBox,
                                                 tanDir, tanFaceSide);
      dropOrder = !uniqueFace;

      if ( !dropOrder )
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              if (a_coarseFineIVS.contains(   a_face.gridIndex(sit())) ||
                 a_coarseFineIVS.contains(otherFace.gridIndex(sit())))
                {
                  dropOrder = true;
                }
            }
        }

      if (dropOrder)
        {
          a_sten.add(a_face, 1.0);
        }
      else
        {
          CH_assert(a_face.isDefined());
          CH_assert(otherFace.isDefined());
          Real dist = Abs(faceCentroid[tanDir]);
          Real thisWeight = 1 - dist;
          Real otherWeight =  dist;
          a_sten.add(a_face, thisWeight);
          a_sten.add(otherFace, otherWeight);
        }
    }
  else
    {
      int tanDirs[2];
      Side::LoHiSide tanFaceSides[2];
      {
        int itan = 0;
        for (int idir = 0; idir < 3; idir++)
          {
            if (idir != faceDir)
              {
                tanDirs[itan] = idir;
                if (faceCentroid[tanDirs[itan]] > 0)
                  tanFaceSides[itan] = Side::Hi;
                else
                  tanFaceSides[itan] = Side::Lo;
                itan++;
              }
          }
      }

      bool dropOrder = false;
      FaceIndex face01, face10, face11;
      int ixDir = tanDirs[0];
      int iyDir = tanDirs[1];
      Side::LoHiSide xSide = tanFaceSides[0];
      Side::LoHiSide ySide = tanFaceSides[1];
      //figure out whether to drop order.
      //if not,  get the faces involved
      VolIndex vofLo = a_face.getVoF(Side::Lo);
      VolIndex vofHi = a_face.getVoF(Side::Hi);

      //first get face10 which is the in
      //the  ixdir direction from the input face
      if (!dropOrder)
        {
          bool uniqueFace10
            = EBArith::getAdjacentFace(face10, a_face, a_ebisBox, a_domainBox, ixDir, xSide);
          dropOrder = !uniqueFace10;
          //now get face11 which is the in
          //the  ixdir,iydir direction from the input face (diagonal)
          if (uniqueFace10)
            {
              bool uniqueFace11
                = EBArith::getAdjacentFace(face11, face10, a_ebisBox, a_domainBox, iyDir, ySide);
              dropOrder = !uniqueFace11;
            }
        }
      //the  ixdir,iydir direction from the input face (diagonal)
      //now get face01 which is the in
      //the  iydir direction from the input face
      if (!dropOrder)
        {
          //first get face01 which is the in
          //the  iydir direction from the input face
          bool uniqueFace01
            = EBArith::getAdjacentFace(face01, a_face, a_ebisBox, a_domainBox, iyDir, ySide);
          dropOrder = !uniqueFace01;
          //now get face11 which is the in
          //the  ixdir,iydir direction from the input face (diagonal)
          //compute temp face and see if it is the same as
          //the one computed from the other direction.
          //if not , drop order
          FaceIndex face11temp;
          if (uniqueFace01)
            {
              bool uniqueFace11
                = EBArith::getAdjacentFace(face11temp, face01, a_ebisBox, a_domainBox, ixDir, xSide);
              dropOrder = !uniqueFace11;
              if ((!dropOrder) && !(face11temp == face11) )
                {
                  dropOrder = true;
                }
            }

          //finally if any of the stencil faces are in the coarse-fine interface, drop order
          if (!dropOrder)
            {
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  if ((a_coarseFineIVS.contains(a_face.gridIndex(sit()))) ||
                     (a_coarseFineIVS.contains(face01.gridIndex(sit()))) ||
                     (a_coarseFineIVS.contains(face10.gridIndex(sit()))) ||
                     (a_coarseFineIVS.contains(face11.gridIndex(sit()))))
                    {
                      dropOrder = true;
                    }
                }
            }
        }
      ///////////////////////
      //construct the stencils
      ///////////////////////
      if (dropOrder)
        {
          a_sten.add(a_face, 1.0);
        }
      else
        {
          FaceIndex face00 = a_face;
          Real xbar = Abs(faceCentroid[ixDir]);
          Real ybar = Abs(faceCentroid[iyDir]);
          Real f00coef = 1.0 - xbar - ybar + xbar*ybar;
          Real f10coef = xbar - xbar*ybar;
          Real f01coef = ybar - xbar*ybar;
          Real f11coef = xbar*ybar;

          for (SideIterator sit; sit.ok(); ++sit)
            {
              if (!face00.isBoundary())
                CH_assert(face00.cellIndex(sit()) >= 0);
              if (!face01.isBoundary())
                CH_assert(face01.cellIndex(sit()) >= 0);
              if (!face10.isBoundary())
                CH_assert(face10.cellIndex(sit()) >= 0);
              if (!face11.isBoundary())
                CH_assert(face11.cellIndex(sit()) >= 0);
            }
          CH_assert(face00.isDefined());
          CH_assert(face01.isDefined());
          CH_assert(face10.isDefined());
          CH_assert(face11.isDefined());

          a_sten.add(face00, f00coef);
          a_sten.add(face01, f01coef);
          a_sten.add(face10, f10coef);
          a_sten.add(face11, f11coef);
        }
    }
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
