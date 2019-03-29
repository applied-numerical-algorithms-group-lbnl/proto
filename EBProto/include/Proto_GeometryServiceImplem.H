namespace Proto
{
/*********************************************/
  template <int order>
  GeometryService<order>::
  GeometryService(const shared_ptr<BaseIF>  &      a_baseIF,
               const RealVect            &      a_origin,
               const Real                &      a_dx,
               const Box                 &      a_domain)
  {
    CH_TIME("GeometryService::GeometryService");
    m_baseIF  = a_baseIF;

    m_domain  = a_domain;
    m_origin  = a_origin;
            
    m_order   = 0;
    m_degreeP = CH_EBIS_ORDER + 1;
  }


/*********************************************/
  template <int order>
  void 
  GeometryService<order>::
  fillGraph(HostBoxData<int>&    a_regIrregCovered,
            vector<inode>&       a_nodes,
            const Box&           a_validRegion,
            const Box&           a_ghostRegion,
            const Box&           a_domain,
            const RealVect&      a_origin,
            const Real&          a_dx,
            const DataIndex&     a_di ) const
  {
    CH_TIME("GeometryService::fillGraph");

    CH_assert(a_domain.contains(a_ghostRegion));

    PointSet ivsirreg = PointSet(DensePointSet(a_ghostRegion, false));
    Real thrshd = 1.0e-10;
    {
      CH_TIME("boxiterator loop");
      for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit)
      {
        const IntVect iv =bit();
        RvgDim cellCenter;
        for (int idir = 0;idir < SpaceDim; ++idir)
        {
          cellCenter[idir] = a_dx*(iv[idir] +0.5) + a_origin[idir];
        }

        //member data: sign(chosen from -1,0,1) of each vertex,
        //location of each edge intersection, cellCenter,normal and gradNormal

        IndexTM<Real, SpaceDim> vectdx;
        vectdx.setAll(a_dx);

        int degreeP = m_order + m_degreeP;
        IFData<SpaceDim> edgeData(*m_baseIF, vectdx, cellCenter,  degreeP);


        //create a CutCellMoment object, in order to detect whether any face coincides with the interface
        CutCellMoments <SpaceDim> cutCell(edgeData);
        if (cutCell.isCovered())
        {
          //set covered cells to -1
          a_regIrregCovered(iv, 0) = -1;
        }
        else if (cutCell.isRegular())
        {
          //set regular cells to 1
          a_regIrregCovered(iv, 0) =  1;
        }
        else
        {
          //set irregular cells to 0
          //irregular if any face coincides with interface and edgeData.m_allVerticesIn = true
          a_regIrregCovered(iv, 0) =  0;
          if (a_validRegion.contains(iv))
          {
            ivsirreg |= iv;
          }
        }
      }
    }
    //if a regular is next to a  covered, change to irregular with correct arcs and so on.
    for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit)
    {
      const IntVect iv =bit();

      if(a_regIrregCovered(iv, 0) == -1)
      {
        fixRegularCellsNextToCovered(a_nodes, a_regIrregCovered, a_validRegion, a_domain.domainBox(), iv, a_dx);
      }
    }

    //now loop through irregular cells and make nodes for each  one.
    PointSet ivsdrop; //volumes too small to keep

    for (PointSetIterator ivsit(ivsirreg); ivsit.ok(); ++ivsit)
    {
      IntVect iv = ivsit();
      CH_TIME("fillGraph::endOfirregularCellLoop");
      inode newNode;

      fillNewNode(newNode,
                  ivsirreg,
                  a_domain,
                  a_origin,
                  a_dx,
                  ivsit());

      Real volFrac = newNode.m_volFrac;
      if (volFrac < thrshd)
      {
        ivsdrop |= iv;
        a_regIrregCovered(iv, 0) = -1;
      }//CP record these nodes to be removed
      else
      {
        //bool fixMoments = true;
        //checkNodeMoments(newNode, a_dx, fixMoments, s_relativeTol);      
        a_nodes.push_back(newNode);
      }
    } //end loop over cells in the box
    // CP: fix sweep that removes cells with volFrac less than a certain threshold
    for(PointSetIterator ivsit(ivsdrop); ivsit.ok(); ++ivsit)
    {
      const IntVect& iv = ivsit();
  
      for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
      {
        for (SideIterator sit; sit.ok(); ++sit)
        {
          int isign = sign(sit());
          IntVect otherIV = iv + isign*BASISV(faceDir);
          if (a_validRegion.contains(otherIV))
          {
            if (a_regIrregCovered(otherIV,0) == 0)
            {
              // i am in the case where the other cell
              // is also irregular.   I just made a previously
              // irregular cell covered so I have to check to
              // see if it had any faces pointed this way.
              int inode = -1;
              bool found = false;
              for (int ivec = 0; ivec < a_nodes.size() && ! found; ivec++)
              {
                if (a_nodes[ivec].m_cell == otherIV)
                {
                  inode = ivec;
                  found = true;
                }
              }
              if (!found && a_validRegion.contains(otherIV))
              {
                MayDay::Abort("something wrong in our logic2");
              }
              if (found)
              {
                int arcindex = a_nodes[inode].index(faceDir, flip(sit()));
                a_nodes[inode].m_arc[         arcindex].resize(0);
                a_nodes[inode].m_areaFrac[    arcindex].resize(0);
                a_nodes[inode].m_faceCentroid[arcindex].resize(0);
                a_nodes[inode].m_faceMoments [arcindex].setToZero();
              }
            }
          }//valid region
        }//sit
      }//facedir
        
      //also need to fix regular cells next to new covered cell
      fixRegularCellsNextToCovered(a_nodes, a_regIrregCovered, a_validRegion, a_domain.domainBox(), iv, a_dx);

    }//ivsdrop
  }
/*************/
  template <int order>
  void
  GeometryService<order>::
  fixRegularCellsNextToCovered(vector<inode>       & a_nodes, 
                               HostBoxData<int>    & a_regIrregCovered,
                               const Box           & a_validRegion,
                               const Box           & a_domain,
                               const IntVect       & a_iv,
                               const Real          & a_dx) const

  {
    Box grownBox(a_iv, a_iv);
    grownBox.grow(1);
    grownBox  &= a_domain;
    PointSet ivstocheck(grownBox);
    ivstocheck -= a_iv;
    Box ghostRegion = a_regIrregCovered.box();
    //first check neighbors in each direction.  
    //If any of these are regular, they are replaced 
    //by irregular cells with a boundary face facing the covered cell.
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      for(SideIterator sit; sit.ok(); ++sit)
      {
        int ishift = sign(sit());
        IntVect ivshift = a_iv + ishift*BASISV(idir);
        ivstocheck -= ivshift;
        int bfvalshift = -1;
        if(ghostRegion.contains(ivshift))
        {
          bfvalshift = a_regIrregCovered(ivshift, 0);
        }
        if(bfvalshift  == 1)
        {
          a_regIrregCovered(ivshift, 0) =  0;

          if(a_validRegion.contains(ivshift))
          {
            inode newNode;
            getFullNodeWithCoveredFace(newNode, 
                                       a_regIrregCovered,
                                       ivshift, 
                                       a_dx,
                                       a_domain);
            a_nodes.push_back(newNode);
          }

        }
      }
    }
    //next we loop through the remaining cells (corner cells in 2d, corner and edge cells in 3D)
    //if any of these are regular, we change them to irregular 
    for(PointSetIterator ivsit(ivstocheck); ivsit.ok(); ++ivsit)
    {
      const IntVect& iv = ivsit();
      if(ghostRegion.contains(iv))
      {
        if(a_regIrregCovered(iv, 0) == 1)
        {
          a_regIrregCovered(iv, 0) = 0;
          inode newNode;
          newNode.makeRegular(iv, a_dx, a_domain);
          a_nodes.push_back(newNode);
        }
      }
    }
  }
/**********************************************/
  template <int order>
  void
  GeometryService<order>::
  getFullNodeWithCoveredFace(inode                   & a_newNode, 
                             const HostBoxData<int>  & a_regIrregCovered,
                             const IntVect           & a_iv,
                             const Real              & a_dx,
                             const Box               & a_domain) const
  {

    a_newNode.m_cell          = a_iv;
    a_newNode.m_volFrac       = 1.0;
    a_newNode.m_cellIndex     = 0;
    a_newNode.m_volCentroid   = RealVect::Zero;
    //set regular cell values then fix up
    a_newNode.m_bndryCentroid = RealVect::Zero;
    //set all moments to regular and then zero out the appropriate face mometns
    a_newNode.setMomentsToRegular(a_dx);
    int coveredDir;
    Side::LoHiSide coveredSide;
    bool found = false;

    for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
    {
      for(SideIterator sit; sit.ok(); ++sit)
      {
        int ishift = sign(sit());
        IntVect ivshift = a_iv + ishift*BASISV(faceDir);
        vector<int> arc;
        vector<Real> areaFrac;
        vector<RealVect> faceCentroid;
        if(!a_domain.contains(ivshift))
        {
          // boundary arcs always -1
          arc.resize(1,-1);
          areaFrac.resize(1, 1.0);
          faceCentroid.resize(1, RealVect::Zero);
        }
        else if (a_regIrregCovered(ivshift, 0) >= 0)
        {
          //irregular cell or regular cell
          //compute vof internals returns something special if 
          //connected to a regular cell but EBGraph treats both the  same.
          //it just  knows that the cell index of a regular cell is 0
          arc.resize(1,0);
          areaFrac.resize(1, 1.0);
          faceCentroid.resize(1, RealVect::Zero);
        }
        else if (a_regIrregCovered(ivshift, 0) < 0)
        {
          found = true;
          coveredDir= faceDir;
          coveredSide = sit();
          // covered face!
          arc.resize(0);
          areaFrac.resize(0);
          faceCentroid.resize(0);
        }
        else
        {
          MayDay::Error("logic error");
        }
          
        int nodeInd = a_newNode.index(faceDir, sit());
        a_newNode.m_arc[nodeInd]          = arc;
        a_newNode.m_areaFrac[nodeInd]     = areaFrac;
        a_newNode.m_faceCentroid[nodeInd] = faceCentroid;
        a_newNode.m_faceMoments[nodeInd].setToZero();
      }
    }
    //fix boundary centroid
    if(found)
    {
      int centsign = sign(coveredSide);
      a_newNode.m_bndryCentroid[coveredDir] =  centsign*0.5;
    }
  }

  /******/
  template <int order>
  void 
  GeometryService<order>::
  fillNewNode(inode               &     a_node,
              const PointSet      &     a_ivsIrreg,
              const Box           &     a_domain,
              const RealVect      &     a_origin,
              const Real          &     a_dx,
              const IntVect       &     a_iv) const
  {
    CH_TIME("GeometryService::ComputeVofInternals");

    //for each CutCellMoments<dim>, we record the cell Center
    //(in physical coordinates at the global dimension)
    RvgDim cellCenter;
    for (int idir = 0;idir < SpaceDim; ++idir)
    {
      cellCenter[idir] = a_dx*(a_iv[idir] +0.5) + a_origin[idir];
    }

    // member data: sign (chosen from -1,0,1) of each vertex,
    // location of each edge intersection, cellCenter, normal and gradNormal
    //int degreeP = m_degreeP;
    int degreeP = m_order + m_degreeP;
    int orderP  = 0;
    IndexTM<Real, SpaceDim> vectdx;
    vectdx.setAll(a_dx);
    IFData<SpaceDim> edgeData(*m_baseIF,vectdx,cellCenter,  degreeP);

    //construct data holders for all moments
    MinimalCCCM<SpaceDim> computeThisVof(edgeData);

    //compute the moments and save answers in thisVof
    computeThisVof.computeMoments(orderP,degreeP);

    CutCellMoments<SpaceDim> thisVof = computeThisVof.m_cutCellMoments;

    vector<int> loArc[SpaceDim];
    vector<int> hiArc[SpaceDim];

    a_node.m_cell        = a_iv;
    a_node.m_cellIndex = 0;

    //now the arcs (0 is lo, 1 is high)
    fillArc(loArc, thisVof, 0, a_ivsIrreg, a_iv);
    fillArc(hiArc, thisVof, 1, a_ivsIrreg, a_iv);
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      int indexlo = inode::index(idir, Side::Lo);
      int indexhi = inode::index(idir, Side::Hi);
      a_node.m_arc[indexlo] = loArc[idir];
      a_node.m_arc[indexhi] = hiArc[idir];
    }

    a_node.m_volumeMoments = thisVof.m_moments;          
    a_node.m_EBMoments     = thisVof.m_EBmoments;

    IndexTM<Real,CH_SPACEDIM> point;
    for (int idir = 0; idir < SpaceDim; ++idir)
    {
      Real cellCent = (a_iv[idir]+ 0.5)*a_dx;
      point[idir] = cellCent;
    }  

    NormalDerivativeNew<SpaceDim> normalDerivative;

    int maxOrder = CH_EBIS_ORDER;

    IndexTM<Real, SpaceDim> itmpoint;
    EBArith::convertToITM(itmpoint, point);
    IFSlicer<SpaceDim> ifSlicer(*m_baseIF);
    NormalDerivativeNew<SpaceDim>::NormalDerivativeMap ndMap = 
      normalDerivative.calculateAll(maxOrder,
                                    itmpoint,
                                    &ifSlicer);

    for (int idir = 0; idir < SpaceDim; ++idir)
    {
      MomItSpaceDim momit;
      for(momit.reset(); momit.ok(); ++momit)
      {
        Real derivVal = ndMap[momit()][idir];
        a_node.m_normalPartialDeriv[idir][momit()] = derivVal;
      }
    }

  
    Iv2 bdId;
    for(SideIterator sit; sit.ok(); ++sit)
    {
      int hilo = 0;
      if(sit() == Side::Hi) hilo = 1;
      
      for (int idir = 0; idir < SpaceDim; ++idir)
      {
        int iindex = inode::index(idir, sit());
        CH_assert((a_node.m_arc[iindex].size() == 1) || (a_node.m_arc[iindex].size() == 0));
        if(a_node.m_arc[iindex].size() == 1)
        {
          bdId[BDID_HILO] = hilo;
          bdId[BDID_DIR] = idir;
          const CutCellMoments<SpaceDim-1>& bdccm = thisVof.getBdCutCellMoments(bdId);
          a_node.m_faceMoments[iindex] = bdccm.m_moments;
        }
        else
        {
          a_node.m_faceMoments[iindex].setToZero();
        }
      }
    }
    //this fills volfrac, areafrac, centroids...
    a_node.setNormalizedStuff(a_dx);
  }

  //records connectivity between vofs
  template <int order>
  void 
  GeometryService<order>::
  fillArc(vector<int>                          a_arc[SpaceDim],
          CutCellMoments<SpaceDim>       &     a_cutCellMoments,
          const int                      &     a_hilo,
          const PointSet                 &     a_ivsIrreg,
          const IntVect                  &     a_curriv) const
  {
    Iv2 bdId;
    //a_hilo is 0 or 1
    bdId[BDID_HILO] = a_hilo;
    for (int idir = 0; idir < SpaceDim; ++idir)
    {
      bdId[BDID_DIR] = idir;
      bool covered =  a_cutCellMoments.getBdCutCellMoments(bdId).isCovered();

      if (covered)
      {
        a_arc[idir].resize(0);
      }
      else
      {
        a_arc[idir].resize(1);

        //otherIV is the iv in the idir direction on the a_hilo side
        IntVect otherIV = a_curriv;
        otherIV[idir] += (a_hilo*2) - 1;

        if (m_domain.contains(otherIV))
        {
          int otherCellIndex;
          if (a_ivsIrreg.contains(otherIV))
          {
            otherCellIndex = 0;
          }
          else
          {
            //arc to regular cell
            otherCellIndex = -2;
          }
          a_arc[idir][0]=otherCellIndex;
        }
        else if (!m_domain.contains(otherIV))
        {
          //boundary arcs always -1
          a_arc[idir][0] = -1;
        }
      }
    }
  }
}
