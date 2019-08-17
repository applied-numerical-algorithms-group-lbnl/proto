#ifndef __Proto_EBGraphIMPLEM__
#define __Proto_EBGraphIMPLEM__

namespace Proto
{

  template <CENTERING cent> 
  bool
  EBGraph::
  getAdjacentFace(EBIndex<cent>       & a_adjFace, 
                  const EBIndex<cent> & a_startFace,
                  int a_tanDir, Side::LoHiSide a_side) const
  {
    PR_assert(a_tanDir != a_startFace.direction());
    bool hasFace = true;
    Point hiiv = a_startFace.m_pt + sign(a_side)*Point::Basis(a_tanDir, 1);
    Point loiv = hiiv - Point::Basis(a_startFace.direction(), 1);
    if((!m_domain.contains(loiv))  || (!m_domain.contains(hiiv)))
    {
      hasFace = false;
      return hasFace;
    }
    EBIndex<CELL> loVoF = a_startFace.getVoF(Side::Lo);
    EBIndex<CELL> hiVoF = a_startFace.getVoF(Side::Hi);

    if(a_tanDir == 0)
    {
      vector<EBIndex<XFACE> > faceslo = getXFaces(loVoF,  a_side);
      vector<EBIndex<XFACE> > faceshi = getXFaces(hiVoF,  a_side);
      EBIndex<CELL> voflo, vofhi;
      bool found = false;
      for(int ilo = 0; ilo < faceslo.size(); ilo++)
      {
        voflo = faceslo[ilo].getVoF(a_side);
        for(int ihi = 0; ihi < faceshi.size(); ihi++)
        {
          vofhi = faceshi[ihi].getVoF(a_side);
          if(isConnected(voflo, vofhi))
          {
            found = true;
            break;
          }
        }
      }
      if(!found)
      {
        hasFace = false;
        return hasFace;
      }
      else
      {
        a_adjFace = EBIndex<cent>(voflo, vofhi);
      }
    }
    else if(a_tanDir == 1)
    {
      vector<EBIndex<YFACE> > faceslo = getYFaces(loVoF,    a_side);
      vector<EBIndex<YFACE> > faceshi = getYFaces(hiVoF,    a_side);
      EBIndex<CELL> voflo, vofhi;
      bool found = false;
      for(int ilo = 0; ilo < faceslo.size(); ilo++)
      {
        voflo = faceslo[ilo].getVoF(a_side);
        for(int ihi = 0; ihi < faceshi.size(); ihi++)
        {
          vofhi = faceshi[ihi].getVoF(a_side);
          if(isConnected(voflo, vofhi))
          {
            found = true;
            break;
          }
        }
      }
      if(!found)
      {
        hasFace = false;
        return hasFace;
      }
      else
      {
        a_adjFace = EBIndex<cent>(voflo, vofhi);
      }
    }
#if DIM==3
    else if(a_tanDir == 2)
    {
      vector<EBIndex<ZFACE> > faceslo = getZFaces(loVoF,  a_side);
      vector<EBIndex<ZFACE> > faceshi = getZFaces(hiVoF,  a_side);
      EBIndex<CELL> voflo, vofhi;
      bool found = false;
      for(int ilo = 0; ilo < faceslo.size(); ilo++)
      {
        voflo = faceslo[ilo].getVoF(a_side);
        for(int ihi = 0; ihi < faceshi.size(); ihi++)
        {
          vofhi = faceshi[ihi].getVoF(a_side);
          if(isConnected(voflo, vofhi))
          {
            found = true;
            break;
          }
        }
      }
      if(!found)
      {
        hasFace = false;
        return hasFace;
      }
      else
      {
        a_adjFace = EBIndex<cent>(voflo, vofhi);
      }
    }
#endif
    else
    {
      PR_error("bogus tandir");
    }
    return hasFace;
  }

  void
  CellNode::
  getFaceInfo(int   & a_vecsize,
              Point & a_mypt,
              int   & a_myVoFID,
              int   & a_otherVoFID,
              bool  & a_isBoundary,
              const int            & a_dir,
              const  EBIndex<CELL> & a_vof,
              const Side::LoHiSide & a_sd,
              const Box            & a_domain) const

  {
    if(isCovered())
    {
      a_vecsize = 0;
    }
    else
    {
      a_vecsize = 1;
      Point inputIV    = a_vof.m_pt;
      Point otherIV    = a_vof.m_pt +sign(a_sd)*Point::Basis(a_dir, 1);

      int inputIndex = a_vof.m_vofIDMe;
      int otherIndex = 0;
      a_isBoundary = false;
      if(!a_domain.contains(otherIV)) 
      {
        otherIndex = -1;
        a_isBoundary = true;
      }

      if(isIrregular())
      {
        vector<IrrGraphNode>& nodeVec = *m_cellList;
        //again, we are dealing with single-valued stuff at this level for now
        //multivalued can only happen by coarsening
        PR_assert(nodeVec.size() == 1);
        int index = ebp_index(a_dir, a_sd);
        vector<int>& arcvec = nodeVec[0].m_arc[index];
        if(arcvec.size() == 0)
        {
          a_vecsize = 0;
        }
        //there would be more complicated stuff here if multivalued was possible.
        //as it is, -1 or 0 are the only possibilities for the other index
      }

      if(a_sd == Side::Hi)
      {
        a_mypt       = otherIV;
        a_myVoFID    = otherIndex;
        a_otherVoFID = inputIndex;
      }
      else 
      {
        a_mypt       = inputIV;
        a_myVoFID    = inputIndex;
        a_otherVoFID = otherIndex;
      }
    }
  }
  /*******************************/
  vector<EBIndex<XFACE> >
  CellNode::
  getXFaces(const  EBIndex<CELL>        & a_vof,
            const Side::LoHiSide & a_sd,
            const Box            & a_domain) const
  {
    PR_assert(a_domain.contains(a_vof.m_pt));
    int    vecsize;
    Point  pt;
    int    myVoFID;
    int    otherVoFID;
    bool   isBoundary;
    getFaceInfo(vecsize, pt, myVoFID, otherVoFID, isBoundary, 0, a_vof, a_sd, a_domain);
    
    vector<EBIndex<XFACE> > retval(vecsize);
    //single valued stuff at this level
    PR_assert(vecsize <= 1);
    if(vecsize > 0)
    {
      retval[0].m_pt         = pt;
      retval[0].m_vofIDMe    = myVoFID;
      retval[0].m_vofIDLo    = otherVoFID;
      retval[0].m_isBoundary = isBoundary;
    }
    return retval;
  }

  /*******************************/
  vector<EBIndex<YFACE> >
  CellNode::
  getYFaces(const  EBIndex<CELL>        & a_vof,
            const Side::LoHiSide & a_sd,
            const Box            & a_domain) const
  {
    PR_assert(a_domain.contains(a_vof.m_pt));
    int    vecsize;
    Point  pt;
    int    myVoFID;
    int    otherVoFID;
    bool   isBoundary;
    getFaceInfo(vecsize, pt, myVoFID, otherVoFID, isBoundary, 1, a_vof, a_sd, a_domain);
    
    vector<EBIndex<YFACE> > retval(vecsize);
    //single valued stuff at this level
    PR_assert(vecsize <= 1);
    if(vecsize > 0)
    {
      retval[0].m_pt         = pt;
      retval[0].m_vofIDMe    = myVoFID;
      retval[0].m_vofIDLo    = otherVoFID;
      retval[0].m_isBoundary = isBoundary;
    }
    return retval;
  }

  /*******************************/
#if DIM==3
  vector<EBIndex<ZFACE> >
  CellNode::
  getZFaces(const  EBIndex<CELL> & a_vof,
            const Side::LoHiSide & a_sd,
            const Box            & a_domain) const
  {
    PR_assert(a_domain.contains(a_vof.m_pt));
    int    vecsize;
    Point  pt;
    int    myVoFID;
    int    otherVoFID;
    bool   isBoundary;
    getFaceInfo(vecsize, pt, myVoFID, otherVoFID, isBoundary, 2, a_vof, a_sd, a_domain);
    
    vector<EBIndex<ZFACE> > retval(vecsize);
    //single valued stuff at this level
    PR_assert(vecsize <= 1);
    if(vecsize > 0)
    {
      retval[0].m_pt         = pt;
      retval[0].m_vofIDMe    = myVoFID;
      retval[0].m_vofIDLo    = otherVoFID;
      retval[0].m_isBoundary = isBoundary;
    }
    return retval;
  }
#endif
  vector< EBIndex<CELL> >
  CellNode::
  getVoFs(const Point& a_pt) const
  {
    vector<EBIndex<CELL> > retval;
    if(isCovered())
    {
      //nothing to return
    }
    else
    {
     if(isRegular())
     {
      retval.resize(1);
     }
     else
     {
      retval.resize(m_cellList->size());
     }
     for(int ivof = 0; ivof < retval.size(); ivof++)
     {
       retval[0].m_pt = a_pt;
       retval[0].m_vofIDMe = ivof;
     }
    }
    return retval;
  }


  vector< EBIndex<BOUNDARY> >
  CellNode::
  getCutFaces(const Point& a_pt) const
  {
    vector<EBIndex<BOUNDARY> > retval;
    //only something to do if irregular
    if(isIrregular())
    {
      retval.resize(m_cellList->size());
      for(int ivof = 0; ivof < retval.size(); ivof++)
      {
        retval[0].m_pt = a_pt;
        retval[0].m_vofIDMe = ivof;
      }
    }
    return retval;
  }
    //this deals with the fine to coarse thing
  void
  EBGraph::
  coarsenVoFs(EBGraph       &  a_fineEBGraph,
              const Box     &  a_vofRegionCoar)
  {
    PROTO_ASSERT(false, "not implemented");
  }
  ///
  void
  EBGraph::
  coarsenFaces(const EBGraph  &  a_fineEBGraph,
               const Box      &  a_faceRegionCoar)
  {
    PROTO_ASSERT(false, "not implemented");

  }
  /*******************************/
  void EBGraph::defineAsAllRegular(const Box& a_region, const Box & a_domain)
  {
    m_isDefined = true;
    m_domain = a_domain;
    m_region = a_region;
    m_tag = AllRegular;
    m_internals = shared_ptr<localData>();
  }

/*******************************/
  void EBGraph::defineAsAllCovered(const Box& a_region, const Box & a_domain)
  {

    m_isDefined = true;
    m_domain = a_domain;
    m_region = a_region;
    m_tag = AllCovered;
    m_internals = shared_ptr<localData>();
  }


/*******************************/
  void EBGraph::defineFromPrimitives(const HostBoxData<int> & a_regIrregCovered,
                                     const vector<IrregGraph>& a_irregGraph,
                                     const Box              & a_region,
                                     const Box              & a_domain,
                                     bool singleValuedCellsOnly)
  {
    PR_TIME("EBGraph::buildGraph");
    //this eliminates all the stuff in Chombo's version referring back to the original vector index and searching though them.
    //this simplifies the construction process enormously.
    //multivalued cells can still happen via coarsening.
    PROTO_ASSERT(singleValuedCellsOnly, "EBGraph::defineFromPrimitives is currently not general enough to handle multivalued inputs");
    m_isDefined = true;
    m_domain = a_domain;
    m_region = a_region;
    m_tag = HasIrregular;
    m_internals = shared_ptr<localData>(new localData());

    m_internals->m_graph.define(m_region);

    // set regular and covered cells
    int num_irreg = 0;
    for (BoxIterator bit = m_region.begin(); bit != m_region.end(); ++bit)
    {
      const Point& iv = *bit;
      if (a_regIrregCovered(iv, 0) == 1) //regular cell
      {
        m_internals->m_graph(iv, 0).defineAsRegular();
      }
      else if (a_regIrregCovered(iv, 0) == -1) // covered cell
      {
        m_internals->m_graph(iv, 0).defineAsCovered();
      }
      else if (a_regIrregCovered(iv, 0) == 0)
      {
        m_internals->m_graph(iv, 0).defineAsIrregular();
        num_irreg++;
      }
      else 
      {
        PR_error("invalid flag");
      }
    }
    //it can be greater if there are multivalued cells
    PROTO_ASSERT(a_irregGraph.size() >= num_irreg, "mismatch in ebgraph::define inputs");

    // now for irregular cells
    for (int ivecIrreg = 0; ivecIrreg < a_irregGraph.size(); ivecIrreg++)
    {
      const IrregGraph & inputNode = a_irregGraph[ivecIrreg];
      const Point     & iv = inputNode.m_cell;
      m_internals->m_irregIVS |= iv;

      CellNode::IrrGraphNode currNode;
      for (int idir = 0; idir < DIM; idir++)
      {
        for (SideIterator sit; sit.ok(); ++sit)
        {
          int irregIndex = ebp_index(idir, sit());
          const vector<int>& irregArcs = inputNode.m_arc[irregIndex];
          //comes from the singleValuedOnly thing
          PR_assert(irregArcs.size() <= 1);
          currNode.m_arc[irregIndex] = irregArcs;
        }
      }
      m_internals->m_graph(iv, 0).m_cellList->push_back(currNode);
    }
  }


/*******************************/
  vector< EBIndex<CELL> > EBGraph::getVoFs(const Point& a_iv) const
  {
    PR_assert(isDefined());

    vector< EBIndex<CELL> > retvec;
    if (m_tag == AllRegular)
    {
      retvec.push_back( EBIndex<CELL> (a_iv, 0));
    }
    else if (m_tag == AllCovered)
    {
      //return an empty vector
    }
    else if (m_tag == HasIrregular)
    {
      const CellNode& node = m_internals->m_graph(a_iv, 0);
      retvec = node.getVoFs(a_iv);
    }
    return retvec;
  }

/*******************************/
  vector< EBIndex<BOUNDARY> > EBGraph::getCutFaces(const Point& a_iv) const
  {
    PR_assert(isDefined());

    vector< EBIndex<BOUNDARY> > retvec;
    if (m_tag == HasIrregular) //otherwise return an empty vector
    {
      const CellNode& node = m_internals->m_graph(a_iv, 0);
      if(node.isIrregular())
      {
        retvec = node.getCutFaces(a_iv);
      }
    }
    return retvec;
  }

/*******************************/
  bool EBGraph::isRegular(const Point& a_iv) const
  {
    PR_assert(isDefined());
    bool retval;
    if (m_tag == AllRegular)
    {
      retval = true;
    }
    else if (m_tag == AllCovered)
    {
      retval = false;
    }
    else if (m_tag == HasIrregular)
    {
      const CellNode& node = m_internals->m_graph(a_iv, 0);
      retval = node.isRegular();
    }
    else
    {
      retval = false;
      PR_error("EBGraph::isRegular:Bogus Tag");
    }
    return retval;
  }


/*******************************/
  bool EBGraph::isIrregular(const Point& a_iv) const
  {
    PR_assert(isDefined());

    bool retval;
    if (m_tag == AllRegular)
    {
      retval = false;
    }
    else if (m_tag == AllCovered)
    {
      retval = false;
    }
    else if (m_tag == HasIrregular)
    {
      const CellNode& node = m_internals->m_graph(a_iv, 0);
      retval = node.isIrregular();
    }
    else
    {
      retval = false;
      PR_error("EBGraph::isIrregular:Bogus Tag");
    }
    return retval;
  }


/*******************************/
  bool EBGraph::isCovered(const Point& a_iv) const
  {
    PR_assert(isDefined());

    bool retval;
    if (m_tag == AllRegular)
    {
      retval = false;
    }
    else if (m_tag == AllCovered)
    {
      retval = true;
    }
    else if (m_tag == HasIrregular)
    {
      const CellNode& node = m_internals->m_graph(a_iv, 0);
      retval = node.isCovered();
    }
    else
    {
      retval = false;
      PR_error("EBGraph::isIrregular:Bogus Tag");
    }
    return retval;
  }

/*******************************/
  bool EBGraph::isCovered(const Box& a_box) const
  {
    PR_assert(isDefined());

    bool retval;
    if (m_tag == AllRegular)
    {
      retval = false;
    }
    else if (m_tag == AllCovered)
    {
      retval = true;
    }
    else if (m_tag == HasIrregular)
    {

      retval = true;
      BoxIterator bit(a_box);
      for (BoxIterator bit = a_box.begin(); bit != a_box.end(); ++bit)
      {
        if (!isCovered(*bit))
        {
          retval = false;
          break;
        }
      }
    }
    else
    {
      retval = false;
      PR_error("EBGraph::isCovered:Bogus Tag");
    }
    return retval;
  }

/*******************************/
  bool EBGraph::isRegular(const Box& a_box) const
  {
    PR_assert(isDefined());

    bool retval;
    if (m_tag == AllRegular)
    {
      retval = true;
    }
    else if (m_tag == AllCovered)
    {
      retval = false;
    }
    else if (m_tag == HasIrregular)
    {
      retval = true;
      for (BoxIterator bit = a_box.begin(); bit != a_box.end(); ++bit)
      {
        if (!isRegular(*bit))
        {
          retval = false;
          break;
        }
      }
    }
    else
    {
      retval = false;
      PR_error("EBGraph::isRegular:Bogus Tag");
    }
    return retval;
  }
/*******************************/
  vector< EBIndex<XFACE> > 
  EBGraph::
  getXFaces(const  EBIndex<CELL>        & a_vof,
            const Side::LoHiSide & a_sd) const
  {
    PR_assert(isDefined());

    vector< EBIndex<XFACE> > retval;
    if (m_tag == AllRegular)
    {
      int idir = 0;
      Point otherIV = a_vof.m_pt
        + sign(a_sd)*Point::Basis(idir);
      int otherCellInd = 0;
      if (!m_domain.contains(otherIV))
      {
        otherCellInd = -1;
      }
      retval.resize(1);
      if(a_sd == Side::Lo)
      {
        retval[0].m_pt      = a_vof.m_pt;
        retval[0].m_vofIDMe = a_vof.m_vofIDMe;
        retval[0].m_vofIDLo = otherCellInd;
      }
      else
      {
        retval[0].m_pt      = otherIV;
        retval[0].m_vofIDMe = otherCellInd;
        retval[0].m_vofIDLo = a_vof.m_vofIDMe;
      }
      retval[0].m_isBoundary = (!m_domain.contains(otherIV));
    }
    else if (m_tag == AllCovered)
    {
      //return empty vector
    }
    else if (m_tag == HasIrregular)
    {
      const Point& iv = a_vof.m_pt;
      const CellNode& node = m_internals->m_graph(iv, 0);
      retval = node.getXFaces(a_vof, a_sd, m_domain);
    }

    return retval;
  }


/*******************************/
  vector< EBIndex<YFACE> > 
  EBGraph::
  getYFaces(const  EBIndex<CELL>        & a_vof,
            const Side::LoHiSide & a_sd) const
  {
    PR_assert(isDefined());

    vector< EBIndex<YFACE> > retval;
    if (m_tag == AllRegular)
    {
      int idir = 1;
      Point otherIV = a_vof.m_pt
        + sign(a_sd)*Point::Basis(idir);
      int otherCellInd = 0;
      if (!m_domain.contains(otherIV))
      {
        otherCellInd = -1;
      }
      retval.resize(1);
      if(a_sd == Side::Lo)
      {
        retval[0].m_pt      = a_vof.m_pt;
        retval[0].m_vofIDMe = a_vof.m_vofIDMe;
        retval[0].m_vofIDLo = otherCellInd;
      }
      else
      {
        retval[0].m_pt      = otherIV;
        retval[0].m_vofIDMe = otherCellInd;
        retval[0].m_vofIDLo = a_vof.m_vofIDMe;
      }
      retval[0].m_isBoundary = (!m_domain.contains(otherIV));
    }
    else if (m_tag == AllCovered)
    {
      //return empty vector
    }
    else if (m_tag == HasIrregular)
    {
      const Point& iv = a_vof.m_pt;
      const CellNode& node = m_internals->m_graph(iv, 0);
      retval = node.getYFaces(a_vof, a_sd, m_domain);
    }

    return retval;
  }
#if DIM==3
/*******************************/
  vector< EBIndex<ZFACE> > 
  EBGraph::
  getZFaces(const  EBIndex<CELL>        & a_vof,
            const Side::LoHiSide & a_sd) const
  {
    PR_assert(isDefined());

    vector< EBIndex<ZFACE> > retval;
    if (m_tag == AllRegular)
    {
      int idir = 2;
      Point otherIV = a_vof.m_pt
        + sign(a_sd)*Point::Basis(idir, 1);
      int otherCellInd = 0;
      if (!m_domain.contains(otherIV))
      {
        otherCellInd = -1;
      }
      retval.resize(1);
      if(a_sd == Side::Lo)
      {
        retval[0].m_pt      = a_vof.m_pt;
        retval[0].m_vofIDMe = a_vof.m_vofIDMe;
        retval[0].m_vofIDLo = otherCellInd;
      }
      else
      {
        retval[0].m_pt      = otherIV;
        retval[0].m_vofIDMe = otherCellInd;
        retval[0].m_vofIDLo = a_vof.m_vofIDMe;
      }
      retval[0].m_isBoundary = (!m_domain.contains(otherIV));
    }
    else if (m_tag == AllCovered)
    {
      //return empty vector
    }
    else if (m_tag == HasIrregular)
    {
      const Point& iv = a_vof.m_pt;
      const CellNode& node = m_internals->m_graph(iv, 0);
      retval = node.getZFaces(a_vof, a_sd, m_domain);
    }

    return retval;
  }
#endif


/*******************************/
  bool EBGraph::isConnected(const  EBIndex<CELL> & a_vof1,
                            const  EBIndex<CELL> & a_vof2) const
  {
    PR_TIME("EBGraph::isConnected");

    PR_assert(isDefined());
    const Point& iv1 = a_vof1.m_pt;
    const Point& iv2 = a_vof2.m_pt;
    if(iv1 == iv2) return false;

    PR_assert(m_region.contains(iv1));
    PR_assert(m_region.contains(iv2));
    EBIndex<CELL>  vofLo, vofHi;
    //use the lexical greaterThan because it will
    //always give the correct answer since only one component
    //is different.
    if (iv1 < (iv2))
    {
      vofLo = a_vof1;
      vofHi = a_vof2;
    }
    else
    {
      vofLo = a_vof2;
      vofHi = a_vof1;
    }
    int direction;
    bool dirfound;
    const Point& ivLo = vofLo.m_pt;
    const Point& ivHi = vofHi.m_pt;

    dirfound = false;
    for (int idir = 0; ((idir<DIM) && !dirfound); idir++)
    {
      if ((ivHi - ivLo) == Point::Basis(idir, 1))
      {
        direction = idir;
        dirfound = true;
      }
    }
    //if not neigboring intvects, no way it can be connected
    if (!dirfound) return false;

    bool voffound = false;
    if(direction == 0)
    {
      vector<EBIndex<XFACE> > faces = getXFaces(vofLo, Side::Hi);
      for (int iface = 0; iface < faces.size(); iface++)
      {
        const EBIndex<XFACE>& face = faces[iface];
        if (face.m_vofIDMe == vofHi.m_vofIDMe)
        {
          voffound = true;
        }
      }
    }
    else if(direction == 1)
    {
      vector<EBIndex<YFACE> > faces = getYFaces(vofLo, Side::Hi);
      for (int iface = 0; iface < faces.size(); iface++)
      {
        const EBIndex<YFACE>& face = faces[iface];
        if (face.m_vofIDMe == vofHi.m_vofIDMe)
        {
          voffound = true;
        }
      }
    }
#if DIM==3
    else if(direction == 2)
    {
      vector<EBIndex<ZFACE> > faces = getZFaces(vofLo, Side::Hi);
      for (int iface = 0; iface < faces.size(); iface++)
      {
        const EBIndex<ZFACE>& face = faces[iface];
        if (face.m_vofIDMe == vofHi.m_vofIDMe)
        {
          voffound = true;
        }
      }
    }
#endif
    else
    {
      PR_error("bogus direction");
    }
    return voffound;
  }

}

#endif
