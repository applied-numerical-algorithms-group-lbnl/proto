#ifndef __Proto_EBGraphIMPLEM__
#define __Proto_EBGraphIMPLEM__

namespace Proto
{
  void
  GraphCell::
  getFaceInfo(int   & a_vecsize,
              Point & a_mypt,
              int   & a_myVoFID,
              int   & a_otherVoFID
              bool  & a_isBoundary,
              const int            & a_dir,
              const VolIndex       & a_vof,
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
      Point otherIV    = a_vof.m_pt +sign(a_sd)*BASISV(a_dir);

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
        vector<IrrCellNode>& nodeVec = *m_cellList;
        //again, we are dealing with single-valued stuff at this level for now
        //multivalued can only happen by coarsening
        PR_assert(nodeVec.size() == 1);
        int index = ebp_index(a_dir, a_sd);
        vector<int>& arcvec = nodeVec[0].m_arcs[index];
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
  GraphCell::
  getXFaces(const VolIndex       & a_vof,
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

  }

 /*******************************/
  vector<EBIndex<YFACE> >
  GraphCell::
  getYFaces(const VolIndex       & a_vof,
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

  }

 /*******************************/
#if DIM==3
  vector<EBIndex<YFACE> >
  GraphCell::
  getZFaces(const VolIndex       & a_vof,
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

  }
#endif
 /*******************************/
  PointSet EBGraph::getIrregCells(const Box& a_subbox) const
  {
    static PointSet retval;
    if(m_internals)
    {
      retval = (m_internals->m_irregIVS) & a_subbox;
    }
    
    return retval;
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
        m_localData->m_graph(iv, 0).defineAsRegular();
      }
      else if (a_regIrregCovered(iv, 0) == -1) // covered cell
      {
        m_localData->m_graph(iv, 0).defineAsCovered();
      }
      else if (a_regIrregCovered(iv, 0) == 0)
      {
        m_localData->m_graph(iv, 0).defineAsIrregular();
        num_irreg++;
      }
      else 
      {
        PR_error("invalid flag");
      }
    }
    //it can be greater if there are multivalued cells
    PR_ASSERT(a_irregGraph.size() >= num_irreg, "mismatch in ebgraph::define inputs");

   // now for irregular cells
   for (int ivecIrreg = 0; ivecIrreg < a_irregGraph.size(); ivecIrreg++)
   {
     const IrregGraph & inputNode = a_irregGraph[ivecIrreg];
     const Point     & iv = inputNode.m_cell;
     m_internals->m_irregIVS |= iv;

     GraphCell::IrrCellNode currNode;
     for (int idir = 0; idir < DIM; idir++)
     {
       for (SideIterator sit; sit.ok(); ++sit)
       {
          int irregIndex = ebp_index(idir, sit());
          const vector<int>& irregArcs = inputNode.m_arc[irregIndex];
         //comes from the singleValuedOnly thing
         PR_assert(irregArc.size() <= 1);
         currNode.m_arc[index] = irregArcs;
       }
     }
     m_localData->m_graph(iv, 0).m_cellList->push_back(currNode);
   }
  }


/*******************************/
  vector<VolIndex> EBGraph::getVoFs(const Point& a_iv) const
  {
    PR_assert(isDefined());

    vector<VolIndex> retvec;
    if (m_tag == AllRegular)
    {
      retvec.push_back(VolIndex(a_iv, 0));
    }
    else if (m_tag == AllCovered)
    {
      //return an empty vector
    }
    else if (m_tag == HasIrregular)
    {
      const CellNode& node = m_localData->m_graph(a_iv, 0);
      retvec = node.getVoFs(a_iv);
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
      const CellNode& node = m_localData->m_graph(a_iv, 0);
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
      const CellNode& node = m_graph(a_iv, 0);
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
      const CellNode& node = m_localData->m_graph(a_iv, 0);
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
        if (!isRegular(bit()))
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
  getXFaces(const VolIndex       & a_vof,
            const Side::LoHiSide & a_sd) const
  {
    PR_assert(isDefined());

    vector< EBIndex<XFACE> > retvec;
    if (m_tag == AllRegular)
    {
      int idir = 0;
      Point otherIV = a_vof.m_pt
        + sign(a_sd)*BASISV(idir);
      int otherCellInd = 0;
      bool isBoundary = false;
      if (!m_domain.contains(otherIV))
      {
        otherCellInd = -1;
        isBoundary = true;
      }
      retvec.resize(1);
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
      retval[0].m_isBoundary = isBoundary()
    }
    else if (m_tag == AllCovered)
    {
      //return empty vector
    }
    else if (m_tag == HasIrregular)
    {
      const Point& iv = a_vof.m_pt;
      const CellNode& node = m_localData->m_graph(iv, 0);
      retvec = node.getXFaces(a_vof, a_sd, m_domain);
    }

    return retvec;
  }


/*******************************/
  vector< EBIndex<YFACE> > 
  EBGraph::
  getYFaces(const VolIndex       & a_vof,
            const Side::LoHiSide & a_sd) const
  {
    PR_assert(isDefined());

    vector< EBIndex<YFACE> > retvec;
    if (m_tag == AllRegular)
    {
      int idir = 1;
      Point otherIV = a_vof.m_pt
        + sign(a_sd)*BASISV(idir);
      int otherCellInd = 0;
      bool isBoundary = false;
      if (!m_domain.contains(otherIV))
      {
        otherCellInd = -1;
        isBoundary = true;
      }
      retvec.resize(1);
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
      retval[0].m_isBoundary = isBoundary()
    }
    else if (m_tag == AllCovered)
    {
      //return empty vector
    }
    else if (m_tag == HasIrregular)
    {
      const Point& iv = a_vof.m_pt;
      const CellNode& node = m_localData->m_graph(iv, 0);
      retvec = node.getXFaces(a_vof, a_sd, m_domain);
    }

    return retvec;
  }



/*******************************/
  bool EBGraph::isConnected(const VolIndex& a_vof1,
                            const VolIndex& a_vof2) const
  {
    PR_TIME("EBGraph::isConnected");

    PR_assert(isDefined());
    const Point& iv1 = a_vof1.m_pt;
    const Point& iv2 = a_vof2.m_pt;
    if(iv1 == iv2) return false;

    PR_assert(m_region.contains(iv1));
    PR_assert(m_region.contains(iv2));
    VolIndex vofLo, vofHi;
    //use the lexical greaterThan because it will
    //always give the correct answer since only one component
    //is different.
    if (iv1 < (iv2))
    {
      vofLo = a_vof2;
      vofHi = a_vof1;
    }
    else
    {
      vofLo = a_vof1;
      vofHi = a_vof2;
    }
    int direction;
    bool dirfound;
    const Point& ivLo = vofLo.m_pt;
    const Point& ivHi = vofHi.m_pt;

    dirfound = false;
    for (int idir = 0; ((idir<DIM) && !dirfound); idir++)
    {
      if ((ivHi - ivLo) == BASISV(idir))
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
        const EBIndex<XFACE>& face = faces[iface];
        if (face.m_vofIDMe == vofHi.m_vofIDMe)
        {
          voffound = true;
        }
      }
    }
    return voffound;
  }

}

#endif
