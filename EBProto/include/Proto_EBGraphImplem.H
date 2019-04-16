#ifndef __Proto_EBGraphIMPLEM__
#define __Proto_EBGraphIMPLEM__

namespace Proto
{

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
                                     const vector<IrregNode>& a_irregGraph,
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
        m_graph(iv, 0).defineAsRegular();
      }
      else if (a_regIrregCovered(iv, 0) == -1) // covered cell
      {
        m_graph(iv, 0).defineAsCovered();
      }
      else if (a_regIrregCovered(iv, 0) == 0)
      {
        m_graph(iv, 0).defineAsIrregular();
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
     const IrregNode & inputNode = a_irregGraph[ivecIrreg];
     const Point     & iv = inputNode.m_cell;
     m_internals->m_irregIVS |= iv;

     GraphCell::IrrGraphNode currNode;
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
     m_graph(iv, 0).m_cellList->push_back(currNode);
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
      const GraphNode& node = m_graph(a_iv, 0);
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
      const GraphNode& node = m_graph(a_iv, 0);
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
      const GraphNode& node = m_graph(a_iv, 0);
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
      const GraphNode& node = m_graph(a_iv, 0);
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
  vector<FaceIndex> EBGraph::getFaces(const VolIndex&       a_vof,
                                      const int&            a_idir,
                                      const Side::LoHiSide& a_sd) const
  {
    PR_assert(isDefined());

    vector<FaceIndex> retvec;
    if (m_tag == AllRegular)
    {
      Point otherIV = a_vof.gridIndex()
        + sign(a_sd)*BASISV(a_idir);
      int otherCellInd = 0;
      if (!m_domain.contains(otherIV))
      {
        otherCellInd = -1;
      }
      VolIndex otherVoF(otherIV, otherCellInd);
      retvec.push_back(FaceIndex(a_vof, otherVoF, a_idir));
    }
    else if (m_tag == AllCovered)
    {
      //return empty vector
    }
    else if (m_tag == HasIrregular)
    {
      const Point& iv = a_vof.gridIndex();
      const GraphNode& node = m_graph(iv, 0);
      retvec = node.getFaces(a_vof, a_idir, a_sd, m_domain);
    }

    return retvec;
  }



/*******************************/
  bool EBGraph::isConnected(const VolIndex& a_vof1,
                            const VolIndex& a_vof2) const
  {
    PR_TIME("EBGraph::isConnected");

    PR_assert(isDefined());
    const Point& iv1 = a_vof1.gridIndex();
    const Point& iv2 = a_vof2.gridIndex();
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
    const Point& ivLo = vofLo.gridIndex();
    const Point& ivHi = vofHi.gridIndex();

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

    vector<FaceIndex> faces = getFaces(vofLo, direction, Side::Hi);
    bool voffound = false;
    for (int iface = 0; iface < faces.size(); iface++)
    {
      const FaceIndex& face = faces[iface];
      if (face.getVoF(Side::Hi) == vofHi)
      {
        voffound = true;
      }
    }
    return voffound;
  }

}

#endif
