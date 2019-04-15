#ifndef __Proto_EBGraphIMPLEM__
#define __Proto_EBGraphIMPLEM__

namespace Proto
{

/*******************************/
  PointSet EBGraph::getIrregCells(const Box& a_subbox) const
  {
    static PointSet emptySet;
    if (m_irregIVS == NULL) return emptySet;
    return (*m_irregIVS) & a_subbox;
  }

/*******************************/
  void EBGraph::defineAsAllRegular(const Box& a_region, const Box & a_domain)
  {
    m_isDefined = true;
    m_domain = a_domain;
    m_region = a_region;
    PR_assert(isDefined());
    PR_assert(isDomainSet());
    m_tag = AllRegular;
    m_internals = shared_ptr<localData>();
  }

/*******************************/
  void EBGraph::defineAsAllCovered()
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
                                     const Box              & a_domain)
  {
    CH_TIME("EBGraph::buildGraph");

    m_domain = a_domain;
    m_region = a_region;
    m_tag = HasIrregular;
/**
   bunch of work to do here
   if (m_irregIVS != NULL) delete m_irregIVS;
   if (m_multiIVS != NULL) delete m_multiIVS;
   m_multiIVS = new PointSet(DensePointSet(m_region, false));
   m_irregIVS = new PointSet(DensePointSet(m_region, false));
   m_graph.define(m_region, 1);

   // set regular and covered cells
   for (BoxIterator bit(m_region); bit.ok(); ++bit)
   {
   const Point& iv = bit();
   if (a_regIrregCovered(iv, 0) == 1) //regular cell
   {
   m_graph(iv, 0).defineAsRegular();
   }
   else if (a_regIrregCovered(iv, 0) == -1) // covered cell
   {
   m_graph(iv, 0).defineAsCovered();
   }
   else if (a_regIrregCovered(iv, 0) != 0)
   {
   MayDay::Error("invalid flag");
   }
   }

   // now for irregular cells
   // add the vofs
   for (int ivecIrreg = 0; ivecIrreg < a_irregGraph.size(); ivecIrreg++)
   {
   const IrregNode& inputNode = a_irregGraph[ivecIrreg];
   const Point& iv =inputNode.m_cell;

   GraphNodeImplem newImplem;
   newImplem.m_nodeInd = ivecIrreg;

   m_graph(iv, 0).addIrregularNode(newImplem, inputNode.m_cellIndex);
   (*m_irregIVS) |= iv;
   if (m_graph(iv, 0).size() > 1)
   {
   (*m_multiIVS) |= iv;
   }
   }

   // add the faces
   for (int ivecIrreg = 0; ivecIrreg < a_irregGraph.size(); ivecIrreg++)
   {
   const IrregNode& inputNode = a_irregGraph[ivecIrreg];
   const Point& iv =inputNode.m_cell;
   vector<GraphNodeImplem>& vecNodes =
   *(m_graph(iv, 0).m_cellList);
   // pick out which node we are talking about
   // by maching its nodeInd with ivecIrreg
   bool foundNode = false;
   GraphNodeImplem* nodePtr = NULL;
   for (int ivecGraph = 0; ivecGraph < vecNodes.size(); ivecGraph++)
   {
   if (vecNodes[ivecGraph].m_nodeInd == ivecIrreg)
   {
   foundNode = true;
   nodePtr = &(vecNodes[ivecGraph]);
   }
   }
   if (!foundNode)
   {
   MayDay::Error("EBGraph: internal error in construction");
   }
   // now add the arcs in the input to the node
   GraphNodeImplem& node = *nodePtr;
   for (int idir = 0; idir < SpaceDim; idir++)
   {
   for (SideIterator sit; sit.ok(); ++sit)
   {
   int irregIndex = IrregNode::index(idir, sit());
   int gNodeIndex = IrregNode::index(idir, sit());
   const vector<int>& irregArcs = inputNode.m_arc[irregIndex];
   vector<int>& nodeArcs = node.m_arc[gNodeIndex];
   nodeArcs.resize(irregArcs.size());
   for (int iarc = 0; iarc < irregArcs.size(); iarc++)
   {
   int otherNodeInd = irregArcs[iarc];
   if (otherNodeInd == -1)
   {
   // if otherNodeInd == -1, boundary arc.
   // just make the arc in our node = -1
   // to signify the same
   nodeArcs[iarc] = -1;
   }
   else if (otherNodeInd == -2)
   {
   // this means that the vof is connected
   // to a regular vof,
   // which always have a cell index of 0
   nodeArcs[iarc] =  0;
   }
   else
   {
   nodeArcs[iarc] =  otherNodeInd;
   }
   }
   }
   }
   }
**/
  }


/*******************************/
  vector<VolIndex> EBGraph::getVoFs(const Point& a_iv) const
  {
    PR_assert(isDefined());
    PR_assert(isDomainSet());
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
      PR_assert(m_region.contains(a_iv));
      PR_assert(m_domain.contains(a_iv));
      const GraphNode& node = m_graph(a_iv, 0);
      retvec = node.getVoFs(a_iv);
    }
    return retvec;
  }

/*******************************/
  bool EBGraph::isRegular(const Point& a_iv) const
  {
    PR_assert(isDefined());
    PR_assert(isDomainSet());
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
      //PR_assert(m_region.contains(a_iv)); //picked up my m_graph already
      //PR_assert(m_domain.contains(a_iv));
      const GraphNode& node = m_graph(a_iv, 0);
      retval = node.isRegular();
    }
    else
    {
      retval = false;
      MayDay::Error("EBGraph::isRegular:Bogus Tag");
    }
    return retval;
  }


/*******************************/
  bool EBGraph::isIrregular(const Point& a_iv) const
  {
    PR_assert(isDefined());
    PR_assert(isDomainSet());
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
      PR_assert(m_region.contains(a_iv));
      PR_assert(m_domain.contains(a_iv));
      const GraphNode& node = m_graph(a_iv, 0);
      retval = node.isIrregular();
    }
    else
    {
      retval = false;
      MayDay::Error("EBGraph::isIrregular:Bogus Tag");
    }
    return retval;
  }


/*******************************/
  bool EBGraph::isCovered(const Point& a_iv) const
  {
    PR_assert(isDefined());
    PR_assert(isDomainSet());
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
      //PR_assert(m_region.contains(a_iv)); this check picked up by m_graph
      //PR_assert(m_domain.contains(a_iv));
      const GraphNode& node = m_graph(a_iv, 0);
      retval = node.isCovered();
    }
    else
    {
      retval = false;
      MayDay::Error("EBGraph::isIrregular:Bogus Tag");
    }
    return retval;
  }

/*******************************/
  bool EBGraph::isCovered(const Box& a_box) const
  {
    PR_assert(isDefined());
    PR_assert(isDomainSet());
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
      PR_assert(m_region.contains(a_box));
      PR_assert(m_domain.contains(a_box));
      retval = true;
      BoxIterator bit(a_box);
      for (bit.reset(); bit.ok() && retval; ++bit)
      {
        if (!isCovered(bit()))
        {
          retval = false;
        }
      }
    }
    else
    {
      retval = false;
      MayDay::Error("EBGraph::isCovered:Bogus Tag");
    }
    return retval;
  }

/*******************************/
  bool EBGraph::isRegular(const Box& a_box) const
  {
    PR_assert(isDefined());
    PR_assert(isDomainSet());
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
      PR_assert(m_region.contains(a_box));
      PR_assert(m_domain.contains(a_box));
      BoxIterator bit(a_box);
      retval = true;
      for (bit.reset(); bit.ok() && retval; ++bit)
      {
        if (!isRegular(bit()))
        {
          retval = false;
        }
      }
    }
    else
    {
      retval = false;
      MayDay::Error("EBGraph::isRegular:Bogus Tag");
    }
    return retval;
  }

/*******************************/
  vector<FaceIndex> EBGraph::getFaces(const VolIndex&       a_vof,
                                      const int&            a_idir,
                                      const Side::LoHiSide& a_sd) const
  {
    PR_assert(isDefined());
    PR_assert(isDomainSet());
    PR_assert(m_region.contains(a_vof.gridIndex()));
    PR_assert(m_domain.contains(a_vof.gridIndex()));
    PR_assert((a_idir >= 0) && (a_idir < SpaceDim));
    PR_assert((a_sd == Side::Lo) || (a_sd == Side::Hi));

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
    CH_TIME("EBGraph::isConnected");

    PR_assert(isDefined());
    const Point& iv1 = a_vof1.gridIndex();
    const Point& iv2 = a_vof2.gridIndex();

    PR_assert(m_region.contains(iv1));
    PR_assert(m_region.contains(iv2));
    VolIndex vofLo, vofHi;
    //use the lexical greaterThan because it will
    //always give the correct answer since only one component
    //is different.
    if (iv1.lexGT(iv2))
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
    for (int idir = 0; ((idir<SpaceDim) && !dirfound); idir++)
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
