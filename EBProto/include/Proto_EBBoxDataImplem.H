
#ifndef _Proto_EBBoxDataImplem_H_
#define _Proto_EBBoxDataImplem_H_

namespace Proto
{
/// define function template specializations
  template <typename T, unsigned int C>
  EBBoxData<CELL, T, C>::
  define(const Box& a_box, const EBGraph& a_graph)
  {
    m_graph = a_graph;
    m_irreg = a_graph.getIrregCells(a_box);
    vector< EBIndex<Cell> > vofs = m_graph.getIrregVoFs(a_box);

    m_irrData.define(vofs);  
    m_regData.define(a_box);  
  }

///
  template <typename T, unsigned int C>
  EBBoxData<BOUNDARY, T, C>::
  define(const Box& a_box, const EBGraph& a_graph)
  {
    m_graph = a_graph;
    m_irreg = a_graph.getIrregCells(a_box);
    vector< EBIndex<Cell> > vofs = m_graph.getIrregVoFs(a_box);

    m_irrData.define(vofs);  
    // leave regular grid data undefined
  }


///
  template <typename T, unsigned int C>
  EBBoxData<XFACE, T, C>::
  define(const Box& a_box, const EBGraph& a_graph)
  {
    m_graph = a_graph;
    m_irreg = a_graph.getIrregCells(a_box);
    vector< EBIndex<Cell> > facs = m_graph.getIrregXFaces(a_box);

    m_irrData.define(facs);  
    Box facbox = a_box.growHi(0, 1);
    m_regFAB.define(facbox);
  }


  template <typename T, unsigned int C>
  EBBoxData<YFACE, T, C>::
  define(const Box& a_box, const EBGraph& a_graph)
  {
    m_graph = a_graph;
    m_irreg = a_graph.getIrregCells(a_box);
    vector< EBIndex<Cell> > facs = m_graph.getIrregYFaces(a_box);

    m_irrData.define(facs);  
    Box facbox = a_box.growHi(1, 1);
    m_regFAB.define(facbox);
  }

#if DIM==3
  template <typename T, unsigned int C>
  EBBoxData<ZFACE, T, C>::
  define(const Box& a_box, const EBGraph& a_graph)
  {
    m_graph = a_graph;
    m_irreg = a_graph.getIrregCells(a_box);
    vector< EBIndex<Cell> > facs = m_graph.getIrregZFaces(a_box);

    m_irrData.define(facs);  
    Box facbox = a_box.growHi(2, 1);
    m_regFAB.define(facbox);
  }
#endif


/// indexInfo function template specializations
  template <typename T, unsigned int C>
  EBDataLoc
  EBBoxData<CELL, T, C>::    
  indexInfo(const EBIndex<CELL>& a_index, int a_ivar) const
  {
    EBDataLoc retval;
    if(m_irreg.contains(a_index.m_pt))
    {
      retval.m_dataType = 1;
      retval.m_offset   = m_irrFAB.offset(a_index, a_ivar);
    }
    else
    {
      retval.m_dataType = 1;
      retval.m_offset   = m_regFAB.offset(a_index, a_ivar);
    }
    return retval;
  }
///
  template <typename T, unsigned int C>
  EBDataLoc
  EBBoxData<BOUNDARY, T, C>::    
  indexInfo(const EBIndex<BOUNDARY>& a_index, int a_ivar) const
  {
    EBDataLoc retval;
    retval.m_dataType = 1;
    retval.m_offset   = m_irrFAB.offset(a_index, a_ivar);
    return retval;
  }
///
  template <typename T, unsigned int C>
  EBDataLoc
  EBBoxData<XFACE, T, C>::    
  indexInfo(const EBIndex<XFACE>& a_index, int a_ivar) const
  {
    Point hipt = a_index.m_pt;
    Point lopt = a_index.m_pt - Point::Basis(0, 1);
    EBDataLoc retval;
    if(m_irreg.contains(hipt) || m_irreg.contains(lopt) )
    {
      retval.m_dataType = 1;
      retval.m_offset   = m_irrFAB.offset(a_index, a_ivar);
    }
    else
    {
      retval.m_dataType = 1;
      retval.m_offset   = m_regFAB.offset(a_index, a_ivar);
    }
    return retval;
  }

///
  template <typename T, unsigned int C>
  EBDataLoc
  EBBoxData<YFACE, T, C>::    
  indexInfo(const EBIndex<YFACE>& a_index, int a_ivar) const
  {
    Point hipt = a_index.m_pt;
    Point lopt = a_index.m_pt - Point::Basis(1, 1);
    EBDataLoc retval;
    if(m_irreg.contains(hipt) || m_irreg.contains(lopt) )
    {
      retval.m_dataType = 1;
      retval.m_offset   = m_irrFAB.offset(a_index, a_ivar);
    }
    else
    {
      retval.m_dataType = 1;
      retval.m_offset   = m_regFAB.offset(a_index, a_ivar);
    }
    return retval;
  }

#if DIM==3
  template <typename T, unsigned int C>
  EBDataLoc
  EBBoxData<ZFACE, T, C>::    
  indexInfo(const EBIndex<ZFACE>& a_index, int a_ivar) const
  {
    Point hipt = a_index.m_pt;
    Point lopt = a_index.m_pt - Point::Basis(1, 1);
    EBDataLoc retval;
    if(m_irreg.contains(hipt) || m_irreg.contains(lopt) )
    {
      retval.m_dataType = 1;
      retval.m_offset   = m_irrFAB.offset(a_index, a_ivar);
    }
    else
    {
      retval.m_dataType = 1;
      retval.m_offset   = m_regFAB.offset(a_index, a_ivar);
    }
    return retval;
  }


#endif



}



#endif
