
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


/// dataType function template specializations
  template <typename T, unsigned int C>
  int 
  EBBoxData<CELL, T, C>::    
  dataType(const EBIndex<CELL>& a_index) const
  {
    int retval = 0;
    if(m_irreg.contains(a_index.m_pt))
    {
      retval = 1;
    }
    return retval;
  }
///
  template <typename T, unsigned int C>
  int 
  EBBoxData<BOUNDARY, T, C>::    
  dataType(const EBIndex<BOUNDARY>& a_index) const
  {
    return 1;//only irreg data defined
  }
///
  template <typename T, unsigned int C>
  int 
  EBBoxData<XFACE, T, C>::    
  dataType(const EBIndex<XFACE>& a_index) const
  {
    Point hipt = a_index.m_pt;
    Point loPt = hipt - Point::Basis(0, 1);
    int retval = 0;
    if(m_irreg.contains(lopt) || m_irreg.contains(hipt))
    {
      retval = 1;
    }
    return retval;
  }

///
  template <typename T, unsigned int C>
  int 
  EBBoxData<YFACE, T, C>::    
  dataType(const EBIndex<YFACE>& a_index) const
  {
    Point hipt = a_index.m_pt;
    Point loPt = hipt - Point::Basis(1, 1);
    int retval = 0;
    if(m_irreg.contains(lopt) || m_irreg.contains(hipt))
    {
      retval = 1;
    }
    return retval;
  }

#if DIM==3
///
  template <typename T, unsigned int C>
  int 
  EBBoxData<ZFACE, T, C>::    
  dataType(const EBIndex<ZFACE>& a_index) const
  {
    Point hipt = a_index.m_pt;
    Point loPt = hipt - Point::Basis(2, 1);
    int retval = 0;
    if(m_irreg.contains(lopt) || m_irreg.contains(hipt))
    {
      retval = 1;
    }
    return retval;
  }

#endif


///template specializations for offset
}



#endif
