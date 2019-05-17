
#ifndef _Proto_EBBoxDATA_H_
#define _Proto_EBBoxDATA_H_

#include <cmath>
#include <cstdlib>
#include <memory>
#include "Proto.H"
#include "EBProto.H"

using std::vector;
using std::shared_ptr;
namespace Proto
{
///
/**
*/
  ///
  struct EBDataLoc
  {
    int    m_dataID;   //which pointer to index
    size_t m_offset;   //distance from start of pointer 
  };

  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  class BaseEBBoxData 
  {
  public:
    ///
    BaseEBBoxData()
    {
    }

    virtual ~BaseEBBoxData()
    {
    }


    ///
    virtual  void setVal(data_t a_val)
    {
      m_regData.setVal(a_val);
      m_irrData.setVal(a_val);
    }
    

    ///for AggStencil
    inline const data_t* dataPtr(int a_dataType, int a_ivar) const
    {
      const data_t*  retval = nullptr;
      if(a_dataType == 0)
      {
        retval = m_regData.dataPtr(a_ivar);
      }
      else
      {
        retval = m_irrData.dataPtr(a_ivar);
      }
    }


    ///for AggStencil
    inline data_t* dataPtr(int a_dataType, int a_ivar) 
    {
      data_t*  retval = nullptr;
      if(a_dataType == 0)
      {
        retval = m_regData.dataPtr(a_ivar);
      }
      else
      {
        retval = m_irrData.dataPtr(a_ivar);
      }
    }

    ///needed to apply stencils
    BoxData<data_t,ncomp>& getRegData() 
    {
      return m_regData;
    }

    ///needed to apply stencils
    const BoxData<data_t,ncomp>& getRegData() const
    {
      return m_regData;
    }

  protected:
    EBGraph                            m_graph;
    PointSet                           m_irreg;
    BoxData<        data_t, ncomp>     m_regData;
    IrregData<cent, data_t, ncomp>     m_irrData;
  };


  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  class EBBoxData: public BaseEBBoxData<cent, data_t, ncomp>
  {
  public:
    ///
    EBBoxData()
    {
    }

    virtual ~EBBoxData()
    {
    }

    /// input box is cell centered
    EBBoxData(const Box& a_box, const EBGraph& a_graph)
    {
      define(a_box, a_graph);
    }
    
    ///needed for leveldata
    EBBoxData(const Box& a_box)
    {  }

  

    /// input box is cell centered---semantic is somewhat different depending on centering
    inline void define(const Box& a_box, const EBGraph& a_graph);

    /// for AggStencil -- needs to be specialized for different centerings
    inline EBDataLoc indexInfo(const EBIndex<cent>& a_baseInd, int a_ivar) const;


  };

  ///template specializations
  template <typename  data_t, unsigned int ncomp>
  class EBBoxData<CELL, data_t, ncomp>: public BaseEBBoxData<CELL, data_t, ncomp>
  {
  public:
    /// input box is cell centered
    EBBoxData(const Box& a_box, const EBGraph& a_graph)
    {
      define(a_box, a_graph);
    }

    void
    define(const Box& a_box, const EBGraph& a_graph)
    {
      m_graph = a_graph;
      m_irreg = a_graph.getIrregCells(a_box);
      vector< EBIndex<CELL> > vofs = m_graph.getIrregVoFs(a_box);

      m_irrData.define(a_box, m_graph);  
      m_regData.define(a_box);  
    }

    ///needed for leveldata
    EBBoxData(const Box& a_box)
    {  }

    ///AggStencil stuff
    EBDataLoc
    indexInfo(const EBIndex<CELL>& a_index, int a_ivar) const
    {
      EBDataLoc retval;
      if(m_irreg.contains(a_index.m_pt))
      {
        retval.m_dataID   = 1;
        retval.m_offset   = m_irrData.index(a_index, a_ivar);
      }
      else
      {
        retval.m_dataID   = 0;
        retval.m_offset   = m_regData.index(a_index.m_pt, a_ivar);
      }
      return retval;
    }
  protected:
    EBGraph                            m_graph;
    PointSet                           m_irreg;
    BoxData<        data_t, ncomp>     m_regData;
    IrregData<CELL, data_t, ncomp>     m_irrData;


  };

  template <typename  data_t, unsigned int ncomp>
  class EBBoxData<BOUNDARY, data_t, ncomp>: public BaseEBBoxData<BOUNDARY, data_t, ncomp>
  {
  public:
    /// input box is cell centered
    EBBoxData(const Box& a_box, const EBGraph& a_graph)
    {
      define(a_box, a_graph);
    }

    ///needed for leveldata
    EBBoxData(const Box& a_box)
    {  }

    void
    define(const Box& a_box, const EBGraph& a_graph)
    {
      m_graph = a_graph;
      m_irreg = a_graph.getIrregCells(a_box);
      vector< EBIndex<BOUNDARY> > vofs = m_graph.getCutFaces(a_box);

      m_irrData.define(vofs);  
      // leave regular grid data undefined
    }


    EBDataLoc
    indexInfo(const EBIndex<BOUNDARY>& a_index, int a_ivar) const
    {
      EBDataLoc retval;
      retval.m_dataID   = 1;
      retval.m_offset   = m_irrData.index(a_index, a_ivar);
      return retval;
    }

  protected:
    EBGraph                            m_graph;
    PointSet                           m_irreg;
    BoxData<        data_t, ncomp>     m_regData;
    IrregData<BOUNDARY, data_t, ncomp>     m_irrData;
  };


  template <typename  data_t, unsigned int ncomp>
  class EBBoxData<XFACE, data_t, ncomp>: public BaseEBBoxData<XFACE, data_t, ncomp>
  {
  public:
    /// input box is cell centered
    EBBoxData(const Box& a_box, const EBGraph& a_graph)
    {
      define(a_box, a_graph);
    }

    ///needed for leveldata
    EBBoxData(const Box& a_box)
    {  }

    void
    define(const Box& a_box, const EBGraph& a_graph)
    {
      m_graph = a_graph;
      m_irreg = a_graph.getIrregCells(a_box);
      vector< EBIndex<XFACE> > facs = m_graph.getIrregXFaces(a_box);

      m_irrData.define(facs);  
      Box facbox = a_box.growHi(0, 1);
      m_regData.define(facbox);
    }

    EBDataLoc
    indexInfo(const EBIndex<XFACE>& a_index, int a_ivar) const
    {
      Point hipt = a_index.m_pt;
      Point lopt = a_index.m_pt - Point::Basis(0, 1);
      EBDataLoc retval;
      if(m_irreg.contains(hipt) || m_irreg.contains(lopt) )
      {
        retval.m_dataID   = 1;
        retval.m_offset   = m_irrData.index(a_index, a_ivar);
      }
      else
      {
        retval.m_dataID   = 0;
        retval.m_offset   = m_regData.index(hipt, a_ivar);
      }
      return retval;
    }
  protected:
    EBGraph                            m_graph;
    PointSet                           m_irreg;
    BoxData<        data_t, ncomp>     m_regData;
    IrregData<XFACE, data_t, ncomp>     m_irrData;

  };


  template <typename  data_t, unsigned int ncomp>
  class EBBoxData<YFACE, data_t, ncomp>: public BaseEBBoxData<YFACE, data_t, ncomp>
  {
  public:
    /// input box is cell centered
    EBBoxData(const Box& a_box, const EBGraph& a_graph)
    {
      define(a_box, a_graph);
    }

    ///needed for leveldata
    EBBoxData(const Box& a_box)
    {  }

    void
    define(const Box& a_box, const EBGraph& a_graph)
    {
      m_graph = a_graph;
      m_irreg = a_graph.getIrregCells(a_box);
      vector< EBIndex<YFACE> > facs = m_graph.getIrregYFaces(a_box);

      m_irrData.define(facs);  
      Box facbox = a_box.growHi(1, 1);
      m_regData.define(facbox);
    }

    EBDataLoc
    indexInfo(const EBIndex<YFACE>& a_index, int a_ivar) const
    {
      Point hipt = a_index.m_pt;
      Point lopt = a_index.m_pt - Point::Basis(1, 1);
      EBDataLoc retval;
      if(m_irreg.contains(hipt) || m_irreg.contains(lopt) )
      {
        retval.m_dataID   = 1;
        retval.m_offset   = m_irrData.index(a_index, a_ivar);
      }
      else
      {
        retval.m_dataID   = 0;
        retval.m_offset   = m_regData.index(hipt, a_ivar);
      }
      return retval;
    }

  protected:
    EBGraph                            m_graph;
    PointSet                           m_irreg;
    BoxData<        data_t, ncomp>     m_regData;
    IrregData<YFACE, data_t, ncomp>     m_irrData;
  };


  template <typename  data_t, unsigned int ncomp>
  class EBBoxData<ZFACE, data_t, ncomp>: public BaseEBBoxData<ZFACE, data_t, ncomp>
  {
  public:
    /// input box is cell centered
    EBBoxData(const Box& a_box, const EBGraph& a_graph)
    {
      define(a_box, a_graph);
    }

    ///needed for leveldata
    EBBoxData(const Box& a_box)
    {  }

    void
    define(const Box& a_box, const EBGraph& a_graph)
    {
      
      m_graph = a_graph;
      m_irreg = a_graph.getIrregCells(a_box);
#if DIM==3
      vector< EBIndex<ZFACE> > facs = m_graph.getIrregZFaces(a_box);

      m_irrData.define(facs);
      Box facbox = a_box.growHi(2, 1);
      m_regData.define(facbox);
#endif      
    }

    EBDataLoc
    indexInfo(const EBIndex<ZFACE>& a_index, int a_ivar) const
    {
      Point hipt = a_index.m_pt;
      Point lopt = a_index.m_pt - Point::Basis(1, 1);
      EBDataLoc retval;
      if(m_irreg.contains(hipt) || m_irreg.contains(lopt) )
      {
        retval.m_dataID = 1;
        retval.m_offset = m_irrData.index(a_index, a_ivar);
      }
      else
      {
        retval.m_dataID   = 0;
        retval.m_offset   = m_regData.indesx(hipt, a_ivar);
      }
      return retval;
    }

  protected:
    EBGraph                             m_graph;
    PointSet                            m_irreg;
    BoxData<         data_t, ncomp>     m_regData;
    IrregData<ZFACE, data_t, ncomp>     m_irrData;
  };
}

#endif
