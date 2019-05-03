
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
    int    m_dataType; //which pointer to index
    size_t m_offset;   //distance from start of pointer
  };

  template <CENTERING cent, typename T=double, unsigned int C=1>
  class EBBoxData 
  {
  public:
    ///
    EBBoxData()
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

    ///
    inline void setVal(const data_t& a_val)
    {
      m_regData.setVal(a_val);
      m_irrData.setVal(a_val);
    }
  

    /// input box is cell centered---semantic is somewhat different depending on centering
    inline void define(const Box& a_box, const EBGraph& a_graph);

    /// for AggStencil -- needs to be specialized for different centerings
    inline EBDataLoc indexInfo(const EBIndex<cent>& a_baseInd, int a_ivar) const;

    ///for AggStencil
    inline const T* dataPtr(int a_dataType, int a_ivar) const
    {
      T*  retval = nullptr;
      if(a_dataType == 0)
      {
        retval = m_regData.dataPtr(a_ivar);
      }
      else
      {
        retval = m_irrData.dataPtr(a_ivar);
      }
    }

    //needed to apply stencils
    BoxData<T,C>& getRegData() 
    {
      return m_regData;
    }

    const BoxData<T,C>& getRegData() const
    {
      return m_regData;
    }

  protected:
    EBGraph                       m_graph;
    PointSet                      m_irreg;

    BoxData<T,  C>                m_regData;
    IrregData<cent, T, C>         m_irrData;
  };

}
#include "Proto_EBBoxDataImplem.H"

#endif
