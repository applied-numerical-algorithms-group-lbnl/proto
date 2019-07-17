
#ifndef _Proto_EBHostBoxDATA_H_
#define _Proto_EBHostBoxDATA_H_

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
  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  class EBHostData 
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
    

    HostBoxData<data_t,ncomp>& getRegData() 
    {
      return m_regData;
    }

    ///needed to apply stencils
    const HostBoxData<data_t,ncomp>& getRegData() const
    {
      return m_regData;
    }



    ///semantics slightly different depending on centering
    virtual void    define(const Box& a_box, const EBGraph& a_graph)
    {
      m_graph = a_graph;
      m_irreg = a_graph.getIrregCells(a_box);  

      //this picks out the right set of stuff depending on centering
      m_irrData.define(a_box, m_graph);  
      m_regData.define(a_box);  
    }
  protected:
    EBGraph                                m_graph;
    PointSet                               m_irreg;
    HostBoxData<        data_t, ncomp>     m_regData;
    HostIrregData<cent, data_t, ncomp>     m_irrData;
  };


}

#endif
