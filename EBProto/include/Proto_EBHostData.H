
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
    EBHostData()
    {
    }

    virtual ~EBHostData()
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

    ///
    HostIrregData<cent, data_t,ncomp>& getIrrData() 
    {
      return m_irrData;
    }

    ///needed to apply stencils
    const HostIrregData<cent, data_t,ncomp>& getIrrData() const
    {
      return m_irrData;
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

    ///for data transfers.  datasize does NOT include sizeof(data_t) (to conform with boxdata::size)
    void dataInfo(data_t*& a_dataPtr, size_t& a_datasize, unsigned int a_whichfab)
    {
      PR_assert((a_whichfab == 0) || (a_whichfab == 1));
      if(a_whichfab == 0)
      {
        dataPtr  =  m_regData.data();
        datasize =  m_regData.size();
      }
      else 
      {
        dataPtr  =  m_irrData.data();
        datasize =  m_irrData.size();
      }
    }
  protected:
    EBGraph                                m_graph;
    PointSet                               m_irreg;
    HostBoxData<        data_t, ncomp>     m_regData;
    HostIrregData<cent, data_t, ncomp>     m_irrData;
  };


}

#endif
