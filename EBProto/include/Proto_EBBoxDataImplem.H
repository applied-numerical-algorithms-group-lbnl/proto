
#ifndef _Proto_IrregDATA_H_
#define _Proto_IrregDATA_H_

#include <cmath>
#include <cstdlib>
#include <memory>
#include "Proto_Box.H"

using std::vector;
using std::shared_ptr;
namespace Proto
{
///
/**
*/
  template <CENTERING cent, typename T, unsigned int C=1>
  class EBBoxData 
  {
  public:
    ///
    EBBoxData()
    {
    }

    ///
    EBBoxData(const Box& a_box)
    {
      define(a_box);
    }

    
    ///needed for leveldata
    EBBoxData(const Box& a_box)
    {  }
    ///



    ///
    inline std::size_t size() const {return m_indicies->size()*C;};

    ///
    inline bool defined() const {return bool(m_data);};
   

    ///
    inline void setVal(const data_t& a_val);
  

    ///
    T* data()
    {
      return m_rawPtr;
    }

    ///
    const T* data() const
    {
      return m_rawPtr;
    }

    ///for AggStencil
    inline size_t offset(const index_t& a_vof, const int& a_ivar) const
    {
      unsigned int i = 0;
      for(;i<m_indicies->size(); i++)
        if(*vof == (*m_indicies)[i]) break;
      i+=a_ivar*m_indicies->size();
      return i;
    }


  protected:

    shared_ptr<vector<index_t> >  m_indicies;
    shared_ptr<data_t>            m_data;
    data_t*                       m_rawPtr;
  };

}


#include "Proto_EBBoxDataImplem.H"

#endif
