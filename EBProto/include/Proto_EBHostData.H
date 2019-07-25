
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

    ///needed for Proto::LevelData.  meaningless here becasue we need the graph.
    EBHostData(const Box& a_box, bool a_stackAlloc)
    {
    }

    EBHostData(const Box& a_box, const EBGraph& a_graph)
    {
      define(a_box, a_graph);
    }

    virtual ~EBHostData()
    {
    }
 
    void copy(const EBHostData<cent, data_t, ncomp>   &  a_src,
              const Box                               &  a_srcbox,
              int                                        a_srccomp,
              const Box                               &  a_dstbox,
              int                                        a_dstcomp,
              int                                        a_numcomp)
    {
      m_regData.copy(a_src.m_regData, a_srcbox, a_srccomp, a_dstbox, a_dstcomp, a_numcomp);
      m_irrData.copy(a_src.m_irrData, a_srcbox, a_srccomp, a_dstbox, a_dstcomp, a_numcomp);
    }

    ///this DOES include sizeof(T)  
    virtual std::size_t charsize(const Box&   a_bx,
                                 unsigned int a_startcomp,
                                 unsigned int a_numcomps) const
    {
      PR_assert((a_startcomp + a_numcomps) < ncomp);
      size_t regsize = m_regData.charsize(a_bx, a_startcomp, a_numcomps);
      size_t irrsize = m_irrData.charsize(a_bx, a_startcomp, a_numcomps);
      return regsize+irrsize;
      
    }
    virtual void  linearIn(void*                  a_buf,
                           const Box       &      a_bx,
                           unsigned int           a_startcomp,
                           unsigned int           a_numcomps)
    {
      PR_assert((a_startcomp + a_numcomps) < ncomp);
      char* charbuf = (char*)a_buf;
      m_regData.linearIn(charbuf,   a_bx, a_startcomp, a_numcomps);
      charbuf += m_regData.charsize(a_bx, a_startcomp, a_numcomps);
      m_irrData.linearIn(charbuf,   a_bx, a_startcomp, a_numcomps);
      charbuf += m_irrData.charsize(a_bx, a_startcomp, a_numcomps);
    }
    virtual void linearOut(void*               a_buf,
                           const Box       &   a_bx,
                           unsigned int        a_startcomp,
                           unsigned int        a_numcomps) const
    {
      PR_assert((a_startcomp + a_numcomps) < ncomp);
      char* charbuf = (char*)a_buf;
      m_regData.linearOut(charbuf,  a_bx, a_startcomp, a_numcomps);
      charbuf += m_regData.charsize(a_bx, a_startcomp, a_numcomps);
      m_irrData.linearOut(charbuf,  a_bx, a_startcomp, a_numcomps);
      charbuf += m_irrData.charsize(a_bx, a_startcomp, a_numcomps);
    }
    ///
    static int preAllocatable()
    {
      return 1; 
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
        a_dataPtr  =  m_regData.data();
        a_datasize =  m_regData.size();
      }
      else 
      {
        a_dataPtr  =  m_irrData.data();
        a_datasize =  m_irrData.size();
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
