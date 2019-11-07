
#ifndef _Proto_EBFluxDATA_H_
#define _Proto_EBFluxDATA_H_

#include <cmath>
#include <cstdlib>
#include <memory>
#include "Proto.H"
#include "EBProto.H"

using std::vector;
using std::shared_ptr;
namespace Proto
{
/// Convenience structure that holds coordinate face data in all directions.
/**
   This is just so I can just declare one thing and have all the fluxes.
 */
  ///
  template <typename  data_t, unsigned int ncomp>
  class EBFluxData
  {
  public:
    ///
    EBFluxData()
    {
    }

    virtual ~EBFluxData()
    {
    }


    /// sets all face data to the input value
    inline void setVal(data_t a_val)
    {
      m_xflux.setVal(a_val);
      m_yflux.setVal(a_val);
#if DIM==3      
      m_zflux.setVal(a_val);
#endif      
    }

    /// 
   /**
      The input box is cell centered.  Each face data holder has a box = surrounding nodes of the box in its direction.  Calls define function.
    */
    EBFluxData(const Box& a_box, const EBGraph& a_graph)
    {
      define(a_box, a_graph);
    }
    
    /// define all data
    /**
       The input box is cell centered.  Each face data holder has a box = surrounding nodes of the box in its direction.
    */
    inline void   define(const Box& a_box, const EBGraph& a_graph)
    {
      Box xbox = a_box.growHi(0, 1);
      Box ybox = a_box.growHi(1, 1);

      auto* ptrx = new EBBoxData<XFACE, data_t, ncomp>(xbox, a_graph);
      auto* ptry = new EBBoxData<YFACE, data_t, ncomp>(ybox, a_graph);
      m_xflux = shared_ptr< EBBoxData<XFACE, data_t, ncomp> >(ptrx); 
      m_yflux = shared_ptr< EBBoxData<YFACE, data_t, ncomp> >(ptry);
#if DIM==3      
      Box zbox = a_box.growHi(2, 1);
      auto* ptrz = new EBBoxData<ZFACE, data_t, ncomp>(zbox, a_graph);
      m_zflux = shared_ptr< EBBoxData<ZFACE, data_t, ncomp> >(ptrz);
#endif      

    }
    ///data  on x coordinate faces
    shared_ptr< EBBoxData<XFACE, data_t, ncomp> >    m_xflux;
    ///data  on y coordinate faces
    shared_ptr< EBBoxData<YFACE, data_t, ncomp> >    m_yflux;
    ///data  on z coordinate faces
    shared_ptr< EBBoxData<ZFACE, data_t, ncomp> >    m_zflux;


  protected:
    EBFluxData(const            EBFluxData<data_t, ncomp>& a_src);
    EBFluxData& operator=(const EBFluxData<data_t, ncomp>& a_src);

    

  };

}

#endif
