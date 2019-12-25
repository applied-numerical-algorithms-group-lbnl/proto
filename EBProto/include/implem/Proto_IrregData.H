
#ifndef _Proto_IrregDATA_H_
#define _Proto_IrregDATA_H_

#include <cmath>
#include <cstdlib>
#include <memory>
#include "Proto_Box.H"
#include "EBProto.H"

using std::vector;
using std::shared_ptr;
namespace Proto
{


///
/**
*/
  template <CENTERING cent, class data_t,  unsigned int ncomp>
  class BaseIrregData
  {
  public:
    typedef EBIndex<cent> index_t;
    ///
    BaseIrregData()
    {
      m_boxdefined = false;
    }

    ///
    virtual ~BaseIrregData()
    {
    }

    ///
    BaseIrregData( const shared_ptr< vector<index_t> >& a_indicies)
    {
      define(a_indicies);
    }

    ///needed for leveldata
    BaseIrregData(const Box& a_box, bool a_stackAllocation=true)
    {  }

    ///
    inline void define( const shared_ptr< vector< index_t > >& a_indicies);

    ///
    inline std::size_t size() const {return m_indicies->size()*ncomp;};

    ///
    inline bool defined() const {return bool(m_data);};
   

    ///
    inline void setVal(const data_t& a_val);
  

    ///
    data_t* dataPtr(unsigned int a_ivar = 0)
    {
      PR_assert(a_ivar < ncomp);
      size_t ioff = a_ivar*m_indicies->size();
      return &(m_rawPtr[ioff]);
    }

    ///
    const data_t* dataPtr(unsigned int a_ivar = 0) const
    {
      PR_assert(a_ivar < ncomp);
      size_t ioff = a_ivar*m_indicies->size();
      return &(m_rawPtr[ioff]);
    }

    data_t* data()
    {
      return m_rawPtr;
    }


    const data_t* data() const
    {
      return m_rawPtr;
    }
    
    ///
    inline size_t  index(const EBIndex<cent>& a_vof, const int& a_ivar) const
    {
      size_t ivec = 0;
      bool found = false;
      for(unsigned int i = 0; i < m_indicies->size(); i++)
      {        
        if(a_vof == (*m_indicies)[i])
        {
          ivec = i;
          found = true;
          break;
        }
      }
      if(!found)
      {
        std::cout << "index<cell> vof not found" << endl;
        abort();
      }
      PROTO_ASSERT(found, "vof not found");

      ivec += a_ivar*m_indicies->size();

      if(ivec >= size())
      {
        std::cout << "index<cell> offset too big" << endl;
        abort();
      }

      return ivec;
    }

    BaseIrregData<cent, data_t, ncomp>& 
    operator+=(const BaseIrregData<cent, data_t, ncomp> & a_input)
    {
#ifdef PROTO_CUDA
      if(m_boxdefined && (m_inputBox == a_input.m_inputBox))
      {
        T* castPtrInp = (T*)(a_input.m_rawPtr);
        thrust::device_ptr<T> devptrLoc = thrust::device_pointer_cast(m_rawPtr);
        thrust::device_ptr<T> devptrInp = thrust::device_pointer_cast(castPtrInp);
        unsigned int nsize = size();
        thrust::transform(thrust::cuda::par.on(DisjointBoxLayout::getCurrentStream()),
                          devptrLoc, devptrLoc+nsize, devptrInp, devptrLoc, thrust::plus<data_t>());
        
      }
      else
      {
        PROTO_ASSERT(false, "operator+= not implemented for irreg datas of different sizes");
      }
#else
      for(int i = 0; i < m_indicies->size(); i++)
      {
        for(int icomp = 0; icomp < ncomp; icomp++)
        {
          size_t dstind =             index((*m_indicies)[i], icomp);
          size_t srcind =     a_input.index((*m_indicies)[i], icomp);
          m_rawPtr[dstind] += a_input.m_rawPtr[srcind];
        }
      }

#endif
      return *this;
    }


    BaseIrregData<cent, data_t, ncomp>& 
    operator-=(const BaseIrregData<cent, data_t, ncomp> & a_input)
    {
#ifdef PROTO_CUDA
      if(m_boxdefined && (m_inputBox == a_input.m_inputBox))
      {
        T* castPtrInp = (T*)(a_input.m_rawPtr);
        thrust::device_ptr<T> devptrLoc = thrust::device_pointer_cast(m_rawPtr);
        thrust::device_ptr<T> devptrInp = thrust::device_pointer_cast(castPtrInp);
        unsigned int nsize = size();
        thrust::transform(thrust::cuda::par.on(DisjointBoxLayout::getCurrentStream()),
                          devptrLoc, devptrLoc+nsize, devptrInp, devptrLoc, thrust::plus<data_t>());
        
      }
      else
      {
        PROTO_ASSERT(false, "operator+= not implemented for irreg datas of different sizes");
      }
#else
      for(int i = 0; i < m_indicies->size(); i++)
      {
        for(int icomp = 0; icomp < ncomp; icomp++)
        {
          size_t dstind =             index((*m_indicies)[i], icomp);
          size_t srcind =     a_input.index((*m_indicies)[i], icomp);
          m_rawPtr[dstind] -= a_input.m_rawPtr[srcind];
        }
      }

#endif
      return *this;
    }


    BaseIrregData<cent, data_t, ncomp>& 
    operator*=(const BaseIrregData<cent, data_t, ncomp> & a_input)
    {
#ifdef PROTO_CUDA
      if(m_boxdefined && (m_inputBox == a_input.m_inputBox))
      {
        T* castPtrInp = (T*)(a_input.m_rawPtr);
        thrust::device_ptr<T> devptrLoc = thrust::device_pointer_cast(m_rawPtr);
        thrust::device_ptr<T> devptrInp = thrust::device_pointer_cast(castPtrInp);
        unsigned int nsize = size();
        thrust::transform(thrust::cuda::par.on(DisjointBoxLayout::getCurrentStream()),
                          devptrLoc, devptrLoc+nsize, devptrInp, devptrLoc, thrust::multiplies<data_t>());
        
      }
      else
      {
        PROTO_ASSERT(false, "operator+= not implemented for irreg datas of different sizes");
      }
#else
      for(int i = 0; i < m_indicies->size(); i++)
      {
        for(int icomp = 0; icomp < ncomp; icomp++)
        {
          size_t dstind =             index((*m_indicies)[i], icomp);
          size_t srcind =     a_input.index((*m_indicies)[i], icomp);
          m_rawPtr[dstind] *= a_input.m_rawPtr[srcind];
        }
      }

#endif
      return *this;
    }
    
    /// debugging only--will not work on the device
    data_t operator() (const index_t& a_vof, int a_comp) const
    {
      size_t ind = index(a_vof, a_comp);
      return m_rawPtr[ind];
    }

    ///might need a different semantic for face-centered data.
    vector<EBIndex<cent> > getIndices(const Box& a_box) const
    {
      vector<EBIndex<cent> > retval;
      for(int i = 0; i < m_indicies->size(); i++)
      {
        if(a_box.contains((*m_indicies)[i].m_pt))
        {
          retval.push_back((*m_indicies)[i]);
        }
      }
      return retval;
    }

    unsigned int vecsize() const
    {
      return m_indicies->size();
    }

    //aliasing define
    template<unsigned int srccomp>
    void define(BaseIrregData<cent, data_t, srccomp>& a_input,
                unsigned int                        & a_comp)
    {
      m_indicies = a_input.getIndicies();
      m_data     = a_input.getData();
      m_rawPtr   = a_input.dataPtr(a_comp);
    }           

    //do not use this
    shared_ptr<vector<  index_t > >  getIndicies() const
    {
      return m_indicies;
    }

    //do not use this
    shared_ptr<data_t> getData() const
    {
      return  m_data;
    }

    
  protected:

    bool m_boxdefined;
    Box m_inputBox; //exists if defined as part of a boxdata
    shared_ptr<vector<  index_t > >        m_indicies;
    shared_ptr<data_t>         m_data;
    data_t*                    m_rawPtr;
    BaseIrregData(const BaseIrregData<cent, data_t,ncomp>& a_src);
    BaseIrregData& operator=(const BaseIrregData<cent, data_t,ncomp>& a_src);


  };


  template <CENTERING cent, class data_t,  unsigned int ncomp>
  class IrregData: public BaseIrregData<cent, data_t, ncomp>
  {
  public:
    typedef EBIndex<cent> index_t;
    ///
    IrregData()
    {
    }

    virtual ~IrregData()
    {
    }

    ///
    IrregData( const shared_ptr< vector<index_t> >& a_indicies)
    {
      define(a_indicies);
    }


    inline void define(const Box& a_box,  const EBGraph& a_graph);

  };

  ///template specialization for cell based data
  template < class data_t,  unsigned int ncomp>
  class IrregData <CELL, data_t, ncomp>: public BaseIrregData<CELL, data_t, ncomp>
  {
  public:

    template<unsigned int srccomp>
    void define(BaseIrregData<CELL, data_t, srccomp>& a_input,
                unsigned int                        & a_comp)
    {
      BaseIrregData<CELL, data_t, ncomp>::define(a_input, a_comp);
    }
    
    inline void define(const Box& a_box,  const EBGraph& a_graph)
    {
      this->m_boxdefined = true;
      this->m_inputBox   = a_box;
      Box interBox = a_box & a_graph.getDomain();
      shared_ptr<vector<EBIndex<CELL> > > indicies(new vector< EBIndex<CELL> >() );
      *indicies = a_graph.getIrregVoFs(interBox);

      BaseIrregData<CELL, data_t, ncomp>::define(indicies);
    }

};


  ///template specialization for boundary data
  template < class data_t,  unsigned int ncomp>
  class IrregData <BOUNDARY, data_t, ncomp>: public BaseIrregData<BOUNDARY, data_t, ncomp>
  {
  public:
    inline void define(const Box& a_box, const EBGraph& a_graph)
    {
      this->m_boxdefined = true;
      this->m_inputBox   = a_box;
      Box interBox = a_box & a_graph.getDomain();

      shared_ptr<vector< EBIndex<BOUNDARY> > > indicies(new vector< EBIndex<BOUNDARY> >() );
      *indicies = a_graph.getCutFaces(interBox);
      BaseIrregData<BOUNDARY, data_t, ncomp>::define(indicies);

    }
  };

  ///template specialization for xface data
  template < class data_t,  unsigned int ncomp>
  class IrregData <XFACE, data_t, ncomp>: public BaseIrregData<XFACE, data_t, ncomp>
  {
  public:
    inline void define(const Box& a_box, const EBGraph& a_graph)
    {
      this->m_boxdefined = true;
      this->m_inputBox   = a_box;
      Box interBox = a_box & a_graph.getDomain();

      shared_ptr<vector< EBIndex<XFACE> > > indicies(new vector< EBIndex<XFACE> >() );
      *indicies = a_graph.getIrregXFaces(interBox);
      BaseIrregData<XFACE, data_t, ncomp>::define(indicies);

    }

  };

  ///template specialization for yface data
  template < class data_t,  unsigned int ncomp>
  class IrregData <YFACE, data_t, ncomp>: public BaseIrregData<YFACE, data_t, ncomp>
  {
  public:
    inline void define(const Box& a_box, const EBGraph& a_graph)
    {
      this->m_boxdefined = true;
      this->m_inputBox   = a_box;
      Box interBox = a_box & a_graph.getDomain();

      shared_ptr<vector< EBIndex<YFACE> > > indicies(new vector< EBIndex<YFACE> >() );
      *indicies = a_graph.getIrregYFaces(interBox);
      BaseIrregData<YFACE, data_t, ncomp>::define(indicies);

    }

  };

#if DIM==3
  ///template specialization for zface data
  template < class data_t,  unsigned int ncomp>
  class IrregData <ZFACE, data_t, ncomp>: public BaseIrregData<ZFACE, data_t, ncomp>
  {
  public:
    inline void define(const Box& a_box, const EBGraph& a_graph)
    {
      this->m_boxdefined = true;
      this->m_inputBox   = a_box;
      Box interBox = a_box & a_graph.getDomain();

      shared_ptr<vector< EBIndex<ZFACE> > > indicies(new vector< EBIndex<ZFACE> >() );
      *indicies = a_graph.getIrregZFaces(interBox);
      BaseIrregData<ZFACE, data_t, ncomp>::define(indicies);

    }
  };
#endif
}


#include "Proto_IrregDataImplem.H"

#endif
