#ifndef _Proto_IrregDataImplem_H_
#define _Proto_IrregDataImplem_H_

namespace Proto
{
  ////
  template <class data_t, class index_t, unsigned int ncomp>
  void 
  IrregData<data_t, index_t, ncomp>::
  define(const shared_ptr<vector<index_t> >& a_indicies)
  {
    m_indicies = a_indicies;
#ifdef PROTO_CUDA
    cudaMalloc(&m_rawPtr, size()*sizeof(T));
    cudaError err = cudaGetLastError();
    if (err != cudaSuccess)
    {
      fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
              __FILE__, __LINE__, cudaGetErrorString(err));
    }

    m_data = ::std::shared_ptr<T>(m_rawPtr, [](T* p){cudaFree(p);});
#else
    std::size_t count = this->size();
    m_data=::std::shared_ptr<T>(new T [count], [](T* p) { delete[] p;});
    m_rawPtr = m_data.get();
#endif
  }

  ////
  template <class data_t, class index_t, unsigned int ncomp>
  void 
  IrregData<data_t, index_t, ncomp>::
  setVal(const data_t& a_val)
  {
#ifdef PROTO_CUDA
    thrust::device_ptr<T> devptr = thrust::device_pointer_cast(m_rawPtr);
    unsigned int nsize = size();

    thrust::fill(thrust::device, devptr, devptr+nsize, a_val);

    cudaError err = cudaGetLastError();
    if (err != cudaSuccess)
    {
      fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
              __FILE__, __LINE__, cudaGetErrorString(err));
    }

#else
    for(int ival = 0; ival < size(); ival++)
    {
      m_rawPtr[ival] = a_val;
    }
#endif
  }

/// define function template specializations
  template <typename T, unsigned int C>
  IrregData<CELL, T, C>::
  define(const Box& a_validBox, const Point& a_grow, const EBGraph& a_graph)
  {
    Box grownBox = a_validBox.grow(a_grow) & a_graph.getDomain();
    shared_ptr<EBIndex<CELL> > indicies(new vector< EBIndex<CELL> >() );
    *indicies = a_graph.getIrregVoFs(grownBox);

    define(indicies);
  }

///
  template <typename T, unsigned int C>
  IrregData<BOUNDARY, T, C>::
  define(const Box& a_validBox, const Point& a_grow, const EBGraph& a_graph)
  {
    Box grownBox = a_validBox.grow(a_grow) & a_graph.getDomain();
    shared_ptr<EBIndex<BOUNDARY> > indicies(new vector< EBIndex<BOUNDARY> >() );
    *indicies = a_graph.getCutFaces(grownBox);
    define(indicies);
  }


///
  template <typename T, unsigned int C>
  IrregData<XFACE, T, C>::
  define(const Box& a_validBox, const Point& a_grow, const EBGraph& a_graph)(const Box& a_box, const EBGraph& a_graph)
  {
    Box grownBox = a_validBox.grow(a_grow) & a_graph.getDomain();
    shared_ptr<EBIndex<XFACE> > indicies(new vector< EBIndex<XFACE> >() );
    *indicies = a_graph.getIrregXFaces(grownBox);
    define(indicies);
  }


  template <typename T, unsigned int C>
  IrregData<YFACE, T, C>::
  define(const Box& a_validBox, const Point& a_grow, const EBGraph& a_graph)(const Box& a_box, const EBGraph& a_graph)
  {
    Box grownBox = a_validBox.grow(a_grow) & a_graph.getDomain();
    shared_ptr<EBIndex<YFACE> > indicies(new vector< EBIndex<YFACE> >() );
    *indicies = a_graph.getIrregYFaces(grownBox);
    define(indicies);
  }

#if DIM==3
  template <typename T, unsigned int C>
  IrregData<ZFACE, T, C>::
  define(const Box& a_validBox, const Point& a_grow, const EBGraph& a_graph)
  {
    Box grownBox = a_validBox.grow(a_grow) & a_graph.getDomain();
    shared_ptr<EBIndex<ZFACE> > indicies(new vector< EBIndex<ZFACE> >() );
    *indicies = a_graph.getIrregZFaces(grownBox);
    define(indicies);
  }
#endif



}
#endif

