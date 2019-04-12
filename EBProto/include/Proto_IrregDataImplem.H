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

}
#endif

