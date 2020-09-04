#pragma once

template<typename T>
class FusedStencil
{
  public:

#ifdef PROTO_CUDA
  void close()
  {
    for(auto& it = m_stencil)
	it->close(); // copy information to the device

  }
#endif

  void copyInfo()
  {
    int nsize = m_stencils.size();
#ifdef PROTO_CUDA
    close();
    protoMalloc(&m_coefs, sizeof(T*),nsize);
    protoMalloc(&m_offsets, sizeof(Point*) * nsize);
    protoMalloc(&m_ssize, sizeof(int) * nsize);
    T hostSsize[nsize];
#else
    m_coefs = new T*[nsize];
    m_offsets = new Point*[nsize];
#endif
    for(int it = 0 ; it < m_stencils.size() ; it++)
    {
      auto $ sten =  m_stencils[it];
      auto &vec_coefs = sten.coefs();
      auto &vec_offsets = sten.offsets();
      const int csize = vec_coefs.size();
      const int osize = vec_offsets.size();
      assert(csize == osize);
#ifdef PROTO_CUDA
      hostSsize[it] = csize;
      protoMalloc(&m_coefs[it], sizeof(T),nsize);
      protoMalloc(&m_offsets[it], sizeof(Point) * nsize);

      T* ptrDevCoefs = sten.devCoefs();
      Point* ptrDevOffsets = sten.devOffsets();

      unsigned int nbBitsCoefs = sizeof(T) * csize;
      unsigned int nbBitsOffsets = sizeof(Point) * osize;

      protoMemcpy(m_coefs[it], ptrDevCoefs, nbBitsCoefs, protoMemcpyDeviceToDevice);
      protoMemcpy(m_offsets[it], ptrDevOffsets, nbBitsOffsets, protoMemcpyDeviceToDevice);
#else
      m_coefs[it] = new T[csize];
      m_offsets[it] = new Point[osize];
      std::copy(vec_coefs.data(), vec_coefs.data() + csize, m_coefs[it]);
      std::copy(vec_offsets.data(), vec_offsets.data() + osize, m_offsets[it]);
#endif
    }
#ifdef PROTO_CUDA
    protoMemcpy(m_ssize, hostSsize,nsize *sizeof(int) protoMemcpyHostToDevice);
#endif
  }

  const int boxInfoSize = 12;

  struct boxInfo{
	int nx[boxInfoSize];
	int ny[boxInfoSize];
	int nz[boxInfoSize];
	Point bx[boxInfoSize];
  };

  template<typename T, size_t C>
  __device__
  void stencil(int a_nx, int a_ny, int a_nz, T* dstPtr, const T* srcPtr, Point srcStride, Point dstStride,
                     Point a_srcLo, Point a_srcHi, Point a_dstLo, Point a_dstHi,
                     Point a_dst, bool a_initToZero,
                     const int ssize, const T* a_coefs,
                     const int* srcIndices, T a_scale)
  {



  const int idx = id % a_nx;
#if DIM==3
  const int nxy = a_nx*a_ny;
  const int idy = (id % nxy) / a_nx;
  const int idz = id / nxy;
  if(idx < a_nx && idy < a_ny && idz < ksize)
#else
  const int idy = id / a_nx;
  if(idx < a_nx && idy < a_ny)
#endif
  {
    a_dst[0]+=idx;
    a_dst[1]+=idy;
#if DIM==3
    a_dst[2]+=idz;
#endif
    unsigned int s_off = Box::sindex(a_dst, a_srcLo, a_srcHi);
    unsigned int d_off = Box::sindex(a_dst, a_dstLo, a_dstHi);

    const double *src= srcPtr + s_off;
    double* dst = dstPtr + d_off;

    for(int c=0; c<C; c++, src +=  srcStride[0], dst += dstStride[0])
    {
      double accum=0;
#pragma unroll
      for(int s = 0; s<ssize; s++)
        accum+=a_coeff[s]*src[srcIndices[s]];

      if(a_initToZero) *dst = accum*a_scale;
      else *dst += accum*a_scale;
    }
  }

  }


  template<typename T, size_t C>
  __global__
  void stencilIndexerFused( boxInfo a_info T* dstPtr, const T* srcPtr, Point srcStride, Point dstStride,
                     Point a_srcLo, Point a_srcHi, Point a_dstLo, Point a_dstHi,
                     bool a_initToZero,
		     const int nStencil,
                     const int* a_ssize, const T** a_coefs,
                     const Point** a_srcIndices, T a_scale)
  {
    const int id =  threadIdx.x+blockIdx.x*blockDim.x;
    __shared__ int s_offset[boxInfoSize][PR_MAX_COEFFS];

    const int tid = threadIdx.x ;
    if( (tid < ssize)&&(tid < PR_MAX_COEFFS) )
    {
      for(int it = 0 ; it < nStencil ; it++)
        s_offset[it][tid]= Box::soffset(a_srcIndices[it][tid], a_srcLo, a_srcHi); 
    }

    __syncthreads();

    for(int it = 0 ; it < nStencil ; it++)
    {
      unsigned int nx = a_info.nx[it];
      unsigned int ny = a_info.ny[it];
      unsigned int nz = a_info.nz[it];
      Point a_dst = a_info.bx[it];
      unsigned int ssize = a_ssize[it];
      const T * coefs = a_coefs[it];
      const int * srcIndices = s_offset[it];
      stencil<T,C>( nx, ny, nz,
	dstPtr, srcPtr,
	srcStride, dstStride,
	a_srcLo, a_srcHi, 
	a_dstLo, a_dstHi, 
	a_dst,
	a_initToZero,
	ssize, coefs,
	srcIndices, a_scale);
	
    }
  }


  template<size_t C, size_t D, size_t E>
  void cudaApplyFused( const BoxData<T,C,D,E> & a_src,
                            BoxData<T,C,D,E>  & a_dst,
                            const Box         * a_box,
                            bool                a_initToZero,
                            T                   a_scale)
  {
    boxInfo<T> t_info;
    unsigned int nbElemMax = 0;
    const int nStencil = m_stencils.size;
    // set value according to size boxes
    for(int it = 0 ; it < nStencil ; it++)
    {
      nbElemMax = std::max(nbElemMax, bx[it].size());
      t_info.nx[it] = a_box[it].size(0);
      t_info.ny[it] = a_box[it].size(1);
      t_info.nz[it] = a_box[it].size(2);
      t_info.bx[it] = a_box[it].low();
    } 

    int blocks;
    int strides;
    v100tuning(nbElemMax, &blocs, &strides);

    for (int ee = 0; ee < E; ee++)
    {
      for (int dd = 0; dd < D; dd++)
      {
        const T* srcPtr = a_src.dataPtr(0, dd, ee);
        T* dstPtr = a_dst.dataPtr(0, dd, ee);
        protoStream_t curstream = DisjointBoxLayout::getCurrentStream();
        // launch my kernel
	ProtoLaunchKernelMemAsync(stencilIndexerFused<T,C>,
			blocks, strides,
			0, curstream,
			t_info, dstPtr, srcPtr,
			a_src.box().low(), a_src.box().high(),
			a_dst.box().low(), a_dst.box().high(),
			a_initToZero,
			nStencil,
			m_ssize, m_coefs,
			m_offsets, a_scale);

	  
      }
    }
  }

  private:
  std::vector <Stencil<T>*> m_stencils; 
  T** m_coefs;
  int* m_ssize;
  Point** m_offsets;

}
