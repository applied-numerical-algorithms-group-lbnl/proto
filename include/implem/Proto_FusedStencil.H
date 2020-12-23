#pragma once
#ifdef PROTO_CUDA


// should be superior to the number of stencil used.
static const int boxInfoSize = 32;

// used to store some data into the constant memory
struct boxInfo
{
  	int nx[boxInfoSize];
	int ny[boxInfoSize];
	int nz[boxInfoSize];
	Point bxLo[boxInfoSize];
        Point bxHi[boxInfoSize];
        Point dstShift[boxInfoSize];
        Point srcRef[boxInfoSize];
        Point dstRef[boxInfoSize];
};

// Aim: apply different pairs of stencil/box on one patch
template<typename T> class FusedStencil
{
  using Stencil = Proto::Stencil<T,MemType::DEVICE>; //This stencil is designed for the device only!
  using Box = Proto::Box;

  public:

  // send data into device memory 
  // m_stencils must be filled before
  void copyInfo();

  void define(std::vector<Stencil> a_stencils[], std::vector<Box> a_boxes[], unsigned int a_nDim);

  std::vector<Stencil*>& getStencils();

  // Apply the stencil kernel on several stencils and boxes on a same in patch and out patch
  // Be carefull that m_stencil and *a_box are filled with the same order
  template <unsigned int C, unsigned char D, unsigned char E>
  void cudaApplyFused( const Proto::BoxData<T,C,D,E,MemType::DEVICE> & a_src,
                             Proto::BoxData<T,C,D,E,MemType::DEVICE> & a_dst,
                             bool                a_initToZero,
                             T                   a_scale) const ;
 
  private:

#ifdef PROTO_CUDA
  // Used to call the close function on each stencil
  void close();
#endif
  // Should have the same size
  std::vector<Stencil*> m_stencils;  ///< to store stencil pointers 
  std::vector<Box> m_boxes;          ///< to store boxes 
  const T** m_coefs = nullptr;       ///< to store m_stencils[i]->devCoefs()
  int* m_ssize = nullptr;            ///< to store sizes of m_coefs and m_offsets
  const Point** m_offsets = nullptr; ///< to store m_sencils[i]->devOffsets
  std::size_t m_nbElemMax = 0;
  boxInfo m_info;
};

#ifdef PROTO_CUDA
  template<typename T>
  void FusedStencil<T>::close()
  {
    for(std::size_t it = 0 ; it < m_stencils.size() ; it++)
      m_stencils[it]->close(); // copy information to the device
  }
#endif

  template<typename T>
  void FusedStencil<T>::copyInfo()
  {
    int nsize = m_stencils.size();

    if(nsize == 0) return; 
    assert(nsize > 0);
#ifdef PROTO_CUDA
    close();
    protoMalloc(m_coefs, sizeof(const T*) * nsize);
    protoMalloc(m_offsets, sizeof(const Point*) * nsize);
    protoMalloc(m_ssize, sizeof(int) * nsize);
    int hostSsize[nsize];
    const T* t_coefs[nsize];
    const Point* t_offsets[nsize];
#else
    m_coefs = new const T*[nsize];
    m_offsets = new const Point*[nsize];
#endif
    for(int it = 0 ; it < nsize ; it++)
    {
      auto& sten =  m_stencils[it];
      auto& vec_coefs = sten->coefs();
      auto& vec_offsets = sten->offsets();
      const int csize = vec_coefs.size();
      const int osize = vec_offsets.size();
      assert(csize == osize);
#ifdef PROTO_CUDA
      hostSsize[it] = csize;
      t_coefs[it] = sten->devCoefs();
      t_offsets[it] = sten->devOffsets(); 
#else
      m_coefs[it] = new T[csize];
      m_offsets[it] = new Point[osize];
      std::copy(vec_coefs.data(), vec_coefs.data() + csize, m_coefs[it]);
      std::copy(vec_offsets.data(), vec_offsets.data() + osize, m_offsets[it]);
#endif
    }
#ifdef PROTO_CUDA
    unsigned int nbBitsCoefs   = sizeof(T*) * nsize;
    unsigned int nbBitsOffsets = sizeof(Point*) * nsize;
    unsigned int nbBitsSsize   = sizeof(int) * nsize;
    protoMemcpy(m_ssize,   hostSsize , nbBitsSsize   , protoMemcpyHostToDevice);
    protoMemcpy(m_coefs,   t_coefs   , nbBitsCoefs   , protoMemcpyHostToDevice);
    protoMemcpy(m_offsets, t_offsets , nbBitsOffsets , protoMemcpyHostToDevice);
#endif
  }
  // Return a reference on the member m_stencil - Used to fill this field
  template<typename T>
  std::vector<Proto::Stencil<T,MemType::DEVICE>*>& FusedStencil<T>::getStencils() 
  { 
	  return m_stencils; 
  }

  
#ifdef PROTO_CUDA
  template<typename T, size_t C>
  __global__
  void stencilIndexerFused( boxInfo a_info,  T* dstPtr, const T* srcPtr, 
		     std::size_t srcStride, std::size_t dstStride,
                     Proto::Point a_srcLo, Proto::Point a_srcHi, Proto::Point a_dstLo, Proto::Point a_dstHi,
                     bool a_initToZero,
		     const int nStencil,
                     int* a_ssize, const T** a_coefs,
                     const Proto::Point** a_srcIndices, T a_scale) 
  {
    const int id =  threadIdx.x+blockIdx.x*blockDim.x;
    __shared__ int s_offset[boxInfoSize][PR_MAX_COEFFS];

    const int tid = threadIdx.x ;
    if(tid < PR_MAX_COEFFS) 
    {
      for(int it = 0 ; it < nStencil ; it++)
        for(int itid = tid ; itid < a_ssize[it] && itid < PR_MAX_COEFFS ; itid+=blockDim.x)
	  s_offset[it][itid]= Proto::Box::soffset(a_srcIndices[it][itid], a_srcLo, a_srcHi); 
    }

    __syncthreads();

    for(int it = 0 ; it < nStencil ; it++)
    {
      unsigned int nx        = a_info.nx[it];
      unsigned int ny        = a_info.ny[it];
      unsigned int nz        = a_info.nz[it];
      Point& a_bxLo          = a_info.bxLo[it];
      Point& a_bxHi          = a_info.bxHi[it];
      Point& a_dstShi        = a_info.dstShift[it];
      Point& a_srcRef        = a_info.srcRef[it];
      Point& a_dstRef        = a_info.dstRef[it];
      unsigned int ssize     = a_ssize[it];
      const T * coefs        = a_coefs[it];
      const int * srcIndices = s_offset[it];

      Proto::stencil<T,C>(id,
			nx, ny, nz,
			dstPtr, srcPtr,
			srcStride, dstStride,
			a_srcLo, a_srcHi, 
			a_dstLo, a_dstHi, 
			a_bxLo, a_bxHi,
			a_initToZero,
			ssize, coefs,
			srcIndices, a_dstShi, 
			a_srcRef, a_dstRef,
			a_scale
		);
    }
  }
#endif

  template<typename T>
  void FusedStencil<T>::define(std::vector<Proto::Stencil<T,MemType::DEVICE>> a_stencils[], 
                               std::vector<Proto::Box>           a_boxes[],
			       unsigned int a_nDim)
  {
    for(int ibc = 0; ibc < a_nDim; ibc++)
    {
        assert(a_stencils[ibc].size() == a_boxes[ibc].size());
        for(int ivec = 0; ivec < a_stencils[ibc].size(); ivec++)
	{
	  m_stencils.push_back(&a_stencils[ibc][ivec]);
	  m_boxes.push_back(a_boxes[ibc][ivec]);
	}
     }

    unsigned int nStencil = m_stencils.size();

    if(nStencil == 0) return;

    assert(nStencil>0);
    assert(nStencil <= boxInfoSize);

    for(int it = 0 ; it < nStencil ; it++)
    {
      m_nbElemMax = std::max(m_nbElemMax, m_boxes[it].size());
      m_info.nx[it] = m_boxes[it].size(0);
      m_info.ny[it] = m_boxes[it].size(1);
#if DIM==3
      m_info.nz[it] = m_boxes[it].size(2);
#else
      m_info.nz[it] = 1;
#endif
      m_info.bxLo[it] = m_boxes[it].low();
      m_info.bxHi[it] = m_boxes[it].high();
      m_info.dstShift[it] = m_stencils[it]->destShift();
      m_info.srcRef[it] = m_stencils[it]->srcRatio();
      m_info.dstRef[it] = m_stencils[it]->destRatio();
    }

    copyInfo();

  }


  template<typename T>
  template <unsigned int C, unsigned char D, unsigned char E>
  void FusedStencil<T>::cudaApplyFused( const Proto::BoxData<T,C,D,E,MemType::DEVICE> & a_src,
					Proto::BoxData<T,C,D,E,MemType::DEVICE>  & a_dst,
					bool                a_initToZero,
					T                   a_scale) const
  {
    PR_TIME("cudaApply(Fused)");
    const int nStencil = m_stencils.size();
    // set value according to size boxes
    if(nStencil == 0) return;  
 
    assert(nStencil>0); 
    assert(nStencil <= boxInfoSize);

    int blocks=0;
    int strides=0;

#ifdef PROTO_CUDA
    v100tuning(m_nbElemMax, blocks, strides);

    assert(m_nbElemMax > 0);
    assert(blocks > 0);
    assert(strides> 0);

    for (int ee = 0; ee < E; ee++)
    {
      for (int dd = 0; dd < D; dd++)
      {
        const T* srcPtr = a_src.dataPtr(0, dd, ee);
        T* dstPtr = a_dst.dataPtr(0, dd, ee);
        protoStream_t curstream = Proto::DisjointBoxLayout::getCurrentStream();
        // launch my kernel
	protoLaunchKernelMemAsync((stencilIndexerFused<T,C>),
			blocks, strides,
			0, curstream,
			m_info, dstPtr, srcPtr,
			a_src.box().size(), a_dst.box().size(),
			a_src.box().low(),  a_src.box().high(),
			a_dst.box().low(),  a_dst.box().high(),
			a_initToZero,
			nStencil,
			m_ssize, m_coefs,
			m_offsets, a_scale);
      }
    }
#else
    std::cout << " the CPU version of cudaApplyFused is not implemented " << std::endl;
    std::abort();
#endif 
    unsigned long long int c = 0;
    for(int it = 0 ; it < nStencil ; it++)
      c += m_stencils[it]->numFlops(m_boxes[it]);
	    
    PR_FLOPS(c);
  }

#else

template<typename T> class FusedStencil
{

};
#endif
