#ifndef __BOXDATAIMPLEM__
#define __BOXDATAIMPLEM__


#pragma once

//=================================================================================================
// Constructors ||
//==============++

#ifdef PROTO_CUDA

//========================================================================
//  CUDA STUFF
//========================================================================
/** @defgroup Macros
    Proto macros
*/
/** @addtogroup Macros */
/*@{*/
#define PROTO_KERNEL_START inline __device__
#define PROTO_KERNEL_END(local_name, app_name)                  \
  __device__ decltype(&local_name) app_name = local_name;
#define PROTO_LAMBDA __device__
#else

#ifdef PROTO_OMP_FINE

#define PROTO_KERNEL_START #pragma omp declare simd notinbranch
#define PROTO_KERNEL_END(local_name, app_name)

#else

#define PROTO_KERNEL_START inline
#define PROTO_KERNEL_END(local_name, app_name) constexpr decltype(&local_name) app_name = local_name;
//#define PROTO_KERNEL_END(local_name, app_name) using app_name = local_name;
#define PROTO_LAMBDA 
/*@}*/

#endif //end if PROTO_OMP_FINE
#endif //end if PROTO_CUDA

// Default Constructor
template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>::BoxData()
  : m_rawPtr(NULL) {}

// Box Constructor
template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>::BoxData(const Box& a_box) 
{
  define(a_box);
}

// Initialization Constructor
template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>::BoxData(const Box& a_box, T a_init) 
{
  define(a_box);
  setVal(a_init);
}

// Deep Copy Constructor
template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>::BoxData(const BoxData<T,C,D,E>& a_src)
{
  PR_TIME("BoxData(BoxData&) (copy constructor)");
  //std::cout << "Calling BoxData copy constructor" << std::endl;
    
    
  if (!a_src.isAlias(*this))
  {    
    define(a_src.m_box);
    a_src.copyTo(*this);
  }
}

// Move Constructor
template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>::BoxData(BoxData<T,C,D,E>&& a_src)
  :   m_box(a_src.m_box),
      m_data(std::move(a_src.m_data)),
      m_rawPtr(a_src.m_rawPtr)
{
  //std::cout << "Calling BoxData move constructor" << std::endl;
}

// Define
template <class T, unsigned int C, unsigned char D, unsigned char E> 
void BoxData<T,C,D,E>::define(const Box& a_box)
{
  PR_TIME("BoxData::define(Box) (memory allocation)");
  m_box=a_box;

#ifdef PROTO_CUDA
  cudaMalloc(&m_rawPtr, size()*sizeof(T));
  m_data = ::std::shared_ptr<T>(m_rawPtr, [](T* p){cudaFree(p);});
#else
  std::size_t count = this->size();
  m_data=::std::shared_ptr<T>(new T [count], [](T* p) { delete[] p;});
  m_rawPtr = m_data.get();
#endif
}

// Destructor
template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>::~BoxData()
{}

// Slice Constructor (Private)
template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>::BoxData(::std::shared_ptr<T> a_data,const  T* a_ptr, const Box& a_box)
{
  m_data = a_data;
  m_box = (a_box);
  m_rawPtr = (T*)a_ptr;
}

// Raw Pointer Constructor
template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>::BoxData(const T* a_ptr, const Box& a_box, int a_ncomp)
{
  define(a_ptr, a_box, a_ncomp);
}

// Raw Pointer Define
template <class T, unsigned int C, unsigned char D, unsigned char E> 
void BoxData<T,C,D,E>::define(const T* a_ptr, const Box& a_box, int a_ncomp)
{
  PROTO_ASSERT((a_ncomp==C && D==1 && E==1), 
               "component mismatch in boxdata alias");
  T* castPtr = const_cast<T*>(a_ptr);
  m_data = ::std::shared_ptr<T>(castPtr, &(null_deleter_boxdata));
  m_rawPtr = castPtr;
  m_box = a_box;
}

template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>::BoxData(LazyStencil<T, C, D, E>&& a_op)
{
  define(a_op.m_range);
  a_op.apply(*this, true);
}
//=================================================================================================
// Data Movement ||
//===============++

// Copy Assignment Operator
template <class T, unsigned int C, unsigned char D, unsigned char E>
BoxData<T,C,D,E>& BoxData<T,C,D,E>::operator=(const BoxData<T,C,D,E>& a_src)
{
  PR_TIME("BoxData::operator=(BoxData&) (copy assign)");
  if (!a_src.isAlias(*this))
  {    
    m_box == a_src.m_box;
    define(a_src.m_box);
    a_src.copyTo(*this);
  }
  return *this;
}

// Move Assignment Operator
template <class T, unsigned int C, unsigned char D, unsigned char E>
BoxData<T,C,D,E>& BoxData<T,C,D,E>::operator=(BoxData<T,C,D,E>&& a_src)
{
  PR_TIME("BoxData::operator=(BoxData&& (move assign)");
    
  if (!a_src.isAlias(*this))
  {    
    std::swap<Box>(m_box,a_src.m_box);
    std::swap<std::shared_ptr<T>>(m_data, a_src.m_data);
    m_rawPtr = a_src.m_rawPtr;
    a_src.m_rawPtr = NULL;
  }
  return *this;
}

// Copy on Intersection
template <class T, unsigned int C, unsigned char D, unsigned char E> 
void BoxData<T,C,D,E>::copyTo(BoxData<T,C,D,E>& a_dest) const
{
  PR_TIME("BoxData::copyTo(BoxData&)"); 
  Box domain = m_box & a_dest.m_box;
  copyTo(a_dest,domain,{{0,C-1},{0,D-1},{0,E-1}},
         Point::Zeros(),{{0,C-1},{0,D-1},{0,E-1}});
}    

// Copy with Box Argument (and optional shift)
template <typename T, unsigned int C, unsigned char D, unsigned char E>
void BoxData<T,C,D,E>::copyTo(BoxData<T,C,D,E>&   a_dest,
                              const Box&              a_srcBox,
                              const Point&            a_destShift) const
{
  PR_TIME("BoxData::copyTo(BoxData&, Box&, Point&)");
  copyTo(a_dest,a_srcBox,{{0,C-1},{0,D-1},{0,E-1}},
         a_destShift,{{0,C-1},{0,D-1},{0,E-1}}); 
}

// Copy From
template <class T,unsigned int C, unsigned char D, unsigned char E>
template<unsigned int Csrc>
void BoxData<T,C,D,E>::copy(const BoxData<T,Csrc,D,E>&  a_src,
                            const Box&     a_srcBox,
                            unsigned int  a_srcComp,
                            const Box&     a_destBox,
                            unsigned int  a_destComp,
                            unsigned int  a_numcomps)
{
  a_src.copyTo(*this, a_srcBox, {a_srcComp,a_srcComp+a_numcomps-1},
               a_destBox.low()-a_srcBox.low(),
               {a_destComp, a_destComp+a_numcomps-1});
}


#ifdef PROTO_CUDA

enum BoxDataOp
{
  Invalid = -1,
  Copy    = 0,
  Add     = 1,
  Subtract, 
  Multiply, 
  Divide, 
  NUM_OP_TYPES
};



template <typename T>
__global__
void boxdataIndexer(unsigned int  a_begin, 
                    unsigned int a_end,
                    T* a_src,  T* a_dst,
                    Box a_box, Box a_srcBox, Box a_dstBox,
                    Point  a_dstShift, BoxDataOp a_op)
{
  unsigned int idx = a_begin + threadIdx.x+blockIdx.x*blockDim.x;

  if((idx >= a_begin) && (idx< a_end))
  {

    Point boxPt = a_box[idx];
    Point srcPt = boxPt;
    Point dstPt = boxPt + a_dstShift;

    unsigned int dstOff = a_dstBox.index(dstPt);
    unsigned int srcOff = a_srcBox.index(srcPt);

    if(a_op == BoxDataOp::Copy)
    {
      a_dst[dstOff] = a_src[srcOff];
    }
    else if(a_op == BoxDataOp::Add)
    {
      a_dst[dstOff] += a_src[srcOff];
    }
    else if(a_op == BoxDataOp::Subtract)
    {
      a_dst[dstOff] -= a_src[srcOff];
    }
    else if(a_op == BoxDataOp::Multiply)
    {
      a_dst[dstOff] *= a_src[srcOff];
    }
    else if(a_op == BoxDataOp::Divide)
    {
      a_dst[dstOff] /= a_src[srcOff];
    }
    else
    {
      printf("boxdataIndexer error: bogus boxdata op input!!!\n");
    }
  }
}


template <typename T>
__global__
void scalarIndexer(int  a_begin, int a_end, T* a_dst, T a_scalar, Box a_box, BoxDataOp a_op)
{
  unsigned int idx = a_begin + threadIdx.x+blockIdx.x*blockDim.x;
  if((idx >= a_begin) && (idx< a_end))
  {

    Point boxPt = a_box[idx];
    unsigned int dstOff = a_box.index(boxPt);

    if(a_op == BoxDataOp::Copy)
    {
      a_dst[dstOff] = a_scalar;
    }
    else if(a_op == BoxDataOp::Add)
    {
      a_dst[dstOff] += a_scalar;
    }
    else if(a_op == BoxDataOp::Subtract)
    {
      a_dst[dstOff] -= a_scalar;
    }
    else if(a_op == BoxDataOp::Multiply)
    {
      a_dst[dstOff] *= a_scalar;
    }
    else if(a_op == BoxDataOp::Divide)
    {
      a_dst[dstOff] /= a_scalar;
    }
    else
    {
      printf("scalarIndexer error: bogus boxdata op input!!!\n");
    }
  }
}
#endif
// General Copy
template< typename T, unsigned int C, unsigned char D, unsigned char E>
template< unsigned int  Cdest,
          unsigned char Ddest,
          unsigned char Edest>
void BoxData<T,C,D,E>::copyTo(BoxData<T,Cdest,Ddest,Edest>& a_dst,
                              const Box&     a_srcBox,
                              CInterval     a_srcComps,
                              const Point&  a_dstShift,
                              CInterval     a_dstComps) const
{
  PR_TIME("boxdata::copyto");
#ifdef PROTO_MEM_CHECK
  memcheck::numcopies++;
#endif
  PROTO_ASSERT(contains(a_srcComps),
               "BoxData::copyTo(...) invalid.\
    Source CInterval must be contained in source data component space.");
  PROTO_ASSERT(a_dst.contains(a_dstComps),
               "BoxData::copyTo(...) invalid.\
    Destination CInterval must be contained in destination data component space.");
  for (int ii = 0; ii < DIM; ii++)
  {
    PROTO_ASSERT(a_srcComps.size(ii) == a_dstComps.size(ii),
                 "BoxData::copyTo(...) invalid.\
      Src interval size %i =/= dst interval size %i on axis %i.",
                 a_srcComps.size(ii), a_dstComps.size(ii),ii);
  }
  PROTO_ASSERT(m_box.contains(a_srcBox),
               "BoxData::copyTo(...) invalid.\
    Source Box input is not a subset of source array domain.");
  PROTO_ASSERT(a_dst.box().contains(a_srcBox.shift(a_dstShift)),
               "BoxData::copyTo(...) invalid.\
    Shifted destination Box input is not a subset of destination array domain.");
     

  if (a_srcBox.size() <= 0){return;}
  int doFastCopy = true;
  doFastCopy &= (m_box == a_srcBox);
  doFastCopy &= (a_dst.box() == a_srcBox.shift(a_dstShift));
  doFastCopy &= (C == Cdest);
  doFastCopy &= (D == Ddest);
  doFastCopy &= (E == Edest);

#ifdef PROTO_CUDA
  if (doFastCopy)
  {
    unsigned int nsize = size();
    T* castSrc = (T*) m_rawPtr;
    thrust::device_ptr<T> devptrSrc = thrust::device_pointer_cast(castSrc);
    T* dstPtr = (T*)(a_dst.data());
    thrust::device_ptr<T> devptrDst = thrust::device_pointer_cast(dstPtr);
    thrust::copy(thrust::device, devptrSrc, devptrSrc+nsize, devptrDst);
  }
  else
  {
    for (unsigned int ee = 0; ee < a_srcComps.size(2); ee++)
    {
      for (unsigned int dd = 0; dd < a_srcComps.size(1); dd++)
      {
        for (unsigned int cc = 0; cc < a_srcComps.size(0); cc++)
        {
          int ccSrc = cc + a_srcComps.low(0);
          int ddSrc = dd + a_srcComps.low(1);
          int eeSrc = ee + a_srcComps.low(2);
          int ccDst = cc + a_dstComps.low(0);     // 
          int ddDst = dd + a_dstComps.low(1);
          int eeDst = ee + a_dstComps.low(2);


          BoxData<T,1,1,1> srcbd = slice(*this, ccSrc, ddSrc, eeSrc);
          BoxData<T,1,1,1> dstbd = slice(a_dst, ccDst, ddDst, eeDst);

          unsigned int begin  = 0;
          unsigned int end    = a_srcBox.size();

          unsigned int stride = a_srcBox.size(0);
          unsigned int blocks = a_srcBox.flatten(0).size();

          T* d_srcPtr = (T*)srcbd.data();
          T* d_dstPtr = (T*)dstbd.data();

          {
            PR_TIME("boxindexer");
            boxdataIndexer<<<blocks, stride>>>(begin, end, d_srcPtr, d_dstPtr, 
                                               a_srcBox, srcbd.box(), dstbd.box(),
                                               a_dstShift, BoxDataOp::Copy);
          }

          cudaError err = cudaGetLastError();
          if (err != cudaSuccess)
          {
            fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                    __FILE__, __LINE__, cudaGetErrorString(err));
          }
        }
      }
    }



  }
#else

  if (doFastCopy)
  {
    for (int ii = 0; ii < size(); ii++)
    {
      (*a_dst[ii]) = m_rawPtr[ii];
    }
  }
  else
  {
    Box cross = a_srcBox.flatten(0);
    int pencilLength = a_srcBox.size(0);
    const T* srcVal;
    T* destVal;
    for (unsigned int ee = 0; ee < a_srcComps.size(2); ee++)
    {
      for (unsigned int dd = 0; dd < a_srcComps.size(1); dd++)
      {
        for (unsigned int cc = 0; cc < a_srcComps.size(0); cc++)
        {
          for (auto iter = cross.begin(); iter != cross.end(); ++iter)
          {
            srcVal =  this->data(*iter,
                                 a_srcComps.low(0) + cc,
                                 a_srcComps.low(1) + dd,
                                 a_srcComps.low(2) + ee);
            destVal = a_dst.data((*iter) + a_dstShift,
                                 a_dstComps.low(0) + cc,
                                 a_dstComps.low(1) + dd,
                                 a_dstComps.low(2) + ee);
            for (int ii = 0; ii < pencilLength; ii++)
            {
              destVal[ii] = srcVal[ii];
            }
          }
        }
      }
    }
  }
#endif
}


// Index Accessor (Non-Const)
template <class T, unsigned int C, unsigned char D, unsigned char E>
T* BoxData<T,C,D,E>::operator[](unsigned int a_index)
{
  PROTO_ASSERT(a_index < size(),
               "BoxData::operator[](uint a_index) invalid for a_index = %i. \
    a_index must be in [0,this->size() = %i)",a_index, (int)size());
    
  return &m_rawPtr[a_index];
}

// Index Accessor (Const)
template <class T, unsigned int C, unsigned char D, unsigned char E>
const T* BoxData<T,C,D,E>::operator[](unsigned int a_index) const
{
  PROTO_ASSERT(a_index < size(),
               "BoxData::operator[](uint a_index) invalid for a_index = %i.\
    a_index must be in [0,this->size() = %i)",a_index, (int)size());
  return &m_rawPtr[a_index];
}


// Create Pointwise Variable (Non-Const)
template <class T, unsigned int C, unsigned char D, unsigned char E>
Var<T,C,D,E> BoxData<T,C,D,E>::var(const Point& a_pt)
{
//  PR_TIME("BoxData::var");
  PROTO_ASSERT(m_box.contains(a_pt),
               "BoxData::var(Point a_pt) invalid for a_pt not in this->box()");
  Var<T,C,D,E> rtn;
  for(unsigned int ee = 0; ee < E; ++ee)
    for(unsigned int dd = 0; dd < D; ++dd)
      for(unsigned int cc = 0; cc < C; ++cc)
      {
        rtn.m_ptrs[cc + C*dd + C*D*ee] = m_rawPtr + index(a_pt,cc,dd,ee);
      }
  return rtn;
}

// Create Pointwise Variable (Const)
template <class T, unsigned int C, unsigned char D, unsigned char E>
Var<T,C,D,E> BoxData<T,C,D,E>::var(const Point& a_pt) const
{
//  PR_TIME("BoxData::var_const");
  PROTO_ASSERT(m_box.contains(a_pt),
               "BoxData::var(Point a_pt) invalid for a_pt not in this->box()");
  Var<T,C,D,E> rtn;
  for(unsigned int ee = 0; ee < E; ++ee)
    for(unsigned int dd = 0; dd < D; ++dd)
      for(unsigned int cc = 0; cc < C; ++cc)
      {
        rtn.m_ptrs[cc + C*dd + C*D*ee] = m_rawPtr + index(a_pt,cc,dd,ee);
      }
  return rtn;
}

//=================================================================================================
// Algebraic Operations ||
//======================++

// Operator |= (BoxData&)
/*
  template <class T, unsigned int C, unsigned char D, unsigned char E> 
  BoxData<T,C,D,E>& BoxData<T,C,D,E>::operator|=(const BoxData<T,C,D,E>& a_rhs)
  {
  a_rhs.copyTo(*this);
  return *this;
  }
*/

// Operator += (BoxData&)
template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>& BoxData<T,C,D,E>::operator+=(const BoxData<T,C,D,E>& a_rhs)
{
  PR_TIME("BoxData::operator+=");
  unsigned long long int numflops = size();
  PR_FLOPS(numflops);
#ifdef PROTO_CUDA

  if(m_box == a_rhs.m_box)
  {
    T* castPtrInp = (T*)(a_rhs.m_rawPtr);
    thrust::device_ptr<T> devptrLoc = thrust::device_pointer_cast(m_rawPtr);
    thrust::device_ptr<T> devptrInp = thrust::device_pointer_cast(castPtrInp);
    unsigned int nsize = size();

    thrust::transform(devptrLoc, devptrLoc+nsize, devptrInp, devptrLoc, thrust::plus<T>());
  }
  else
  {

    Box domain = m_box & a_rhs.m_box;
    if (domain.size() <= 0){return *this;}

    for (unsigned int ee = 0; ee < E; ee++)
      for (unsigned int dd = 0; dd < D; dd++)
        for (unsigned int cc = 0; cc < C ; cc++)
        {
          BoxData<T,1,1,1> srcbd = slice(a_rhs, cc, dd, ee);
          BoxData<T,1,1,1> dstbd = slice(*this, cc, dd, ee);

          unsigned int begin  = 0;
          unsigned int end    = domain.size();

          unsigned int stride = domain.size(0);
          unsigned int blocks = domain.flatten(0).size();

          T* d_srcPtr = (T*)srcbd.data();
          T* d_dstPtr = (T*)dstbd.data();

          boxdataIndexer<<<blocks, stride>>>(begin, end, d_srcPtr, d_dstPtr, 
                                             domain, srcbd.box(), dstbd.box(),
                                             Point::Zeros(), BoxDataOp::Add);

          cudaError err = cudaGetLastError();
          if (err != cudaSuccess)
          {
            fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                    __FILE__, __LINE__, cudaGetErrorString(err));
          }
        }



  }
#else
  if(m_box == a_rhs.m_box)
  {
    const unsigned int n = this->size();
    PROTO_ASSERT(n == a_rhs.size(),
                 "Something has gone horribly wrong.\
        Two BoxData<T,C,D,E> with the same m_box have different sizes.");
    for (int ii = 0; ii < n; ii++)
    {
      m_rawPtr[ii] += *a_rhs[ii];
    }
  }
  else
  {
    Box domain = m_box & a_rhs.m_box;
    if (domain.size() <= 0){return *this;}

    Box cross = domain.flatten(0);
    unsigned int npencil = domain.size(0);

    T* thisVal;
    const T* rhsVal;
    for(unsigned char e=0; e<E; e++)
      for(unsigned char d=0; d<D; d++)
        for(unsigned int c=0; c<C; c++)
        {
          for (auto iter = cross.begin(); iter != cross.end(); ++iter)
          {
            thisVal = this->data(*iter,c,d,e);
            rhsVal = a_rhs.data(*iter,c,d,e);
            for (int ii = 0; ii < npencil; ii++)
            {
              thisVal[ii] += rhsVal[ii];
            }
          }
        }
  }
#endif    
  return *this;
}

// Operator -= (BoxData&)
template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>& BoxData<T,C,D,E>::operator-=(const BoxData<T,C,D,E>& a_rhs)
{
  PR_TIME("BoxData::operator-=(BoxData&)");
  unsigned long long int numflops = size();
  PR_FLOPS(numflops);
#ifdef PROTO_CUDA
  if(m_box == a_rhs.m_box)
  {
    T* castPtrInp = (T*)(a_rhs.m_rawPtr);
    thrust::device_ptr<T> devptrLoc = thrust::device_pointer_cast(m_rawPtr);
    thrust::device_ptr<T> devptrInp = thrust::device_pointer_cast(castPtrInp);
    unsigned int nsize = size();

    thrust::transform(devptrLoc, devptrLoc+nsize, devptrInp, devptrLoc, thrust::minus<T>());
  }
  else
  {
    Box domain = m_box & a_rhs.m_box;
    if (domain.size() <= 0){return *this;}

    for (unsigned int ee = 0; ee < E; ee++)
      for (unsigned int dd = 0; dd < D ; dd++)
        for (unsigned int cc = 0; cc < C; cc++)
        {

          BoxData<T,1,1,1> srcbd = slice(a_rhs, cc, dd, ee);
          BoxData<T,1,1,1> dstbd = slice(*this, cc, dd, ee);

          unsigned int begin  = 0;
          unsigned int end    = domain.size();
          unsigned int stride = domain.size(0);
          unsigned int blocks = domain.flatten(0).size();

          T* d_srcPtr = (T*)srcbd.data();
          T* d_dstPtr = (T*)dstbd.data();


          boxdataIndexer<<<blocks, stride>>>(begin, end, d_srcPtr, d_dstPtr, 
                                             domain, srcbd.box(), dstbd.box(),
                                             Point::Zeros(), BoxDataOp::Subtract);

          cudaError err = cudaGetLastError();
          if (err != cudaSuccess)
          {
            fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                    __FILE__, __LINE__, cudaGetErrorString(err));
          }
        }


    
  }
#else
    
  if(m_box == a_rhs.m_box)
  {
    const unsigned int n = this->size();
    PROTO_ASSERT(n == a_rhs.size(),
                 "Something has gone horribly wrong.\
        Two BoxData<T,C,D,E> with the same m_box have different sizes.");
    for (int ii = 0; ii < n; ii++)
    {
      m_rawPtr[ii] -= *a_rhs[ii];
    }
  }
  else
  {
    Box domain = m_box & a_rhs.m_box;
    if (domain.size() <= 0){return *this;}

    Box cross = domain.flatten(0);
    unsigned int npencil = domain.size(0);

    T* thisVal;
    const T* rhsVal;
    for(unsigned char e=0; e<E; e++)
      for(unsigned char d=0; d<D; d++)
        for(unsigned int c=0; c<C; c++)
        {
          for (auto iter = cross.begin(); iter != cross.end(); ++iter)
          {
            thisVal = this->data(*iter,c,d,e);
            rhsVal =  a_rhs.data(*iter,c,d,e);
            for (int ii = 0; ii < npencil; ii++)
            {
              thisVal[ii] -= rhsVal[ii];
            }
          }
        }
  }
#endif
  return *this;
}

// Operator *= (BoxData&)
template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>& BoxData<T,C,D,E>::operator*=(const BoxData<T,C,D,E>& a_rhs)
{
    
  PR_TIME("BoxData::operator*=");
  unsigned long long int numflops = size();
  PR_FLOPS(numflops);
#ifdef PROTO_CUDA
  if(m_box == a_rhs.m_box)
  {
    T* castPtrInp = (T*)(a_rhs.m_rawPtr);
    thrust::device_ptr<T> devptrLoc = thrust::device_pointer_cast(m_rawPtr);
    thrust::device_ptr<T> devptrInp = thrust::device_pointer_cast(castPtrInp);
    unsigned int nsize = size();

    thrust::transform(devptrLoc, devptrLoc+nsize, devptrInp, devptrLoc, thrust::multiplies<T>());
  }
  else
  {
    Box domain = m_box & a_rhs.m_box;
    if (domain.size() <= 0){return *this;}

    for (unsigned int ee = 0; ee < E; ee++)
      for (unsigned int dd = 0; dd < D ; dd++)
        for (unsigned int cc = 0; cc < C; cc++)
        {

          BoxData<T,1,1,1> srcbd = slice(a_rhs, cc, dd, ee);
          BoxData<T,1,1,1> dstbd = slice(*this, cc, dd, ee);

          unsigned int begin  = 0;
          unsigned int end    = domain.size();
          unsigned int stride = domain.size(0);
          unsigned int blocks = domain.flatten(0).size();

          T* d_srcPtr = (T*)srcbd.data();
          T* d_dstPtr = (T*)dstbd.data();

          boxdataIndexer<<<blocks, stride>>>(begin, end, d_srcPtr, d_dstPtr, 
                                             domain, srcbd.box(), dstbd.box(),
                                             Point::Zeros(), BoxDataOp::Multiply);

          cudaError err = cudaGetLastError();
          if (err != cudaSuccess)
          {
            fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                    __FILE__, __LINE__, cudaGetErrorString(err));
          }
        }


  }
#else
  if(m_box == a_rhs.m_box)
  {
    const unsigned int n = this->size();
    PROTO_ASSERT(n == a_rhs.size(),
                 "Something has gone horribly wrong.\
        Two BoxData<T,C,D,E> with the same m_box have different sizes.");
    for (int ii = 0; ii < n; ii++)
    {
      m_rawPtr[ii] *= *a_rhs[ii];
    }
  }
  else
  {
    Box domain = m_box & a_rhs.m_box;
    if (domain.size() <= 0){return *this;}

    Box cross = domain.flatten(0);
    unsigned int npencil = domain.size(0);

    T* thisVal;
    const T* rhsVal;
    for(unsigned char e=0; e<E; e++)
      for(unsigned char d=0; d<D; d++)
        for(unsigned int c=0; c<C; c++)
        {
          for (auto iter = cross.begin(); iter != cross.end(); ++iter)
          {
            thisVal = this->data(*iter,c,d,e);
            rhsVal =  a_rhs.data(*iter,c,d,e);
            for (int ii = 0; ii < npencil; ii++)
            {
              thisVal[ii] *= rhsVal[ii];
            }
          }
        }
  }
#endif
  return *this;
}

// Operator /= (BoxData&)
template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>& BoxData<T,C,D,E>::operator/=(const BoxData<T,C,D,E>& a_rhs)
{
    
  PR_TIME("BoxData::operator/=");
  unsigned long long int numflops = size();
  PR_FLOPS(numflops);
#ifdef PROTO_CUDA
  if(m_box == a_rhs.m_box)
  {
    T* castPtrInp = (T*)(a_rhs.m_rawPtr);
    thrust::device_ptr<T> devptrLoc = thrust::device_pointer_cast(m_rawPtr);
    thrust::device_ptr<T> devptrInp = thrust::device_pointer_cast(castPtrInp);
    unsigned int nsize = size();

    thrust::transform(devptrLoc, devptrLoc+nsize, devptrInp, devptrLoc, thrust::divides<T>());
  }
  else
  {
    Box domain = m_box & a_rhs.m_box;
    if (domain.size() <= 0){return *this;}

    for (unsigned int ee = 0; ee < E; ee++)
      for (unsigned int dd = 0; dd < D ; dd++)
        for (unsigned int cc = 0; cc < C; cc++)
        {

          BoxData<T,1,1,1> srcbd = slice(a_rhs, cc, dd, ee);
          BoxData<T,1,1,1> dstbd = slice(*this, cc, dd, ee);

          unsigned int begin  = 0;
          unsigned int end    = domain.size();
          unsigned int stride = domain.size(0);
          unsigned int blocks = domain.flatten(0).size();

          T* d_srcPtr = (T*)srcbd.data();
          T* d_dstPtr = (T*)dstbd.data();

          boxdataIndexer<<<blocks, stride>>>(begin, end, d_srcPtr, d_dstPtr, 
                                             domain, srcbd.box(), dstbd.box(),
                                             Point::Zeros(), BoxDataOp::Divide);

          cudaError err = cudaGetLastError();
          if (err != cudaSuccess)
          {
            fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                    __FILE__, __LINE__, cudaGetErrorString(err));
          }
        }


  }
#else
  if(m_box == a_rhs.m_box)
  {
    const unsigned int n = this->size();
    PROTO_ASSERT(n == a_rhs.size(),
                 "Something has gone horribly wrong.\
        Two BoxData<T,C,D,E> with the same m_box have different sizes.");
    for (int ii = 0; ii < n; ii++)
    {
      PROTO_ASSERT(a_rhs[ii] != 0,
                   "BoxData::operator/=(BoxData&) invalid. Divide by zero at index %i.",ii);
      m_rawPtr[ii] /= *a_rhs[ii];
    }
  }
  else
  {
    Box domain = m_box & a_rhs.m_box;
    if (domain.size() <= 0){return *this;}

    Box cross = domain.flatten(0);
    unsigned int npencil = domain.size(0);

    T* thisVal;
    const T* rhsVal;
    for(unsigned char e=0; e<E; e++)
      for(unsigned char d=0; d<D; d++)
        for(unsigned int c=0; c<C; c++)
        {
          for (auto iter = cross.begin(); iter != cross.end(); ++iter)
          {
            thisVal = this->data(*iter,c,d,e);
            rhsVal =  a_rhs.data(*iter,c,d,e);
            for (int ii = 0; ii < npencil; ii++)
            {
              PROTO_ASSERT(a_rhs[ii] != 0,
                           "BoxData::operator/=(BoxData&) invalid.\
                    Divide by zero at index %i.",a_rhs.index(*iter,c,d,e));
              thisVal[ii] /= rhsVal[ii];
            }
          }
        }
  }
#endif
  return *this;
}


template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>& BoxData<T,C,D,E>::operator+=(T a_scale)
{
  PR_TIME("boxdatta += scalar");
  unsigned long long int numflops = size();
  PR_FLOPS(numflops);
#ifdef PROTO_CUDA

  for (unsigned int ee = 0; ee < E; ee++)
    for (unsigned int dd = 0; dd < D; dd++)
      for (unsigned int cc = 0; cc < C; cc++)
      {

        BoxData<T,1,1,1> dstbd = slice(*this, cc, dd, ee);

        unsigned int begin  = 0;
        unsigned int end    = m_box.size();

        unsigned int stride = m_box.size(0);
        unsigned int blocks = m_box.flatten(0).size();

        T* d_dstPtr = (T*)dstbd.data();

        scalarIndexer<<<blocks, stride>>>(begin, end, d_dstPtr, a_scale,
                                          dstbd.box(), BoxDataOp::Add); 

        cudaError err = cudaGetLastError();
        if (err != cudaSuccess)
        {
          fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                  __FILE__, __LINE__, cudaGetErrorString(err));
        }
      }



#else  
  for (int ii = 0; ii < this->size(); ii++)
  {
    m_rawPtr[ii] += a_scale;
  }
#endif
  return *this;
}

template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>& BoxData<T,C,D,E>::operator-=(T a_scale)
{
  PR_TIME("boxdatta -= scalar");
  unsigned long long int numflops = size();
  PR_FLOPS(numflops);
#ifdef PROTO_CUDA
  for (unsigned int ee = 0; ee < E; ee++)
    for (unsigned int dd = 0; dd < D; dd++)
      for (unsigned int cc = 0; cc < C; cc++)
      {

        BoxData<T,1,1,1> dstbd = slice(*this, cc, dd, ee);

        unsigned int begin  = 0;
        unsigned int end    = m_box.size();
        unsigned int stride = m_box.size(0);
        unsigned int blocks = m_box.flatten(0).size();

        T* d_dstPtr = (T*)dstbd.data();

        scalarIndexer<<<blocks, stride>>>(begin, end, d_dstPtr, a_scale,
                                          dstbd.box(), BoxDataOp::Subtract); 

        cudaError err = cudaGetLastError();
        if (err != cudaSuccess)
        {
          fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                  __FILE__, __LINE__, cudaGetErrorString(err));

        }
      }


//

#else
  for (int ii = 0; ii < this->size(); ii++)
  {
    m_rawPtr[ii] -= a_scale;
  }
#endif
  return *this;
}

template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>& BoxData<T,C,D,E>::operator*=(T a_scale)
{
  PR_TIME("BoxData::operator*=(T)");
  unsigned long long int numflops = size();
  PR_FLOPS(numflops);
#ifdef PROTO_CUDA

  for (unsigned int ee = 0; ee < E; ee++)
    for (unsigned int dd = 0; dd < D; dd++)
      for (unsigned int cc = 0; cc < C; cc++)
      {

        BoxData<T,1,1,1> dstbd = slice(*this, cc, dd, ee);

        unsigned int begin  = 0;
        unsigned int end    = m_box.size();
        unsigned int stride = m_box.size(0);
        unsigned int blocks = m_box.flatten(0).size();

        T* d_dstPtr = (T*)dstbd.data();

        scalarIndexer<<<blocks, stride>>>(begin, end, d_dstPtr, a_scale,
                                          dstbd.box(), BoxDataOp::Multiply); 

        cudaError err = cudaGetLastError();
        if (err != cudaSuccess)
        {
          fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                  __FILE__, __LINE__, cudaGetErrorString(err));
        }
      }

//
#else  
    
  for (int ii = 0; ii < this->size(); ii++)
  {
    m_rawPtr[ii] *= a_scale;
  }
#endif
  return *this;
}

template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>& BoxData<T,C,D,E>::operator/=(T a_scale)
{
  PR_TIME("BoxData::operator*=(T)");
  unsigned long long int numflops = size();
  PR_FLOPS(numflops);

  PROTO_ASSERT(a_scale != 0,
               "BoxData::operator/=(T a_scale) invalid for a_scale = 0. Division by 0");
#ifdef PROTO_CUDA

  for (unsigned int ee = 0; ee < E; ee++)
    for (unsigned int dd = 0; dd < D; dd++)
      for (unsigned int cc = 0; cc < C; cc++)
      {
        BoxData<T,1,1,1> dstbd = slice(*this, cc, dd, ee);

        unsigned int begin  = 0;
        unsigned int end    = m_box.size();
        unsigned int stride = m_box.size(0);
        unsigned int blocks = m_box.flatten(0).size();

        T* d_dstPtr = (T*)dstbd.data();
        Point dstLo = dstbd.box().low();
        Point dstHi = dstbd.box().high();

        scalarIndexer<<<blocks, stride>>>(begin, end, d_dstPtr, a_scale,
                                          dstbd.box(), BoxDataOp::Divide); 

        cudaError err = cudaGetLastError();
        if (err != cudaSuccess)
        {
          fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                  __FILE__, __LINE__, cudaGetErrorString(err));
        }
      }




#else  
  for (int ii = 0; ii < this->size(); ii++)
  {
    m_rawPtr[ii] /= a_scale;
  }
#endif
  return *this;
}

//=================================================================================================
// Utility Functions ||
//===================++



template <class T, class B>
inline void setBoxDataVal(BoxData<T>& a_data, const B& a_box, T a_val)
{
#ifdef PROTO_CUDA

  std::cout << "should not get here with cuda" << std::endl;
  abort();
  
#else
//  std::cout << "setBoxDataVal2" << std::endl;
  Box intersect = a_box & a_data.box();
  Box face = intersect.flatten(0);
  int npencil = intersect.size(0);
        
  for (auto iter = face.begin(); iter != face.end(); ++iter)
  {
    T* ptr = a_data.data(*iter);
    for (int ii = 0; ii < npencil; ii++)
    {
      ptr[ii] = a_val;
    }
  }
#endif
}


template <class T, unsigned int C, unsigned char D, unsigned char E>
void BoxData<T,C,D,E>::setVal(const T& a_val)
{

  PR_TIME("BoxData::setval");
#ifdef PROTO_CUDA
//  std::cout << "setValCuda" << std::endl;

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
//  std::cout << "setValNonCuda" << std::endl;
  int iset = 0;
  for (unsigned int ee = 0; ee < E; ee++)
  {
    for (unsigned int dd = 0; dd < D; dd++)
    {
      for (unsigned int cc = 0; cc < C; cc++)
      {
        auto slicebd = slice(*this, cc, dd, ee);
        setBoxDataVal<T,Box>(slicebd, box(), a_val);
        iset++;
      }
    }
  }
#endif
}

// Set Values Within Box
template <class T, unsigned int C, unsigned char D, unsigned char E>
void BoxData<T,C,D,E>::setVal(const T&   a_val,
                              const Box&  a_box)
{
  for (unsigned int ee = 0; ee < E; ee++)
  {
    for (unsigned int dd = 0; dd < D; dd++)
    {
      for (unsigned int cc = 0; cc < C; cc++)
      {
        auto slicebd = slice(*this, cc, dd, ee);
        setBoxDataVal<T>(slicebd, a_box, a_val);
      }
    }
  }
}

// Set Component Values Within Box
template <class T, unsigned int C, unsigned char D, unsigned char E>
void BoxData<T,C,D,E>::setVal(const T& a_val,
                              const Box& a_box,
                              int a_c,
                              int a_d,
                              int a_e)
{
    
  PR_TIME("BoxData::setVal(T&, Box, int, int, int)");
  PROTO_ASSERT((a_c < C),
               "BoxData::setVal(..., int a_c, int a_d, int a_e) invalid for a_c = %i.\
     a_c must be in [0,%i=C).",
               a_c, DIM);
  PROTO_ASSERT((a_d < D),
               "BoxData::setVal(..., int a_c, int a_d, int a_e) invalid for a_d = %i. \
    a_d must be in [0,%i=D).",
               a_d, DIM);
  PROTO_ASSERT((a_e < E),
               "BoxData::setVal(..., int a_c, int a_d, int a_e) invalid for a_e = %i. \
    a_e must be in [0,%i=E).",
               a_e, DIM);
    
  auto slicebd = slice(*this, a_c, a_d, a_e);
  setBoxDataVal<T>(slicebd, a_box, a_val);
}
#ifdef PROTO_CUDA
template<typename T>
struct absolute_value 
{
  __host__ __device__ T operator()(const T &x) const
  {
    return x < T(0) ? -x : x;
  }
};
#endif


/// Maximum Absolute Value (Global)
template <class T, unsigned int C, unsigned char D, unsigned char E>
T BoxData<T,C,D,E>::absMax() const
{

  T minval = this->min();
  T maxval = this->max();
  T absmax = std::max(std::abs(minval), std::abs(maxval));

  return absmax;
}

/// Maximum Absolute Value (Componentwise)
template <class T, unsigned int C, unsigned char D, unsigned char E>
T BoxData<T,C,D,E>::absMax(int a_c, int a_d, int a_e) const
{
  auto slicebd = slice(*this, a_c, a_d, a_e);
  return slicebd.absMax();
}


/// Minimum Value (Global)
template <class T, unsigned int C, unsigned char D, unsigned char E>
T BoxData<T,C,D,E>::min() const
{
    
  PR_TIME("BoxData::min");
  T min = std::numeric_limits<T>::max();
#ifdef PROTO_CUDA

  size_t nsize = size();

  thrust::device_ptr<T> devptr(m_rawPtr);
  cudaError err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }

  T min1 = thrust::reduce(thrust::device, devptr, devptr + nsize, min, thrust::minimum<T>());

  min = min1;
  err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }

#else

  T m = std::numeric_limits<T>::max();
  for (int ii = 0; ii < this->size(); ii++)
  {
    if (m_rawPtr[ii] < m){m = m_rawPtr[ii];}
  }
  min = m;
#endif
  return min;
}

/// Minimum Value (Componentwise)
template <class T, unsigned int C, unsigned char D, unsigned char E>
T BoxData<T,C,D,E>::min(int a_c, int a_d, int a_e) const
{
  auto slicebd = slice(*this, a_c, a_d, a_e);
  return slicebd.min();
}



template <class T, unsigned int C, unsigned char D, unsigned char E>
T BoxData<T,C,D,E>::sum(int a_c, int a_d, int a_e) const
{
  auto slicebd = slice(*this, a_c, a_d, a_e);
  return slicebd.sum();
}

template <class T, unsigned int C, unsigned char D, unsigned char E>
T BoxData<T,C,D,E>::sum() const
{
  PR_TIME("BoxData::sum");
  T sum = 0;
#ifdef PROTO_CUDA

  size_t nsize = size();

  thrust::device_ptr<T> devptr(m_rawPtr);
  cudaError err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }

  T sum1 = thrust::reduce(thrust::device, devptr, devptr + nsize, sum, thrust::plus<T>());

  sum = sum1;
  err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }

#else

  T m = 0;
  for (int ii = 0; ii < this->size(); ii++)
  {
    m += m_rawPtr[ii];
  }
  sum =  m;
#endif
  return sum;
}

/// Maximum Value (Global)
template <class T, unsigned int C, unsigned char D, unsigned char E>
T BoxData<T,C,D,E>::max() const
{
  PR_TIME("BoxData::max");
  T max = -std::numeric_limits<T>::max();
#ifdef PROTO_CUDA

  size_t nsize = size();

  thrust::device_ptr<T> devptr(m_rawPtr);
  cudaError err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }

  T max1 = thrust::reduce(thrust::device, devptr, devptr + nsize, max, thrust::maximum<T>());

  max = max1;
  err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }

#else

  T m = -(std::numeric_limits<T>::max());
  for (int ii = 0; ii < this->size(); ii++)
  {
    if (m_rawPtr[ii] > m){m = m_rawPtr[ii];}
  }
  max =  m;
#endif
  return max;
}

/// Maximum Value (Componentwise)
template <class T, unsigned int C, unsigned char D, unsigned char E>
T BoxData<T,C,D,E>::max(int a_c, int a_d, int a_e) const
{
    
  auto slicebd = slice(*this, a_c, a_d, a_e);
  return slicebd.max();
}

template <class T, unsigned int C, unsigned char D, unsigned char E>
void BoxData<T,C,D,E>::linearOut(void*     a_buf,
                                 const Box& a_box,
                                 CInterval a_comps) const
{
    
  PR_TIME("BoxData::linearout");
  PROTO_ASSERT(m_box.contains(a_box),
               "BoxData::linearOut(void* a_buf, Box& a_box, CInterval a_comps) invalid. \
    a_box must be a subset of this->box().");
  PROTO_ASSERT(this->contains(a_comps),
               "BoxData::linearOut(void* a_buf, Box& a_box, CInterval a_comps) invalid. \
    a_comps must be contained in the component space of *this."); 
    
   
  Box cross = a_box.flatten(0);
  int npencil = a_box.size(0);
    
  T* buffer = (T*)a_buf;
  const T* val;
    
  for (int ee = a_comps.low(2); ee <= a_comps.high(2); ee++)
  {  
    for (int dd = a_comps.low(1); dd <= a_comps.high(1); dd++)
    {
      for (int cc = a_comps.low(0); cc <= a_comps.high(0); cc++)
      { 
        for (auto iter = cross.begin(); iter != cross.end(); ++iter)
        {
          val = this->data(*iter,cc,dd,ee);
          for (int ii = 0; ii < npencil; ii++)
          {
            *buffer = val[ii];
            ++buffer;
          }
        }
      }
    }
  }
}

template <class T, unsigned int C, unsigned char D, unsigned char E>
void BoxData<T,C,D,E>::linearIn(void*     a_buf,
                                const Box& a_box,
                                CInterval a_comps)
{
    
  PR_TIME("BoxData::linearin");
  PROTO_ASSERT(m_box.contains(a_box),
               "BoxData::linearOut(void* a_buf, Box& a_box, CInterval a_comps) invalid. \
    a_box must be a subset of this->box().");
  PROTO_ASSERT(this->contains(a_comps),
               "BoxData::linearOut(void* a_buf, Box& a_box, CInterval a_comps) invalid. \
    a_comps must be contained in the component space of *this."); 
    
  Box cross = a_box.flatten(0);
  int npencil = a_box.size(0);
    
  T* buffer = (T*)a_buf;
  T* val;
    
  for (int ee = a_comps.low(2); ee <= a_comps.high(2); ee++)
  {  
    for (int dd = a_comps.low(1); dd <= a_comps.high(1); dd++)
    {  
      for (int cc = a_comps.low(0); cc <= a_comps.high(0); cc++)
      { 
        for (auto iter = cross.begin(); iter != cross.end(); ++iter)
        {
          val = this->data(*iter,cc,dd,ee);
          for (int ii = 0; ii < npencil; ii++)
          {
            val[ii] = *buffer;
            ++buffer;
          }
        }
      }
    }
  }
}

template <typename T, unsigned int C, unsigned char D, unsigned char E>
void BoxData<T,C,D,E>::print() const
{
    
  using std::cout;
  using std::endl;

  cout << "BoxData<" << C << ", " << (int)D << ", " << (int)E << ">" << endl;
  cout << "Domain: " << m_box << " Size: " << this->size() << " Values" << endl;
  if (this->size() <= 0){return;} 
  cout << "Extremum (Max / Min): (" << this->max() << " / " << this->min() << ")" << endl;
  for (int ee = 0; ee < E; ee++)
    for (int dd = 0; dd < D; dd++)
      for (int cc = 0; cc < C; cc++)
      {
        cout << "\tComponent (" << cc << ", " << dd << ", " << ee << "): ";
        cout << "(" << this->max(cc,dd,ee) << " / " << this->min(cc,dd,ee) << ")" << endl;
      }
}


template <typename T, unsigned int C, unsigned char D, unsigned char E>
void BoxData<T,C,D,E>::printData(int a_prec ) const
{
    
  printData(m_box,-1,-1,-1,a_prec);
}

template <typename T, unsigned int C, unsigned char D, unsigned char E>
void BoxData<T,C,D,E>::printData(const Box& a_box, int a_prec) const
{
    
  printData(a_box,-1,-1,-1, a_prec);
}

template <typename T, unsigned int C, unsigned char D, unsigned char E>
void BoxData<T,C,D,E>::printData(const Box& a_box, int a_c, int a_d, int a_e, int a_prec) const
{
  using std::cout;
  using std::endl;
  if (!defined())
  {
    cout << "Undefined BoxData: No data to print." << endl;
    return;
  }
  if (m_box.size() <= 0)
  {
    cout << "BoxData has empty domain: No data to print." << endl;
    return;
  }
  PROTO_ASSERT(m_box.contains(a_box),
               "BoxData::printData(a_box,...) invalid. Input Box is not a subset of the domain of *this.");
  cout << "Domain: " << m_box << " Subdomain: " << a_box << endl;
  if (a_c != -1)
  {
    PROTO_ASSERT((a_c < C),
                 "BoxData::printData(a_box,a_c,a_d,a_e) invalid. a_c = %i is out of bounds of [0,C=%i).",
                 a_c,C);
  }
  if (a_d != -1)
  {
    PROTO_ASSERT((a_d < D),
                 "BoxData::printData(a_box,a_c,a_d,a_e) invalid. a_d = %i is out of bounds of [0,D=%i).",
                 a_d,D);
  }
  if (a_e != -1)
  {
    PROTO_ASSERT((a_e < E),
                 "BoxData::printData(a_box,a_c,a_d,a_e) invalid. a_e = %i is out of bounds of [0,E=%i).",
                 a_e,E);
  }

  int width = log10(this->absMax());
  int prec = a_prec; 
  width = std::min(width, 6);
  width += (prec+3);

  for (int ee = 0; ee < E; ee++)
    for (int dd = 0; dd < D; dd++)
      for (int cc = 0; cc < C; cc++)
      {
        bool proceed = true;
        proceed &= ((ee == a_e) || (a_e < 0));
        proceed &= ((dd == a_d) || (a_d < 0));
        proceed &= ((cc == a_c) || (a_c < 0));
        if (!proceed){continue;}

        cout << endl << "Component (" << cc << ", " << dd << ", " << ee << "): " << endl;
        cout << std::setfill('=') << std::setw(std::min((int)(a_box.size(0)*(width+1)+8),100)) << " " << endl;
#if DIM < 2
        cout << "\t";
        for (auto iter = a_box.begin(); iter != a_box.end(); ++iter)
        {
          cout << std::fixed << std::setprecision(prec) << std::setfill(' ');
          cout << std::setw(width) << (*this)(*iter,cc,dd,ee) << " ";
        }
        cout << endl;
#elif DIM < 3
        int nx = a_box.size(0);
        int ny = a_box.size(1);
        T* ptr = m_rawPtr + index(a_box.flatten(0).high(),cc,dd,ee);
        for (int jj = ny-1; jj >= 0; jj--)
        {
          cout << "\t";
          for (int ii = 0; ii < nx; ii++)
          {
            cout << std::fixed << std::setprecision(prec) << std::setfill(' ');
            cout << std::setw(width) << ptr[ii] << " ";
          }
          cout << endl;
          ptr -= m_box.size(0);
        }
#else
        Box layerBox = a_box.flatten(0).flatten(1);
        for (auto iter = layerBox.begin(); iter != layerBox.end(); ++iter)
        {
          cout << "Layer (--,--,";
          for (int kk = 2; kk < DIM; kk++)
          {
            cout << (*iter)[kk];
            if ((DIM-kk) > 1){cout << ", ";}
          }
          cout << "):" << endl;
          cout << std::setfill('-') << std::setw(std::min((int)(a_box.size(0)*6+8),100)) << " " << endl;

          Point start = (*iter);
          start[1] += (a_box.size(1)-1);
          T* ptr = m_rawPtr + index(start,cc,dd,ee);
          int nx = a_box.size(0);
          int ny = a_box.size(1);
          for (int jj = ny-1; jj >= 0; jj--)
          {
            cout << "\t";
            for (int ii = 0; ii < nx; ii++)
            {
              cout << std::fixed << std::setprecision(prec) << std::setfill(' ');
              cout << std::setw(width) << ptr[ii] << " ";
            }
            cout << endl;
            ptr -= m_box.size(0);
          }
        }
#endif

      }
  cout << endl;
//    cout << defaultfloat;
}
/// Alias (Non-Const)
template<class T, unsigned int C, unsigned char D, unsigned char E>
BoxData<T,C,D,E> alias(BoxData<T,C,D,E>& a_original, const Point& shift)
{ 
  PR_TIME("alias(BoxData<T,C,D,E>&)");
  const Box& b=a_original.box();
  BoxData<T,C,D,E> rtn(a_original.m_data, a_original.m_rawPtr, b);
  rtn.shift(shift);
  return rtn; // will invoke the alias constructor, thus carry the pointer values, not a deep copy of the rtn object
}

/// Alias (Const)
template<class T, unsigned int C, unsigned char D, unsigned char E>
const BoxData<T,C,D,E> alias(const BoxData<T,C,D,E>& a_original, const Point& shift)
{
  PR_TIME("alias(const BoxData<T,C,D,E>&)");
  const Box& b=a_original.box();
  BoxData<T,C,D,E>* src = const_cast<BoxData<T,C,D,E>*>(&a_original);
  const BoxData<T,C,D,E> rtn(src->m_data, src->m_rawPtr, b.shift(shift));
  return rtn; // will invoke the move constructor, thus carry the pointer values, not a deep copy of the rtn object
}



template<typename T, unsigned int C, unsigned char D, unsigned char E>
BoxData<T> slice(const BoxData<T,C,D,E>& a_src,
                 unsigned int a_c,
                 unsigned int a_d,
                 unsigned int a_e)
{
  PR_TIME("slice(BoxData<T,C,D,E>&, int, int, int)");
  PROTO_ASSERT((a_c < C),
               "slice(BoxData& a_src, uint a_c, uint a_d, uint a_e) invalid.\
    Component index a_c = %i should be positive and less than C = %i.",
               a_c, C);
  PROTO_ASSERT((a_d < D),
               "slice(BoxData& a_src, uint a_c, uint a_d, uint a_e) invalid.\
    Component index a_d = %i should be positive and less than D = %i.",
               a_d, D);
  PROTO_ASSERT((a_e < E),
               "slice(BoxData& a_src, uint a_c, uint a_d, uint a_e) invalid.\
    Component index a_e = %i should be positive and less than E = %i.",
               a_e, E);
  const Box& b = a_src.box();
  T* ptr = (T*)(a_src.data(b.low(),a_c,a_d,a_e));
  BoxData<T> rtn(a_src.getData(), ptr, b);
  return rtn; 
}


template<typename T>
inline T p_ref(const T& a_s, const Point& a_p) {return a_s;}

template<typename T, unsigned int C, unsigned char D, unsigned char E>
inline Var<T,C,D,E> p_ref(BoxData<T,C,D,E>& a_data, const Point& a_p){return --(a_data.var(a_p));}

template<typename T, unsigned int C, unsigned char D, unsigned char E>
inline Var<T,C,D,E> p_ref(const BoxData<T,C,D,E>& a_data, const Point& a_p){return --(a_data.var(a_p));}

template<typename T>
inline T var_incr(const T& t){return t;}


template<typename T, unsigned int C, unsigned char D, unsigned char E>
Var<T, C, D, E>& var_incr(Var<T, C, D, E>& t) {return ++t;}




template<typename Func, typename... T>
inline void pencilFunc(const Func& F, int count, T... vars)
{
  for(int i=0; i < count; ++i)
  {
    F(var_incr(vars)...);
  }
}

template<typename Func, typename... T>
inline void pencilFunc_p(const Func& F, int count, Point a_pt, T... vars)
{
  for(int i=0; i < count; ++i, ++a_pt[0])
  {
    F(a_pt, var_incr(vars)...);
  }
}

////
template <typename T, unsigned int C, unsigned char D, unsigned char E>
void forall_parse(Box& a_domain, const BoxData<T,C,D,E>& a_data)
{
  // Can do anything you want here later
  if(a_domain.empty()) 
  {
    a_domain=a_data.box();
  }
  else 
  {
    a_domain = a_domain & a_data.box();
  }
}

////
template <typename T>
void forall_parse(Box& a_domain, T a_scalar)
{
  // This is a dummy function; it's meant to be empty
}

////
template <typename Last>
void forall_parse_args(Box& a_domain, Last& a_last)
{
  forall_parse(a_domain, a_last);
}

////
template <typename T, unsigned int C, unsigned char D, unsigned char E, typename... Rest>
void forall_parse_args(Box& a_domain, const BoxData<T,C,D,E>& a_first, Rest&... a_rest)
{
  forall_parse(a_domain, a_first);
  forall_parse_args(a_domain, a_rest...);
}

////
template <typename T, unsigned int C, unsigned char D, unsigned char E>
void forall_parse_args(Box& a_domain, const BoxData<T,C,D,E>& a_first)
{
  forall_parse(a_domain, a_first);
//  a_domain=a_first.box();
}

////
template <typename First , typename... Rest>
void forall_parse_args(Box& a_domain, First& a_first, Rest&... a_rest)
{
  forall_parse(a_domain, a_first);
  forall_parse_args(a_domain, a_rest...);
}

//========================================================================
// Base Forall functionality
//========================================================================

// Main version
template<typename Func, typename... Srcs>
void forallInPlaceBase(const Func& a_F, Box a_box, Srcs&... a_srcs)
{
  PR_TIME("forallInPlaceBase");
#ifdef PROTO_CUDA
  cudaForall(a_F, a_box, a_srcs...);
#else
  Box cross = a_box.flatten(0);
  int npencil = a_box.size(0);

  for (auto iter = cross.begin(); iter != cross.end(); ++iter)
  {
    pencilFunc(a_F, npencil, (p_ref(a_srcs,*iter))...);
  }
#endif
}

// Point input version
template<typename Func, typename... Srcs>
inline void forallInPlaceBase_p(const Func& a_F,  Box a_box, Srcs&... a_srcs)
{
  PR_TIME("forallInPlaceBase_p");
#ifdef PROTO_CUDA

  return cudaForall_p(a_F, a_box, a_srcs...) ;

#else
  Box cross = a_box.flatten(0);
  int npencil = a_box.size(0);

  for (auto iter = cross.begin(); iter != cross.end(); ++iter)
  {
    pencilFunc_p(a_F, npencil, *iter, p_ref(a_srcs,*iter)...);
  }
#endif
}

//========================================================================
//  Frontend Forall Implementation
//========================================================================

// Main Forall
template<typename T, unsigned int C, unsigned char D, unsigned char E,
         typename Func, typename... Srcs>
BoxData<T,C,D,E> forall(const Func& a_F, Srcs&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect,a_srcs...);
 
  PROTO_ASSERT(intersect.size() > 0,
               "forall(...) inferred the domain Box to be empty.");

  BoxData<T,C,D,E> ret(intersect);
  forallInPlaceBase(a_F, intersect, ret, a_srcs...);
  return ret;
}

// Overload with Optional Box Argument
template<typename T, unsigned int C, unsigned char D, unsigned char E,
         typename Func, typename... Srcs>
BoxData<T,C,D,E> forall(const Func& a_F, Box a_box, Srcs&... a_srcs)
{
  BoxData<T,C,D,E> ret(a_box);
  forallInPlaceBase(a_F, a_box, ret, a_srcs...);
  return ret;
}

// In Place Version
template<typename Func, typename... Srcs>
void forallInPlace(const Func& a_F, Srcs&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect, a_srcs...);

  PROTO_ASSERT(intersect.size() > 0,
               "forall(...) inferred the domain Box to be empty.");

  forallInPlaceBase(a_F, intersect, a_srcs...);
}

// In Place Version with Box Argument
template<typename Func, typename... Srcs>
void forallInPlace(const Func& a_F, Box a_box, Srcs&... a_srcs)
{
  forallInPlaceBase(a_F, a_box, a_srcs...);
}

// Main version with Point Access
template<typename T, unsigned int C, unsigned char D, unsigned char E,
         typename Func, typename... Srcs>
BoxData<T,C,D,E> forall_p(const Func& a_F, Srcs&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect,a_srcs...);

  PROTO_ASSERT(intersect.size() > 0, "forall(...) inferred the domain Box to be empty.");

  BoxData<T,C,D,E> ret(intersect);
  forallInPlaceBase_p(a_F, intersect, ret, a_srcs...);
  return ret;
}

// Main version with Point Access and Box Argument
template<typename T, unsigned int C, unsigned char D, unsigned char E,
         typename Func, typename... Srcs>
BoxData<T,C,D,E> forall_p(const Func& a_F,  Box a_box, Srcs&... a_srcs)
{
  BoxData<T,C,D,E> ret(a_box);
  forallInPlaceBase_p(a_F, a_box, ret, a_srcs...);
  return ret;
}

// In Place version with Point Access
template<typename Func, typename... Srcs>
inline void forallInPlace_p(const Func& a_F, Srcs&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect, a_srcs...);
    
  PROTO_ASSERT(intersect.size() > 0,
               "forall(...) inferred the domain Box to be empty.");
  forallInPlaceBase_p(a_F, intersect, a_srcs...);
}

// In Place version with Point Access and Box Argument
template<typename Func, typename... Srcs>
inline void forallInPlace_p(const Func& a_F,  Box a_box, Srcs&... a_srcs)
{
  forallInPlaceBase_p(a_F, a_box, a_srcs...);
}

#ifdef PROTO_CUDA

//template<typename T, unsigned int C, unsigned char D, unsigned char E>
//inline 
//template<typename Func, typename... Vars>
//__global__
//void indexer_grid(int a_size[3], Func a_body, Vars... a_vars) 
//{
//    uint i = threadIdx.x + blockIdx.x*blockDim.x;
//    uint j = threadIdx.y + blockIdx.y*blockDim.y;
//    uint k = threadIdx.z + blockIdx.z*blockDim.z;
//    if (i<a_size[0] && j<a_size[1] && k<a_size[2])
//    {
//        a_body(a_vars...);
//    }
//}

template<typename Func, typename... Srcs>
__global__
void indexer(int a_begin, int a_end, Func a_body, Srcs... a_srcs) 
{
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  if (idx >= a_begin && idx < a_end)
  {
    a_body(a_srcs...);
  }
}

template<typename Func, typename... Srcs>
__global__
void indexer_p(int a_begin, int a_end, Box a_box, 
               Func a_body, Srcs... a_srcs) 
{
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  if (idx >= a_begin && idx < a_end)
  {
    Point p = a_box[idx];
    a_body(p,a_srcs...);
  }
}


template<typename Func>
inline Func mapper(const Func& device_f)
{
  PR_TIME("mapper");

  static std::unordered_map<Func, Func> mapCache;
  auto  got = mapCache.find(device_f);
  if(got != mapCache.end())
   {
    return mapCache[device_f];
  }

  Func rtn(device_f);
  cudaError_t err = cudaSuccess;
  {
    PR_TIME("cudaMemcpyFromSymbol");
    err = cudaMemcpyFromSymbol(&rtn, device_f, sizeof(Func), 0, cudaMemcpyDeviceToHost);
    
  }

  if (err != cudaSuccess)
  {
    printf("FAILED to get SYMBOL\n");
    fprintf(stderr, "cudaMemcpyFromSymbol() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }
  mapCache[device_f] = rtn;
  return rtn;
}

template<typename Func, typename... Srcs>
inline void makeVars(const Func& a_F, Box a_box, Srcs... a_srcs)
{
  PR_TIME("BoxData::makevars");
  int N = a_box.size();
  int stride = a_box.size(0);
  int blocks = a_box.flatten(0).size();
  {
    PR_TIME("indexer");
    indexer<<<blocks, stride>>>(0, N, mapper(a_F), a_srcs...);
  }
  cudaError err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }
}
template<typename Func, typename... Srcs>
inline void makeVarsStream(cudaStream_t& a_stream, const Func& a_F, Box a_box, Srcs... a_srcs)
{
  PR_TIME("BoxData::makevarsstream");
  int N = a_box.size();
  int stride = a_box.size(0);
  int blocks = a_box.flatten(0).size();
  size_t smem = 0;
  {
    PR_TIME("indexer");
    indexer<<<blocks, stride, smem, a_stream>>>(0, N, mapper(a_F), a_srcs...);
  }
  cudaError err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }
}

template<typename FuncStruct, typename... Srcs>
__global__
void structIndexer(int a_begin, int a_end, FuncStruct a_body, Srcs... a_srcs) 
{
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  if (idx >= a_begin && idx < a_end)
  {
    a_body.op(a_srcs...);
  }
}

template<typename FuncStruct, typename... Srcs>
inline void makeVarsStruct(cudaStream_t& a_stream, const FuncStruct& a_F, Box a_box, Srcs... a_srcs)
{
  PR_TIME("BoxData::makevarsstruct");
  int N = a_box.size();
  int stride = a_box.size(0);
  int blocks = a_box.flatten(0).size();
  size_t smem = 0;
  {
    PR_TIME("structIndexer");
    structIndexer<<<blocks, stride, smem, a_stream>>>(0, N, a_F, a_srcs...);
  }
  cudaError err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }
}


template<typename... Srcs>
__device__
inline void dummyBlah(Srcs... a_srcs)
{
}

template<typename T>
__device__
inline T var_incr_z(unsigned int zinc, const T& t ){return t;}

template<typename T, unsigned int C, unsigned char D, unsigned char E>
__device__
inline Var<T, C, D, E>& var_incr_z(unsigned int zinc, Var<T, C, D, E>& t) {return (t += zinc);}

template<typename FuncStruct, typename... Srcs>
__global__
void zincStructIndexer(unsigned int a_begin, unsigned int a_Nplane, unsigned int a_Nz, unsigned int a_zinc, FuncStruct a_body, Srcs... a_srcs) 
{
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  if (idx >= a_begin && idx < a_Nplane)
  {
    for(unsigned int zloc = 0; zloc < a_Nz; zloc++)
    {
      a_body.op(a_srcs...);
      dummyBlah(var_incr_z(a_zinc, a_srcs)...);
    }
      
  }
}


template<typename FuncStruct>
__global__
void emptyIndexer(FuncStruct a_body)
{
}

template<typename FuncStruct, typename... Srcs>
inline void makeVarsEmptyIndexer(cudaStream_t& a_stream, const FuncStruct& a_F, Box a_box, Box a_srcBox)
{
  //a_srcs is coming into here as vars
  PR_TIME("BoxData::makevarsEmptyIndexer");

  int stride = a_box.size(0);
  int blocks = a_box.size(1);
#if DIM==3
  int Nplane    = a_box.flatten(2).size();
  int Nz        = a_box.size(2);
  int zinc      = a_srcBox.flatten(2).size();
#else
  int Nplane    = a_box.size();
  int Nz        = 1;
  int zinc      = 1; //has to be non-zero or we have an infinite loop
#endif

  size_t smem = 0;
  {
    PR_TIME("zincStructIndexer");
    emptyIndexer<<<blocks, stride, smem, a_stream>>>(a_F);
  }
  cudaError err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }
}

template<typename FuncStruct, typename... Srcs>
inline void makeVarsZincStruct(cudaStream_t& a_stream, const FuncStruct& a_F, Box a_box, Box a_srcBox, Srcs... a_srcs)
{
  //a_srcs is coming into here as vars
  PR_TIME("BoxData::makevarszincstruct");

  int stride = a_box.size(0);
  int blocks = a_box.size(1);
#if DIM==3
  int Nplane    = a_box.flatten(2).size();
  int Nz        = a_box.size(2);
  int zinc      = a_srcBox.flatten(2).size();
#else
  int Nplane    = a_box.size();
  int Nz        = 1;
  int zinc      = 1; //has to be non-zero or we have an infinite loop
#endif

  size_t smem = 0;
  {
    PR_TIME("zincStructIndexer");
    zincStructIndexer<<<blocks, stride, smem, a_stream>>>(0, Nplane, Nz, zinc, a_F, a_srcs...);
  }
  cudaError err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }
}

template<typename Func, typename... Srcs>
inline void makeVars_p(const Func& a_F, Box a_box, Srcs... a_srcs)
{
  PR_TIME("BoxData::makevars_p");
  int N = a_box.size();
  int stride = a_box.size(0);
  Box cross = a_box.flatten(0);
  int blocks = cross.size();
  {
    PR_TIME("indexer_p");
    indexer_p<<<blocks, stride>>>(0, N, cross,
                                  mapper(a_F), a_srcs...);
  }
  cudaError err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }
}



template<typename T>
inline T p_ref_cuda(const T& a_s, const Point& a_p) {return a_s;}


template<typename T, unsigned int C, unsigned char D, unsigned char E>
inline Var<T,C,D,E> p_ref_cuda(BoxData<T,C,D,E>& a_data, const Point& a_p)
{
  return a_data.var(a_p);
}

template<typename T, unsigned int C, unsigned char D, unsigned char E>
inline Var<T,C,D,E> p_ref_cuda(const BoxData<T,C,D,E>& a_data, const Point& a_p)
{
  return a_data.var(a_p);
}

template<typename Func, typename... Srcs>
inline void cudaForall(const Func& a_F, Box a_box, Srcs&... a_srcs)
{
  PR_TIME("cudaForall");
  makeVars(a_F, a_box, p_ref_cuda(a_srcs, a_box.low())...);
}

template<typename Func, typename... Srcs>
inline void cudaForallStream(cudaStream_t& a_stream, const Func& a_F, Box a_box, Srcs&... a_srcs)
{
  PR_TIME("cudaForallStream");
  makeVarsStream(a_stream, a_F, a_box, p_ref_cuda(a_srcs, a_box.low())...);
}

template<typename FuncStruct, typename... Srcs>
inline void cudaForallStruct(cudaStream_t& a_stream, const FuncStruct& a_F, Box a_box, Srcs&... a_srcs)
{
  PR_TIME("cudaForallStruct");
  makeVarsStruct(a_stream, a_F, a_box, p_ref_cuda(a_srcs, a_box.low())...);
}

template<typename FuncStruct, typename... Srcs>
inline void cudaForallZincStructInt(cudaStream_t& a_stream, const FuncStruct& a_F, Box a_box, Box a_srcBox, Srcs&... a_srcs)
{
  PR_TIME("cudaForallZincStructInt");
  makeVarsZincStruct(a_stream, a_F, a_box, a_srcBox, p_ref_cuda(a_srcs, a_box.low())...);
}


template<typename FuncStruct, typename T, unsigned int C, unsigned char D, unsigned char E, typename... Srcs>
inline void cudaForallZincStruct(cudaStream_t& a_stream, const FuncStruct& a_F, Box a_box, BoxData<T, C, D, E>& a_firstBd, Srcs&... a_srcs)
{
  PR_TIME("cudaForallZincStruct");
  cudaForallZincStructInt(a_stream, a_F, a_box, a_firstBd.box(), a_firstBd, a_srcs...);
}

template<typename FuncStruct, typename... Srcs>
inline void cudaForallEmptyIndexer(cudaStream_t& a_stream, const FuncStruct& a_F, Box a_box, Box a_srcBox)
{
  PR_TIME("cudaForallEmptyIndexer");
  makeVarsEmptyIndexer(a_stream, a_F, a_box, a_srcBox);
}


template<typename Func, typename... Srcs>
inline void cudaForall_p(const Func& a_F, Box a_box, Srcs&... a_srcs) 
{
  PR_TIME("cudaForall_p");
  makeVars_p(a_F, a_box, p_ref_cuda(a_srcs, a_box.low())...);
}

#if 0

template<typename FuncStruct, typename... Srcs>
__global__
void deviceMVIndexer(int a_begin, int a_end, FuncStruct a_body, Srcs... a_srcs) 
{
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  if (idx >= a_begin && idx < a_end)
  {
    a_body.op(a_srcs...);
  }
}

template<typename FuncStruct, typename... Srcs>
inline void makeVarsDeviceMV(cudaStream_t& a_stream, const FuncStruct& a_F, Box a_box, Box a_srcBox, Srcs... a_srcs)
{
  PR_TIME("BoxData::makevarsstruct");
  //srcs are boxdatas and such coming in (Not Vars)
  int N = a_box.size();
  int stride = a_box.size(0);
  int blocks = a_box.flatten(0).size();
  size_t smem = 0;
  {
    PR_TIME("structIndexer");
    deviceMVIndexer<<<blocks, stride, smem, a_stream>>>(0, N, a_F, a_box, a_srcBox, a_srcs...);
  }
  cudaError err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }
}


template<typename FuncStruct, typename... Srcs>
inline void cudaForallDeviceMVInt(cudaStream_t& a_stream, const FuncStruct& a_F, Box a_box, Box a_srcBox, Srcs&... a_srcs)
{
  PR_TIME("cudaForallStruct");
  makeVarsDeviceMV(a_stream, a_F, a_box, a_srcBox, a_srcs...);
}
template<typename FuncStruct, typename T, unsigned int C, unsigned char D, unsigned char E, typename... Srcs>
inline void cudaForallDeviceMV(cudaStream_t& a_stream, const FuncStruct& a_F, Box a_box, BoxData<T, C, D, E>& a_firstBd, Srcs&... a_srcs)
{
  PR_TIME("cudaForallZincStruct");
  cudaForallDeviceMVInt(a_stream, a_F, a_box, a_firstBd.box(), a_firstBd, a_srcs...);
}

#endif

#endif

#endif
