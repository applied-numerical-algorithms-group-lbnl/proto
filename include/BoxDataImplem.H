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
BoxData<T,C,D,E>::BoxData(const Bx& a_box) 
{
  define(a_box);
}

// Initialization Constructor
template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>::BoxData(const Bx& a_box, T a_init) 
{
  define(a_box);
  setVal(a_init);
}

// Deep Copy Constructor
template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>::BoxData(const BoxData<T,C,D,E>& a_src)
{
  //PR_TIMERS("BoxData(BoxData&) (copy constructor)");
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
void BoxData<T,C,D,E>::define(const Bx& a_box)
{
  //PR_TIMERS("BoxData::define(Box) (memory allocation)");
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
BoxData<T,C,D,E>::BoxData(::std::shared_ptr<T> a_data,const  T* a_ptr, const Bx& a_box)
{
  T* castPtr = (T*)a_ptr;
  m_data = ::std::shared_ptr<T>(castPtr, &(null_deleter_boxdata));
  m_box = (a_box);
  m_rawPtr = (T*)a_ptr;
}

// Raw Pointer Constructor
template <class T, unsigned int C, unsigned char D, unsigned char E> 
BoxData<T,C,D,E>::BoxData(const T* a_ptr, const Bx& a_box, int a_ncomp)
{
  define(a_ptr, a_box, a_ncomp);
}

// Raw Pointer Define
template <class T, unsigned int C, unsigned char D, unsigned char E> 
void BoxData<T,C,D,E>::define(const T* a_ptr, const Bx& a_box, int a_ncomp)
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
  for (int ii = 0; ii < a_op.size(); ii++)
  {
    PROTO_ASSERT(!((a_op.m_src[ii])->isAlias(*this)),
                 "Source and destination data are aliased.");
  }
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
  //PR_TIMERS("BoxData::operator=(BoxData&) (copy assign)");
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
  //PR_TIMERS("BoxData::operator=(BoxData&& (move assign)");
    
  if (!a_src.isAlias(*this))
  {    
    std::swap<Bx>(m_box,a_src.m_box);
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
  //PR_TIMERS("BoxData::copyTo(BoxData&)"); 
  Bx domain = m_box & a_dest.m_box;
  copyTo(a_dest,domain,{{0,C-1},{0,D-1},{0,E-1}},
         Point::Zeros(),{{0,C-1},{0,D-1},{0,E-1}});
}    

// Copy with Box Argument (and optional shift)
template <typename T, unsigned int C, unsigned char D, unsigned char E>
void BoxData<T,C,D,E>::copyTo(BoxData<T,C,D,E>&   a_dest,
                              const Bx&              a_srcBox,
                              const Point&            a_destShift) const
{
  //PR_TIMERS("BoxData::copyTo(BoxData&, Box&, Point&)");
  copyTo(a_dest,a_srcBox,{{0,C-1},{0,D-1},{0,E-1}},
         a_destShift,{{0,C-1},{0,D-1},{0,E-1}}); 
}

// Copy From
template <class T,unsigned int C, unsigned char D, unsigned char E>
template<unsigned int Csrc>
void BoxData<T,C,D,E>::copy(const BoxData<T,Csrc,D,E>&  a_src,
          const Bx&     a_srcBox,
          unsigned int  a_srcComp,
          const Bx&     a_destBox,
          unsigned int  a_destComp,
          unsigned int  a_numcomps)
{
  a_src.copyTo(*this, a_srcBox, {a_srcComp,a_srcComp+a_numcomps-1},
               a_destBox.low()-a_srcBox.low(),
               {a_destComp, a_destComp+a_numcomps-1});
}
// General Copy
template< typename T, unsigned int C, unsigned char D, unsigned char E>
template< unsigned int  Cdest,
          unsigned char Ddest,
          unsigned char Edest>
void BoxData<T,C,D,E>::copyTo(BoxData<T,Cdest,Ddest,Edest>& a_dest,
                              const Bx&     a_srcBox,
                              CInterval     a_srcComps,
                              const Point&  a_destShift,
                              CInterval     a_destComps) const
{
#ifdef PROTO_MEM_CHECK
  memcheck::numcopies++;
#endif
  PROTO_ASSERT(contains(a_srcComps),
    "BoxData::copyTo(...) invalid.\
    Source CInterval must be contained in source data component space.");
  PROTO_ASSERT(a_dest.contains(a_destComps),
    "BoxData::copyTo(...) invalid.\
    Destination CInterval must be contained in destination data component space.");
  for (int ii = 0; ii < DIM; ii++)
  {
    PROTO_ASSERT(a_srcComps.size(ii) == a_destComps.size(ii),
      "BoxData::copyTo(...) invalid.\
      Src interval size %i =/= dst interval size %i on axis %i.",
      a_srcComps.size(ii), a_destComps.size(ii),ii);
  }
  PROTO_ASSERT(m_box.contains(a_srcBox),
    "BoxData::copyTo(...) invalid.\
    Source Box input is not a subset of source array domain.");
  PROTO_ASSERT(a_dest.box().contains(a_srcBox.shift(a_destShift)),
    "BoxData::copyTo(...) invalid.\
    Shifted destination Box input is not a subset of destination array domain.");
     
  if (m_data == a_dest.getData()) {return;}
  if (a_srcBox.size() <= 0){return;}
  int doFastCopy = true;
  doFastCopy &= (m_box == a_srcBox);
  doFastCopy &= (a_dest.box() == a_srcBox.shift(a_destShift));
  doFastCopy &= (C == Cdest);
  doFastCopy &= (D == Ddest);
  doFastCopy &= (E == Edest);

  if (doFastCopy)
  {
    unsigned int nsize = size();

#ifdef PROTO_CUDA
    T* castSrc = (T*) m_rawPtr;
    thrust::device_ptr<T> devptrSrc = thrust::device_pointer_cast(castSrc);
    T* dstPtr = (T*)(a_dest.data());
    thrust::device_ptr<T> devptrDst = thrust::device_pointer_cast(dstPtr);
    thrust::copy(thrust::device, devptrSrc, devptrSrc+nsize, devptrDst);
#else
    for (int ii = 0; ii < size(); ii++)
    {
      (*a_dest[ii]) = m_rawPtr[ii];
    }
#endif    
    return;
  }
  Bx cross = a_srcBox.flatten(0);
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
          //TODO: this will need reworking when operator() doesn't return a ref -clg
          srcVal = this->data(*iter,
                              a_srcComps.low(0) + cc,
                              a_srcComps.low(1) + dd,
                              a_srcComps.low(2) + ee);
          destVal = a_dest.data((*iter) + a_destShift,
                              a_destComps.low(0) + cc,
                              a_destComps.low(1) + dd,
                              a_destComps.low(2) + ee);
#ifdef PROTO_CUDA
          T* castSrc = (T*) srcVal;
          thrust::device_ptr<T> devptrSrc = thrust::device_pointer_cast(castSrc);
          thrust::device_ptr<T> devptrDst = thrust::device_pointer_cast(destVal);
          thrust::copy(thrust::device, devptrSrc, devptrSrc+pencilLength, devptrDst);
#else
          for (int ii = 0; ii < pencilLength; ii++)
          {
            destVal[ii] = srcVal[ii];
          }
#endif
        }
      }
    }
  }
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

#ifdef PROTO_CUDA
// Create Pointwise Variable for GPU (Non-Const)
/* TODO: Flagging for removal
template <class T, unsigned int C, unsigned char D, unsigned char E>
Var<T,C,D,E> BoxData<T,C,D,E>::cudaVar(const Point& a_pt)
{
    
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

// Create Pointwise Variable for GPU (Const)
template <class T, unsigned int C, unsigned char D, unsigned char E>
Var<T,C,D,E> BoxData<T,C,D,E>::cudaVar(const Point& a_pt) const
{
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
*/
#endif
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
  //PR_TIMERS("BoxData::operator+=(BoxData&)");
#ifdef PROTO_CUDA
 PROTO_ASSERT(m_box == a_rhs.m_box, "cuda version of arithmetic only codes for equal boxes");

 T* castPtrInp = (T*)(a_rhs.m_rawPtr);
 thrust::device_ptr<T> devptrLoc = thrust::device_pointer_cast(m_rawPtr);
 thrust::device_ptr<T> devptrInp = thrust::device_pointer_cast(castPtrInp);
 unsigned int nsize = size();

 thrust::transform(devptrLoc, devptrLoc+nsize, devptrInp, devptrLoc, thrust::plus<T>());
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
    Bx domain = m_box & a_rhs.m_box;
    if (domain.size() <= 0){return *this;}

    Bx cross = domain.flatten(0);
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
#ifdef PROTO_CUDA
 PROTO_ASSERT(m_box == a_rhs.m_box, "cuda version of arithmetic only codes for equal boxes");

 T* castPtrInp = (T*)(a_rhs.m_rawPtr);
 thrust::device_ptr<T> devptrLoc = thrust::device_pointer_cast(m_rawPtr);
 thrust::device_ptr<T> devptrInp = thrust::device_pointer_cast(castPtrInp);
 unsigned int nsize = size();

 thrust::transform(devptrLoc, devptrLoc+nsize, devptrInp, devptrLoc, thrust::minus<T>());
#else
  //PR_TIMERS("BoxData::operator-=(BoxData&)");
    
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
    Bx domain = m_box & a_rhs.m_box;
    if (domain.size() <= 0){return *this;}

    Bx cross = domain.flatten(0);
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
  //PR_TIMERS("BoxData::operator*=(BoxData&)");
    
#ifdef PROTO_CUDA
 PROTO_ASSERT(m_box == a_rhs.m_box, "cuda version of arithmetic only codes for equal boxes");

 T* castPtrInp = (T*)(a_rhs.m_rawPtr);
 thrust::device_ptr<T> devptrLoc = thrust::device_pointer_cast(m_rawPtr);
 thrust::device_ptr<T> devptrInp = thrust::device_pointer_cast(castPtrInp);
 unsigned int nsize = size();

 thrust::transform(devptrLoc, devptrLoc+nsize, devptrInp, devptrLoc, thrust::multiplies<T>());
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
    Bx domain = m_box & a_rhs.m_box;
    if (domain.size() <= 0){return *this;}

    Bx cross = domain.flatten(0);
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
  //PR_TIMERS("BoxData::operator/=(BoxData&)");
    
#ifdef PROTO_CUDA
 PROTO_ASSERT(m_box == a_rhs.m_box, "cuda version of arithmetic only codes for equal boxes");

 T* castPtrInp = (T*)(a_rhs.m_rawPtr);
 thrust::device_ptr<T> devptrLoc = thrust::device_pointer_cast(m_rawPtr);
 thrust::device_ptr<T> devptrInp = thrust::device_pointer_cast(castPtrInp);
 unsigned int nsize = size();

 thrust::transform(devptrLoc, devptrLoc+nsize, devptrInp, devptrLoc, thrust::divides<T>());
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
    Bx domain = m_box & a_rhs.m_box;
    if (domain.size() <= 0){return *this;}

    Bx cross = domain.flatten(0);
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

#ifdef PROTO_CUDA

  BoxData<T,C,D,E> val(m_box);
  val.setVal(a_scale);
  (*this) += val;

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
#ifdef PROTO_CUDA
  BoxData<T,C,D,E> val(m_box);
  val.setVal(a_scale);
  (*this) -= val;
    
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
  //PR_TIMERS("BoxData::operator*=(T)");
#ifdef PROTO_CUDA

  BoxData<T,C,D,E> val(m_box);
  val.setVal(a_scale);
  (*this) *= val;

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
  PROTO_ASSERT(a_scale != 0,
               "BoxData::operator/=(T a_scale) invalid for a_scale = 0. Division by 0");
#ifdef PROTO_CUDA

  BoxData<T,C,D,E> val(m_box);
  val.setVal(a_scale);
  (*this) /= val;

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
  Bx intersect = a_box & a_data.box();
  Bx face = intersect.flatten(0);
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
        setBoxDataVal<T,Bx>(slicebd, box(), a_val);
        iset++;
      }
    }
  }
#endif
}

// Set Values Within Box
template <class T, unsigned int C, unsigned char D, unsigned char E>
void BoxData<T,C,D,E>::setVal(const T&   a_val,
                              const Bx&  a_box)
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
                              const Bx& a_box,
                              int a_c,
                              int a_d,
                              int a_e)
{
    
  PROTO_ASSERT((a_c >= 0) && (a_c < C),
               "BoxData::setVal(..., int a_c, int a_d, int a_e) invalid for a_c = %i.\
     a_c must be in [0,%i=C).",
               a_c, DIM);
  PROTO_ASSERT((a_d >= 0) && (a_d < D),
               "BoxData::setVal(..., int a_c, int a_d, int a_e) invalid for a_d = %i. \
    a_d must be in [0,%i=D).",
               a_d, DIM);
  PROTO_ASSERT((a_e >= 0) && (a_e < E),
               "BoxData::setVal(..., int a_c, int a_d, int a_e) invalid for a_e = %i. \
    a_e must be in [0,%i=E).",
               a_e, DIM);
  //PR_TIMERS("BoxData::setVal(T&, Box, int, int, int)");
    
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
//#ifdef PROTO_CUDA
//
// size_t nsize = size();
//
//
// thrust::device_ptr<T> devptr(m_rawPtr);
// cudaError err = cudaGetLastError();
// if (err != cudaSuccess)
// {
//   fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
//           __FILE__, __LINE__, cudaGetErrorString(err));
// }
//
// T absmax1 = thrust::transform_reduce(devptr, devptr + nsize, 
//                                   absolute_value<T>(),
//                                   0,
//                                   thrust::maximum<T>());
// absmax = absmax1;
// err = cudaGetLastError();
// if (err != cudaSuccess)
// {
//   fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
//           __FILE__, __LINE__, cudaGetErrorString(err));
// }
//
//#else
//  for (int ii = 0; ii < size(); ii++)
//  {
//    T val = std::abs(m_rawPtr[ii]);
//    if (val > absmax)
//    {
//      absmax = val;
//    }
//  }
//#endif
//  return absmax;
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

/// Maximum Value (Global)
template <class T, unsigned int C, unsigned char D, unsigned char E>
T BoxData<T,C,D,E>::max() const
{
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
                                 const Bx& a_box,
                                 CInterval a_comps) const
{
    
  PROTO_ASSERT(m_box.contains(a_box),
    "BoxData::linearOut(void* a_buf, Box& a_box, CInterval a_comps) invalid. \
    a_box must be a subset of this->box().");
  PROTO_ASSERT(this->contains(a_comps),
    "BoxData::linearOut(void* a_buf, Box& a_box, CInterval a_comps) invalid. \
    a_comps must be contained in the component space of *this."); 
    
   
  Bx cross = a_box.flatten(0);
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
                                const Bx& a_box,
                                CInterval a_comps)
{
    
  PROTO_ASSERT(m_box.contains(a_box),
    "BoxData::linearOut(void* a_buf, Box& a_box, CInterval a_comps) invalid. \
    a_box must be a subset of this->box().");
  PROTO_ASSERT(this->contains(a_comps),
    "BoxData::linearOut(void* a_buf, Box& a_box, CInterval a_comps) invalid. \
    a_comps must be contained in the component space of *this."); 
    
  Bx cross = a_box.flatten(0);
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
void BoxData<T,C,D,E>::printData(const Bx& a_box, int a_prec) const
{
    
  printData(a_box,-1,-1,-1, a_prec);
}

template <typename T, unsigned int C, unsigned char D, unsigned char E>
void BoxData<T,C,D,E>::printData(const Bx& a_box, int a_c, int a_d, int a_e, int a_prec) const
{
    
  using namespace std;
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
    PROTO_ASSERT((a_c >= 0) && (a_c < C),
                 "BoxData::printData(a_box,a_c,a_d,a_e) invalid. a_c = %i is out of bounds of [0,C=%i).",
                 a_c,C);
  }
  if (a_d != -1)
  {
    PROTO_ASSERT((a_d >= 0) && (a_d < D),
                 "BoxData::printData(a_box,a_c,a_d,a_e) invalid. a_d = %i is out of bounds of [0,D=%i).",
                 a_d,D);
  }
  if (a_e != -1)
  {
    PROTO_ASSERT((a_e >= 0) && (a_e < E),
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
        cout << setfill('=') << setw(std::min((int)(a_box.size(0)*(width+1)+8),100)) << " " << endl;
#if DIM < 2
        cout << "\t";
        for (auto iter = a_box.begin(); iter != a_box.end(); ++iter)
        {
          cout << fixed << setprecision(prec) << setfill(' ');
          cout << setw(width) << (*this)(*iter,cc,dd,ee) << " ";
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
            cout << fixed << setprecision(prec) << setfill(' ');
            cout << setw(width) << ptr[ii] << " ";
          }
          cout << endl;
          ptr -= m_box.size(0);
        }
#else
        Bx layerBox = a_box.flatten(0).flatten(1);
        for (auto iter = layerBox.begin(); iter != layerBox.end(); ++iter)
        {
          cout << "Layer (--,--,";
          for (int kk = 2; kk < DIM; kk++)
          {
            cout << (*iter)[kk];
            if ((DIM-kk) > 1){cout << ", ";}
          }
          cout << "):" << endl;
          cout << setfill('-') << setw(std::min((int)(a_box.size(0)*6+8),100)) << " " << endl;

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
              cout << fixed << setprecision(prec) << setfill(' ');
              cout << setw(width) << ptr[ii] << " ";
            }
            cout << endl;
            ptr -= m_box.size(0);
          }
        }
#endif

      }
//    cout << defaultfloat;
}
/// Alias (Non-Const)
template<class T, unsigned int C, unsigned char D, unsigned char E>
BoxData<T,C,D,E> alias(BoxData<T,C,D,E>& a_original, const Point& shift)
{ 
  //PR_TIMERS("alias(BoxData<T,C,D,E>&)");
  const Bx& b=a_original.box();
  BoxData<T,C,D,E> rtn(a_original.m_data, a_original.m_rawPtr, b);
  rtn.shift(shift);
  return rtn; // will invoke the alias constructor, thus carry the pointer values, not a deep copy of the rtn object
}

/// Alias (Const)
template<class T, unsigned int C, unsigned char D, unsigned char E>
const BoxData<T,C,D,E> alias(const BoxData<T,C,D,E>& a_original, const Point& shift)
{
  //PR_TIMERS("alias(const BoxData<T,C,D,E>&)");
  const Bx& b=a_original.box();
  BoxData<T,C,D,E>* src = const_cast<BoxData<T,C,D,E>*>(&a_original);
  const BoxData<T,C,D,E> rtn(src->m_data, src->m_rawPtr, b.shift(shift));
  return rtn; // will invoke the move constructor, thus carry the pointer values, not a deep copy of the rtn object
}

template<typename T, unsigned int C, unsigned char D, unsigned char E>
BoxData<T> slice(BoxData<T,C,D,E>& a_src,
                 unsigned int a_c,
                 unsigned int a_d,
                 unsigned int a_e)
{
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
  //PR_TIMERS("slice(BoxData<T,C,D,E>&, int, int, int)");
  const Bx& b = a_src.box();
  T* ptr = a_src.m_rawPtr + a_src.index(b.low(),a_c,a_d,a_e);
  BoxData<T> rtn(a_src.m_data, ptr, b);
  return rtn; 
}

template<typename T, unsigned int C, unsigned char D, unsigned char E>
const BoxData<T> slice(const BoxData<T,C,D,E>& a_src, int a_c, int a_d, int a_e)
{
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
  //PR_TIMERS("slice(const BoxData<T,C,D,E>&, int, int, int)");
  const Bx& b = a_src.box();
  
  
  const T* srcPtr =  a_src[0];
  shared_ptr<T> data = a_src.getData();
  const T* ptr = srcPtr + a_src.index(b.low(),a_c,a_d,a_e);

  BoxData<T> rtn(data, ptr, b);
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
void forall_parse(Bx& a_domain, const BoxData<T,C,D,E>& a_data)
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
void forall_parse(Bx& a_domain, T a_scalar)
{
  // This is a dummy function; it's meant to be empty
}

////
template <typename Last>
void forall_parse_args(Bx& a_domain, Last& a_last)
{
  forall_parse(a_domain, a_last);
}

////
template <typename T, unsigned int C, unsigned char D, unsigned char E, typename... Rest>
void forall_parse_args(Bx& a_domain, const BoxData<T,C,D,E>& a_first, Rest&... a_rest)
{
  forall_parse(a_domain, a_first);
  forall_parse_args(a_domain, a_rest...);
}

////
template <typename T, unsigned int C, unsigned char D, unsigned char E>
void forall_parse_args(Bx& a_domain, const BoxData<T,C,D,E>& a_first)
{
  forall_parse(a_domain, a_first);
//  a_domain=a_first.box();
}

////
template <typename First , typename... Rest>
void forall_parse_args(Bx& a_domain, First& a_first, Rest&... a_rest)
{
  forall_parse(a_domain, a_first);
  forall_parse_args(a_domain, a_rest...);
}

//========================================================================
// Base Forall functionality
//========================================================================

// Main version
template<typename Func, typename... Srcs>
void forallInPlaceBase(const Func& a_F, Bx a_box, Srcs&... a_srcs)
{
#ifdef PROTO_CUDA
  cudaForall(a_F, a_box, a_srcs...);
#else
  Bx cross = a_box.flatten(0);
  int npencil = a_box.size(0);

  for (auto iter = cross.begin(); iter != cross.end(); ++iter)
  {
    pencilFunc(a_F, npencil, (p_ref(a_srcs,*iter))...);
  }
#endif
}

// Point input version
template<typename Func, typename... Srcs>
inline void forallInPlaceBase_p(const Func& a_F,  Bx a_box, Srcs&... a_srcs)
{
#ifdef PROTO_CUDA

  return cudaForall_p(a_F, a_box, a_srcs...) ;

#else
  Bx cross = a_box.flatten(0);
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
  Bx intersect;
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
BoxData<T,C,D,E> forall(const Func& a_F, Bx a_box, Srcs&... a_srcs)
{
  BoxData<T,C,D,E> ret(a_box);
  forallInPlaceBase(a_F, a_box, ret, a_srcs...);
  return ret;
}

// In Place Version
template<typename Func, typename... Srcs>
void forallInPlace(const Func& a_F, Srcs&... a_srcs)
{
  Bx intersect;
  forall_parse_args(intersect, a_srcs...);

  PROTO_ASSERT(intersect.size() > 0,
               "forall(...) inferred the domain Box to be empty.");

  forallInPlaceBase(a_F, intersect, a_srcs...);
}

// In Place Version with Box Argument
template<typename Func, typename... Srcs>
void forallInPlace(const Func& a_F, Bx a_box, Srcs&... a_srcs)
{
  forallInPlaceBase(a_F, a_box, a_srcs...);
}

// Main version with Point Access
template<typename T, unsigned int C, unsigned char D, unsigned char E,
         typename Func, typename... Srcs>
BoxData<T,C,D,E> forall_p(const Func& a_F, Srcs&... a_srcs)
{
  Bx intersect;
  forall_parse_args(intersect,a_srcs...);

  PROTO_ASSERT(intersect.size() > 0, "forall(...) inferred the domain Box to be empty.");

  BoxData<T,C,D,E> ret(intersect);
  forallInPlaceBase_p(a_F, intersect, ret, a_srcs...);
  return ret;
}

// Main version with Point Access and Box Argument
template<typename T, unsigned int C, unsigned char D, unsigned char E,
         typename Func, typename... Srcs>
BoxData<T,C,D,E> forall_p(const Func& a_F,  Bx a_box, Srcs&... a_srcs)
{
  BoxData<T,C,D,E> ret(a_box);
  forallInPlaceBase_p(a_F, a_box, ret, a_srcs...);
  return ret;
}

// In Place version with Point Access
template<typename Func, typename... Srcs>
inline void forallInPlace_p(const Func& a_F, Srcs&... a_srcs)
{
  Bx intersect;
  forall_parse_args(intersect, a_srcs...);
    
  PROTO_ASSERT(intersect.size() > 0,
               "forall(...) inferred the domain Box to be empty.");
  forallInPlaceBase_p(a_F, intersect, a_srcs...);
}

// In Place version with Point Access and Box Argument
template<typename Func, typename... Srcs>
inline void forallInPlace_p(const Func& a_F,  Bx a_box, Srcs&... a_srcs)
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

__device__ __host__
inline Point compute_point(Point a_low, Point a_high, unsigned int a_index)
{
  int tuple[DIM];
  for (int ii = 0; ii < DIM; ii++)
  {
    int factor = (a_high[ii] - a_low[ii] + 1);
    int kred = a_index % factor;
    tuple[ii] = kred + a_low[ii];
    a_index = (a_index - kred)/factor;
  }
  return Point(tuple);
}

template<typename Func, typename... Srcs>
__global__
void indexer_p(int a_begin, int a_end, Point a_low, Point a_high,
               Func a_body, Srcs... a_srcs) 
{
  int idx = threadIdx.x + blockIdx.x*blockDim.x;
  if (idx >= a_begin && idx < a_end)
  {
    Point p = compute_point(a_low, a_high, blockIdx.x);
    p[0] += threadIdx.x;
    a_body(p,a_srcs...);
  }
}


template<typename Func>
inline Func mapper(const Func& device_f)
{
  Func rtn(device_f);

  cudaError_t err = cudaMemcpyFromSymbol(&rtn, device_f, sizeof(Func), 0, cudaMemcpyDeviceToHost);
  if (err != cudaSuccess)
  {
    printf("FAILED to get SYMBOL\n");
    fprintf(stderr, "cudaMemcpyFromSymbol() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }
  return rtn;
}

template<typename Func, typename... Srcs>
inline void makeVars(const Func& a_F, Bx a_box, Srcs... a_srcs)
{
  int N = a_box.size();
  int stride = a_box.size(0);
  int blocks = a_box.flatten(0).size();
  indexer<<<blocks, stride>>>(0, N, mapper(a_F), a_srcs...);
  cudaError err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }
}

template<typename Func, typename... Srcs>
inline void makeVars_p(const Func& a_F, Bx a_box, Srcs... a_srcs)
{
  int N = a_box.size();
  int stride = a_box.size(0);
  Bx cross = a_box.flatten(0);
  int blocks = cross.size();
  indexer_p<<<blocks, stride>>>(0, N, cross.low(), cross.high(),
                                mapper(a_F), a_srcs...);
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
inline void cudaForall(const Func& a_F, Bx a_box, Srcs&... a_srcs)
{
  makeVars(a_F, a_box, p_ref_cuda(a_srcs, a_box.low())...);
}

template<typename Func, typename... Srcs>
        
inline void cudaForall_p(const Func& a_F, Bx a_box, Srcs&... a_srcs) 
{
  makeVars_p(a_F, a_box, p_ref_cuda(a_srcs, a_box.low())...);
}


#endif

#endif
