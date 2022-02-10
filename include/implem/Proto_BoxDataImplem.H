#ifndef __BOXDATAIMPLEM__
#define __BOXDATAIMPLEM__

#pragma once

//=================================================================================================
// Constructors ||
//==============++

#define PROTO_KERNEL_START inline CUDA_DECORATION
#define PROTO_KERNEL_END(local_name, app_name)                  \
struct struct_##local_name {   \
    template <typename... T> \
    inline CUDA_DECORATION void operator()(T... args) const { local_name(args...);} \
    const char* myname = #app_name; \
};        \
static struct_##local_name app_name;
#define PROTO_LAMBDA CUDA_DECORATION

//=================================================================================================
// Constructors ||
//==============++

// Default Constructor
template <class T, unsigned int C,MemType MEMTYPE, unsigned char D, unsigned char E> 
BoxData<T,C,MEMTYPE,D,E>::BoxData() : m_rawPtr(NULL) {}

// Box Constructor
template <class T, unsigned int C,MemType MEMTYPE, unsigned char D, unsigned char E> 
BoxData<T,C,MEMTYPE,D,E>::BoxData(const Box& a_box)  //stack allocation defaults to true
{
    define(a_box);
}

// Initialization Constructor
template <class T, unsigned int C,MemType MEMTYPE, unsigned char D, unsigned char E> 
BoxData<T,C,MEMTYPE,D,E>::BoxData(const Box& a_box, T a_init) 
{
    define(a_box);
    setVal(a_init);
}

// Move Constructor
// TODO: *NEED* to make sure that move semantics are functioning properly
template <class T, unsigned int C,MemType MEMTYPE, unsigned char D, unsigned char E> 
BoxData<T,C,MEMTYPE,D,E>::BoxData(BoxData<T,C,MEMTYPE,D,E>&& a_src)
  :   m_box(a_src.m_box),
      m_data(std::move(a_src.m_data)),
      m_rawPtr(a_src.m_rawPtr)
{
    std::cout << "Calling BoxData move constructor" << std::endl;
}

// Define
template <class T, unsigned int C,MemType MEMTYPE, unsigned char D, unsigned char E> 
void BoxData<T,C,MEMTYPE,D,E>::define(const Box& a_box)
{
    m_box=a_box;
    auto& stack = Stack<MEMTYPE>::getStack();
    if (stack.enabled())
    {
        m_rawPtr = (T*)(stack.alloc(size()*sizeof(T)));
        PROTO_ASSERT(m_rawPtr != NULL,
            "BoxData::define | Error: Stack returned a NULL pointer when it shouldn't have.");
        // Deletion is managed by Stack
        m_data = std::shared_ptr<T>(m_rawPtr, &(null_deleter_boxdata<MEMTYPE>));
        m_stackAlloc = true;
    } else {
        protoMalloc(MEMTYPE, m_rawPtr, size()*sizeof(T));
        m_data = std::shared_ptr<T>(m_rawPtr, [](T* p){ protoFree(MEMTYPE,p);}); 
        m_stackAlloc = false;
    }
}

// Destructor
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
BoxData<T,C,MEMTYPE,D,E>::~BoxData()
{}

// Slice Constructor (Private)
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
BoxData<T,C,MEMTYPE,D,E>::BoxData(::std::shared_ptr<T> a_data,const  T* a_ptr, const Box& a_box)
{
    m_data = a_data;
    m_box = (a_box);
    m_rawPtr = (T*)a_ptr;
}

// Raw Pointer Constructor
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
BoxData<T,C,MEMTYPE,D,E>::BoxData(const T* a_ptr, const Box& a_box, int a_ncomp)
{
    define(a_ptr, a_box, a_ncomp);
}

// Raw Pointer Define
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::define(const T* a_ptr, const Box& a_box, int a_ncomp)
{
    PROTO_ASSERT((a_ncomp==C && D==1 && E==1),
        "BoxData::define(T*, Box&, int) | Error: Component mismatch in boxdata alias");
    T* castPtr = const_cast<T*>(a_ptr);
    m_data = ::std::shared_ptr<T>(castPtr, &(null_deleter_boxdata<MEMTYPE>));
    m_rawPtr = castPtr;
    m_box = a_box;
}

template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
BoxData<T,C,MEMTYPE,D,E>::BoxData(LazyStencil<T, C, MEMTYPE, D, E>&& a_op)
{
    define(a_op.m_range);
    a_op.apply(*this, true);
}

//=================================================================================================
// Data Movement ||
//===============++

// Move Assignment Operator
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
BoxData<T,C,MEMTYPE,D,E>& BoxData<T,C,MEMTYPE,D,E>::operator=(BoxData<T,C,MEMTYPE,D,E>&& a_src)
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
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::copyTo(BoxData<T,C,MEMTYPE,D,E>& a_dest) const
{
    PR_TIME("BoxData::copyTo(BoxData&)");
    Box domain = m_box & a_dest.m_box;
    copyTo(a_dest,domain,{{0,C-1},{0,D-1},{0,E-1}},
            Point::Zeros(),{{0,C-1},{0,D-1},{0,E-1}});
}

// Copy with Box Argument (and optional shift)
template <typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::copyTo(BoxData<T,C,MEMTYPE,D,E>&   a_dest,
        const Box&              a_srcBox,
        const Point&            a_destShift) const
{
    PR_TIME("BoxData::copyTo(BoxData&, Box&, Point&)");
    copyTo(a_dest,a_srcBox,{{0,C-1},{0,D-1},{0,E-1}},
            a_destShift,{{0,C-1},{0,D-1},{0,E-1}});
}

// Copy From
// TODO: Get rid of this if we can
template <class T,unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
template<unsigned int Csrc>
void BoxData<T,C,MEMTYPE,D,E>::copy(const BoxData<T,Csrc,MEMTYPE,D,E>&  a_src,
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


template <BoxDataOp Op, typename T, size_t C>
struct boxdataIndexer {
    static void cpu(unsigned int  a_begin,
            unsigned int a_end,
            T* a_src,  T* a_dst,
            Box a_box, Box a_srcBox, Box a_dstBox,
            Point  a_dstShift, int srcSize, int dstSize)
    {
        for(int idx = a_begin ; idx < a_end ; idx++)
        {
            Point boxPt = a_box[idx];
            Point srcPt = boxPt;
            Point dstPt = boxPt + a_dstShift;

            unsigned int dstOff = a_dstBox.index(dstPt);
            unsigned int srcOff = a_srcBox.index(srcPt);

            for(int c = 0 ; c < C; c++)
            {
                T* tmp_src = a_src + c * srcSize ; //(a_end == nb_elem);
                T* tmp_dst = a_dst + c * dstSize ;

                if(Op == BoxDataOp::Copy)
                {
                    tmp_dst[dstOff] = tmp_src[srcOff];
                }
                else if(Op == BoxDataOp::Add)
                {
                    tmp_dst[dstOff] += tmp_src[srcOff];
                }
                else if(Op == BoxDataOp::Subtract)
                {
                    tmp_dst[dstOff] -= tmp_src[srcOff];
                }
                else if(Op == BoxDataOp::Multiply)
                {
                    tmp_dst[dstOff] *= tmp_src[srcOff];
                }
                else if(Op == BoxDataOp::Divide)
                {
                    tmp_dst[dstOff] /= tmp_src[srcOff];
                }
                else if(Op == BoxDataOp::Assign)
                {
                    tmp_dst[dstOff] = tmp_src[srcOff];
                }
                else
                {
                    printf("boxdataIndexer error: bogus boxdata op input!!!\n");
                }
            }
        } // end for idx
    } // end cpu
#ifdef PROTO_CUDA
    __device__ static void gpu(unsigned int  a_begin,
            unsigned int a_end,
            T* a_src,  T* a_dst,
            const Box& a_box, const Box& a_srcBox, const Box& a_dstBox,
            const Point&  a_dstShift, int srcSize, int dstSize)
    {
        unsigned int idx = a_begin + threadIdx.x+blockIdx.x*blockDim.x;

        if((idx >= a_begin) && (idx< a_end))
        {
            Point boxPt = a_box[idx];
            Point srcPt = boxPt;
            Point dstPt = boxPt + a_dstShift;

            unsigned int dstOff = a_dstBox.index(dstPt);
            unsigned int srcOff = a_srcBox.index(srcPt);

            for(int c = 0 ; c < C; c++)
            {
                T* tmp_src = a_src + c * srcSize ; //(a_end == nb_elem);
                T* tmp_dst = a_dst + c * dstSize ;

                if(Op == BoxDataOp::Copy)
                {
                    tmp_dst[dstOff] = tmp_src[srcOff];
                }
                else if(Op == BoxDataOp::Add)
                {
                    tmp_dst[dstOff] += tmp_src[srcOff];
                }
                else if(Op == BoxDataOp::Subtract)
                {
                    tmp_dst[dstOff] -= tmp_src[srcOff];
                }
                else if(Op == BoxDataOp::Multiply)
                {
                    tmp_dst[dstOff] *= tmp_src[srcOff];
                }
                else if(Op == BoxDataOp::Divide)
                {
                    tmp_dst[dstOff] /= tmp_src[srcOff];
                }
                else if(Op == BoxDataOp::Assign)
                {
                    tmp_dst[dstOff] = tmp_src[srcOff];
                }
                else
                {
                    printf("boxdataIndexer error: bogus boxdata op input!!!\n");
                }
            }
        }
    }
#endif
}; // end boxData indexer

template <BoxDataOp Op, typename T, size_t C>
struct scalarIndexer
{
    static void cpu(int a_begin, int a_end, T* a_dst, T a_scalar, Box a_box, unsigned int dstSize)
    {
        for(int idx = a_begin ; idx < a_end ; idx++)
        {
            Point boxPt = a_box[idx];
            unsigned int dstOff = a_box.index(boxPt);

            T* tmp_dst;
            for(int c = 0 ; c < C ; c++)
            {
                tmp_dst = a_dst + c * dstSize;
                if(Op == BoxDataOp::Copy)
                {
                    tmp_dst[dstOff] = a_scalar;
                }
                else if(Op == BoxDataOp::Add)
                {
                    tmp_dst[dstOff] += a_scalar;
                }
                else if(Op == BoxDataOp::Subtract)
                {
                    tmp_dst[dstOff] -= a_scalar;
                }
                else if(Op == BoxDataOp::Multiply)
                {
                    tmp_dst[dstOff] *= a_scalar;
                }
                else if(Op == BoxDataOp::Divide)
                {
                    tmp_dst[dstOff] /= a_scalar;
                }
                else if(Op == BoxDataOp::Assign)
                {
                    tmp_dst[dstOff] = a_scalar;
                }
                else
                {
                    printf("scalarIndexer error: bogus boxdata op input!!!\n");
                }
            }
        }
    }
 #ifdef PROTO_CUDA
    __device__ static  void gpu(int a_begin, int a_end, T* a_dst, T a_scalar, const Box& a_box, unsigned int dstSize)
    {
        unsigned int idx = a_begin + threadIdx.x+blockIdx.x*blockDim.x;
        if((idx >= a_begin) && (idx< a_end))
        {
            Point boxPt = a_box[idx];
            unsigned int dstOff = a_box.index(boxPt);

            T* tmp_dst;
            for(int c = 0 ; c < C ; c++)
            {
                tmp_dst = a_dst + c * dstSize;
                if(Op == BoxDataOp::Copy)
                {
                    tmp_dst[dstOff] = a_scalar;
                }
                else if(Op == BoxDataOp::Add)
                {
                    tmp_dst[dstOff] += a_scalar;
                }
                else if(Op == BoxDataOp::Subtract)
                {
                    tmp_dst[dstOff] -= a_scalar;
                }
                else if(Op == BoxDataOp::Multiply)
                {
                    tmp_dst[dstOff] *= a_scalar;
                }
                else if(Op == BoxDataOp::Divide)
                {
                    tmp_dst[dstOff] /= a_scalar;
                }
                else if(Op == BoxDataOp::Assign)
                {
                    tmp_dst[dstOff] = a_scalar;
                }
                else
                {
                    printf("scalarIndexer error: bogus boxdata op input!!!\n");
                }
            }
        }
    }
#endif
}; // end scalar indexer

template<typename T,BoxDataOp op>
struct opKernel
{
    static void cpu(T* x, T y, const unsigned int size)
    {
        switch(op)
        {
            case BoxDataOp::Assign:
                for(int id = 0 ; id < size ; id++)
                    x[id] = y;
                break;
        }
    }
    static void cpu(T* x, T* y, const unsigned int size)
    {
        switch(op)
        {
            case BoxDataOp::Add:
                for(int id = 0 ; id < size ; id++)
                    x[id] += y[id];
                break;
            case BoxDataOp::Subtract:
                for(int id = 0 ; id < size ; id++)
                    x[id] -= y[id];
                break;
            case BoxDataOp::Divide:
                for(int id = 0 ; id < size ; id++)
                    x[id] /= y[id];
                break;
            case BoxDataOp::Multiply:
                for(int id = 0 ; id < size ; id++)
                    x[id] *= y[id];
                break;
        }
    }
#ifdef PROTO_CUDA
    __device__ static void gpu(T* x, T* y, const unsigned int size)
    {
        int id = blockIdx.x*blockDim.x+threadIdx.x;
        if (id < size)
        {
            switch(op)
            {
                case BoxDataOp::Add:
                    x[id] += y[id];
                    break;
                case BoxDataOp::Subtract:
                    x[id] -= y[id];
                    break;
                case BoxDataOp::Divide:
                    x[id] /= y[id];
                    break;
                case BoxDataOp::Multiply:
                    x[id] *= y[id];
                    break;
                case BoxDataOp::Assign:
                    x[id] = y[id];
                    break;
            }
        }
    }

    __device__ static void gpu(T* x, T y, const unsigned int size)
    {
        int id = blockIdx.x*blockDim.x+threadIdx.x;
        if (id < size)
        {
            switch(op)
            {
                case BoxDataOp::Assign:
                    x[id] = y;
                    break;
            }
        }
    }
#endif
};

//Specialized copy in which a_src and a_dst have the same number of components.
//See general copy below
template< typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::copyTo(
        BoxData<T,C,MEMTYPE,D,E>& a_dst,
        const Box&  a_srcBox,
        const Box&  a_dstBox) const
{
    CInterval cint(0,C-1);
    //Box srcBox = a_srcBox&a_dstBox;
    Point shift = Point::Zeros();
    this->copyTo(a_dst,a_srcBox,cint,shift,cint);
}

// General Copy
template< typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
template< unsigned int  Cdest,
          unsigned char Ddest,
          unsigned char Edest>
void BoxData<T,C,MEMTYPE,D,E>::copyTo(
        BoxData<T,Cdest,MEMTYPE,Ddest,Edest>& a_dst,
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
    doFastCopy &= (a_srcComps == a_dstComps);

    if (doFastCopy)
    {
        unsigned int nsize = size();
        T* castSrc = (T*) m_rawPtr;
        T* dstPtr = (T*)(a_dst.data());
        protoMemcpyAsync(MEMTYPE,dstPtr, castSrc, nsize*sizeof(T),
                protoMemcpyDeviceToDevice,
                protoGetCurrentStream);
    }
    else
    {
        for (unsigned int ee = 0; ee < a_srcComps.size(2); ee++)
        {
            for (unsigned int dd = 0; dd < a_srcComps.size(1); dd++)
            {
                int ccSrc = a_srcComps.low(0);
                int ddSrc = dd + a_srcComps.low(1);
                int eeSrc = ee + a_srcComps.low(2);
                int ccDst = a_dstComps.low(0);
                int ddDst = dd + a_dstComps.low(1);
                int eeDst = ee + a_dstComps.low(2);

                BoxData<T,1,MEMTYPE> srcbd = slice(*this, ccSrc, ddSrc, eeSrc);
                BoxData<T,1,MEMTYPE> dstbd = slice(a_dst, ccDst, ddDst, eeDst);

                unsigned int begin  = 0;
                unsigned int end    = a_srcBox.size();

                unsigned int stride = a_srcBox.size(0);
                unsigned int blocks = a_srcBox.flatten(0).size();

                T* d_srcPtr = (T*)srcbd.data();
                T* d_dstPtr = (T*)dstbd.data();
                {
                    PR_TIME("boxindexer");

                    protoLaunchKernelMemAsyncT<MEMTYPE,
                        boxdataIndexer<BoxDataOp::Copy,T,C>>(
                                blocks, stride,0,
                                protoGetCurrentStream, begin, end, d_srcPtr, d_dstPtr,
                                a_srcBox, srcbd.box(), dstbd.box(),
                                a_dstShift, srcbd.box().size(), dstbd.box().size());
                } // end indexer
            } // end for dd
        } // end for ee
    } 
}


// Index Accessor (Non-Const)
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
T* BoxData<T,C,MEMTYPE,D,E>::operator[](unsigned int a_index)
{
    PROTO_ASSERT(a_index < size(),
            "BoxData::operator[](uint a_index) invalid for a_index = %i. \
            a_index must be in [0,this->size() = %i)",a_index, (int)size());

    return &m_rawPtr[a_index];
}

// Index Accessor (Const)
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
const T* BoxData<T,C,MEMTYPE,D,E>::operator[](unsigned int a_index) const
{
    PROTO_ASSERT(a_index < size(),
            "BoxData::operator[](uint a_index) invalid for a_index = %i.\
            a_index must be in [0,this->size() = %i)",a_index, (int)size());
    return &m_rawPtr[a_index];
}


// Create Pointwise Variable (Non-Const)
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
Var<T,C,MEMTYPE,D,E> BoxData<T,C,MEMTYPE,D,E>::var(const Point& a_pt)
{
//  PR_TIME("BoxData::var");
  PROTO_ASSERT(m_box.contains(a_pt),
               "BoxData::var(Point a_pt) invalid for a_pt not in this->box()");
  Var<T,C,MEMTYPE,D,E> rtn;
  for(unsigned int ee = 0; ee < E; ++ee)
    for(unsigned int dd = 0; dd < D; ++dd)
      for(unsigned int cc = 0; cc < C; ++cc)
      {
        rtn.m_ptrs[cc + C*dd + C*D*ee] = m_rawPtr + index(a_pt,cc,dd,ee);
      }
  rtn.boxDimX=this->box().size(0);
  rtn.boxDimY=this->box().size(1);

  return rtn;
}

// Create Pointwise Variable (Const)
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
Var<T,C,MEMTYPE,D,E> BoxData<T,C,MEMTYPE,D,E>::var(const Point& a_pt) const
{
//  PR_TIME("BoxData::var_const");
  PROTO_ASSERT(m_box.contains(a_pt),
               "BoxData::var(Point a_pt) invalid for a_pt not in this->box()");
  Var<T,C,MEMTYPE,D,E> rtn;
  for(unsigned int ee = 0; ee < E; ++ee)
    for(unsigned int dd = 0; dd < D; ++dd)
      for(unsigned int cc = 0; cc < C; ++cc)
      {
        rtn.m_ptrs[cc + C*dd + C*D*ee] = m_rawPtr + index(a_pt,cc,dd,ee);
      }
  rtn.boxDimX=this->box().size(0);
  rtn.boxDimY=this->box().size(1);
  return rtn;
}

//=================================================================================================
// Algebraic Operations ||
//======================++

// Operator |= (BoxData&)
/*
  template <class T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE>
  BoxData<T,C,MEMTYPE,D,E>& BoxData<T,C,MEMTYPE,D,E>::operator|=(const BoxData<T,C,MEMTYPE,D,E>& a_rhs)
  {
  a_rhs.copyTo(*this);
  return *this;
  }
*/


// Operator += (BoxData&)
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
template <BoxDataOp op>
void BoxData<T,C,MEMTYPE,D,E>::operatorT(const BoxData<T,C,MEMTYPE,D,E>& a_rhs)
{
  // string message;
  // if(op == BoxDataOp::Add)
  //   message ="BoxData::operator+=";
  // else if(op == BoxDataOp::Subtract)
  //   message = "BoxData::operator-=(BoxData&)";
  // else if(op == BoxDataOp::Divide)
  //   message = "BoxData::operator/=";
  // else if(op == BoxDataOp::Multiply)
  //   message = "BoxData::operator*=";

  PR_TIME("BoxData::operator(BoxData&)");


  if(m_box == a_rhs.m_box)
  {
    unsigned int nsize = size();
    int blockSize = 512;
    int nBlocks = (nsize + blockSize - 1) / blockSize;
    protoLaunchKernelMemAsyncT<MEMTYPE,
                               opKernel<T,op>>( 
		    nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_rhs.m_rawPtr, nsize);
    PR_FLOPS(size());
  }
  else
  {

    Box domain = m_box & a_rhs.m_box;
    if (domain.size() <= 0){return;}

    for (unsigned int ee = 0; ee < E; ee++)
      for (unsigned int dd = 0; dd < D; dd++)
        {
          BoxData<T,1> srcbd = slice(a_rhs, 0, dd, ee);
          BoxData<T,1> dstbd = slice(*this, 0, dd, ee);

          unsigned int begin  = 0;
          unsigned int end    = domain.size();

          unsigned int stride = domain.size(0);
          unsigned int blocks = domain.flatten(0).size();

          T* d_srcPtr = (T*)srcbd.data();
          T* d_dstPtr = (T*)dstbd.data();

          protoLaunchKernelT<MEMTYPE,
                             boxdataIndexer<op, T,C>>( 
                            blocks, stride, begin, end, d_srcPtr, d_dstPtr,
                            domain, srcbd.box(), dstbd.box(),
                            Point::Zeros(), srcbd.box().size(), dstbd.box().size());
 

        }

    PR_FLOPS(domain.size()*E*D);
  }
}

// Operator += (BoxData&)
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
BoxData<T,C,MEMTYPE,D,E>& BoxData<T,C,MEMTYPE,D,E>::operator+=(const BoxData<T,C,MEMTYPE,D,E>& a_rhs)
{
  PR_TIME("BoxData::operator+=");
  this->operatorT<BoxDataOp::Add>(a_rhs);
  return *this;
}

// Operator -= (BoxData&)
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
BoxData<T,C,MEMTYPE,D,E>& BoxData<T,C,MEMTYPE,D,E>::operator-=(const BoxData<T,C,MEMTYPE,D,E>& a_rhs)
{
  PR_TIME("BoxData::operator-=");
  this->operatorT<BoxDataOp::Subtract>(a_rhs);
  return *this;
}

// Operator *= (BoxData&)
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
BoxData<T,C,MEMTYPE,D,E>& BoxData<T,C,MEMTYPE,D,E>::operator*=(const BoxData<T,C,MEMTYPE,D,E>& a_rhs)
{
  PR_TIME("BoxData::operator*=");
  this->operatorT<BoxDataOp::Multiply>(a_rhs);
  return *this;
}

// Operator /= (BoxData&)
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
BoxData<T,C,MEMTYPE,D,E>& BoxData<T,C,MEMTYPE,D,E>::operator/=(const BoxData<T,C,MEMTYPE,D,E>& a_rhs)
{
  PR_TIME("BoxData::operator/=");
  this->operatorT<BoxDataOp::Divide>(a_rhs);
  return *this;
}

template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
template <BoxDataOp op>
void BoxData<T,C,MEMTYPE,D,E>::operatorT(T a_scale)
{
  // string message;
  // if(op == BoxDataOp::Add)
  //   message ="BoxData::operator+=";
  // else if(op == BoxDataOp::Subtract)
  //   message = "BoxData::operator-=(BoxData&)";
  // else if(op == BoxDataOp::Divide)
  //   message = "BoxData::operator/=";
  // else if(op == BoxDataOp::Multiply)
  //   message = "BoxData::operator*=";

  PR_TIME("BoxData::operator");
  PR_FLOPS(size());


  for (unsigned int ee = 0; ee < E; ee++)
    for (unsigned int dd = 0; dd < D; dd++)
      {
        BoxData<T,1> dstbd = slice(*this, 0, dd, ee);

        unsigned int begin  = 0;
        unsigned int end    = m_box.size();

        unsigned int stride = m_box.size(0);
        unsigned int blocks = m_box.flatten(0).size();

        T* d_dstPtr = (T*)dstbd.data();
        size_t smem=0;
        protoLaunchKernelMemAsyncT<MEMTYPE,
                                   scalarIndexer<op, T,C>>( 
                     blocks, stride, smem, protoGetCurrentStream, begin, end, d_dstPtr, a_scale,
                     dstbd.box(), dstbd.box().size());
      }
}

template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
BoxData<T,C,MEMTYPE,D,E>& BoxData<T,C,MEMTYPE,D,E>::operator+=(T a_scale)
{
  PR_TIME("BoxData::operator+=(scale)");
  this->operatorT<BoxDataOp::Add>(a_scale);
  return *this;
}

template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
BoxData<T,C,MEMTYPE,D,E>& BoxData<T,C,MEMTYPE,D,E>::operator-=(T a_scale)
{
  PR_TIME("BoxData::operator-=(scale)");
  this->operatorT<BoxDataOp::Subtract>(a_scale);
  return *this;
}

template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
BoxData<T,C,MEMTYPE,D,E>& BoxData<T,C,MEMTYPE,D,E>::operator*=(T a_scale)
{
  PR_TIME("BoxData::operator*=(scale)");
  this->operatorT<BoxDataOp::Multiply>(a_scale);
  return *this;
}

template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
BoxData<T,C,MEMTYPE,D,E>& BoxData<T,C,MEMTYPE,D,E>::operator/=(T a_scale)
{
  PR_TIME("BoxData::operator/=(scale)");
  this->operatorT<BoxDataOp::Divide>(a_scale);
  return *this;
}

//=================================================================================================
// Utility Functions ||
//===================++



template <class T, MemType MEMTYPE, class B>
inline void setBoxDataVal(BoxData<T,1,MEMTYPE,1,1>& a_data, const B& a_box, T a_val)
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

template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::setVal(const T& a_val)
{
  PR_TIME("BoxData::setval");
  unsigned int nsize = size();
  int blockSize = 512;
  int nBlocks = (nsize + blockSize - 1) / blockSize;
  protoLaunchKernelMemAsyncT<MEMTYPE,
                             opKernel<T,BoxDataOp::Assign>>(
		  nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_val, nsize);
}

// Set Values Within Box
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::setVal(const T&   a_val,
                                      const Box&  a_box)
{
  for (unsigned int ee = 0; ee < E; ee++)
  {
    for (unsigned int dd = 0; dd < D; dd++)
    {
      for (unsigned int cc = 0; cc < C; cc++)
      {
        auto slicebd = slice(*this, cc, dd, ee);
        setBoxDataVal<T,MEMTYPE>(slicebd, a_box, a_val);
      }
    }
  }
}

// Set Component Values Within Box
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::setVal(const T& a_val,
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
  setBoxDataVal<T,MEMTYPE>(slicebd, a_box, a_val);
}

template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::setToZero()
{
    setVal(0);
}

/// Compute Reduction (Global)
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
template<Proto::Operation OP>
void BoxData<T,C,MEMTYPE,D,E>::reduce(Reduction<T,OP>& a_rxn) const
{
    PR_TIME("BoxData::reduce");
    a_rxn.reduce(m_rawPtr,size());
}

/// Compute Reduction (Componentwise)
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
template<Proto::Operation OP>
void BoxData<T,C,MEMTYPE,D,E>::reduce(Reduction<T, OP>& a_rxn, int a_c, int a_d, int a_e) const
{
    auto temp = slice(*this, a_c, a_d, a_e);
    temp.reduce(a_rxn);
}

/// Maximum Absolute Value (Global)
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::absMax(Reduction<T,Abs>& a_rxn) const
{
    reduce(a_rxn);
}

/// Maximum Absolute Value (Componentwise)
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::absMax(Reduction<T, Abs>& a_rxn, int a_c, int a_d, int a_e) const
{
    reduce(a_rxn, a_c, a_d, a_e);
}

/// Maximum Absolute Value (Global)
template <class T, unsigned int C,  MemType MEMTYPE,unsigned char D, unsigned char E>
T BoxData<T,C,MEMTYPE,D,E>::absMax() const
{
    Reduction<T,Abs> rxn(true); 
    reduce(rxn);
    return rxn.fetch();
}

/// Maximum Absolute Value (Componentwise)
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
T BoxData<T,C,MEMTYPE,D,E>::absMax(int a_c, int a_d, int a_e) const
{
    Reduction<T,Abs> rxn(true); 
    reduce(rxn, a_c, a_d, a_e);
    return rxn.fetch();
}

/// Minimum Value (Global)
template <class T, unsigned int C,  MemType MEMTYPE,unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::min(Reduction<T, Min>& a_rxn) const
{
    reduce(a_rxn);
}

/// Minimum Value (Componentwise)
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::min(Reduction<T, Min>& a_rxn,
    int a_c, int a_d, int a_e) const
{
    reduce(a_rxn, a_c, a_d, a_e);
}

/// Minimum Value (Global)
template <class T, unsigned int C,  MemType MEMTYPE,unsigned char D, unsigned char E>
T BoxData<T,C,MEMTYPE,D,E>::min() const
{
    Reduction<T,Min> rxn(true);
    reduce(rxn);
    return rxn.fetch();
}

/// Minimum Value (Componentwise)
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
T BoxData<T,C,MEMTYPE,D,E>::min(int a_c, int a_d, int a_e) const
{
    Reduction<T,Min> rxn(true);
    reduce(rxn, a_c, a_d, a_e);
    return rxn.fetch();
}

template <class T, unsigned int C, MemType MEMTYPE,unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::sum(Reduction<T, Sum>& a_rxn) const
{
    reduce(a_rxn);
}

template <class T, unsigned int C, MemType MEMTYPE,unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::sum(
    Reduction<T, Sum>& a_rxn,
    int a_c, int a_d, int a_e) const
{
    reduce(a_rxn, a_c, a_d, a_e);
}


/// Global Sum
template <class T, unsigned int C, MemType MEMTYPE,unsigned char D, unsigned char E>
T BoxData<T,C,MEMTYPE,D,E>::sum() const
{
    Reduction<T,Sum> rxn(true);
    reduce(rxn);
    return rxn.fetch();
}

template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
T BoxData<T,C,MEMTYPE,D,E>::sum(int a_c, int a_d, int a_e) const
{
    Reduction<T,Sum> rxn(true);
    reduce(rxn, a_c, a_d, a_e);
    return rxn.fetch();
}

template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
T BoxData<T,C,MEMTYPE,D,E>::integrate(T a_dx, int a_c, int a_d, int a_e) const
{
    return sum(a_c, a_c, a_d, a_e) * pow(a_dx, DIM);
}

/// Maximum Value (Global)
template <class T, unsigned int C,MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::max(Reduction<T, Max>& a_rxn) const
{
    reduce(a_rxn);
}

/// Maximum Value (Componentwise)
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::max(Reduction<T, Max>& a_rxn,
    int a_c, int a_d, int a_e) const
{
    reduce(a_rxn, a_c, a_d, a_e);
}



/// Maximum Value (Global)
template <class T, unsigned int C,MemType MEMTYPE, unsigned char D, unsigned char E>
T BoxData<T,C,MEMTYPE,D,E>::max() const
{
    Reduction<T,Max> rxn(true);
    reduce(rxn);
    return rxn.fetch();
}

/// Maximum Value (Componentwise)
template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
T BoxData<T,C,MEMTYPE,D,E>::max(int a_c, int a_d, int a_e) const
{
    Reduction<T,Max> rxn(true);
    reduce(rxn, a_c, a_d, a_e);
    return rxn.fetch();
}

template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::linearOut(void*                  a_buf,
                                         const ::Proto::Box&      a_bx,
                                         unsigned int a_startcomp,
                                         unsigned int a_numcomps) const
{
    CInterval comps(a_startcomp, a_startcomp + a_numcomps-1);
    linearOut(a_buf, a_bx, comps);
}

template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::linearOut(void*     a_buf,
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

    BoxData<T,C,MEMTYPE,D,E> tmp((T*)a_buf, a_box, C);
    copyTo(tmp, a_box, a_comps,  Point::Zeros(), CInterval( 0, a_comps.size(0)-1));
}


template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::linearIn(void*                  a_buf,
                                        const ::Proto::Box&      a_bx,
                                        unsigned int a_startcomp,
                                        unsigned int a_numcomps)
{
    CInterval comps(a_startcomp, a_startcomp + a_numcomps-1);
    linearIn(a_buf, a_bx, comps);
}

template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::linearIn(void*     a_buf,
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

   BoxData<T,C,MEMTYPE,D,E> tmp((T*) a_buf, a_box, C);
   tmp.copyTo(*this, a_box, CInterval(0, a_comps.size(0)-1),  Point::Zeros(), a_comps);
}

template <typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::print() const
{
  using std::cout;
  using std::endl;

  cout << "BoxData<" << C << ", " << (int)D << ", " << (int)E << ">" << endl;
  cout << "Domain: " << m_box << " Size: " << this->size() << " Values" << endl;
  if (this->size() <= 0){return;}
  for (int ee = 0; ee < E; ee++)
    for (int dd = 0; dd < D; dd++)
      for (int cc = 0; cc < C; cc++)
      {
        cout << "\tComponent (" << cc << ", " << dd << ", " << ee << "): ";
        cout << "(" << this->max(cc,dd,ee) << " / " << this->min(cc,dd,ee) << ")" << endl;
      }
}


template <typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::printData(int a_prec ) const
{
  printData(m_box,-1,-1,-1,a_prec);
}

template <typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::printData(const Box& a_box, int a_prec) const
{
  printData(a_box,-1,-1,-1, a_prec);
}

template <typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void BoxData<T,C,MEMTYPE,D,E>::printData(const Box& a_box, int a_c, int a_d, int a_e, int a_prec) const
{
  using std::cout;
  using std::endl;
#ifdef PR_MPI
    std::ofstream out;
    std::string fname = "pr_out_" + to_string(Proto::procID()) + ".out";
    out.open(fname);
#else
    std::ostream& out = std::cout;
#endif
  if (!defined())
  {
    out << "Undefined BoxData: No data to print." << endl;
    return;
  }
  if (m_box.size() <= 0)
  {
    out << "BoxData has empty domain: No data to print." << endl;
    return;
  }
  if (this->absMax() > 1.0e5 || this->absMax() < pow(10.0, -a_prec))
  {
      this->print();
      return;
  }

  PROTO_ASSERT(m_box.contains(a_box),
               "BoxData::printData(a_box,...) invalid. Input Box is not a subset of the domain of *this.");
  out << "Domain: " << m_box << " Subdomain: " << a_box << endl;
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

        out << endl << "Component (" << cc << ", " << dd << ", " << ee << "): " << endl;
        out << std::setfill('=') << std::setw(std::min((int)(a_box.size(0)*(width+1)+8),100)) << " " << endl;
#if DIM < 2
        out << "\t";
        for (auto iter = a_box.begin(); iter != a_box.end(); ++iter)
        {
          out << std::fixed << std::setprecision(prec) << std::setfill(' ');
          out << std::setw(width) << (*this)(*iter,cc,dd,ee) << " ";
        }
        out << endl;
#elif DIM < 3
        int nx = a_box.size(0);
        int ny = a_box.size(1);
        T* ptr = m_rawPtr + index(a_box.flatten(0).high(),cc,dd,ee);
        for (int jj = ny-1; jj >= 0; jj--)
        {
          out << "\t";
          for (int ii = 0; ii < nx; ii++)
          {
            out << std::fixed << std::setprecision(prec) << std::setfill(' ');
            out << std::setw(width) << ptr[ii] << " ";
          }
          out << endl;
          ptr -= m_box.size(0);
        }
#else
        Box layerBox = a_box.flatten(0).flatten(1);
        for (auto iter = layerBox.begin(); iter != layerBox.end(); ++iter)
        {
          out << "Layer (--,--,";
          for (int kk = 2; kk < DIM; kk++)
          {
            out << (*iter)[kk];
            if ((DIM-kk) > 1){cout << ", ";}
          }
          out << "):" << endl;
          out << std::setfill('-') << std::setw(std::min((int)(a_box.size(0)*6+8),100)) << " " << endl;

          Point start = (*iter);
          start[1] += (a_box.size(1)-1);
          T* ptr = m_rawPtr + index(start,cc,dd,ee);
          int nx = a_box.size(0);
          int ny = a_box.size(1);
          for (int jj = ny-1; jj >= 0; jj--)
          {
            out << "\t";
            for (int ii = 0; ii < nx; ii++)
            {
              out << std::fixed << std::setprecision(prec) << std::setfill(' ');
              out << std::setw(width) << ptr[ii] << " ";
            }
            out << endl;
            ptr -= m_box.size(0);
          }
        }
#endif
      }
  out << endl;
#ifdef PR_MPI
  out.close();
#endif
}
/// Alias (Non-Const)
template<class T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE>
BoxData<T,C,MEMTYPE,D,E> alias(BoxData<T,C,MEMTYPE,D,E>& a_original, const Point& shift)
{
  PR_TIME("alias(BoxData<T,C,MEMTYPE,D,E>&)");
  const Box& b=a_original.box();
  BoxData<T,C,MEMTYPE,D,E> rtn(a_original.m_data, a_original.m_rawPtr, b);
  rtn.shift(shift);
  return rtn; // will invoke the alias constructor, thus carry the pointer values, not a deep copy of the rtn object
}

/// Alias (Const)
template<class T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE>
const BoxData<T,C,MEMTYPE,D,E> alias(const BoxData<T,C,MEMTYPE,D,E>& a_original, const Point& shift)
{
  PR_TIME("alias(const BoxData<T,C,MEMTYPE,D,E>&)");
  const Box& b=a_original.box();
  BoxData<T,C,MEMTYPE,D,E>* src = const_cast<BoxData<T,C,MEMTYPE,D,E>*>(&a_original);
  const BoxData<T,C,MEMTYPE,D,E> rtn(src->m_data, src->m_rawPtr, b.shift(shift));
  return rtn; // will invoke the move constructor, thus carry the pointer values, not a deep copy of the rtn object
}

template<typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
BoxData<T,1,MEMTYPE> slice(const BoxData<T,C,MEMTYPE,D,E>& a_src,
                               unsigned int a_c,
                               unsigned int a_d,
                               unsigned int a_e)
{
  PR_TIME("slice(BoxData<T,C,MEMTYPE,D,E>&, int, int, int)");
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
  BoxData<T,1,MEMTYPE> rtn(a_src.getData(), ptr, b);
  return rtn;
}

template<typename T, unsigned int C, unsigned char CC, MemType MEMTYPE>
BoxData<T,CC,MEMTYPE> slice(const BoxData<T,C,MEMTYPE>& a_src,
                                unsigned int a_nstart)
{
  PR_TIME("slice(BoxData<T,C,1,1>&, int)");
  PROTO_ASSERT((a_nstart + CC <= C),
    "Error in slice(BoxData<T,C,1,1>, int) | Invalid slicing range")
  const Box& b = a_src.box();
  T* ptr = (T*)(a_src.data(b.low(),a_nstart));
  BoxData<T,CC,MEMTYPE> rtn(a_src.getData(), ptr, b);
  return rtn;
}

////
template <typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void forall_parse(Box& a_domain, const BoxData<T,C,MEMTYPE,D,E>& a_data)
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
template <typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE, typename... Rest>
void forall_parse_args(Box& a_domain, const BoxData<T,C,MEMTYPE,D,E>& a_first, Rest&&... a_rest)
{
  forall_parse(a_domain, a_first);
  forall_parse_args(a_domain, std::forward<Rest>(a_rest)...);
}

////
template <typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void forall_parse_args(Box& a_domain, const BoxData<T,C,MEMTYPE,D,E>& a_first)
{
  forall_parse(a_domain, a_first);
}

////
template <typename First , typename... Rest>
void forall_parse_args(Box& a_domain, First& a_first, Rest&&... a_rest)
{
  forall_parse(a_domain, a_first);
  forall_parse_args(a_domain, std::forward<Rest>(a_rest)...);
}

//========================================================================
// Base Forall functionality
//========================================================================

// Main version
template<typename Func, typename... Srcs>
void forallInPlaceBaseOp(unsigned long long int a_num_flops_point,
                         const char* a_timername,
                         const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME(a_timername);
  forallInPlaceBase(a_F, a_box, std::forward<Srcs>(a_srcs)...);
  unsigned long long int boxfloops = a_num_flops_point*a_box.size();
  PR_FLOPS(boxfloops);
}

// Main version
template<typename Func, typename... Srcs>
void forallInPlaceBase(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("forallInPlaceBase");
  protoForall(a_F, a_box, std::forward<Srcs>(a_srcs)...);
}

template<typename Func, typename... Srcs>
void forallInPlaceBaseOp_p(unsigned long long int a_num_flops_point,
                           const char* a_timername,
                           const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME(a_timername);
  protoForall_p(a_F, a_box, std::forward<Srcs>(a_srcs)...) ;

  unsigned long long int  boxfloops = a_num_flops_point*a_box.size();
  PR_FLOPS(boxfloops);
}


template<typename Func, typename... Srcs>
void forallInPlace_i(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  protoForall_i(a_F, a_box, std::forward<Srcs>(a_srcs)...) ;
}

template<typename Func, typename... Srcs>
void forallInPlaceBaseOp_i(unsigned long long int a_num_flops_point,
                           const char* a_timername,
                           const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME(a_timername);

  forallInPlace_i(a_F, a_box, std::forward<Srcs>(a_srcs)...);

  unsigned long long int boxfloops = a_num_flops_point*a_box.size();
  PR_FLOPS(boxfloops);
}



// Point input version
template<typename Func, typename... Srcs>
inline void forallInPlaceBase_p(const Func& a_F,  Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("forallInPlaceBase_p");
  return protoForall_p(a_F, a_box, std::forward<Srcs>(a_srcs)...) ;
}

//========================================================================
//  Frontend Forall Implementation
//========================================================================

template<typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE,
         typename Func, typename... Srcs>
BoxData<T,C,MEMTYPE,D,E> forallOp(unsigned long long int a_num_flops_point,   // Main ForallOp
                                  const char* a_timername,
                                  const Func& a_F, Srcs&&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect,std::forward<Srcs>(a_srcs)...);

  PROTO_ASSERT(intersect.size() > 0,
               "forall(...) inferred the domain Box to be empty.");

  BoxData<T,C,MEMTYPE,D,E> ret(intersect);
  forallInPlaceBaseOp(a_num_flops_point, a_timername, a_F, intersect, ret, std::forward<Srcs>(a_srcs)...);
  return ret;
}

// Main Forall
template<typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE,
         typename Func, typename... Srcs>
BoxData<T,C,MEMTYPE,D,E> forall(const Func& a_F, Srcs&&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect,std::forward<Srcs>(a_srcs)...);

  PROTO_ASSERT(intersect.size() > 0,
               "forall(...) inferred the domain Box to be empty.");

  BoxData<T,C,MEMTYPE,D,E> ret(intersect);
  forallInPlaceBase(a_F, intersect, ret, std::forward<Srcs>(a_srcs)...);
  return ret;
}

// Overload forallOp with Optional Box Argument
template<typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE,
         typename Func, typename... Srcs>
BoxData<T,C,MEMTYPE,D,E> forallOp(unsigned long long int a_num_flops_point,
                          const char* a_timername,
                          const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  BoxData<T,C,MEMTYPE,D,E> ret(a_box);
  forallInPlaceBaseOp(a_num_flops_point, a_timername, a_F, a_box, ret, std::forward<Srcs>(a_srcs)...);
  return ret;
}


// Overload forall with Optional Box Argument
template<typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE,
         typename Func, typename... Srcs>
BoxData<T,C,MEMTYPE,D,E> forall(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  BoxData<T,C,MEMTYPE,D,E> ret(a_box);
  forallInPlaceBase(a_F, a_box, ret, std::forward<Srcs>(a_srcs)...);
  return ret;
}

// In Place forallOp Version
template<typename Func, typename... Srcs>
void forallInPlaceOp(unsigned long long int a_num_flops_point,
                     const char* a_timername,
                     const Func& a_F, Srcs&&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect, std::forward<Srcs>(a_srcs)...);

  PROTO_ASSERT(intersect.size() > 0,
               "forall(...) inferred the domain Box to be empty.");

  forallInPlaceBaseOp(a_num_flops_point, a_timername, a_F, intersect, std::forward<Srcs>(a_srcs)...);
}

// In Place forall Version
template<typename Func, typename... Srcs>
void forallInPlace(const Func& a_F, Srcs&&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect, std::forward<Srcs>(a_srcs)...);

  PROTO_ASSERT(intersect.size() > 0,
               "forall(...) inferred the domain Box to be empty.");

  forallInPlaceBase(a_F, intersect, std::forward<Srcs>(a_srcs)...);
}

// In Place forallop Version with Box Argument
template<typename Func, typename... Srcs>
void forallInPlaceOp(unsigned long long int a_num_flops_point,
                     const char* a_timername,
                     const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  forallInPlaceBaseOp(a_num_flops_point, a_timername, a_F, a_box, std::forward<Srcs>(a_srcs)...);
}

// In Place forall Version with Box Argument
template<typename Func, typename... Srcs>
void forallInPlace(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  forallInPlaceBase(a_F, a_box, std::forward<Srcs>(a_srcs)...);
}

// Main forallOp version with Point Access
template<typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE,
         typename Func, typename... Srcs>
BoxData<T,C,MEMTYPE,D,E> forallOp_p(unsigned long long int a_num_flops_point,
                                    std::string a_timername,
                                    const Func& a_F, Srcs&&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect, std::forward<Srcs>(a_srcs)...);

  PROTO_ASSERT(intersect.size() > 0, "forall(...) inferred the domain Box to be empty.");

  BoxData<T,C,MEMTYPE,D,E> ret(intersect, true);
  forallInPlaceBaseOp_p(a_num_flops_point, a_timername, a_F, intersect, ret, std::forward<Srcs>(a_srcs)...);
  return ret;
}

// Main forall version with Point Access
template<typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE,
         typename Func, typename... Srcs>
BoxData<T,C,MEMTYPE,D,E> forall_p(const Func& a_F, Srcs&&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect, std::forward<Srcs>(a_srcs)...);

  PROTO_ASSERT(intersect.size() > 0, "forall(...) inferred the domain Box to be empty.");

  BoxData<T,C,MEMTYPE,D,E> ret(intersect, true);
  forallInPlaceBase_p(a_F, intersect, ret, std::forward<Srcs>(a_srcs)...);
  return ret;
}

// Main forallOp version with Point Access and Box Argument
template<typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE,
         typename Func, typename... Srcs>
BoxData<T,C,MEMTYPE,D,E> forallOp_p(unsigned long long int a_num_flops_point,
                                    const char* a_timername,
                                    const Func& a_F,  Box a_box, Srcs&&... a_srcs)
{
  BoxData<T,C,MEMTYPE,D,E> ret(a_box, true);
  forallInPlaceBaseOp_p(a_num_flops_point, a_timername, a_F, a_box, ret, std::forward<Srcs>(a_srcs)...);
  return ret;
}

// Main forall version with Point Access and Box Argument
template<typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE,
         typename Func, typename... Srcs>
BoxData<T,C,MEMTYPE,D,E> forall_p(const Func& a_F,  Box a_box, Srcs&&... a_srcs)
{
  BoxData<T,C,MEMTYPE,D,E> ret(a_box, true);
  forallInPlaceBase_p(a_F, a_box, ret, std::forward<Srcs>(a_srcs)...);
  return ret;
}

// In Place forallOp version with Point Access
template<typename Func, typename... Srcs>
inline void forallInPlaceOp_p(unsigned long long int a_num_flops_point,
                              const char* a_timername,
                              const Func& a_F, Srcs&&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect, std::forward<Srcs>(a_srcs)...);

  PROTO_ASSERT(intersect.size() > 0,
               "forall(...) inferred the domain Box to be empty.");
  forallInPlaceBaseOp_p(a_num_flops_point, a_timername, a_F, intersect, std::forward<Srcs>(a_srcs)...);
}

// In Place forall version with Point Access
template<typename Func, typename... Srcs>
inline void forallInPlace_p(const Func& a_F, Srcs&&... a_srcs)
{
  Box intersect;
  forall_parse_args(intersect, std::forward<Srcs>(a_srcs)...);

  PROTO_ASSERT(intersect.size() > 0,
               "forall(...) inferred the domain Box to be empty.");
  forallInPlaceBase_p(a_F, intersect, std::forward<Srcs>(a_srcs)...);
}

// In Place forallOp version with Point Access and Box Argument
template<typename Func, typename... Srcs>
inline void forallInPlace_p(unsigned long long int a_num_flops_point,
                            const char* a_timername,
                            const Func& a_F,  Box a_box, Srcs&&... a_srcs)
{
  forallInPlaceBase_p(a_num_flops_point, a_timername, a_F, a_box, std::forward<Srcs>(a_srcs)...);
}


// In Place forall version with Point Access and Box Argument
template<typename Func, typename... Srcs>
inline void forallInPlace_p(const Func& a_F,  Box a_box, Srcs&&... a_srcs)
{
  forallInPlaceBase_p(a_F, a_box, std::forward<Srcs>(a_srcs)...);
}


template<typename FuncStruct>
struct emptyIndexer
{
#ifdef PROTO_CUDA
  __device__ static 
  void gpu(FuncStruct a_body)
  {
  }
#endif
  static inline void cpu(FuncStruct a_body){}
};

template<typename... T>
void emptyT(T... in)
{
}


template<typename FuncStruct, typename... Srcs>
struct structIndexer
{
  #ifdef PROTO_CUDA
  __device__ static 
  void gpu(int a_begin, int a_end, FuncStruct a_body, Srcs... a_srcs)
  {
    int idx = threadIdx.x + blockIdx.x*blockDim.x;
    if (idx >= a_begin && idx < a_end)
      {
        a_body.op(a_srcs...);
      }
  }
  #endif
  
  static void cpu(int begin, int end, FuncStruct a_body, Srcs... a_srcs)
  {
    for(int i=begin;i<end; i++)
      {
        a_body.op(a_srcs...);
        emptyT(var_incr(a_srcs)...);
      }
  }
};

template<typename T>
inline T p_ref(const T& a_s, const Point& a_p) {return a_s;}

template<typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
inline Var<T,C,MEMTYPE,D,E> p_ref(Var<T,C,MEMTYPE,D,E>&  a_data, const Point& a_p){  a_data += a_p; return a_data;}

template<typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
inline Var<T,C,MEMTYPE,D,E> p_ref(const Var<T,C,MEMTYPE,D,E>& a_data, const Point& a_p){a_data += a_p; return a_data;}

template<typename T>
inline T var_incr(const T& t){return t;}

template<typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
Var<T, C, MEMTYPE, D, E>& var_incr(Var<T, C, MEMTYPE, D, E>& t) {return ++t;}


template<typename Func, typename... T>
inline void pencilFunc(const Func& F, int count, T... vars)
{
  for(int i=0; i < count; ++i)
  {
    F(vars...);
    emptyT(var_incr(vars)...);
  }
}

template<typename Func, typename... T>
inline void pencilFunc_p(const Func& F, int count, Point a_pt, T... vars)
{
  for(int i=0; i < count; ++i, ++a_pt[0])
  {
    F(a_pt,vars...);
    emptyT(var_incr(vars)...);
  }
}

template<typename Func, typename... T>
inline void pencilFunc_i(const Func& F, int count, Point a_pt, T... vars)
{
  for(int i=0; i < count; ++i, ++a_pt[0])
  {
    //F(a_pt.m_tuple,vars...);
    F(a_pt,vars...);
    emptyT(var_incr(vars)...);
  }
}

template<typename Func, typename... Srcs>
struct indexer {

  static void cpu(int a_begin, int a_end, const Box& a_box, Func a_body, Srcs... a_srcs)
  {
    Box cross = a_box.flatten(0);
    int npencil = a_box.size(0);
    
    auto last = a_box.low();
    
    for (auto iter = cross.begin(); iter != cross.end(); ++iter)
      {
        pencilFunc(a_body, npencil, p_ref(a_srcs,*iter-last)...);
        last = *iter;
      }
  }
#ifdef PROTO_CUDA
  __device__ static
  void gpu(int a_begin, int a_end, const Box& a_box, Func a_body, Srcs... a_srcs)
  {
    int idx = threadIdx.x + blockIdx.x*blockDim.x;
    if (idx >= a_begin && idx < a_end)
      {
        a_body(a_srcs...);
      }
  }
#endif
};

template<typename Func, typename... Srcs>
struct indexer_p {
  
  static void cpu(int a_begin, int a_end, const Box& a_box,
                  Func a_body, Srcs... a_srcs)
  {
    
    Box cross = a_box.flatten(0);
    int npencil = a_box.size(0);

    auto last = a_box.low();
    for (auto iter = cross.begin(); iter != cross.end(); ++iter)
      {
        pencilFunc_p(a_body, npencil, *iter, (p_ref(a_srcs,*iter-last))...);
        last = *iter;
      }
  }

  #ifdef PROTO_CUDA
  __device__ static
  void gpu(int a_begin, int a_end, const Box& a_box,
           Func a_body, Srcs... a_srcs)
  {
    int idx = threadIdx.x + blockIdx.x*blockDim.x;
    int idy = threadIdx.y + blockIdx.y*blockDim.y;
    int idz = threadIdx.z + blockIdx.z*blockDim.z;
    if (idx >= a_begin && idx < a_end)
      {
        Point p = a_box(idx,idy,idz);
        a_body(p,a_srcs...);
      }
  }
#endif
  
};

template<typename Func, typename... Srcs>
struct indexer_i {
  static void cpu(int a_begin, int a_end, const Box& a_box,
                  Func a_body, Srcs... a_srcs)
  {
    Box cross = a_box.flatten(0);
    int npencil = a_box.size(0);
    auto last = a_box.low();
    
    for (auto iter = cross.begin(); iter != cross.end(); ++iter)
      {
        pencilFunc_i(a_body, npencil, *iter, p_ref(a_srcs,*iter-last)...);
        last = *iter;
      }
  }

#ifdef PROTO_CUDA
  __device__ static
  void gpu(int a_begin, int a_end, const Box& a_box,
           Func a_body, Srcs... a_srcs)
  {
    int idx = threadIdx.x + blockIdx.x*blockDim.x;
    int idy = threadIdx.y + blockIdx.y*blockDim.y;
    int idz = threadIdx.z + blockIdx.z*blockDim.z;
    
    if (idx >= a_begin && idx < a_end)
      {
        Point p = a_box(idx,idy,idz);
        a_body(p.m_tuple,a_srcs...);
      }
  }
#endif
};

#endif // END INCLUDE GUARD


template<typename T>
struct getMemType
{
  //static constexpr MemType type() {return Proto::MemType::BOTH;}
  //static MemType type_eval()
  //{
  //	  return Proto::MemType::BOTH;
  // }
  static constexpr MemType type() {return MEMTYPE_DEFAULT;}
  static MemType type_eval()
  {
	  return MEMTYPE_DEFAULT;
  }
};

template<typename T, unsigned int C, MemType MEMTYPE,unsigned int D, unsigned int E>
struct getMemType<BoxData<T,C,MEMTYPE,D,E>>
{
  static constexpr MemType type()
  {
    return MEMTYPE;
  }
  static MemType type_eval()
  {
	  return MEMTYPE;
  }
};

template<typename T, unsigned int C, MemType MEMTYPE,unsigned int D, unsigned int E>
struct getMemType<Var<T,C,MEMTYPE,D,E>>
{
  static constexpr MemType type()
  {
    return MEMTYPE;
  }
};

template<typename First>
constexpr MemType getMemTypeFromSrcs()
{
  return getMemType<First>::type();
 
}

template<typename First, typename Second, typename... Rest>
constexpr MemType getMemTypeFromSrcs()
{
  return (MemType)(getMemType<First>::type() & getMemTypeFromSrcs<Second, Rest...>());
}


inline void v100tuningBox(Box a_box, dim3& blocks, int& blockSize)
 {
 #if DIM == 3
   if(a_box.size(0)<=512)
   {
    blocks = dim3(1,a_box.size(1),a_box.size(2));
    blockSize = a_box.size(0);
  }
  else
  {
    blockSize = 512;
    blocks = dim3((a_box.size(0)+blockSize-1)/blockSize, a_box.size(1),a_box.size(2));
  }

#else
  if(a_box.size(0)<=512)
  {
    blocks = dim3(1,a_box.size(1),1);
    blockSize = a_box.size(0);
  }
  else
  {
    blockSize = 512;
    blocks = dim3((a_box.size(0)+blockSize-1)/blockSize, a_box.size(1),1);
  }
#endif
}

inline void v100tuningBox(int a_stride, Box a_cross, dim3& blocks, int& blockSize)
{

#if DIM == 3
  if(a_stride<=512)
  {
    blocks = dim3(1,a_cross.size(1),a_cross.size(2));
    blockSize = a_stride;
  }
  else
  {
    blockSize = 512;
    blocks = dim3((a_stride+blockSize-1)/blockSize, a_cross.size(1),a_cross.size(2));
  }

#else
  if(a_stride<=512)
  {
    blocks = dim3(1,a_cross.size(1),1);
    blockSize = a_stride;
  }
  else
  {
    blockSize = 512;
    blocks = dim3((a_stride+blockSize-1)/blockSize, a_cross.size(1),1);
  }
#endif
}



template<typename Func, typename... Srcs>
inline void makeVars(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("BoxData::makevars");
  int N = a_box.size();
  int stride = a_box.size(0);

#if DIM == 3
  dim3 blocks(1,a_box.size(1),a_box.size(2));
#else
  dim3 blocks(1,a_box.size(1),1);
#endif
  {
    PR_TIME("indexer");
    size_t smem = 0;

#ifdef superDebug
    std::cout << "Try Kernel Name is: " << a_F.myname << std::endl;
#endif

    assert(stride<1024);

    protoLaunchKernelMemAsyncT<getMemTypeFromSrcs<Srcs...>(),indexer<Func, Srcs...>>( blocks, stride, smem, protoGetCurrentStream, 0, N, a_box, a_F, std::forward<Srcs>(a_srcs)...);


#ifdef superDebug
    std::cout << "Success Kernel Name is: " << a_F.myname << std::endl;
#endif
  }
}

template<typename Func, typename... Srcs>
inline void makeVarsStream(protoStream_t& a_stream, const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("BoxData::makevarsstream");
  int N = a_box.size();
  int stride = a_box.size(0);
  size_t smem = 0;
#if DIM == 3
  dim3 blocks(1,a_box.size(1),a_box.size(2));
#else
  dim3 blocks(1,a_box.size(1),1);
#endif
  {
    PR_TIME("indexer");

    assert(stride<1024);
    protoLaunchKernelMemAsyncT<getMemTypeFromSrcs<Srcs...>(), indexer<Func, Srcs...>>(blocks, stride, smem, a_stream, 0, N, a_box, a_F, std::forward<Srcs>(a_srcs)...);
  }
}


template<typename FuncStruct, typename... Srcs>
inline void makeVarsStruct(protoStream_t& a_stream, const FuncStruct& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("BoxData::makevarsstruct");
  int N = a_box.size();
  int stride = a_box.size(0);
#if DIM == 3
  dim3 blocks(1,a_box.size(1),a_box.size(2));
#else
  dim3 blocks(1,a_box.size(1),1);
#endif
  size_t smem = 0;
  {
    PR_TIME("structIndexer");

    assert(stride<1024);
    protoLaunchKernelMemAsyncT<getMemTypeFromSrcs<Srcs...>,structIndexer<FuncStruct, Srcs...>>( blocks, stride, smem, a_stream, 0, N, a_F, std::forward<Srcs>(a_srcs)...);
  }
}

template<typename FuncStruct, typename... Srcs>
inline void makeVarsEmptyIndexer(protoStream_t& a_stream, const FuncStruct& a_F, Box a_box, Box a_srcBox)
{
  //a_srcs is coming into here as vars
  PR_TIME("BoxData::makevarsEmptyIndexer");

  int stride = a_box.size(0);
  int blocks = a_box.size(1);
  /*
#if DIM==3
  int Nplane    = a_box.flatten(2).size();
  int Nz        = a_box.size(2);
  int zinc      = a_srcBox.flatten(2).size();
#else
  int Nplane    = a_box.size();
  int Nz        = 1;
  int zinc      = 1; //has to be non-zero or we have an infinite loop
#endif
  */
  size_t smem = 0;
  {
    PR_TIME("zincStructIndexer");
    //Doesn't matter what memtype we pass in because this is an empyt kernel
    assert(stride<1024);
    protoLaunchKernelMemAsyncT<getMemTypeFromSrcs<Srcs...>,emptyIndexer<FuncStruct>>( blocks, stride, smem, a_stream, a_F);
  }
}

template<typename Func, typename... Srcs>
inline void makeVars_p(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("BoxData::makevars_p");
  int N = a_box.size();
  int stride = a_box.size(0);
  Box cross = a_box.flatten(0);
//  int blocks =  cross.size();
#if DIM == 3
  dim3 blocks(1,cross.size(1),cross.size(2));
#else
  dim3 blocks(1,cross.size(1),1);
#endif

  {
    PR_TIME("indexer_p");
    // MemType memtype=getMemTypeFromSrcs(std::forward<Srcs>(a_srcs)...);
    protoLaunchKernelT<getMemTypeFromSrcs<Srcs...>(),indexer_p<Func, Srcs...>>(blocks, stride, 0, N, a_box,
                      a_F, std::forward<Srcs>(a_srcs)...);
  }
}



template<typename Func, typename... Srcs>
inline void makeVars_i(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("BoxData::makevars_i");
  int N = a_box.size();
  int stride = a_box.size(0);
  Box cross = a_box.flatten(0);
#if DIM == 3
  dim3 blocks(1,cross.size(1),cross.size(2));
#else
  dim3 blocks(1,cross.size(1),1);
#endif
  {
    PR_TIME("indexer_i");
    protoLaunchKernelT<getMemTypeFromSrcs<Srcs...>(),indexer_i<Func, Srcs...>>(blocks, stride, 0, N, a_box,
                      a_F, std::forward<Srcs>(a_srcs)...);
  }
}


template<typename T>
inline T p_ref_cuda(const T& a_s, const Point& a_p) {return a_s;}


template<typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
inline Var<T,C,MEMTYPE,D,E> p_ref_cuda(BoxData<T,C,MEMTYPE,D,E>& a_data, const Point& a_p)
{
  return a_data.var(a_p);
}

template<typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
inline Var<T,C,MEMTYPE,D,E> p_ref_cuda(const BoxData<T,C,MEMTYPE,D,E>& a_data, const Point& a_p)
{
  return a_data.var(a_p);
}



template<typename Func, typename... Srcs>
inline void protoForall(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("protoForall");
  makeVars(a_F, a_box, p_ref_cuda(a_srcs, a_box.low())...);
}

template<typename Func, typename... Srcs>
inline void protoForallStream(protoStream_t& a_stream, const Func& a_F, Box a_box, Srcs&... a_srcs)
{
  PR_TIME("protoForallStream");
  makeVarsStream(a_stream, a_F, a_box, p_ref_cuda(a_srcs, a_box.low())...);
}

template<typename FuncStruct, typename... Srcs>
inline void protoForallStruct(protoStream_t& a_stream, const FuncStruct& a_F, Box a_box, Srcs&... a_srcs)
{
  PR_TIME("protoForallStruct");
  makeVarsStruct(a_stream, a_F, a_box, p_ref_cuda(a_srcs, a_box.low())...);
}

template<typename FuncStruct, typename... Srcs>
inline void protoForallEmptyIndexer(protoStream_t& a_stream, const FuncStruct& a_F, Box a_box, Box a_srcBox)
{
  PR_TIME("protoForallEmptyIndexer");
  makeVarsEmptyIndexer(a_stream, a_F, a_box, a_srcBox);
}


template<typename Func, typename... Srcs>
inline void protoForall_p(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("protoForall_p");
  makeVars_p(a_F, a_box, p_ref_cuda(a_srcs, a_box.low())...);
}


template<typename Func, typename... Srcs>
inline void protoForall_i(const Func& a_F, Box a_box, Srcs&&... a_srcs)
{
  PR_TIME("protoForall_i");
  makeVars_i(a_F, a_box, p_ref_cuda(a_srcs, a_box.low())...);
}

