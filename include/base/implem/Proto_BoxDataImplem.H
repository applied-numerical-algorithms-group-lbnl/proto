
#define PROTO_KERNEL_START inline ACCEL_DECORATION
#define PROTO_KERNEL_END(local_name, app_name)                  \
struct struct_##local_name {   \
    template <typename... T> \
    inline ACCEL_DECORATION void operator()(T... args) const { local_name(args...);} \
    const char* myname = #app_name; \
};        \
static struct_##local_name app_name;
#define PROTO_LAMBDA ACCEL_DECORATION

//=================================================================================================
// Constructors ||
//==============++

// Default Constructor
template <class T, unsigned int C,MemType MEM, unsigned int D, unsigned int E> 
BoxData<T,C,MEM,D,E>::BoxData() : m_rawPtr(nullptr)
{
    m_alias = false;
}

// Box Constructor
template <class T, unsigned int C,MemType MEM, unsigned int D, unsigned int E> 
BoxData<T,C,MEM,D,E>::BoxData(const Box& a_box)
{
    define(a_box);
}

// Initialization Constructor
template <class T, unsigned int C,MemType MEM, unsigned int D, unsigned int E> 
BoxData<T,C,MEM,D,E>::BoxData(const Box& a_box, T a_init) 
{
    define(a_box);
    setVal(a_init);
}

// Move Constructor
// TODO: *NEED* to make sure that move semantics are functioning properly
template <class T, unsigned int C,MemType MEM, unsigned int D, unsigned int E> 
BoxData<T,C,MEM,D,E>::BoxData(BoxData<T,C,MEM,D,E>&& a_src)
  :   m_box(a_src.m_box),
      m_data(std::move(a_src.m_data)),
      m_rawPtr(a_src.m_rawPtr)
{
    m_alias = false;
}

// Define
template <class T, unsigned int C,MemType MEM, unsigned int D, unsigned int E> 
void BoxData<T,C,MEM,D,E>::define(const Box& a_box)
{
    PR_TIME("BoxData::define(Box)");
    m_box=a_box;
    auto& stack = Stack<MEM>::getStack();
    if (stack.enabled())
    {
        m_rawPtr = (T*)(stack.alloc(size()*sizeof(T)));
        PROTO_ASSERT(m_rawPtr != nullptr,
            "BoxData::define | Error: Stack returned a nullptr pointer when it shouldn't have.");
        // Deletion is managed by Stack
        m_data = std::shared_ptr<T>(m_rawPtr, &(null_deleter_boxdata<MEM>));
        m_stackAlloc = true;
    } else {
        m_rawPtr = (T*)proto_malloc<MEM>(size()*sizeof(T));
        m_data = std::shared_ptr<T>(m_rawPtr, [](T* p){ proto_free<MEM>(p);}); 
        m_stackAlloc = false;
    }
    m_alias = false;
}

// Alias Define
template <class T, unsigned int C,MemType MEM, unsigned int D, unsigned int E> 
template <unsigned int CC>
void BoxData<T,C,MEM,D,E>::define(
        BoxData<T, CC, MEM, 1, 1>& a_src,
        unsigned int&                  a_comp)
{
    PR_TIME("BoxData::define(Alias)");
    m_box    = a_src.box();
    m_rawPtr = a_src.data(a_comp);
    m_data   = a_src.aliasData();
    m_alias  = true;
}           

// Raw Pointer Constructor
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>::BoxData(const T* a_ptr, const Box& a_box, int a_ncomp)
{
    define(a_ptr, a_box, a_ncomp);
}

// Raw Pointer Define
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::define(const T* a_ptr, const Box& a_box, int a_ncomp)
{
    PR_TIME("BoxData::define(RawPtr)");
    PROTO_ASSERT((a_ncomp==C && D==1 && E==1),
        "BoxData::define(T*, Box&, int) | Error: Component mismatch in boxdata alias");
    T* castPtr = const_cast<T*>(a_ptr);
    m_data = ::std::shared_ptr<T>(castPtr, &(null_deleter_boxdata<MEM>));
    m_rawPtr = castPtr;
    m_box = a_box;
    m_alias = true;
}

// Stencil Constructor
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>::BoxData(LazyStencil<T, C, MEM, D, E>&& a_op)
{
    Box range = a_op.m_range;
    if (range.empty()) { range = a_op.inferredRange(); }
    define(range);
    a_op.apply(*this, true);
    PROTO_ASSERT(!containsInfOrNAN(),
        "A Stencil output generated values of inf or nan");
}

template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>::BoxData(LazyInterpStencil<T, C, MEM, D, E>&& a_op)
{
    Box range = a_op.range;
    if (range.empty()) { range = a_op.inferredRange(); }
    define(range);
    a_op.apply(*this, true);
    PROTO_ASSERT(!containsInfOrNAN(),
        "An InterpStencil output generated values of inf or nan");
}


// Slice Constructor (Private)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>::BoxData(::std::shared_ptr<T> a_data,const  T* a_ptr, const Box& a_box)
{
    m_data = a_data;
    m_box = (a_box);
    m_rawPtr = (T*)a_ptr;
    m_alias = true;
}

// Destructor
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>::~BoxData()
{}


//=============================================================================================
// Data Movement ||
//===============++

// Move Assignment Operator
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>& BoxData<T,C,MEM,D,E>::operator=(BoxData<T,C,MEM,D,E>&& a_src)
{
    PR_TIME("BoxData::operator=(BoxData&& (move assign)");
    if (!a_src.isAlias(*this))
    {
        if (m_alias)
        {
            a_src.copyTo(*this);
        } else {
            std::swap<Box>(m_box,a_src.m_box);
            std::swap<std::shared_ptr<T>>(m_data, a_src.m_data);
            m_rawPtr = a_src.m_rawPtr;
            a_src.m_rawPtr = nullptr;
#ifdef PROTO_MEM_CHECK
            BoxDataMemCheck::numMoveAssign()++;
#endif
        }
    }
    return *this;
}

// Copy on Intersection
template <class T, unsigned int C, MemType MEM_SRC, unsigned int D, unsigned int E>
template <MemType MEM_DST>
void BoxData<T,C,MEM_SRC,D,E>::copyTo(BoxData<T,C,MEM_DST,D,E>& a_dest) const
{
    PR_TIME("BoxData::copyTo(BoxData&)");
    Box domain = m_box & a_dest.box();
    if (domain.empty()) { return; }

    copyTo(a_dest, domain, {{0,C-1},{0,D-1},{0,E-1}},
            Point::Zeros(),{{0,C-1},{0,D-1},{0,E-1}});
}

// Copy Region
template <typename T, unsigned int C, MemType MEM_SRC, unsigned int D, unsigned int E>
template <MemType MEM_DST>
void BoxData<T,C,MEM_SRC,D,E>::copyTo(BoxData<T,C,MEM_DST,D,E>&   a_dest,
        const Box&              a_srcBox,
        const Point&            a_destShift) const
{
    PR_TIME("BoxData::copyTo(BoxData&, Box&, Point&)");
    if (a_srcBox.empty()) { return; }
    copyTo(a_dest,a_srcBox,{{0,C-1},{0,D-1},{0,E-1}},
            a_destShift,{{0,C-1},{0,D-1},{0,E-1}});
}

// Copy From
// TODO: Get rid of this if we can
template <class T,unsigned int C, MemType MEM, unsigned int D, unsigned int E>
template<unsigned int Csrc>
void BoxData<T,C,MEM,D,E>::copy(const BoxData<T,Csrc,MEM,D,E>&  a_src,
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

// Two Box Copy
template< typename T, unsigned int C, MemType MEM_SRC, unsigned int D, unsigned int E>
template <MemType MEM_DST>
void BoxData<T,C,MEM_SRC,D,E>::copyTo(
        BoxData<T,C,MEM_DST,D,E>& a_dst,
        const Box&  a_srcBox,
        const Box&  a_dstBox) const
{
    CInterval cint(0,C-1);
    //Box srcBox = a_srcBox&a_dstBox;
    Point shift = Point::Zeros();
    this->copyTo(a_dst,a_srcBox,cint,shift,cint);
}

// BoxData Indexer
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
/*
                else if (Op == BoxDataOp::ThresholdMax)
                {
                    tmp_dst[dstOff] = (tmp_dst[dstOff] > tmp_src[srcOff]) ? tmp_src[srcOff] : tmp_dst[dstOff];
                } else if (Op == BoxDataOp::ThresholdMin)
                {
                    tmp_dst[dstOff] = (tmp_dst[dstOff] < tmp_src[srcOff]) ? tmp_src[srcOff] : tmp_dst[dstOff];
                }
  */
                else
                {
                    printf("boxdataIndexer error: bogus boxdata op input!!!\n");
                }
            }
        } // end for idx
    } // end cpu
#ifdef PROTO_ACCEL
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
                /*
                else if (Op == BoxDataOp::ThresholdMax)
                {
                    tmp_dst[dstOff] = (tmp_dst[dstOff] > tmp_src[srcOff]) ? tmp_src[srcOff] : tmp_dst[dstOff];
                } else if (Op == BoxDataOp::ThresholdMin)
                {
                    tmp_dst[dstOff] = (tmp_dst[dstOff] < tmp_src[srcOff]) ? tmp_src[srcOff] : tmp_dst[dstOff];
                }
                */
                else
                {
                    printf("boxdataIndexer error: bogus boxdata op input!!!\n");
                }
            }
        }
    }
#endif
}; // end boxData indexer

// Scalar Indexer
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
                /*
                else if (Op == BoxDataOp::ThresholdMax)
                {
                    tmp_dst[dstOff] = (tmp_dst[dstOff] > a_scalar) ? a_scalar : tmp_dst[dstOff];
                } 
                else if (Op == BoxDataOp::ThresholdMin)
                {
                    tmp_dst[dstOff] = (tmp_dst[dstOff] < a_scalar) ? a_scalar : tmp_dst[dstOff];
                }
                */
                else
                {
                    printf("scalarIndexer error: bogus boxdata op input!!!\n");
                }
            }
        }
    }
 #ifdef PROTO_ACCEL
    __device__ static  void gpu(
            int a_begin, int a_end, T* a_dst, T a_scalar,
            const Box& a_box, unsigned int dstSize)
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
                /*
                else if (Op == BoxDataOp::ThresholdMax)
                {
                    tmp_dst[dstOff] = (tmp_dst[dstOff] > a_scalar) ? a_scalar : tmp_dst[dstOff];
                } 
                else if (Op == BoxDataOp::ThresholdMin)
                {
                    tmp_dst[dstOff] = (tmp_dst[dstOff] < a_scalar) ? a_scalar : tmp_dst[dstOff];
                }
                */
                else
                {
                    printf("scalarIndexer error: bogus boxdata op input!!!\n");
                }
            }
        }
    }
#endif
}; // end scalar indexer


// General Copy
template< typename T, unsigned int C, MemType MEM_SRC, unsigned int D, unsigned int E>
template< unsigned int  Cdest, MemType MEM_DST, unsigned int Ddest, unsigned int Edest>
void BoxData<T,C,MEM_SRC,D,E>::copyTo(
        BoxData<T,Cdest,MEM_DST,Ddest,Edest>& a_dst,
        const Box&    a_srcBox,
        CInterval     a_srcComps,
        const Point&  a_dstShift,
        CInterval     a_dstComps) const
{
    PROTO_ASSERT(containsInfOrNAN(),
        "Detected Inf or NAN values in source data of copyTo operation");
    PROTO_ASSERT(a_srcComps == a_dstComps,
            "BoxData::copyTo | Error: source and destination components are restricted to \
            be the same pending a bug fix in boxDataIndexer.");
    PROTO_ASSERT(C == Cdest,
            "BoxData::copyTo | Error: source and destination components are restricted to \
            be the same pending a bug fix in boxDataIndexer.");
    PROTO_ASSERT(D == Ddest,
            "BoxData::copyTo | Error: source and destination components are restricted to \
            be the same pending a bug fix in boxDataIndexer.");
    PROTO_ASSERT(E == Edest,
            "BoxData::copyTo | Error: source and destination components are restricted to \
            be the same pending a bug fix in boxDataIndexer.");
    
    PR_TIME("BoxData::CopyTo");
    PROTO_ASSERT(contains(a_srcComps),
        "BoxData::copyTo | Error: \
        Source CInterval must be contained in source data component space.");
    PROTO_ASSERT(a_dst.contains(a_dstComps),
        "BoxData::copyTo | Error: \
        Destination CInterval must be contained in destination data component space.");
    for (int ii = 0; ii < 3; ii++)
    {
        PROTO_ASSERT(a_srcComps.size(ii) == a_dstComps.size(ii),
            "BoxData::copyTo | Error: \
            Src interval size %i =/= dst interval size %i on axis %i.",
            a_srcComps.size(ii), a_dstComps.size(ii),ii);
    }
    PROTO_ASSERT(m_box.containsBox(a_srcBox),
        "BoxData::copyTo | Error: \
        Source Box input is not a subset of source array domain.");

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
        if (MEM_SRC == DEVICE && MEM_DST == DEVICE)
        {
            protoMemcpyAsync(MEM_SRC,dstPtr, castSrc, nsize*sizeof(T),
                    protoMemcpyDeviceToDevice,
                    protoGetCurrentStream);
        }
        proto_memcpy<MEM_SRC, MEM_DST>(castSrc, dstPtr, nsize*sizeof(T));
#ifdef PROTO_MEM_CHECK
        BoxDataMemCheck::numCopies()++;
        BoxDataMemCheck::numFastCopies()++;
#endif
    }
    else
    {
        Box dstBox = a_srcBox.shift(a_dstShift) & a_dst.box();
        if (dstBox.empty()) { return; }
        Box srcBox = dstBox.shift(-a_dstShift);
        for (unsigned int ee = 0; ee < a_srcComps.size(2); ee++)
        {
            for (unsigned int dd = 0; dd < a_srcComps.size(1); dd++)
            {
                // DEVICE <--> DEVICE (INDEXER KERNEL LAUNCH)
                if ((MEM_SRC == DEVICE) && (MEM_DST == DEVICE))
                {
                    // No loop over C component needed because each kernel launch
                    // executes copy over all components on the C axis
                    int ccSrc = a_srcComps.low(0);
                    int ddSrc = dd + a_srcComps.low(1);
                    int eeSrc = ee + a_srcComps.low(2);
                    int ccDst = a_dstComps.low(0);
                    int ddDst = dd + a_dstComps.low(1);
                    int eeDst = ee + a_dstComps.low(2);

                    BoxData<T,1,MEM_SRC> srcbd = slice(*this, ccSrc, ddSrc, eeSrc);
                    BoxData<T,1,MEM_DST> dstbd = slice(a_dst, ccDst, ddDst, eeDst);

                    T* d_srcPtr = (T*)(srcbd.data());
                    T* d_dstPtr = (T*)(dstbd.data());

                    unsigned int begin  = 0;
                    unsigned int end    = dstBox.size();
                    unsigned int stride = dstBox.size(0);
                    unsigned int blocks = dstBox.flatten(0).size();

                    protoLaunchKernelMemAsyncT<MEM_SRC,
                        boxdataIndexer<BoxDataOp::Copy,T,C>>(
                                blocks, stride,0,
                                protoGetCurrentStream, begin, end,
                                d_srcPtr, d_dstPtr,
                                srcBox, srcbd.box(), dstbd.box(),
                                a_dstShift, srcbd.box().size(), dstbd.box().size());
                // ALL OTHER CASES (PENCILING)
                } else {
                    // Unlike the DEVICE <--> DEVICE case, here we need the loop over C
                    for (unsigned int cc = 0; cc < a_srcComps.size(0); cc++)
                    {
                        int ccSrc = cc + a_srcComps.low(0);
                        int ddSrc = dd + a_srcComps.low(1);
                        int eeSrc = ee + a_srcComps.low(2);
                        int ccDst = cc + a_dstComps.low(0);
                        int ddDst = dd + a_dstComps.low(1);
                        int eeDst = ee + a_dstComps.low(2);

                        BoxData<T,1,MEM_SRC> srcbd = slice(*this, ccSrc, ddSrc, eeSrc);
                        BoxData<T,1,MEM_DST> dstbd = slice(a_dst, ccDst, ddDst, eeDst);

                        //int srcPencilOffset = srcbd.box().size(0); //unused --dtg
                        //int dstPencilOffset = dstbd.box().size(0); //unused --dtg
                        int pencilLength = dstBox.size(0);
                        //int srcStartOffset = srcbd.box().index(srcBox.low());
                        //int dstStartOffset = dstbd.box().index(dstBox.low());
                        //int numPencils = dstBox.flatten(0).size();
                        
                        Box srcPencilBox = srcBox.flatten(0);
                        for (auto pi : srcPencilBox)
                        {
                            T* d_srcPtr = (T*)(srcbd.data(pi));
                            T* d_dstPtr = (T*)(dstbd.data(pi + a_dstShift));
                            proto_memcpy<MEM_SRC, MEM_DST>(d_srcPtr, d_dstPtr, pencilLength*sizeof(T));
                        }
                    }
                }
            } // end for dd
        } // end for ee
#ifdef PROTO_MEM_CHECK
        BoxDataMemCheck::numCopies()++;
        BoxDataMemCheck::numSlowCopies()++;
#endif
    } // end slow copy 
}

template< typename T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T, C, MEM, D, E>::copyTo(
        BoxData<T, C, MEM, D, E>&   a_dst,
        const CoordPermutation&           a_rotation) const
{
    PR_TIME("BoxData::copyTo (with rotation)");
    unsigned int L = linearSize();
    T* srcData = m_rawPtr;
    T* dstData = a_dst.m_rawPtr;
    if (MEM != HOST)
    {
        srcData = (T*)proto_malloc<HOST>(L);
        proto_memcpy<MEM, HOST>(m_rawPtr, srcData, L);
        dstData = (T*)proto_malloc<HOST>(L);
    }
    a_rotation.rotateBuffer<T,C,D,E>(srcData, dstData, m_box, a_dst.m_box); 
    if (MEM != HOST)
    {
        proto_memcpy<HOST, MEM>(dstData, a_dst.m_rawPtr, L);
        proto_free(srcData);
        proto_free(dstData);
    }
}

template< typename T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T, C, MEM, D, E>::rotate(
        Box                 a_box,
        const CoordPermutation&   a_rotation)
{
    PR_TIME("BoxData::rotate");

    unsigned int L = linearSize();
    T* oldData = (T*)proto_malloc<HOST>(L);
    proto_memcpy<MEM, HOST>(m_rawPtr, oldData, L);
    T* newData;
    if (MEM == HOST)
    {
        newData = m_rawPtr;
    } else {
        newData = (T*)proto_malloc<HOST>(L);
    }
    a_rotation.rotateBuffer<T,C,D,E>(oldData, newData, m_box, a_box); 
    if (MEM != HOST)
    {
        proto_memcpy<HOST, MEM>(newData, m_rawPtr, L);
        proto_free(newData);
    }
    proto_free<HOST>(oldData);
    m_box = a_box;
}

//=============================================================================================
// Accessors ||
//============++

// Compute Index
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
ACCEL_DECORATION
size_t BoxData<T,C,MEM,D,E>::index(
        const Point a_pt,
        unsigned int a_c,
        unsigned int a_d,
        unsigned int a_e) const
{
    PROTO_ASSERT(m_box.containsPoint(a_pt),
            "BoxData::data | Error: Point argument must be in this->box().");   
    PROTO_ASSERT((a_c < C),
            "BoxData::data | Error: Invalid component a_c = %i. a_c must be in [0,%i)",a_c,C);
    PROTO_ASSERT((a_d < D),
            "BoxData::data | Error: Invalid component a_d = %i. a_d must be in [0,%i)",a_d,D);
    PROTO_ASSERT((a_e < E),
            "BoxData::data | Error: Invalid component a_e = %i. a_e must be in [0,%i)",a_e,E);
    size_t m = m_box.size();
    size_t k = m_box.index(a_pt);
    return k+m*a_c+a_d*m*C+a_e*m*C*D;
}

// Point Accessor
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
ACCEL_DECORATION
T& BoxData<T,C,MEM,D,E>::operator()(
        const Point&    a_pt,
        unsigned int    a_c,
        unsigned int   a_d, 
        unsigned int   a_e) const
{
    return m_rawPtr[index(a_pt,a_c,a_d,a_e)];
}

// Array Accessor
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
ACCEL_DECORATION
Array<T, C*D*E> BoxData<T,C,MEM,D,E>::array(const Point& a_pt) const
{
    Array<T, C*D*E> a;
    for (int ee = 0; ee < E; ee++)
    for (int dd = 0; dd < D; dd++)
    for (int cc = 0; cc < C; cc++)
    {
        int index = C*D*ee + C*dd + cc;
        a[index] = (*this)(a_pt, cc, dd, ee);
    }
    return a;
}

// Index Accessor (Const)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
const T* BoxData<T,C,MEM,D,E>::operator[](unsigned int a_index) const
{
    PROTO_ASSERT(a_index < size(),
        "BoxData::operator[] | Error: invalid a_index = %i. \
        a_index must be in [0,this->size() = %i)",a_index, (int)size());
    return &m_rawPtr[a_index];
}


// Index Accessor (Non-Const)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
T* BoxData<T,C,MEM,D,E>::operator[](unsigned int a_index)
{
    PROTO_ASSERT(a_index < size(),
        "BoxData::operator[] | Error: invalid a_index = %i. \
        a_index must be in [0,this->size() = %i)",a_index, (int)size());
    return &m_rawPtr[a_index];
}

// Pointer Accessor (Const)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
const T* BoxData<T,C,MEM,D,E>::data(
        const Point&    a_pt,
        unsigned int    a_c,
        unsigned int   a_d, 
        unsigned int   a_e) const
{
    size_t offset = index(a_pt, a_c, a_d, a_e);
    return &(m_rawPtr[offset]);
}

// Pointer Accessor (Non-Const)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
T* BoxData<T,C,MEM,D,E>::data(
        const Point&    a_pt,
        unsigned int    a_c,
        unsigned int   a_d, 
        unsigned int   a_e)
{
    size_t offset = index(a_pt, a_c, a_d, a_e);
    return &(m_rawPtr[offset]);
}

// Pointer Accessor (Const)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
const T* BoxData<T,C,MEM,D,E>::data(
        unsigned int    a_c,
        unsigned int   a_d, 
        unsigned int   a_e) const
{
    size_t offset = index(m_box.low(), a_c, a_d, a_e);
    return &(m_rawPtr[offset]);
}

// Pointer Accessor (Non-Const)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
T* BoxData<T,C,MEM,D,E>::data(
        unsigned int    a_c,
        unsigned int   a_d, 
        unsigned int   a_e)
{
    size_t offset = index(m_box.low(), a_c, a_d, a_e);
    return &(m_rawPtr[offset]);
}

// Create Pointwise Variable (Non-Const)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
Var<T,C,MEM,D,E> BoxData<T,C,MEM,D,E>::var(const Point& a_pt)
{
    PR_TIME("BoxData::var");
    PROTO_ASSERT(m_box.containsPoint(a_pt),
        "BoxData::var | Error: a_pt not in this->box()");
    Var<T,C,MEM,D,E> rtn;
    for(unsigned int ee = 0; ee < E; ++ee)
        for(unsigned int dd = 0; dd < D; ++dd)
            for(unsigned int cc = 0; cc < C; ++cc)
            {
              //PR_TIME("set pointers");
                rtn.m_ptrs[cc + C*dd + C*D*ee] = m_rawPtr + index(a_pt,cc,dd,ee);
            }
    rtn.boxDimX=this->box().size(0);
    rtn.boxDimY=this->box().size(1);

    return rtn;
}

// Create Pointwise Variable (Const)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
Var<T,C,MEM,D,E> BoxData<T,C,MEM,D,E>::var(const Point& a_pt) const
{
    PROTO_ASSERT(m_box.containsPoint(a_pt),
        "BoxData::var | Error: a_pt not in this->box()");
    Var<T,C,MEM,D,E> rtn;
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

// Kernel Structure
//  Used by the operatorT functions to execute algebraic ops 
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
#ifdef PROTO_ACCEL
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
}; // end struct opKernel

// Algebraic Op Helper (BoxData)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
template <BoxDataOp op>
void BoxData<T,C,MEM,D,E>::operatorT(const BoxData<T,C,MEM,D,E>& a_rhs)
{
    PR_TIME("BoxData::operator(BoxData&)");

    if(m_box == a_rhs.m_box)
    {
        unsigned int nsize = size();
        int blockSize = 512;
        int nBlocks = (nsize + blockSize - 1) / blockSize;
        protoLaunchKernelMemAsyncT<MEM, opKernel<T,op>>( 
            nBlocks, blockSize, 0, protoGetCurrentStream, m_rawPtr, a_rhs.m_rawPtr, nsize
        );
        PR_FLOPS(size());
    }
    else
    {
        Box domain = m_box & a_rhs.m_box;
        if (domain.size() <= 0){return;}

        for (unsigned int ee = 0; ee < E; ee++)
            for (unsigned int dd = 0; dd < D; dd++)
            {
                BoxData<T,1,MEM> srcbd = slice(a_rhs, 0, dd, ee);
                BoxData<T,1,MEM> dstbd = slice(*this, 0, dd, ee);

                unsigned int begin  = 0;
                unsigned int end    = domain.size();

                unsigned int stride = domain.size(0);
                unsigned int blocks = domain.flatten(0).size();

                T* d_srcPtr = (T*)srcbd.data();
                T* d_dstPtr = (T*)dstbd.data();

                protoLaunchKernelT<MEM, boxdataIndexer<op, T,C>>(
                    blocks, stride, begin, end, d_srcPtr, d_dstPtr,
                    domain, srcbd.box(), dstbd.box(),
                    Point::Zeros(), srcbd.box().size(), dstbd.box().size()
                );
            }

        PR_FLOPS(domain.size()*E*D);
    }
}

// Operator + (BoxData&)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
BoxData<T,C,MEM,D,E>::operator+(const BoxData<T,C,MEM,D,E>& a_rhs) const
{
    PR_TIME("BoxData::operator+");
    BoxData<T,C,MEM,D,E> ret(m_box & a_rhs.box());
    this->copyTo(ret);
    ret.template operatorT<BoxDataOp::Add>(a_rhs);
    return ret;
}

// Operator - (BoxData&)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
BoxData<T,C,MEM,D,E>::operator-(const BoxData<T,C,MEM,D,E>& a_rhs) const
{
    PR_TIME("BoxData::operator-");
    BoxData<T,C,MEM,D,E> ret(m_box & a_rhs.box());
    this->copyTo(ret);
    ret.template operatorT<BoxDataOp::Subtract>(a_rhs);
    return ret;
}

// Operator * (BoxData&)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
BoxData<T,C,MEM,D,E>::operator*(const BoxData<T,C,MEM,D,E>& a_rhs) const
{
    PR_TIME("BoxData::operator*");
    BoxData<T,C,MEM,D,E> ret(m_box & a_rhs.box());
    this->copyTo(ret);
    ret.template operatorT<BoxDataOp::Multiply>(a_rhs);
    return ret;
}

// Operator / (BoxData&)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
BoxData<T,C,MEM,D,E>::operator/(const BoxData<T,C,MEM,D,E>& a_rhs) const
{
    PR_TIME("BoxData::operator/");
    BoxData<T,C,MEM,D,E> ret(m_box & a_rhs.box());
    this->copyTo(ret);
    ret.template operatorT<BoxDataOp::Divide>(a_rhs);
    PROTO_ASSERT(!containsInfOrNAN(),
        "A division generated values of inf or nan");
    return ret;
}

// Operator += (BoxData&)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>&
BoxData<T,C,MEM,D,E>::operator+=(const BoxData<T,C,MEM,D,E>& a_rhs)
{
    PR_TIME("BoxData::operator+=");
    this->operatorT<BoxDataOp::Add>(a_rhs);
    return *this;
}

// Operator -= (BoxData&)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>&
BoxData<T,C,MEM,D,E>::operator-=(const BoxData<T,C,MEM,D,E>& a_rhs)
{
    PR_TIME("BoxData::operator-=");
    this->operatorT<BoxDataOp::Subtract>(a_rhs);
    return *this;
}

// Operator *= (BoxData&)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>&
BoxData<T,C,MEM,D,E>::operator*=(const BoxData<T,C,MEM,D,E>& a_rhs)
{
    PR_TIME("BoxData::operator*=");
    // PC: Replacing with forallInPlace
    //this->operatorT<BoxDataOp::Multiply>(a_rhs);
    forallInPlace([] PROTO_LAMBDA (Var<T,C,MEM,D,E>& a_out,Var<T,C,MEM,D,E>& a_in)
                  {
                    for (int ee = 0; ee < E; ee++)
                      for (int dd = 0; dd < D; dd++)
                        for (int cc = 0; cc < C; cc++)
                          a_out(cc,dd,ee) *= a_in(cc,dd,ee);
                  },*this,a_rhs);
    return *this;
}

// Operator /= (BoxData&)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>&
BoxData<T,C,MEM,D,E>::operator/=(const BoxData<T,C,MEM,D,E>& a_rhs)
{
    PR_TIME("BoxData::operator/=");
    this->operatorT<BoxDataOp::Divide>(a_rhs);
    PROTO_ASSERT(!containsInfOrNAN(),
    "A division generated values of inf or nan");
    return *this;
}

// Algebraic Op Helper (Scalar)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
template <BoxDataOp op>
void BoxData<T,C,MEM,D,E>::operatorT(T a_scale)
{
    PR_TIME("BoxData::operator");
    PR_FLOPS(size());

    for (unsigned int ee = 0; ee < E; ee++)
        for (unsigned int dd = 0; dd < D; dd++)
        {
            BoxData<T,1,MEM> dstbd = slice(*this, 0, dd, ee);

            unsigned int begin  = 0;
            unsigned int end    = m_box.size();

            unsigned int stride = m_box.size(0);
            unsigned int blocks = m_box.flatten(0).size();

            T* d_dstPtr = (T*)dstbd.data();
            size_t smem=0;
            protoLaunchKernelMemAsyncT<MEM, scalarIndexer<op, T,C>>( 
                blocks, stride, smem, protoGetCurrentStream, begin, end, d_dstPtr, a_scale,
                dstbd.box(), dstbd.box().size());
        }
}

template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>& BoxData<T,C,MEM,D,E>::operator+=(T a_scale)
{
    PR_TIME("BoxData::operator+=(scale)");
    this->operatorT<BoxDataOp::Add>(a_scale);
    return *this;
}

template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>& BoxData<T,C,MEM,D,E>::operator-=(T a_scale)
{
    PR_TIME("BoxData::operator-=(scale)");
    this->operatorT<BoxDataOp::Subtract>(a_scale);
    return *this;
}

template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>& BoxData<T,C,MEM,D,E>::operator*=(T a_scale)
{
    PR_TIME("BoxData::operator*=(scale)");
    // PC: Replacing with forallInPlace
    //this->operatorT<BoxDataOp::Multiply>(a_scale);
    forallInPlace([] PROTO_LAMBDA (Var<T,C,MEM,D,E>& a_out,T& a_in)
                  {
                    for (int ee = 0; ee < E; ee++)
                      for (int dd = 0; dd < D; dd++)
                        for (int cc = 0; cc < C; cc++)
                          a_out(cc,dd,ee) *= a_in;
                  },*this,a_scale);
    return *this;
}

template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>& BoxData<T,C,MEM,D,E>::operator/=(T a_scale)
{
    PR_TIME("BoxData::operator/=(scale)");
    this->operatorT<BoxDataOp::Divide>(a_scale);
    PROTO_ASSERT(!containsInfOrNAN(),
    "A division generated values of inf or nan");
    return *this;
}

template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>& BoxData<T,C,MEM,D,E>::operator+=(const Array<T, C*D*E>& a_scale)
{
    PR_TIME("BoxData::operator+=(scale)");
    for (int ee = 0; ee < E; ee++)
    for (int dd = 0; dd < D; dd++)
    for (int cc = 0; cc < C; cc++)
    {
        T scale = a_scale[ee*C*D + dd*C + cc];
        auto data = slice(*this,cc,dd,ee);
        data.template operatorT<BoxDataOp::Add>(scale);
    }
    return *this;
}

template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>& BoxData<T,C,MEM,D,E>::operator-=(const Array<T, C*D*E>& a_scale)
{
    PR_TIME("BoxData::operator+=(scale)");
    for (int ee = 0; ee < E; ee++)
    for (int dd = 0; dd < D; dd++)
    for (int cc = 0; cc < C; cc++)
    {
        T scale = a_scale[ee*C*D + dd*C + cc];
        auto data = slice(*this,cc,dd,ee);
        data.template operatorT<BoxDataOp::Subtract>(scale);
    }
    return *this;
}

template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>& BoxData<T,C,MEM,D,E>::operator*=(const Array<T, C*D*E>& a_scale)
{
    PR_TIME("BoxData::operator+=(scale)");
    for (int ee = 0; ee < E; ee++)
    for (int dd = 0; dd < D; dd++)
    for (int cc = 0; cc < C; cc++)
    {
        T scale = a_scale[ee*C*D + dd*C + cc];
        auto data = slice(*this,cc,dd,ee);
        data.template operatorT<BoxDataOp::Multiply>(scale);
    }
    return *this;
}

template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>& BoxData<T,C,MEM,D,E>::operator/=(const Array<T, C*D*E>& a_scale)
{
    PR_TIME("BoxData::operator+=(scale)");
    for (int ee = 0; ee < E; ee++)
    for (int dd = 0; dd < D; dd++)
    for (int cc = 0; cc < C; cc++)
    {
        T scale = a_scale[ee*C*D + dd*C + cc];
        PROTO_ASSERT(scale != 0, "BoxData::operator/= Error: Division by 0");
        auto data = slice(*this,cc,dd,ee);
        data.template operatorT<BoxDataOp::Divide>(scale);
    }
    PROTO_ASSERT(!containsInfOrNAN(),
    "A division generated values of inf or nan");
    return *this;
}


//=============================================================================================
// Utility Functions ||
//===================++

// Global Set Value
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::setVal(const T& a_val)
{
    if (box().empty()) { return; }
    PR_TIME("BoxData::setVal");
    forallInPlace([] PROTO_LAMBDA (Var<T,C,MEM,D,E>& a_out,T& a_in)
                  {
                    for (int ee = 0; ee < E; ee++)
                      for (int dd = 0; dd < D; dd++)
                        for (int cc = 0; cc < C; cc++)
                          a_out(cc,dd,ee) = a_in;
                  },*this,a_val);
}
// Set To Zero
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::setToZero()
{
    setVal(0);
}

template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::setRandom(T a_low, T a_high)
{
    std::random_device rd;  // A random seed generator
    std::mt19937 gen(rd()); // Mersenne Twister PRNG engine
    std::uniform_real_distribution<T> dis(a_low, a_high); // Uniform dist between 0 and 1
  
    if (MEM == HOST)
    {
        for (int ii = 0; ii < size(); ii++)
        {
            m_rawPtr[ii] = dis(gen);
        }
    } else {
        T* buffer = (T*)proto_malloc<HOST>(linearSize());
        for (int ii = 0; ii < size(); ii++)
        {
            buffer[ii] = dis(gen);
        }
        proto_memcpy<HOST,MEM>(buffer, m_rawPtr, linearSize());
        proto_free<HOST>(buffer);
    }
}

    // Helper function for BoxData::setVal with a Box argument
template <class T, MemType MEM, class B>
inline void setBoxDataVal(BoxData<T,1,MEM,1,1>& a_data, const B& a_box, T a_val)
{
    if (a_box.empty()) { return; }
#ifdef PROTO_ACCEL
    MayDay<void>::Abort(
        "::setBoxDataVal | Error: This function should never be called when PROTO_ACCEL=TRUE");
#else
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

// Set Values Within Box
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::setVal(const T& a_val, const Box&  a_box)
{
    Box intersect = m_box & a_box;
    if (intersect.empty()) { return; }
    for (unsigned int ee = 0; ee < E; ee++)
    {
        for (unsigned int dd = 0; dd < D; dd++)
        {
            for (unsigned int cc = 0; cc < C; cc++)
            {
                auto slicebd = slice(*this, cc, dd, ee);
                setBoxDataVal<T,MEM>(slicebd, intersect, a_val);
            }
        }
    }
}

// Set Component Values Within Box
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::setVal(
        const T&    a_val,
        const Box&  a_box,
        int         a_c,
        int         a_d,
        int         a_e)
{

    PR_TIME("BoxData::setVal(T&, Box, int, int, int)");
    PROTO_ASSERT((a_c < C),
            "BoxData::data | Error: Invalid component a_c = %i. a_c must be in [0,%i)",a_c,C);
    PROTO_ASSERT((a_d < D),
            "BoxData::data | Error: Invalid component a_d = %i. a_d must be in [0,%i)",a_d,D);
    PROTO_ASSERT((a_e < E),
            "BoxData::data | Error: Invalid component a_e = %i. a_e must be in [0,%i)",a_e,E);

    auto slicebd = slice(*this, a_c, a_d, a_e);
    setBoxDataVal<T,MEM>(slicebd, a_box, a_val);
}

/// Compute Reduction (Global)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
template<Proto::Operation OP>
void BoxData<T,C,MEM,D,E>::reduce(Reduction<T,OP,MEM>& a_rxn) const
{
    PR_TIME("BoxData::reduce");
    a_rxn.reduce(m_rawPtr,size());
}

/// Compute Reduction (Componentwise)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
template<Proto::Operation OP>
void BoxData<T,C,MEM,D,E>::reduce(Reduction<T,OP,MEM>& a_rxn, int a_c, int a_d, int a_e) const
{
    auto temp = slice(*this, a_c, a_d, a_e);
    temp.reduce(a_rxn);
}

/// Maximum Absolute Value (Global)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::absMax(Reduction<T,Abs,MEM>& a_rxn) const
{
    reduce(a_rxn);
}

/// Maximum Absolute Value (Componentwise)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::absMax(Reduction<T,Abs,MEM>& a_rxn, int a_c, int a_d, int a_e) const
{
    reduce(a_rxn, a_c, a_d, a_e);
}

/// Maximum Absolute Value (Global)
template <class T, unsigned int C,  MemType MEM,unsigned int D, unsigned int E>
T BoxData<T,C,MEM,D,E>::absMax() const
{
    Reduction<T,Abs,MEM> rxn(true); 
    reduce(rxn);
    return rxn.fetchLocal();
}

/// Maximum Absolute Value (Componentwise)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
T BoxData<T,C,MEM,D,E>::absMax(int a_c, int a_d, int a_e) const
{
    Reduction<T,Abs,MEM> rxn(true); 
    reduce(rxn, a_c, a_d, a_e);
    return rxn.fetchLocal();
}

/// Minimum Value (Global)
template <class T, unsigned int C,  MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::min(Reduction<T,Min,MEM>& a_rxn) const
{
    reduce(a_rxn);
}

/// Minimum Value (Componentwise)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::min(Reduction<T, Min, MEM>& a_rxn,
    int a_c, int a_d, int a_e) const
{
    reduce(a_rxn, a_c, a_d, a_e);
}

/// Minimum Value (Global)
template <class T, unsigned int C,  MemType MEM, unsigned int D, unsigned int E>
T BoxData<T,C,MEM,D,E>::min() const
{
    Reduction<T,Min,MEM> rxn(true);
    reduce(rxn);
    return rxn.fetchLocal();
}

/// Minimum Value (Componentwise)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
T BoxData<T,C,MEM,D,E>::min(int a_c, int a_d, int a_e) const
{
    Reduction<T,Min,MEM> rxn(true);
    reduce(rxn, a_c, a_d, a_e);
    return rxn.fetchLocal();
}

/// Maximum Value (Global)
template <class T, unsigned int C,MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::max(Reduction<T,Max,MEM>& a_rxn) const
{
    reduce(a_rxn);
}

/// Maximum Value (Componentwise)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::max(Reduction<T, Max, MEM>& a_rxn, int a_c, int a_d, int a_e) const
{
    reduce(a_rxn, a_c, a_d, a_e);
}

/// Maximum Value (Global)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
T BoxData<T,C,MEM,D,E>::max() const
{
    Reduction<T,Max,MEM> rxn(true);
    reduce(rxn);
    return rxn.fetchLocal();
}

/// Maximum Value (Componentwise)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
T BoxData<T,C,MEM,D,E>::max(int a_c, int a_d, int a_e) const
{
    Reduction<T,Max,MEM> rxn(true);
    reduce(rxn, a_c, a_d, a_e);
    return rxn.fetchLocal();
}

// Sum of Values (Global)
template <class T, unsigned int C, MemType MEM,unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::sum(Reduction<T, Sum, MEM>& a_rxn) const
{
    reduce(a_rxn);
}

// Sum of Values (Componentwise)
template <class T, unsigned int C, MemType MEM,unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::sum( Reduction<T, Sum, MEM>& a_rxn, int a_c, int a_d, int a_e) const
{
    reduce(a_rxn, a_c, a_d, a_e);
}

// Sum of Values (Global)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
T BoxData<T,C,MEM,D,E>::sum() const
{
    Reduction<T,Sum,MEM> rxn(true);
    reduce(rxn);
    return rxn.fetchLocal();
}

// Sum of Values (Componentwise)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
T BoxData<T,C,MEM,D,E>::sum(int a_c, int a_d, int a_e) const
{
    Reduction<T,Sum,MEM> rxn(true);
    reduce(rxn, a_c, a_d, a_e);
    return rxn.fetchLocal();
}

// Sum of Abs Values (Global)
template <class T, unsigned int C, MemType MEM,unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::sumAbs(Reduction<T, SumAbs, MEM>& a_rxn) const
{
    reduce(a_rxn);
}

// Sum of Abs Values (Componentwise)
template <class T, unsigned int C, MemType MEM,unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::sumAbs( Reduction<T, SumAbs, MEM>& a_rxn, int a_c, int a_d, int a_e) const
{
    reduce(a_rxn, a_c, a_d, a_e);
}

// Sum of Abs Values (Global)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
T BoxData<T,C,MEM,D,E>::sumAbs() const
{
    Reduction<T,SumAbs,MEM> rxn(true);
    reduce(rxn);
    return rxn.fetchLocal();
}

// Sum of Values (Componentwise)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
T BoxData<T,C,MEM,D,E>::sumAbs(int a_c, int a_d, int a_e) const
{
    Reduction<T,SumAbs,MEM> rxn(true);
    reduce(rxn, a_c, a_d, a_e);
    return rxn.fetchLocal();
}

// Sum of Abs Values (Global)
template <class T, unsigned int C, MemType MEM,unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::sumSquare(Reduction<T, SumSquare, MEM>& a_rxn) const
{
    reduce(a_rxn);
}

// Sum of Abs Values (Componentwise)
template <class T, unsigned int C, MemType MEM,unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::sumSquare( Reduction<T, SumSquare, MEM>& a_rxn, int a_c, int a_d, int a_e) const
{
    reduce(a_rxn, a_c, a_d, a_e);
}

// Sum of Abs Values (Global)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
T BoxData<T,C,MEM,D,E>::sumSquare() const
{
    Reduction<T,SumSquare,MEM> rxn(true);
    reduce(rxn);
    return rxn.fetchLocal();
}

// Sum of Values (Componentwise)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
T BoxData<T,C,MEM,D,E>::sumSquare(int a_c, int a_d, int a_e) const
{
    Reduction<T,SumSquare,MEM> rxn(true);
    reduce(rxn, a_c, a_d, a_e);
    return rxn.fetchLocal();
}


// Integrate (Isotropic Grid)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
T BoxData<T,C,MEM,D,E>::integrate(T a_dx, int a_c, int a_d, int a_e) const
{
    return sum(a_c, a_d, a_e) * pow(a_dx, DIM);
}

// Integrate (Anisotropic Grid)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
T BoxData<T,C,MEM,D,E>::integrate(Array<T, DIM> a_dx, int a_c, int a_d, int a_e) const
{
    T factor = 1.0;
    for (int dir = 0; dir < DIM; dir++) { factor *= a_dx[dir]; }
    return sum(a_c, a_d, a_e) * factor;
}

// Shift Domain
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::shift(const Point& a_shift) const
{
    m_box = m_box.shift(a_shift);
}

// Buffer Size
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
size_t BoxData<T,C,MEM,D,E>::linearSize() const
{
    return size()*sizeof(T);
}

// Buffer Write
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::linearOut(
        void*      a_buf,
        const Box& a_box,
        CInterval  a_comps) const
{

    PR_TIME("BoxData::linearOut");
    PROTO_ASSERT(m_box.containsBox(a_box),
            "BoxData::linearOut | Error: \
            a_box must be a subset of this->box().");
    PROTO_ASSERT(this->contains(a_comps),
            "BoxData::linearOut | Error: \
            a_comps must be contained in the component space of *this.");
    PROTO_ASSERT(pointerMemType(a_buf)==MEM,
            "BoxData::linearOut | Error: \
            Buffer and BoxData have different MemType.");
    T *ptr = static_cast<T*>(a_buf);
    unsigned int c_lo = a_comps.low(0), d_lo = a_comps.low(1), e_lo = a_comps.low(2);
    unsigned int c_hi = a_comps.high(0), d_hi = a_comps.high(1), e_hi = a_comps.high(2);
    for (int i=c_lo; i<=c_hi; i++)
        for (int j=d_lo; j<=d_hi; j++)
            for (int k=e_lo; k<=e_hi; k++) {
                BoxData<T,1,MEM> dest(ptr, a_box, 1);
                BoxData<T,1,MEM> src = slice(*this,i,j,k);
                src.copyTo(dest, a_box);
                ptr += a_box.size();
            }
}

// Buffer Write (Primitive)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::linearOut(
        void*               a_buf,
        const ::Proto::Box& a_bx,
        unsigned int        a_startcomp,
        unsigned int        a_numcomps) const
{
    CInterval comps(a_startcomp, a_startcomp + a_numcomps-1);
    linearOut(a_buf, a_bx, comps);
}

// Buffer Read
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::linearIn(
        void*      a_buf,
        const Box& a_box,
        CInterval  a_comps)
{

    PR_TIME("BoxData::linearIn");
    PROTO_ASSERT(m_box.containsBox(a_box),
            "BoxData::linearIn | Error: \
            a_box must be a subset of this->box().");
    PROTO_ASSERT(this->contains(a_comps),
            "BoxData::linearIn | Error: \
            a_comps must be contained in the component space of *this.");
    PROTO_ASSERT(pointerMemType(a_buf)==MEM,
            "BoxData::linearIn | Error: \
            Buffer and BoxData have different MemType.");

    T *ptr = static_cast<T*>(a_buf);
    unsigned int c_lo = a_comps.low(0), d_lo = a_comps.low(1), e_lo = a_comps.low(2);
    unsigned int c_hi = a_comps.high(0), d_hi = a_comps.high(1), e_hi = a_comps.high(2);
    for (int i=c_lo; i<=c_hi; i++)
        for (int j=d_lo; j<=d_hi; j++)
            for (int k=e_lo; k<=e_hi; k++) {
                BoxData<T,1,MEM> src(ptr, a_box);
                BoxData<T,1,MEM> dest = slice(*this,i,j,k);
                src.copyTo(dest, a_box);
                ptr += a_box.size();
            }
}

// Buffer Read (Primitive)
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::linearIn(
        void*               a_buf,
        const ::Proto::Box& a_bx,
        unsigned int        a_startcomp,
        unsigned int        a_numcomps)
{
    CInterval comps(a_startcomp, a_startcomp + a_numcomps-1);
    linearIn(a_buf, a_bx, comps);
}

// Contains Component Interval Query
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
bool BoxData<T,C,MEM,D,E>::contains(CInterval a_interval) const
{
    bool ret = true;
    ret &= (a_interval.high(0) < C);
    ret &= (a_interval.high(1) < D);
    ret &= (a_interval.high(2) < E);
    return ret;
}

// Is Alias Query
template <class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
template<unsigned int CC, unsigned int DD, unsigned int EE>
bool BoxData<T,C,MEM,D,E>::isAlias(const BoxData<T, CC, MEM, DD, EE>& a_src) const
{
    return (this->aliasData() == a_src.aliasData());
}

// Simple Print
template <typename T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::print() const
{
    using Proto::pr_out;

    pr_out() << "BoxData<" << C << ", " << (int)D << ", " << (int)E << ">" << std::endl;
    pr_out() << "Domain: " << m_box << " | Size: " << this->size() << std::endl;
    if (this->size() <= 0){return;}
    for (int ee = 0; ee < E; ee++)
        for (int dd = 0; dd < D; dd++)
            for (int cc = 0; cc < C; cc++)
            {
                pr_out() << "\tComponent (" << cc << ", " << dd << ", " << ee << "): ";
                pr_out() << " | Max: " << this->max(cc,dd,ee) << " | Min: " << this->min(cc,dd,ee);
                pr_out() << std::endl;
            }
}

// Print Data
template <typename T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::printData(int a_prec ) const
{
    printData(m_box,-1,-1,-1,a_prec);
}

// Print Data
template <typename T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::printData(const Box& a_box, int a_prec) const
{
    printData(a_box,-1,-1,-1, a_prec);
}

// Print Data
template <typename T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoxData<T,C,MEM,D,E>::printData(
        const Box& a_box,
        int a_c, int a_d, int a_e,
        int a_prec) const
{
    using std::cout;
    using std::endl;
    std::ostream& out = Proto::pr_out();
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
    if (this->absMax() > 1.0e8 || this->absMax() < pow(10.0, -a_prec))
    {
        this->print();
        return;
    }

    PROTO_ASSERT(m_box.containsBox(a_box),
        "BoxData::printData(a_box,...) | Error: \
        Box is not a subset of the domain of *this.");
    
    out << "Domain: " << m_box << " | Subdomain: " << a_box << endl;
    if (a_c != -1)
    {
        PROTO_ASSERT((a_c < C),
            "BoxData::printData | Error: \
            a_c = %i is out of bounds of [0,C=%i).", a_c,C);
    }
    if (a_d != -1)
    {
        PROTO_ASSERT((a_d < D),
            "BoxData::printData | Error: \
            a_d = %i is out of bounds of [0,D=%i).", a_d,D);
    }
    if (a_e != -1)
    {
        PROTO_ASSERT((a_e < E),
            "BoxData::printData | Error: \
            a_e = %i is out of bounds of [0,E=%i).", a_e,E);
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
                out << std::setfill('=');
                out << std::setw(std::min((int)(a_box.size(0)*(width+1)+8),100)) << " " << endl;
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
                    out << std::setfill('-');
                    out << std::setw(std::min((int)(a_box.size(0)*6+8),100)) << " " << endl;

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
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
bool BoxData<T,C,MEM,D,E>::containsInfOrNAN() const
{
    // NAN only applies to floating point data
    if (!std::is_floating_point<T>::value) { return false; }
    for (int ii = 0; ii < numValues(); ii++)
    {
        if (std::isnan(m_rawPtr[ii]) || std::isinf(m_rawPtr[ii])) { return true; }
    }
    return false;
}

//=================================================================================================
// Alias and Slice ||
//=================++

/// Alias (Non-Const)
    template<class T, unsigned int C, unsigned int D, unsigned int E, MemType MEM>
BoxData<T,C,MEM,D,E> alias(BoxData<T,C,MEM,D,E>& a_original, const Point& shift)
{
    PR_TIME("alias(BoxData<T,C,MEM,D,E>&)");
    const Box& b=a_original.box();
    BoxData<T,C,MEM,D,E> rtn(a_original.m_data, a_original.m_rawPtr, b);
    rtn.shift(shift);
    return rtn;
}

/// Alias (Const)
template<class T, unsigned int C, unsigned int D, unsigned int E, MemType MEM>
const BoxData<T,C,MEM,D,E> alias(const BoxData<T,C,MEM,D,E>& a_original, const Point& shift)
{
    PR_TIME("alias(const BoxData<T,C,MEM,D,E>&)");
    const Box& b=a_original.box();
    BoxData<T,C,MEM,D,E>* src = const_cast<BoxData<T,C,MEM,D,E>*>(&a_original);
    const BoxData<T,C,MEM,D,E> rtn(src->m_data, src->m_rawPtr, b.shift(shift));
    return rtn;
}

// Scalar Slice
template<typename T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,1,MEM> slice(
        const BoxData<T,C,MEM,D,E>& a_src,
        unsigned int a_c,
        unsigned int a_d,
        unsigned int a_e)
{
    PR_TIME("slice(BoxData<T,C,MEM,D,E>&, int, int, int)");
    PROTO_ASSERT((a_c < C),
        "slice(BoxData&, uint, uint, uint) | Error: \
        Component index a_c = %i should be positive and less than C = %i.", a_c, C);
    PROTO_ASSERT((a_d < D),
        "slice(BoxData&, uint, uint, uint) | Error: \
        Component index a_d = %i should be positive and less than D = %i.", a_d, D);
    PROTO_ASSERT((a_e < E),
        "slice(BoxData&, uint, uint, uint) | Error: \
        Component index a_e = %i should be positive and less than E = %i.", a_e, E);
    const Box& b = a_src.box();
    T* ptr = (T*)(a_src.data(b.low(),a_c,a_d,a_e));
    BoxData<T,1,MEM> rtn(a_src.m_data, ptr, b);
    return rtn;
}

template<typename T, unsigned int C, unsigned int CC, MemType MEM>
BoxData<T,CC,MEM> slice(
        const BoxData<T,C,MEM>& a_src,
        unsigned int a_nstart)
{
    PR_TIME("slice(BoxData<T,C,1,1>&, int)");
    PROTO_ASSERT((a_nstart + CC <= C),
        "slice(const BoxData&, uint) | Error: Invalid slicing range")
    const Box& b = a_src.box();
    T* ptr = (T*)(a_src.data(b.low(),a_nstart));
    BoxData<T,CC,MEM> rtn(a_src.m_data, ptr, b);
    return rtn;
}

template<typename T, unsigned int C, unsigned int D, MemType MEM>
BoxData<T,C,MEM,1> plane(
        const BoxData<T,C,MEM,D>& a_src,
        unsigned int a_d)
{
    PR_TIME("slice(BoxData<T,C,1,1>&, int)");
    PROTO_ASSERT((a_d < D),
        "slice(const BoxData&, uint) | Error: Invalid slicing range")
    const Box& b = a_src.box();
    T* ptr = (T*)(a_src.data(b.low(),0,a_d));
    BoxData<T,C,MEM,1> rtn(a_src.m_data, ptr, b);
    return rtn;
}
