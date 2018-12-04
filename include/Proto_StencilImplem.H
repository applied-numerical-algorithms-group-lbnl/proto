//=======================================================================================
// Constructors ||
//==============++ 
#include "Proto_Timer.H"
// Default Constructor
template <typename T>
Stencil<T>::Stencil()
{
    m_destRefratio=Point::Ones();
    m_destShift=Point::Zeros();
    m_srcRefratio=Point::Ones();
    m_isClosed = false;
}

// General Constructor
template <typename T>
Stencil<T>::Stencil(Shift a_shift,
                    T     a_coef,
                    Point a_destRefratio,
                    Point a_destShift,
                    Point a_srcRefratio)
{
    m_isClosed = false;
    m_destRefratio = a_destRefratio;
    m_srcRefratio = a_srcRefratio;
    m_destShift = a_destShift;
    augment(a_coef, a_shift.shift());
}

//=======================================================================================
// Operators ||
//===========++ 

// Stencil Composition
template <typename T>
Stencil<T> Stencil<T>::operator* (const Stencil<T>& a_stencil) const
{
    Stencil<T> rtn;
    rtn.m_destRefratio = m_destRefratio;
    rtn.m_srcRefratio  = m_srcRefratio;
    rtn.m_destShift    = m_destShift;
    for(int i=0; i<m_coefs.size(); ++i)
    for(int j=0; j<a_stencil.m_coefs.size(); ++j)
    {
        rtn+=Stencil<T>(Shift(m_offsets[i]+a_stencil.m_offsets[j]),
                        m_coefs[i]*a_stencil.m_coefs[j]);
                            
    }
    return rtn;
}

// Scalar Multiplication
template <typename T>
Stencil<T> Stencil<T>::operator*(const T a_coef) const
{
    Stencil<T> rtn = *this;
    for (int ii = 0; ii < this->size(); ii++)
    {
        rtn.m_coefs[ii] *= a_coef;
    }
    return rtn;
}

// In Place Stencil Composition
template <typename T>
void Stencil<T>::operator*=(const Stencil<T>& a_stencil)
{
    //std::move to avoid copying the product
    (*this) = std::move((*this)*a_stencil);
}   

// In Place Scalar Multiplication
template <typename T>
void Stencil<T>::operator*=(const T a_coef)
{
    for (int l = 0; l < m_coefs.size(); l++)
    {
        m_coefs[l]*=a_coef;
    }
}

// Stencil Addition
template <class T>
Stencil<T> Stencil<T>::operator+(const Stencil<T>& a_stencil) const
{
    PROTO_ASSERT(m_srcRefratio == a_stencil.m_srcRefratio,
    "Stencil::operator+(Stencil a_stencil) invalid.\
    Cannot add stencils with differing srcRefratios.");
    PROTO_ASSERT(m_destRefratio == a_stencil.m_destRefratio,
    "Stencil::operator+(Stencil a_stencil) invalid.\
    Cannot add stencils with differing destRefratios.");
    PROTO_ASSERT(m_destShift == a_stencil.m_destShift,
    "Stencil::operator+(Stencil a_stencil) invalid.\
    Cannot add stencils with differing destShifts.");
    Stencil<T> rtn = (*this);
    rtn += a_stencil;
    return rtn; // should make use of return value optimization
}

// Stencil Subtraction
template <class T>
Stencil<T> Stencil<T>::operator-(const Stencil<T>& a_stencil) const
{
    PROTO_ASSERT(m_srcRefratio == a_stencil.m_srcRefratio,
    "Stencil::operator+(Stencil a_stencil) invalid.\
    Cannot add stencils with differing srcRefratios.");
    PROTO_ASSERT(m_destRefratio == a_stencil.m_destRefratio,
    "Stencil::operator+(Stencil a_stencil) invalid.\
    Cannot add stencils with differing destRefratios.");
    PROTO_ASSERT(m_destShift == a_stencil.m_destShift,
    "Stencil::operator+(Stencil a_stencil) invalid.\
    Cannot add stencils with differing destShifts.");
    Stencil<T> rtn = (*this);
    rtn -= a_stencil;
    return rtn; // should make use of return value optimization
}

// In Place Addition
template <class T>
void Stencil<T>::operator+=(const Stencil<T>& a_stencil)
{
    PROTO_ASSERT(m_srcRefratio == a_stencil.m_srcRefratio,
    "Stencil::operator+=(Stencil a_stencil) invalid.\
    Cannot add stencils with differing srcRefratios.");
    PROTO_ASSERT(m_destRefratio == a_stencil.m_destRefratio,
    "Stencil::operator+=(Stencil a_stencil) invalid.\
    Cannot add stencils with differing destRefratios.");
    PROTO_ASSERT(m_destShift == a_stencil.m_destShift,
    "Stencil::operator+=(Stencil a_stencil) invalid.\
    Cannot add stencils with differing destShifts.");
    for (int ii = 0 ; ii < a_stencil.m_coefs.size();ii++)
    {
        const Point& offset = a_stencil.m_offsets[ii];
        T coef = a_stencil.m_coefs[ii];
        augment(coef, offset);
    }
}

// In Place Subtraction
template <class T>
void Stencil<T>::operator-=(const Stencil<T>& a_stencil)
{
    PROTO_ASSERT(m_srcRefratio == a_stencil.m_srcRefratio,
    "Stencil::operator+=(Stencil a_stencil) invalid.\
    Cannot add stencils with differing srcRefratios.");
    PROTO_ASSERT(m_destRefratio == a_stencil.m_destRefratio,
    "Stencil::operator+=(Stencil a_stencil) invalid.\
    Cannot add stencils with differing destRefratios.");
    PROTO_ASSERT(m_destShift == a_stencil.m_destShift,
    "Stencil::operator+=(Stencil a_stencil) invalid.\
    Cannot add stencils with differing destShifts.");
    for (int ii = 0 ; ii < a_stencil.m_coefs.size();ii++)
    {
        const Point& offset = a_stencil.m_offsets[ii];
        T coef = -a_stencil.m_coefs[ii];
        augment(coef, offset);
    }
}

// Equality Operator
template <typename T>
bool Stencil<T>::operator==(Stencil<T>& a_stencil) const
{
    if ((size() != a_stencil.size()) ||
         (m_srcRefratio != a_stencil.m_srcRefratio) ||
         (m_destRefratio != a_stencil.m_destRefratio) ||
         (m_destShift != a_stencil.m_destShift))
    {
        return false;
    }
    for (int ii = 0; ii < m_offsets.size(); ii++)
    {
        bool found = false;
        for (int jj = 0; jj < a_stencil.m_offsets.size(); jj++)
        {
            if (m_offsets[ii] == a_stencil.m_offsets[jj])
            {
                found = true;
                if (m_coefs[ii] != a_stencil.m_coefs[jj]){return false;}
            }
        }
        if (!found){return false;}
    }
    return true;
}

// Augment Stencil
template <typename T>
void Stencil<T>::augment(T a_coef, Point a_offset)
{
    bool isThere = false;
    int jj = 0;
    for (;jj < m_coefs.size();jj++)
    {
        if (m_offsets[jj] == a_offset)
        {
            isThere = true;
            break;
        }
    }
    if (isThere) 
    {
        m_coefs[jj] += a_coef;
    }
    else
    {
        m_coefs.push_back(a_coef);
        m_offsets.push_back(a_offset);
        if (!m_span.contains(a_offset))
        {
            Point low = m_span.low();
            Point high = m_span.high();
            for (int d = 0; d < DIM; d++)
            {
                low[d] = (a_offset[d] < low[d]) ? a_offset[d] : low[d];
                high[d] = (a_offset[d] > high[d]) ? a_offset[d] : high[d];
            }
            m_span = Box(low,high);
        }
    }
}

//=======================================================================================
// Apply ||
//=======++ 

// Apply Stencil
template <typename T>
template <unsigned int C, unsigned char D, unsigned char E>
void Stencil<T>::apply(const BoxData<T,C,D,E>&  a_src,
                       BoxData<T,C,D,E>&        a_dest,
                       const Box&                a_box,
                       bool                     a_initToZero,
                       T                        a_scale) const                      
{
  PR_TIME("Stencil::apply");
#ifdef PROTO_CUDA
  
    BoxData<T, C, D, E>& castsrc = const_cast<BoxData<T,C,D,E> &>(a_src);
    Box & castbox = const_cast<Box &>(a_box);
    Stencil<T>* castthis = const_cast<Stencil<T>* >(this);
    castthis->cudaApply(castsrc, a_dest, castbox, a_initToZero, a_scale);

#else


/**/
    if (m_coefs.size() == 0){return;}
    // Check input domain Boxes for proper nesting
    PROTO_ASSERT(!a_src.isAlias(a_dest),
    "Stencil::apply(...) invalid. a_src and a_dest alias the same data.");
   
    //CH_TIMERS("Stencil::apply");
    // Initialize a_dest to zero within a_box if flag is set
    for (int ii = 0; ii < DIM; ii++)
    {
        PROTO_ASSERT((m_destShift[ii] >= 0) && (m_destShift[ii] < m_destRefratio[ii]), 
        "m_destShift[ii] must be non-negative and less than m_destRefratio[ii] for all ii in [0,DIM)");
    }
    Box cross = a_box.flatten(0);
    int npencil = a_box.size(0);
    if (a_initToZero)
    {
        for (int ee = 0; ee < E; ee++)
        for (int dd = 0; dd < D; dd++)
        for (int cc = 0; cc < C; cc++)
        for (auto iter = cross.begin(); iter != cross.end(); ++iter)
        {
            Point pt = (*iter)*m_destRefratio + m_destShift;
            T* val = a_dest.data(pt,cc,dd,ee);
            for (int ii = 0; ii < npencil; ii++)
            {
                val[ii*m_destRefratio[0]] = 0;
            }
        }
    }
    
    // if a_scale is 0, nothing left to do...
    if (a_scale == 0){return;}
    
    // ...otherwise make a copy of m_coefs and multiply by a_scale
    std::vector<T> coefs = m_coefs;
    if (a_scale != 1)
    {
        for (int ii = 0; ii < this->size(); ii++)
        {
            coefs[ii] *= a_scale;
        }
    }

    // convert Points in m_offsets to linear offsets in a_src's data
    int offsetFactor[DIM];
    offsetFactor[0] = 1;
    for (int ii = 1; ii < DIM; ii++)
    {
        offsetFactor[ii] = offsetFactor[ii-1]*a_src.box().size(ii-1);
    }
    std::vector<int> offsets;
    offsets.resize(this->size());
    for (int ii = 0; ii < this->size(); ii++)
    {
        offsets[ii] = 0;
        for (int dir = 0; dir < DIM; dir++)
        {
            offsets[ii] += m_offsets[ii][dir]*offsetFactor[dir];
        }
    }

    // apply the stencil
    for (int ee = 0; ee < E; ee++)
    for (int dd = 0; dd < D; dd++)
    for (int cc = 0; cc < C; cc++)
    {    
        for (auto iter = cross.begin(); iter != cross.end(); ++iter)
        {
            for (int jj = 0; jj < this->size(); jj++)
            {
                Point dpt = (*iter)*m_destRefratio + m_destShift;
                if (!a_dest.box().contains(dpt)){continue;} //can happen when destShift is non-trivial
                Point spt = (*iter)*m_srcRefratio;
                T* destVal = a_dest.data(dpt,cc,dd,ee);
                const T* srcVal = a_src.data(spt,cc,dd,ee) + offsets[jj];
                for (int ii = 0; ii < npencil; ii++)
                {
                    destVal[(ii*m_destRefratio[0])] += coefs[jj]*srcVal[(ii*m_srcRefratio[0])];
                }
            }
        }
    }
/**/
#endif //end non-cuda version
    long long int numpts = a_box.size();
    long long int numflops = 3*size(); //multiplies by scale and  by weight and add
    PR_FLOPS(numpts*numflops);
}

//=======================================================================================
// LazyStencil ||
//=============++ 

template <typename T, unsigned int C, unsigned char D, unsigned char E>
LazyStencil<T,C,D,E>::LazyStencil(Stencil<T> a_stencil,
                                  const BoxData<T,C,D,E>* a_src,
                                  Box a_box,
                                  T a_scale)
{
    m_range = a_stencil.range(a_src->box()) & a_box.refine(a_stencil.destRatio());
    m_stencil.push_back(a_stencil);
    m_src.push_back(const_cast<BoxData<T,C,D,E>*>(a_src));
    m_box.push_back(a_box);
    m_scale.push_back(a_scale);
    PROTO_ASSERT(m_stencil.size() == m_src.size(),
    "LazyStencil data corrupted; m_stencil and m_src have different sizes.");
    PROTO_ASSERT(m_stencil.size() == m_src.size(),
    "LazyStencil data corrupted; m_stencil and m_box have different sizes.");
    PROTO_ASSERT(m_stencil.size() == m_src.size(),
    "LazyStencil data corrupted; m_stencil and m_scale have different sizes.");
}

template <typename T, unsigned int C, unsigned char D, unsigned char E>
void LazyStencil<T,C,D,E>::apply(BoxData<T,C,D,E>& a_dest,
                                 bool a_overwrite)
{
    // all stencils in the LazyStencil have the same src- / dest-ratio
    Point destRatio = m_stencil.back().destRatio();
    for (int ii = 0; ii < m_src.size(); ii++)
    {
        Box b = m_box[ii];
        Point shift = m_stencil[ii].destShift();
        b = b & m_stencil[ii].indexDomain(a_dest.box());
        /*
        if ((m_stencil[ii].destRatio() == Point::Ones()))
        {
            b = b & a_dest.box();
        }
        else if ((m_stencil[ii].destRatio() != Point::Ones())
         && (m_stencil[ii].srcRatio() == Point::Ones())) //fix interpolation ranges
        {
            b = b & m_stencil[ii].domain(a_dest.box());
        }
        */
        m_stencil[ii].apply(*(m_src[ii]),a_dest,
                            b, a_overwrite, m_scale[ii]);
        
    }
}

//=======================================================================================
// InterpStencil ||
//===============++ 

template<class T>
template<unsigned int C, unsigned char D, unsigned char E>
LazyStencil<T,C,D,E> InterpStencil<T>::operator()(const BoxData<T,C,D,E>& a_src,
                                                  Box                      a_box,
                                                  T                       a_scale) const
{
    //if (!m_closed){this->close();}

    LazyStencil<T,C,D,E> ret;
    if (m_s.size() <= 0){return ret;}
    Box bsrc = a_src.box();
    for (int ii = 0; ii < m_s.size(); ii++)
    {
        Box b = Box(bsrc.low() - m_s[ii].span().low(), bsrc.high() - m_s[ii].span().high());
        b = b.taperCoarsen(m_s[ii].srcRatio());
        ret.m_stencil.push_back(m_s[ii]);
        ret.m_src.push_back(const_cast<BoxData<T,C,D,E>*>(&a_src));
        ret.m_box.push_back(a_box & b);
        ret.m_scale.push_back(a_scale);
    }
    return ret;
}

template<class T>
template<unsigned int C, unsigned char D, unsigned char E>
inline LazyStencil<T,C,D,E> InterpStencil<T>::operator()(const BoxData<T,C,D,E>& a_src,
                                                         T a_scale) const
{
    //if (!m_closed){this->close();}
    
    LazyStencil<T,C,D,E> ret;
    if (m_s.size() <= 0){return ret;}
    Box bsrc = a_src.box();

    for (int ii = 0; ii < m_s.size(); ii++)
    {
        //compute range for building the output BoxData, if necessary
        ret.m_range = ret.m_range & m_s[ii].range(bsrc).low();
        ret.m_range = ret.m_range & m_s[ii].range(bsrc).high();

        //compute the largest possible computation box (different from range) 
        Box b = Box(bsrc.low() - m_s[ii].span().low(), bsrc.high() - m_s[ii].span().high());
        b = b.taperCoarsen(m_s[ii].srcRatio());
        ret.m_stencil.push_back(m_s[ii]);
        ret.m_src.push_back(const_cast<BoxData<T,C,D,E>*>(&a_src));
        ret.m_box.push_back(b);
        ret.m_scale.push_back(a_scale);
    }
    return ret;
}


//======================================================================
//                            CUDA STUFF                                
//======================================================================
#ifdef PROTO_CUDA

__device__ __host__
inline Point compute_box_point(Point a_low, Point a_high, unsigned int a_index)
{
  Point retval;
  int index = a_index;
  for (int ii = 0; ii < DIM; ii++)
  {
    int factor = (a_high[ii] - a_low[ii] + 1);
    int kred = index % factor;
    retval[ii] = kred + a_low[ii];
    index = (index - kred)/factor;
  }
  return retval;
}
/**/
template <typename T>
__global__
void stencilIndexer(int  a_begin, int a_end, int a_ncoef,
                    thrust::device_ptr<T>     a_coeff,
                    thrust::device_ptr<Point> a_offset, 
                    bool a_initToZero,  T* a_src,  T* a_dst,
                    Point a_boxLo, Point a_boxHi, 
                    Point a_srcLo, Point a_srcHi, 
                    Point a_dstLo, Point a_dstHi,
                    Point a_dstShift, Point a_srcRefRat, Point a_dstRefRat, 
                    T a_scale)
{
  int idx = a_begin + threadIdx.x+blockIdx.x*blockDim.x;

  if((idx >= a_begin) && (idx< a_end))
  {

    Point boxPt = compute_box_point(a_boxLo, a_boxHi, idx);
    Point dstPt = boxPt;
    dstPt *= a_dstRefRat;
    dstPt += a_dstShift;
    unsigned int dstOff = Box::sindex(dstPt, a_dstLo, a_dstHi);
    if(a_initToZero)
    {
      a_dst[dstOff] = 0.;
    }
    for(int ivec = 0; ivec < a_ncoef; ivec++)
    {
      const     T & coef   = a_coeff[ivec];
      const Point & offset = a_offset[ivec];
      Point srcSt = boxPt*a_srcRefRat;
      Point srcPt = srcSt + offset;

      unsigned int srcOff = Box::sindex(srcPt, a_srcLo, a_srcHi);
      
      a_dst[dstOff] += a_scale*coef*a_src[srcOff];
    }
  }
}
/**/
// Apply Stencil
template <typename T>
template <unsigned int C, unsigned char D, unsigned char E>
void Stencil<T>::cudaApply( BoxData<T,C,D,E>&  a_src,
                            BoxData<T,C,D,E>&  a_dst,
                            Box&                a_box,
                            bool               a_initToZero,
                            T                  a_scale) 
{
  // Check input domain Boxes for proper nesting
  PROTO_ASSERT(!a_src.isAlias(a_dst),
               "Stencil::cudaApply(...) invalid. a_src and a_dst alias the same data.");
   
  //CH_TIMERS("Stencil::apply");
  // Initialize a_dst to zero within a_box if flag is set
  for (int ii = 0; ii < DIM; ii++)
  {
    PROTO_ASSERT((m_destShift[ii] >= 0) && (m_destShift[ii] < m_destRefratio[ii]), 
                 "m_destShift[ii] must be non-negative and less than m_destRefratio[ii] for all ii in [0,DIM)");
  }

  //this copies to the device
  if(!m_isClosed)
  {
    d_coeff  = m_coefs;
    d_offset = m_offsets;
    m_isClosed = true;
  }
    

  // apply the stencil
  for (int ee = 0; ee < E; ee++)
  {    
    for (int dd = 0; dd < D; dd++)
    {    
      for (int cc = 0; cc < C; cc++)
      {    
        BoxData<T,1,1,1> srcbd = slice(a_src, cc, dd, ee);
        BoxData<T,1,1,1> dstbd = slice(a_dst, cc, dd, ee);
        int begin  = 0;
        int end    = a_box.size();
        int stride = a_box.size(0);
        int blocks = a_box.flatten(0).size();

        T* d_srcPtr = (T*)srcbd.data();
        T* d_dstPtr = (T*)dstbd.data();
        Point srcLo = srcbd.box().low();
        Point srcHi = srcbd.box().high();
        Point dstLo = dstbd.box().low();
        Point dstHi = dstbd.box().high();
        Point boxLo = a_box.low();
        Point boxHi = a_box.high();
        unsigned int ncoef = m_coefs.size();
        const thrust::device_ptr<T>       coefptr = d_coeff.data();
        const thrust::device_ptr<Point> offsetptr = d_offset.data();
        stencilIndexer<<<blocks, stride>>>(begin, end,  ncoef, coefptr, offsetptr,
                                           a_initToZero, d_srcPtr, d_dstPtr, 
                                           boxLo, boxHi, srcLo, srcHi, dstLo, dstHi, 
                                           m_destShift, m_srcRefratio, m_destRefratio, a_scale);
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

/**/
#endif
