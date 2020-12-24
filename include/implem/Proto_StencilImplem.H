#pragma once
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

template <typename T>
unsigned long long int Stencil<T>::numFlops(const Box& a_box) const
{
  unsigned long long int numpts = a_box.size();
  unsigned long long int numflops = 3*size(); //multiplies by scale and  by weight and add
  return numpts*numflops;
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

template <typename T>
Point Stencil<T>::spanPoint() const
{
    int v[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        v[ii] = std::max(std::abs(m_span.low()[ii]), std::abs(m_span.high()[ii]));
    }
    return Point(v);
}

template <typename T>
Point Stencil<T>::ghost() const
{
    /* I apologize if I broke this, but I don't think this works.
    Point ghost;
    for (int ii = 0; ii < DIM; ii++)
    {
        ghost[ii] = std::max(-m_span.low()[ii], m_span.high()[ii]);
    }
    return ghost;
    */
    return spanPoint();
}

// Augment Stencil
template <typename T>
void
Stencil<T>::augment(T a_coef, Point a_offset)
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

template <typename T>
void
Stencil<T>::invert(int a_dir)
{
    PROTO_ASSERT(a_dir >= 0 && a_dir < DIM,
            "Stencil::transpose(int a, int b) invalid. Both a and b must be in [0,DIM = %i)",DIM);
    for (int ii = 0; ii < size(); ii++)
    {
        int tuple[DIM];
        for (int d = 0; d < DIM; d++)
        {
            if (d == a_dir)
            {
                tuple[d] = -m_offsets[ii][d];
            } else {
                tuple[d] = m_offsets[ii][d];
            }
        }
        Point newOffset(tuple);
        m_offsets[ii] = newOffset;
    }
}

template <typename T>
void
Stencil<T>::transpose(unsigned char a, unsigned char b)
{
    PROTO_ASSERT((a < DIM) && (b < DIM),
            "Stencil::transpose(int a, int b) invalid. Both a and b must be in [0,DIM = %i)",DIM);
#if DIM < 2
    return;
#else
    if (a == b) {return;}
    for (int l = 0; l < size(); l++)
    {
        std::swap(m_offsets[l][a], m_offsets[l][b]);
    }
#endif
}

template <typename T>
Box
Stencil<T>::indexRange(Box a_domain) const
{
    Point newlo, newhi;
    for(int idir = 0; idir < DIM; idir++)
    {
        newlo[idir] = a_domain.low() [idir] + std::abs(m_span.low() [idir]);
        newhi[idir] = a_domain.high()[idir] - std::abs(m_span.high()[idir]);
    }
    Box indexRange(newlo, newhi);
    indexRange = indexRange.taperCoarsen(m_srcRefratio);
    return indexRange;
}

template <typename T>    
Box
Stencil<T>::indexDomain(Box a_range) const
{
    Box indexDomain = a_range.shift(-m_destShift);
    indexDomain = indexDomain.taperCoarsen(m_destRefratio);
    return indexDomain;
}

template <typename T>
Box
Stencil<T>::range(Box a_domain) const
{
    Box range;
    range = indexRange(a_domain);
    range = Box(range.low()*m_destRefratio, range.high()*m_destRefratio);
    range = range.shift(m_destShift);
    return range;
}

template <typename T>
Box
Stencil<T>::domain(Box a_range) const
{
    Box domain = indexDomain(a_range);
    domain = Box(domain.low()*m_srcRefratio, domain.high()*m_srcRefratio);
    domain = Box(domain.low() + m_span.low(), domain.high()+m_span.high());
    return domain;
}

template <typename T>
T
Stencil<T>::diagonalValue() const
{
    T retval = 0;
    for(int ipt = 0; ipt < m_offsets.size(); ipt++)
    {
        if(m_offsets[ipt] == Point::Zeroes())
        {
            retval += m_coefs[ipt];
        }
    }
    return retval;
}

template <typename T>
void
Stencil<T>::print() const
{
    using ::std::cout;
    using ::std::endl;

    cout << "Coeficients and Shifts: { ";
    for (int ii = 0; ii < m_coefs.size(); ii++) {
        if (ii == 0){cout << endl;}
        cout << std::setfill(' ') << std::setw(8) << m_coefs[ii];
        cout << ": " << m_offsets[ii] << endl;
    }
    cout << "}" << endl;
    cout << "srcRefratio:\t" << m_srcRefratio << endl;
    cout << "destRefratio:\t" << m_destRefratio << endl;
    cout << "destShift:\t" << m_destShift << endl;
    cout << "span:\t\t" << m_span << endl;
    cout << endl;
}
//=======================================================================================
// Apply ||
//=======++ 

// Apply Stencil
template <typename T>
template <unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE>
void Stencil<T>::apply(const BoxData<T,C,D,E,MEMTYPE>&  a_src,
                       BoxData<T,C,D,E,MEMTYPE>&        a_dest,
                       const Box&               a_box,
                       bool                     a_initToZero,
                       T                        a_scale) const
{
    PR_TIME("Stencil::apply");
    BoxData<T, C, D, E, MEMTYPE>& castsrc = const_cast<BoxData<T,C,D,E,MEMTYPE> &>(a_src);
    Box & castbox = const_cast<Box &>(a_box);
    Stencil<T>* castthis = const_cast<Stencil<T>* >(this);
    castthis->protoApply(castsrc, a_dest, castbox, a_initToZero, a_scale);
    PR_FLOPS(this->numFlops(a_box));
}

template <typename T>
template<unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE>
LazyStencil<T,C,D,E,MEMTYPE>
Stencil<T>::operator()(const BoxData<T,C,D,E,MEMTYPE>& a_src, T a_scale) const
{
    Box b = this->indexRange(a_src.box());
    return LazyStencil<T,C,D,E,MEMTYPE>(this, &a_src, b, a_scale);
}
template <typename T>
template<unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE>
LazyStencil <T,C,D,E,MEMTYPE>
Stencil<T>::operator()(const BoxData<T,C,D,E,MEMTYPE>&  a_src, Box a_box, T a_scale) const
{
#ifndef NDEBUG
    Box b = this->indexRange(a_src.box());
    PROTO_ASSERT(b.contains(a_box),
            "Stencil::operator() invalid. User provided a Box larger than the valid range of the input data.");
#endif
    return LazyStencil<T, C, D, E, MEMTYPE>(this, &a_src, a_box, a_scale);
}

template <typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE>
BoxData<T,C,D,E,MEMTYPE>&
operator|=(BoxData<T,C,D,E,MEMTYPE>& a_dest, LazyStencil<T,C,D,E,MEMTYPE>&& a_op)
{
    if (a_op.size() > 0)
    {
        for (auto iter = a_op.m_src.begin(); iter != a_op.m_src.end(); ++iter)
        {
            PROTO_ASSERT(!(*iter)->isAlias(a_dest),
            "BoxData::operator|=(BoxData) invalid. Stencil source and destination data are aliased");
        }
        a_op.apply(a_dest, true);
    }
    return a_dest; 
}

template <class T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE>
BoxData<T,C,D,E,MEMTYPE>&
operator+=(BoxData<T,C,D,E,MEMTYPE>& a_dest, LazyStencil<T,C,D,E,MEMTYPE>&& a_op)
{
    if (a_op.size() > 0)
    {
        for (auto iter = a_op.m_src.begin(); iter != a_op.m_src.end(); ++iter)
        {
            PROTO_ASSERT(!(*iter)->isAlias(a_dest),
            "BoxData::operator|=(BoxData) invalid. Stencil source and destination data are aliased");
        }
        a_op.apply(a_dest,false); 
    }
    return a_dest;
}

//=======================================================================================
// LazyStencil ||
//=============++ 

template <typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE>
LazyStencil<T,C,D,E,MEMTYPE>::LazyStencil(const Stencil<T>* a_stencil,
                                  const BoxData<T,C,D,E,MEMTYPE>* a_src,
                                  Box a_box,
                                  T a_scale)
{
    m_range = a_stencil->range(a_src->box()) & a_box.refine(a_stencil->destRatio());
    m_stencil.push_back(a_stencil);
    m_src.push_back(const_cast<BoxData<T,C,D,E,MEMTYPE>*>(a_src));
    m_box.push_back(a_box);
    m_scale.push_back(a_scale);
    PROTO_ASSERT(m_stencil.size() == m_src.size(),
    "LazyStencil data corrupted; m_stencil and m_src have different sizes.");
    PROTO_ASSERT(m_stencil.size() == m_src.size(),
    "LazyStencil data corrupted; m_stencil and m_box have different sizes.");
    PROTO_ASSERT(m_stencil.size() == m_src.size(),
    "LazyStencil data corrupted; m_stencil and m_scale have different sizes.");
}

template <typename T, unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE>
void LazyStencil<T,C,D,E,MEMTYPE>::apply(BoxData<T,C,D,E,MEMTYPE>& a_dest,
					 bool a_overwrite)
{
    // all stencils in the LazyStencil have the same src- / dest-ratio
  //    Point destRatio = m_stencil.back()->destRatio();
    for (int ii = 0; ii < m_src.size(); ii++)
    {
        Box b = m_box[ii];
        //    Point shift = m_stencil[ii]->destShift();
        b = b & m_stencil[ii]->indexDomain(a_dest.box());
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
        m_stencil[ii]->apply(*(m_src[ii]),a_dest,
                            b, a_overwrite, m_scale[ii]);
        
    }
}

//=======================================================================================
// InterpStencil ||
//===============++ 

template<class T>
template<unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE>
LazyStencil<T,C,D,E,MEMTYPE> InterpStencil<T>::operator()(const BoxData<T,C,D,E,MEMTYPE>& a_src,
							  Box                      a_box,
							  T                       a_scale) const
{
    //if (!m_closed){this->close();}

    LazyStencil<T,C,D,E,MEMTYPE> ret;
    if (m_s.size() <= 0){return ret;}
    Box bsrc = a_src.box();
    for (int ii = 0; ii < m_s.size(); ii++)
    {
        Box b = Box(bsrc.low() - m_s[ii].span().low(), bsrc.high() - m_s[ii].span().high());
        b = b.taperCoarsen(m_s[ii].srcRatio());
        ret.m_stencil.push_back(&(m_s[ii]));
        ret.m_src.push_back(const_cast<BoxData<T,C,D,E,MEMTYPE>*>(&a_src));
        ret.m_box.push_back(a_box & b);
        ret.m_scale.push_back(a_scale);
    }
    return ret;
}

template<class T>
template<unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE>
inline LazyStencil<T,C,D,E,MEMTYPE> InterpStencil<T>::operator()(const BoxData<T,C,D,E,MEMTYPE>& a_src,
								 T a_scale) const
{
    //if (!m_closed){this->close();}
    
    LazyStencil<T,C,D,E,MEMTYPE> ret;
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
        ret.m_stencil.push_back(&(m_s[ii]));
        ret.m_src.push_back(const_cast<BoxData<T,C,D,E,MEMTYPE>*>(&a_src));
        ret.m_box.push_back(b);
        ret.m_scale.push_back(a_scale);
    }
    return ret;
}

template <typename T>
inline void  
Stencil<T>::
closeForDevice()
{
  if(!m_isClosed)
  {
    size_t coeffsize   = m_coefs.size()*sizeof(T);
    size_t offsetsize  = m_offsets.size()*sizeof(Point);
    protoMalloc(d_coeff , coeffsize);
    protoMalloc(d_offset, offsetsize);
    T*     h_coeff  = m_coefs.data();
    Point* h_offset = m_offsets.data();
    protoMemcpy(d_coeff , h_coeff,   coeffsize,protoMemcpyHostToDevice);
    protoMemcpy(d_offset, h_offset, offsetsize,protoMemcpyHostToDevice);
    m_isClosed = true;
    PROTO_ASSERT(m_coefs.size()<64, "m_coefs.size()>= 64");
    for(int i=0; i<m_coefs.size(); i++)
    {
          c_coeff.elements[i]=m_coefs[i];
    }
  }
}

//======================================================================
//                            CUDA STUFF                                
//======================================================================
#ifdef PROTO_CUDA
#define PR_MAX_COEFFS 27

/**/
// Apply Stencil

template<typename T,size_t C>
__device__
inline void stencil(const int id, unsigned int a_nx, unsigned int a_ny, unsigned int a_nz, 
		     T* dstPtr, const T* srcPtr, 
		     std::size_t srcStride, std::size_t dstStride,
                     Point a_srcLo, Point a_srcHi, Point a_dstLo, Point a_dstHi,
                     Point a_boxLo, Point a_boxHi, bool a_initToZero,
                     unsigned int ssize, const T* a_coefs,
                     const int* a_offset, Point a_dstShift,
		     Point a_srcRefratio, Point a_destRefratio,
		     T a_scale)
{
  const int idx = id % a_nx;
  
#if DIM==3
  const int nxy = a_nx*a_ny;
  const int idy = (id % nxy) / a_nx;
  const int idz = id / nxy;
  if(idx < a_nx && idy < a_ny && idz < a_nx)
#else
  const int idy = id / a_nx;
  if(idx < a_nx && idy < a_ny)
#endif
  {
    Point a_dst=a_boxLo;
    a_dst[0]+=idx;
    a_dst[1]+=idy;
    #if DIM==3
    a_dst[2]+=idz;
    #endif

    Point dstPt =  a_dst;
    dstPt *= a_destRefratio;
    dstPt += a_dstShift;
    unsigned int dstOff = Box::sindex(dstPt, a_dstLo, a_dstHi);
    Point srcSt = a_dst*a_srcRefratio;
    unsigned int srcOff = Box::sindex(srcSt, a_srcLo, a_srcHi); 

    for(int c=0; c<C; c++, srcPtr +=  srcStride, dstPtr += dstStride)
    {
      T accum = 0;
#pragma unroll
      for(int ivec = 0; ivec < ssize; ivec++)
        accum += a_coefs[ivec]*srcPtr[srcOff+a_offset[ivec]];

      if(a_initToZero) dstPtr[dstOff] = a_scale*accum;
      else dstPtr[dstOff] += a_scale*accum;
    }
  }
}

template<typename T, size_t C>
__global__
void stencilIndexer3(int a_nx, int a_ny, T* dstPtr, const T* srcPtr, int srcStride, int dstStride,
                     Point a_srcLo, Point a_srcHi, Point a_dstLo, Point a_dstHi,
                     Point a_boxLo, Point a_boxHi, bool a_initToZero,
                     const int ssize, const T* a_coefs,
                     const Point* srcIndices, int ksize,
		     Point a_dstShift , Point a_srcRefratio, Point a_destRefratio,
		     T a_scale)
{
  int id = threadIdx.x+blockIdx.x*blockDim.x;
  __shared__ int s_offset[PR_MAX_COEFFS];

  for(int itid = threadIdx.x; (itid < ssize)&&(itid < PR_MAX_COEFFS); itid+=blockDim.x )
  {
    s_offset[itid]= Box::sindex(a_srcLo+srcIndices[itid], a_srcLo, a_srcHi);
  }

  __syncthreads();

  stencil<T,C>(id, a_nx, a_ny, ksize,
		  dstPtr, srcPtr,
		  srcStride, dstStride,
		  a_srcLo, a_srcHi, a_dstLo, a_dstHi,
		  a_boxLo, a_boxHi, a_initToZero,
		  ssize, a_coefs, s_offset,
		  a_dstShift, a_srcRefratio, a_destRefratio,
		  a_scale
	);
}

#else

template<typename T,size_t C>
inline void stencil(const int id, unsigned int a_nx, unsigned int a_ny, unsigned int a_nz, 
		     T* dstPtr, const T* srcPtr, 
		     std::size_t srcStride, std::size_t dstStride,
                     Point a_srcLo, Point a_srcHi, Point a_dstLo, Point a_dstHi,
                     Point a_boxLo, Point a_boxHi, bool a_initToZero,
                     unsigned int ssize, const T* a_coefs,
                     const int* a_offset, Point a_dstShift,
		     Point a_srcRefratio, Point a_destRefratio,
		     T a_scale)
{
  for(int idz = 0; idz < a_nz ; idz++)
    for(int idy = 0; idy < a_ny ; idy++)
      for(int idx = 0; idx < a_nx ; idx++)
      {
        const T* tmp_srcPtr= srcPtr;
        T* tmp_dstPtr= dstPtr;
        Point a_dst=a_boxLo;
        a_dst[0]+=idx;
        a_dst[1]+=idy;
#if DIM==3
        a_dst[2]+=idz;
#endif
        Point dstPt =  a_dst;
        dstPt *= a_destRefratio;
        dstPt += a_dstShift;
        unsigned int dstOff = Box::sindex(dstPt, a_dstLo, a_dstHi);
        Point srcSt = a_dst*a_srcRefratio;
        unsigned int srcOff = Box::sindex(srcSt, a_srcLo, a_srcHi); 

        for(int c=0; c<C; c++, tmp_srcPtr +=  srcStride, tmp_dstPtr += dstStride)
        {
          T accum = 0;
#pragma unroll
          for(int ivec = 0; ivec < ssize; ivec++)
            accum += a_coefs[ivec]*tmp_srcPtr[srcOff+a_offset[ivec]];

          if(a_initToZero) tmp_dstPtr[dstOff] = a_scale*accum;
          else tmp_dstPtr[dstOff] += a_scale*accum;
        }
      }
}

template<typename T, size_t C>
void stencilIndexer3(int a_nx, int a_ny, T* dstPtr, const T* srcPtr, int srcStride, int dstStride,
                     Point a_srcLo, Point a_srcHi, Point a_dstLo, Point a_dstHi,
                     Point a_boxLo, Point a_boxHi, bool a_initToZero,
                     const int ssize, const T* a_coefs,
                     const Point* srcIndices, int ksize,
		     Point a_dstShift , Point a_srcRefratio, Point a_destRefratio,
		     T a_scale)
{
  int s_offset[ssize];
  for(int itid = 0; itid < ssize; itid++)
  {
    s_offset[itid]= Box::sindex(a_srcLo+srcIndices[itid], a_srcLo, a_srcHi);
  }

  unsigned int id = 0; // useless parameter on cpu

  stencil<T,C>(id, a_nx, a_ny, ksize,
		  dstPtr, srcPtr,
		  srcStride, dstStride,
		  a_srcLo, a_srcHi, a_dstLo, a_dstHi,
		  a_boxLo, a_boxHi, a_initToZero,
		  ssize, a_coefs, s_offset,
		  a_dstShift, a_srcRefratio, a_destRefratio,
		  a_scale
	);
}
#endif

template <typename T>
template <unsigned int C, unsigned char D, unsigned char E, MemType MEMTYPE>
void Stencil<T>::protoApply( const BoxData<T,C,D,E,MEMTYPE> &  a_src,
                            BoxData<T,C,D,E,MEMTYPE> &  a_dst,
                            const Box              &  a_box,
                            bool                a_initToZero,
                            T                   a_scale)
{
  PR_TIME("protoApply");

  if(!m_isClosed && MEMTYPE==MemType::DEVICE) closeForDevice();
  if(a_box.size() == 0) return;

#if DIM==3
  int k= a_box.size(2);
#else
  int k=1;
#endif

  for (int ee = 0; ee < E; ee++)
    {
      for (int dd = 0; dd < D; dd++)
        {
          const T* srcPtr = a_src.dataPtr(0, dd, ee);
          T* dstPtr = a_dst.dataPtr(0, dd, ee);
	  T* coeffPtr = (MEMTYPE==MemType::DEVICE) ? d_coeff : m_coefs.data();
	  Point* offsetPtr = (MEMTYPE==MemType::DEVICE) ? d_offset : m_offsets.data();
          protoLaunchKernelMemAsync((stencilIndexer3<T,C>),
	        a_box.flatten(0).size(),	//blocks,//blocks,
		a_box.size(0),//strides,
		0,
            	protoGetCurrentStream, a_box.size(0), a_box.size(1), 
		dstPtr, srcPtr, 
		a_src.box().size(), a_dst.box().size(),
		a_src.box().low(), a_src.box().high(), 
		a_dst.box().low(), a_dst.box().high(), 
		a_box.low(), a_box.high(), a_initToZero,
		m_coefs.size(), coeffPtr, offsetPtr, k, m_destShift, 
		m_srcRefratio, m_destRefratio, a_scale);
        }
    }
  //PR_FLOPS(this->numFlops(a_box));
}

/**/
