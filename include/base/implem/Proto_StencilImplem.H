#pragma once
//=======================================================================================
// Constructors ||
//==============++ 

// Default Constructor
template <typename T>
Stencil<T>::Stencil()
{
    PR_TIME("Stencil::constructor(default)");
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
    PR_TIME("Stencil::constructor(general)");
    m_isClosed = false;
    m_destRefratio = a_destRefratio;
    m_srcRefratio = a_srcRefratio;
    m_destShift = a_destShift;
    addCoef(a_coef, a_shift.point());
}

// Destructor
template <typename T>
Stencil<T>::~Stencil()
{
    PR_TIME("Stencil::destructor");
    if(m_isClosed)
    {
#ifdef PROTO_ACCEL
        proto_free<DEVICE>(d_coeff);
        proto_free<DEVICE>(d_offset);
        //protoFreeGPU(d_coeff);
        //protoFreeGPU(d_offset);
#endif
        m_isClosed = false;
    }
}


template <typename T>
unsigned long long int Stencil<T>::numFlops(const Box& a_box) const
{
  unsigned long long int numpts = a_box.size();
  unsigned long long int numflops = 2*size()+1; //multiplies by scale and  by weight and add
  return numpts*numflops;
}
//=======================================================================================
// Operators ||
//===========++ 

// Stencil Composition
template <typename T>
Stencil<T> Stencil<T>::operator* (const Stencil<T>& a_stencil) const
{
    PR_TIME("Stencil::operator*(Stencil)");
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
    PR_TIME("Stencil::operator*(Scalar)");
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
    PR_TIME("Stencil::operator*=(Stencil)");
    //std::move to avoid copying the product
    (*this) = std::move((*this)*a_stencil);
}   

// In Place Scalar Multiplication
template <typename T>
void Stencil<T>::operator*=(const T a_coef)
{
    PR_TIME("Stencil::operator*=(Scalar)");
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
    PR_TIME("Stencil::operator+=(Stencil)");
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
        addCoef(coef, offset);
    }
}

// In Place Subtraction
template <class T>
void Stencil<T>::operator-=(const Stencil<T>& a_stencil)
{
    PR_TIME("Stencil::operator-=(Stencil)");
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
        addCoef(coef, offset);
    }
}

// Equality Operator
template <typename T>
bool Stencil<T>::operator==(const Stencil<T>& a_stencil) const
{
    PR_TIME("Stencil::operator==(Stencil)");
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
Point Stencil<T>::ghost() const
{
    int v[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        v[ii] = std::max(std::abs(m_span.low()[ii]), std::abs(m_span.high()[ii]));
    }
    return Point(v);
}

// Augment Stencil
template <typename T>
void
Stencil<T>::addCoef(T a_coef, Point a_offset)
{
    PR_TIME("Stencil::addCoef");
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
    PR_TIME("Stencil::invert");
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
    PR_TIME("Stencil::transpose");
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
    PR_TIME("Stencil::indexRange");
    if (a_domain.empty()) { return Box(); }
    Point L = a_domain.low() - m_span.low();
    Point H = a_domain.high() - m_span.high();
    Box indexRange(L, H);
    indexRange = indexRange.taperCoarsen(m_srcRefratio);
    return indexRange;
}

template <typename T>    
Box
Stencil<T>::indexDomain(Box a_range) const
{
    PR_TIME("Stencil::indexDomain");
    if (a_range.empty()) { return Box(); }
    Box indexDomain = a_range.shift(-m_destShift);
    indexDomain = indexDomain.taperCoarsen(m_destRefratio);
    return indexDomain;
}

template <typename T>
Box
Stencil<T>::range(Box a_domain) const
{
    PR_TIME("Stencil::range");
    if (a_domain.empty()) { return Box(); }
    Box range;
    range = indexRange(a_domain);
    // this is intentionally not a refine
    range = Box(range.low()*m_destRefratio, range.high()*m_destRefratio);
    range = range.shift(m_destShift);
    return range;
}

template <typename T>
Box
Stencil<T>::domain(Box a_range) const
{
    PR_TIME("Stencil::domain");
    if (a_range.empty()) { return Box(); }
    Box domain = indexDomain(a_range);
    // this is intentionally not a refine
    domain = Box(domain.low()*m_srcRefratio, domain.high()*m_srcRefratio);
    domain = Box(domain.low() + m_span.low(), domain.high()+m_span.high());
    return domain;
}

template <typename T>
T
Stencil<T>::diagonal() const
{
    PR_TIME("Stencil::diagonal");
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
    PR_TIME("Stencil::print");
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
template <unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void Stencil<T>::apply(const BoxData<T,C,MEMTYPE,D,E>&  a_src,
                       BoxData<T,C,MEMTYPE,D,E>&        a_dest,
                       const Box&               a_box,
                       bool                     a_initToZero,
                       T                        a_scale) const
{
    PR_TIME("Stencil::apply");
    BoxData<T, C, MEMTYPE, D, E>& castsrc = const_cast<BoxData<T,C,MEMTYPE,D,E> &>(a_src);
    Box & castbox = const_cast<Box &>(a_box);
    Stencil<T>* castthis = const_cast<Stencil<T>* >(this);
    // Box inference happens here
    if (a_box.empty())
    {
        castbox = indexDomain(a_dest.box()) & indexRange(a_src.box());
    }
    castthis->protoApply(castsrc, a_dest, castbox, a_initToZero, a_scale);
}

template <typename T>
template<unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
LazyStencil<T,C,MEMTYPE,D,E>
Stencil<T>::operator()(const BoxData<T,C,MEMTYPE,D,E>& a_src, T a_scale) const
{
    //Box b = this->range(a_src.box());
    return LazyStencil<T,C,MEMTYPE,D,E>(this, &a_src, Box(), a_scale);
}
template <typename T>
template<unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
LazyStencil <T,C,MEMTYPE,D,E>
Stencil<T>::operator()(const BoxData<T,C,MEMTYPE,D,E>&  a_src, Box a_box, T a_scale) const
{
#ifndef NDEBUG
    Box b = this->range(a_src.box());
    PROTO_ASSERT(b.contains(a_box),
            "Stencil::operator() | Error: User provided Box is larger than the maximum range.");
    PROTO_ASSERT(!a_box.empty(),
            "Stencil::operator() | Error: User explicitly passed an empty Box.");
#endif
    return LazyStencil<T, C, MEMTYPE, D, E>(this, &a_src, a_box, a_scale);
}

template <typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
BoxData<T,C,MEMTYPE,D,E>&
operator|=(BoxData<T,C,MEMTYPE,D,E>& a_dest, LazyStencil<T,C,MEMTYPE,D,E>&& a_op)
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

template <class T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
BoxData<T,C,MEMTYPE,D,E>&
operator+=(BoxData<T,C,MEMTYPE,D,E>& a_dest, LazyStencil<T,C,MEMTYPE,D,E>&& a_op)
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

template <typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
LazyStencil<T,C,MEMTYPE,D,E>::LazyStencil(const Stencil<T>* a_stencil,
                                  const BoxData<T,C,MEMTYPE,D,E>* a_src,
                                  Box a_box,
                                  T a_scale)
{
    m_range = a_box; // user supplies the range box
    m_stencil.push_back(a_stencil);
    m_src.push_back(const_cast<BoxData<T,C,MEMTYPE,D,E>*>(a_src));
    m_box.push_back(a_stencil->indexDomain(a_box)); // apply box is in index space
    m_scale.push_back(a_scale);
    PROTO_ASSERT(m_stencil.size() == m_src.size(),
    "LazyStencil data corrupted; m_stencil and m_src have different sizes.");
    PROTO_ASSERT(m_stencil.size() == m_src.size(),
    "LazyStencil data corrupted; m_stencil and m_box have different sizes.");
    PROTO_ASSERT(m_stencil.size() == m_src.size(),
    "LazyStencil data corrupted; m_stencil and m_scale have different sizes.");
}

template <typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
void LazyStencil<T,C,MEMTYPE,D,E>::apply(BoxData<T,C,MEMTYPE,D,E>& a_dest,
					 bool a_overwrite)
{
    for (int ii = 0; ii < m_src.size(); ii++)
    {
        m_stencil[ii]->apply(*(m_src[ii]),a_dest,
                            m_box[ii], a_overwrite, m_scale[ii]);
    }
}

template <typename T, unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
Box LazyStencil<T,C,MEMTYPE,D,E>::inferredRange() const
{
    if (size() == 0) {return Box(); }
    Box range;
    for (int ii = 0; ii < size(); ii++)
    {
        Box ri = m_stencil[0]->range(m_src[0]->box());
        if (ii == 0) {range = ri;}
        else {
            range = range & ri.low();
            range = range & ri.high();
        }
    }
    return range;
}

template <typename T>
inline void Stencil<T>::closeForDevice()
{
    if(!m_isClosed)
    {
#ifdef PROTO_ACCEL
        size_t coeffsize   = m_coefs.size()*sizeof(T);
        size_t offsetsize  = m_offsets.size()*sizeof(Point);
        d_coeff = (T*)proto_malloc<DEVICE>(coeffsize);
        d_offset = (Point*)proto_malloc<DEVICE>(offsetsize);
        T*     h_coeff  = m_coefs.data();
        Point* h_offset = m_offsets.data();
        proto_memcpy<HOST, DEVICE>(h_coeff, d_coeff, coeffsize);
        proto_memcpy<HOST, DEVICE>(h_offset, d_offset, offsetsize);
#endif
        m_isClosed = true;
    }
}

//======================================================================
//                            CUDA STUFF                                
//======================================================================
/**/
// Apply Stencil

#ifdef PROTO_ACCEL
template<typename T,size_t C>
__device__
inline void stencil_gpu(const int id, unsigned int a_nx, unsigned int a_ny, unsigned int a_nz, 
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

#endif

template<typename T,size_t C>
inline void stencil_cpu(const int id, unsigned int a_nx, unsigned int a_ny, unsigned int a_nz, 
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

template<typename T, unsigned int C>
struct sIndexer3 {
  static void cpu(int a_nx, int a_ny, T* dstPtr, const T* srcPtr, int srcStride, int dstStride,
                     Point a_srcLo, Point a_srcHi, Point a_dstLo, Point a_dstHi,
                     Point a_boxLo, Point a_boxHi, bool a_initToZero,
                     const int ssize, const T* a_coefs,
                     const Point* srcIndices, int ksize,
		     Point a_dstShift , Point a_srcRefratio, Point a_destRefratio,
		     T a_scale)
  {
    int s_offset[PR_MAX_COEFFS];
    for(int itid = 0; itid < ssize; itid++)
      {
        s_offset[itid]= Box::sindex(a_srcLo+srcIndices[itid], a_srcLo, a_srcHi);
      }

    unsigned int id = 0; // useless parameter on cpu

    stencil_cpu<T,C>(id, a_nx, a_ny, ksize,
		  dstPtr, srcPtr,
		  srcStride, dstStride,
		  a_srcLo, a_srcHi, a_dstLo, a_dstHi,
		  a_boxLo, a_boxHi, a_initToZero,
		  ssize, a_coefs, s_offset,
		  a_dstShift, a_srcRefratio, a_destRefratio,
		  a_scale
                 );
  }
#ifdef PROTO_ACCEL
  __device__ static
  void gpu(int a_nx, int a_ny, T* dstPtr, const T* srcPtr, int srcStride, int dstStride,
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

    stencil_gpu<T,C>(id, a_nx, a_ny, ksize,
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


};


template <typename T>
template <unsigned int C, unsigned char D, unsigned char E>
void Stencil<T>::hostApply(const BoxData<T,C,MemType::HOST,D,E>&  a_src,
        BoxData<T,C,MemType::HOST,D,E>&        a_dest,
        const Box&               a_box,
        bool                     a_initToZero,
        T                        a_scale) const                      
{
    /**/
    PR_TIME("Stencil::hostApply");
    if (m_coefs.size() == 0){return;}

    Box cross = a_box.flatten(0);
    int npencil = a_box.size(0);
    if (a_initToZero)
    {
        PR_TIME("Stencil::hostApply::initToZero");
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
        PR_TIME("Stencil::hostApply::initToZero");
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
    {
        PR_TIME("Stencil::hostApply::linearizeOffsets");
        offsets.resize(this->size());
        for (int ii = 0; ii < this->size(); ii++)
        {
            offsets[ii] = 0;
            for (int dir = 0; dir < DIM; dir++)
            {
                offsets[ii] += m_offsets[ii][dir]*offsetFactor[dir];
            }
        }
    }

    // apply the stencil

    {
        PR_TIME("Stencil::hostApply::applyStencil");
        for (int ee = 0; ee < E; ee++)
        for (int dd = 0; dd < D; dd++)
        for (int cc = 0; cc < C; cc++)
        {    
            for (auto iter = cross.begin(); iter != cross.end(); ++iter)
            {
                Point dpt = (*iter)*m_destRefratio + m_destShift;
                if (!a_dest.box().contains(dpt)){continue;} //can happen when destShift is non-trivial
                for (int jj = 0; jj < this->size(); jj++)
                {
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
    }
}

template <typename T>
template <unsigned int C, unsigned char D, unsigned char E>
void Stencil<T>::protoApply( const BoxData<T,C,HOST,D,E> &  a_src,
                            BoxData<T,C,HOST,D,E> &  a_dst,
                            const Box              &  a_box,
                            bool                a_initToZero,
                            T                   a_scale)
{
    PR_TIME("Stencil::protoApply");
    if (a_box.size() == 0) { return; }
    PR_FLOPS(this->numFlops(a_box));
    hostApply(a_src, a_dst, a_box, a_initToZero, a_scale);
}

template <typename T>
    template <unsigned int C, unsigned char D, unsigned char E>
void Stencil<T>::protoApply( const BoxData<T,C,DEVICE,D,E> &  a_src,
        BoxData<T,C,DEVICE,D,E> &  a_dst,
        const Box              &  a_box,
        bool                a_initToZero,
        T                   a_scale)
{
    PR_TIME("Stencil::protoApply");
    if (a_box.size() == 0) { return; }
    PR_FLOPS(this->numFlops(a_box));
#ifdef PROTO_ACCEL
    // copy coefficient data to DEVICE if needed
    closeForDevice();
#if DIM==3
    int k= a_box.size(2);
#else
    int k=1;
#endif
    int strides = a_box.size(0);
    int blocks  = a_box.flatten(0).size();
    if(strides > 512)
    {
        int nbElem = blocks * strides;
        v100tuning(nbElem,blocks,strides);
    }
    assert(strides < 1024);

    for (int ee = 0; ee < E; ee++)
    {
        for (int dd = 0; dd < D; dd++)
        {
            const T* srcPtr = a_src.data(0, dd, ee);
            T* dstPtr = a_dst.data(0, dd, ee);
            protoLaunchKernelMemAsyncT<DEVICE, sIndexer3<T, C>>(
                    a_box.flatten(0).size(),	//blocks,//blocks,
                    a_box.size(0),//strides,

                    0,
                    protoGetCurrentStream,
                    a_box.size(0), a_box.size(1), 
                    dstPtr, srcPtr, 
                    a_src.box().size(), a_dst.box().size(),
                    a_src.box().low(), a_src.box().high(), 
                    a_dst.box().low(), a_dst.box().high(), 
                    a_box.low(), a_box.high(), a_initToZero,
                    m_coefs.size(), d_coeff, d_offset, k, m_destShift, 
                    m_srcRefratio, m_destRefratio, a_scale);
        }
    }
#endif
}

/**/
