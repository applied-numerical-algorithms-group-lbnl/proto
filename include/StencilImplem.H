//=======================================================================================
// Constructors ||
//==============++ 

// Default Constructor
template <typename T>
Stencil<T>::Stencil()
{
    m_destRefratio=Point::Ones();
    m_destShift=Point::Zeros();
    m_srcRefratio=Point::Ones();
}

// General Constructor
template <typename T>
Stencil<T>::Stencil(Shift a_shift,
                    T     a_coef,
                    Point a_destRefratio,
                    Point a_destShift,
                    Point a_srcRefratio)
{
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
            m_span = Bx(low,high);
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
                       const Bx&                a_box,
                       bool                     a_initToZero,
                       T                        a_scale) const                      
{
#ifdef PROTO_CUDA
  
    BoxData<T, C, D, E>& castsrc = const_cast<BoxData<T,C,D,E> &>(a_src);
    Bx & castbox = const_cast<Bx &>(a_box);
    Stencil<T>* castthis = const_cast<Stencil<T>* >(this);
    castthis->cudaApply(castsrc, a_dest, castbox, a_initToZero, a_scale);

#else

    //debugging cuda on cpu
    /**
    BoxData<T, C, D, E>& castsrc = const_cast<BoxData<T,C,D,E> &>(a_src);
    Bx & castbox = const_cast<Bx &>(a_box);
    Stencil<T>* castthis = const_cast<Stencil<T>* >(this);
    castthis->fauxCudaApply(castsrc, a_dest, castbox, a_initToZero, a_scale);
    **/
    //end debugging

/**/
    if (m_coefs.size() == 0){return;}
    // Check input domain Boxes for proper nesting
    //TODO: Checks on inputs
    PROTO_ASSERT(!a_src.isAlias(a_dest),
    "Stencil::apply(...) invalid. a_src and a_dest alias the same data.");
   
    //CH_TIMERS("Stencil::apply");
    // Initialize a_dest to zero within a_box if flag is set
    for (int ii = 0; ii < DIM; ii++)
    {
        PROTO_ASSERT((m_destShift[ii] >= 0) && (m_destShift[ii] < m_destRefratio[ii]), 
        "m_destShift[ii] must be non-negative and less than m_destRefratio[ii] for all ii in [0,DIM)");
    }
    Bx cross = a_box.flatten(0);
    int npencil = a_box.size(0);
    if (a_initToZero)
    {
        for (int ee = 0; ee < E; ee++)
        for (int dd = 0; dd < D; dd++)
        for (int cc = 0; cc < C; cc++)
        for (auto iter = cross.begin(); iter != cross.end(); ++iter)
        {
            Point pt = (*iter)*m_destRefratio + m_destShift;
            T* val = &a_dest(pt,cc,dd,ee);
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
                T* destVal = &a_dest(dpt,cc,dd,ee);
                const T* srcVal = &a_src(spt,cc,dd,ee) + offsets[jj];
                for (int ii = 0; ii < npencil; ii++)
                {
                    destVal[(ii*m_destRefratio[0])] += coefs[jj]*srcVal[(ii*m_srcRefratio[0])];
                }
            }
        }
    }
/**/
#endif //end non-cuda version
}

//=======================================================================================
// LazyStencil ||
//=============++ 

template <typename T, unsigned int C, unsigned char D, unsigned char E>
LazyStencil<T,C,D,E>::LazyStencil(Stencil<T> a_stencil,
                                  const BoxData<T,C,D,E>* a_src,
                                  Bx a_box,
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
        Bx b = m_box[ii];
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
                                                  Bx                      a_box,
                                                  T                       a_scale) const
{
    LazyStencil<T,C,D,E> ret;
    if (m_s.size() <= 0){return ret;}
    Bx bsrc = a_src.box();
    for (int ii = 0; ii < m_s.size(); ii++)
    {
        Bx b = Bx(bsrc.low() - m_s[ii].span().low(), bsrc.high() - m_s[ii].span().high());
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
    LazyStencil<T,C,D,E> ret;
    if (m_s.size() <= 0){return ret;}
    Bx bsrc = a_src.box();

    for (int ii = 0; ii < m_s.size(); ii++)
    {
        //compute range for building the output BoxData, if necessary
        ret.m_range = ret.m_range & m_s[ii].range(bsrc).low();
        ret.m_range = ret.m_range & m_s[ii].range(bsrc).high();

        //compute the largest possible computation box (different from range) 
        Bx b = Bx(bsrc.low() - m_s[ii].span().low(), bsrc.high() - m_s[ii].span().high());
        b = b.taperCoarsen(m_s[ii].srcRatio());
        ret.m_stencil.push_back(m_s[ii]);
        ret.m_src.push_back(const_cast<BoxData<T,C,D,E>*>(&a_src));
        ret.m_box.push_back(b);
        ret.m_scale.push_back(a_scale);
    }
    return ret;
}

template<class T>
InterpStencil<T> InterpStencil<T>::PiecewiseConstant(Point a_ratio)
{
    InterpStencil<T> ret(a_ratio);
    Bx K = Bx(a_ratio - Point::Ones());
    for (auto iter = K.begin(); iter != K.end(); ++iter)
    {
        Stencil<T>& S = ret(*iter);
        S = 1.0*Shift::Zeros();
        S.destRatio() = a_ratio;
        S.destShift() = (*iter);
    }
    return ret;
}

template<class T>
InterpStencil<T> InterpStencil<T>::PiecewiseLinear(Point a_ratio)
{
    InterpStencil<T> ret(a_ratio);
    Bx K = Bx(a_ratio - Point::Ones());
    Bx U = Bx::Cube(2);
    for (auto kiter = K.begin(); kiter != K.end(); ++kiter)
    {
        Point p = *kiter;
        Point q = a_ratio - p;
        Stencil<double>& S = ret(*kiter);
        for (auto uiter = U.begin(); uiter != U.end(); ++ uiter)
        {
            Point e = *uiter;
            double coef = 1.0;
            for (int ii = 0; ii < DIM; ii++)
            {
                if (e[ii] == 0)
                {
                    coef *= q[ii]/((double)a_ratio[ii]);
                } else {
                    coef *= p[ii]/((double)a_ratio[ii]);
                }
            }
            if (coef != 0)
            {
                S += coef*Shift(e);
            }
        }
        S.destRatio() = a_ratio;
        S.destShift() = p;
    }
    return ret;    
}

// helper function for Quadratic InterpStencil
inline double qcoef(int a_ni, int a_r, int a_side)
{
    PROTO_ASSERT(std::abs(a_side) <= 1,
    "qcoef(ni, r, side) invalide for side = %i. Value must be in [-1,1]", a_side);
    double n = a_ni; double r = a_r;
    double coef;
    if (a_side == -1)
    {
        coef = 1.0 + 3.0*n + 3.0*n*n - 3.0*r - 6.0*n*r + 2.0*r*r;
    } 
    else if (a_side == 0)
    {
        coef = -2.0 - 6.0*n - 6.0*n*n + 3.0*r + 6.0*n*r + 5.0*r*r;
    }
    else if (a_side == 1)
    {
        coef = -(-1.0 - 3.0*n - 3.0*n*n + r*r);
    }
    coef /= (6.0*r*r);
    return coef;
}

template<class T>
InterpStencil<T> InterpStencil<T>::Quadratic(int a_ratio)
{
    //can make it anisotropic later...
    InterpStencil<T> ret(Point::Ones(a_ratio));
    Bx K = Bx(Point::Ones(a_ratio-1));
    Bx U = Bx::Cube(3).shift(Point::Ones(-1));

    double r = a_ratio;
     
    for (auto kiter = K.begin(); kiter != K.end(); ++kiter)
    {
        Point n = *kiter;

        Stencil<double>& S = ret(*kiter);
        for (auto iter = U.begin(); iter != U.end(); ++iter)
        {
            T coef = 1;
            for (int ii = 0; ii < DIM; ii++)
            {
                coef *= (T)qcoef(n[ii],r,(*iter)[ii]);
            }
            S += coef*Shift(*iter);
        }
        S.destRatio() = Point::Ones(a_ratio);
        S.destShift() = n;
    }
    return ret;    
}
template <typename T>
template<unsigned int C, unsigned char D, unsigned char E>   
void 
Stencil<T>::
getIntegerOffsets(std::vector<int>      &  a_srcOffsets,
                  std::vector<int>      &  a_dstOffsets,
                  int                   &  a_bxpts,
                  Bx                    &  a_bx,
                  BoxData<T,C,D,E>      &  a_src,
                  BoxData<T,C,D,E>      &  a_dst) const
{
  a_bxpts = a_bx.size();
  int stnpts = m_coefs.size();
  int dstSize = a_bxpts;
  int srcSize = a_bxpts*stnpts;
  a_srcOffsets.resize(srcSize);
  a_dstOffsets.resize(dstSize);

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

  for(int ibox = 0; ibox < a_bx.size(); ibox++)
  {
    Point boxpt = a_bx[ibox];
    Point dstpt = m_destRefratio*boxpt;
    dstpt += m_destShift;
    a_dstOffsets[ibox] = a_dst.index(dstpt);
    for(int icoef = 0; icoef < m_coefs.size(); icoef++)
    {
      int srcindex = icoef*a_bxpts + ibox;
      Point srcPt = boxpt*m_srcRefratio;
      a_srcOffsets[srcindex] = a_src.index(srcPt) + offsets[icoef];

    }
  }
}

template <typename T>
void 
fauxStencilIndexer(int a_stride, int a_blocks, int  a_begin, int a_end, int a_n, 
                   int a_srcVarOffset, int a_dstVarOffset, bool a_initToZero,
                   T* a_src,  T* a_dst, std::vector<T>& a_coeff,
                   std::vector<int>& a_srcOffset, std::vector<int>& a_dstOffset)
{
  for(int iblock = 0; iblock < a_blocks; iblock++)
  {
    for(int istride = 0; istride < a_stride; istride++)
    {
      int blockDim = a_stride;
      int blockIdx = iblock;
      int threadId = istride;
      int idx1 = iblock*a_stride +  istride;
      int idx = a_begin + threadId+blockIdx*blockDim;
      int ihuh = 0;
      if(idx1 != idx)
      {
        ihuh = 1;
      }
      int n = a_n;
      int bxpts = a_end;
      if((idx >= a_begin) && (idx< a_end))
      {
        if(a_initToZero)
        {
          a_dst[a_dstOffset[idx] + a_dstVarOffset] = 0.;
        }
        for(int isten = 0; isten < n; isten++)
        {
          int isx = isten*bxpts + idx;
          int dstoff = a_dstOffset[idx] + a_dstVarOffset;
          int srcoff = a_srcOffset[isx] + a_srcVarOffset;
          a_dst[dstoff] += a_src[srcoff]*a_coeff[isten];
        }
      }
    }
  }
}

// Apply Stencil
template <typename T>
template <unsigned int C, unsigned char D, unsigned char E>
void Stencil<T>::fauxCudaApply( BoxData<T,C,D,E>&  a_src,
                                BoxData<T,C,D,E>&  a_dst,
                                Bx&                a_box,
                                bool               a_initToZero,
                                T                  a_scale) 
{
  // Check input domain Boxes for proper nesting
  //TODO: Checks on inputs
  PROTO_ASSERT(!a_src.isAlias(a_dst),
               "Stencil::cudaApply(...) invalid. a_src and a_dst alias the same data.");
   
  //CH_TIMERS("Stencil::apply");
  // Initialize a_dst to zero within a_box if flag is set
  for (int ii = 0; ii < DIM; ii++)
  {
    PROTO_ASSERT((m_destShift[ii] >= 0) && (m_destShift[ii] < m_destRefratio[ii]), 
                 "m_destShift[ii] must be non-negative and less than m_destRefratio[ii] for all ii in [0,DIM)");
  }

  //make a copy of m_coefs and multiply by a_scale
  std::vector<T> coefs = m_coefs;
  if (a_scale != 1)
  {
    for (int ii = 0; ii < this->size(); ii++)
    {
      coefs[ii] *= a_scale;
    }
  }

  //m_coeff, m_offsets-lives on host
  //but m_offsets is in points and we need to make that 
  std::vector<T> & h_coeff  = coefs;
  std::vector<int>   h_offsetSrc, h_offsetDst;
  int boxpts;

  //first compute integer offsets
  getIntegerOffsets(h_offsetSrc, h_offsetDst, boxpts, a_box, a_src, a_dst);

  // apply the stencil
  for (int ee = 0; ee < E; ee++)
  {    
    for (int dd = 0; dd < D; dd++)
    {    
      for (int cc = 0; cc < C; cc++)
      {    
        Point pt = a_box.low();
        int varoffsetSrc = a_src.index(pt, cc, dd, ee) - a_src.index(pt, 0, 0, 0);
        int varoffsetDst = a_dst.index(pt, cc, dd, ee) - a_dst.index(pt, 0, 0, 0);
//        cout << "cc dd ee varoffsetsrc varoffsetdst " << cc << " " << dd << " " << ee << " " << varoffsetSrc << " " << varoffsetDst << endl;
        int begin = 0;
        int end= boxpts;
        int n = size();
        T* h_srcPtr = a_src.hostRawPtr();
        T* h_dstPtr = a_dst.hostRawPtr();
        int stride = a_box.size(0);
        int blocks = a_box.flatten(0).size();

        fauxStencilIndexer(stride, blocks, begin, end, n, varoffsetSrc, varoffsetDst, 
                           a_initToZero, 
                           h_srcPtr, h_dstPtr, h_coeff,
                           h_offsetSrc, h_offsetDst);
      }
    }
  }
}
//======================================================================
//                            CUDA STUFF                                
//======================================================================
#ifdef PROTO_CUDA

/**/
template <typename T>
__global__
void 
stencilIndexer(int  a_begin, int a_end, int a_n, 
               int a_srcVarOffset, int a_dstVarOffset, bool a_initToZero,
               T* a_src,  T* a_dst, T* a_coeff,
               int* a_srcOffset, int* a_dstOffset)
{
  int idx = a_begin + threadIdx.x+blockIdx.x*blockDim.x;
  int n = a_n;
  int bxpts = a_end;
//  printf("inside indexer n = %d, begin = %d, end = %d, idx = %d\n", n, a_begin, a_end, idx);
  if((idx >= a_begin) && (idx< a_end))
  {
//    printf(" idx = %d, n = %d \n", idx, n);
    if(a_initToZero)
    {
//      a_dst[a_dstOffset[idx] + a_dstVarOffset] = 45;
      a_dst[a_dstOffset[idx] + a_dstVarOffset] = 0.;
    }
    for(int isten = 0; isten < n; isten++)
    {
      int isx = isten*bxpts + idx;
      int dstoff = a_dstOffset[idx] + a_dstVarOffset;
      int srcoff = a_srcOffset[isx] + a_srcVarOffset;
//      printf("isten = %3d, idx = %3d,  isx = %3d, srcoff = %3d, dstoff = %3d, coeff = %f\n", isten, idx, isx, srcoff, dstoff, a_coeff[isten]);
//      a_src[a_srcOffset[isx]+a_srcVarOffset] = 23;
      a_dst[dstoff] += a_src[srcoff]*a_coeff[isten];
    }
  }
}


/**/
// Apply Stencil
template <typename T>
template <unsigned int C, unsigned char D, unsigned char E>
void Stencil<T>::cudaApply( BoxData<T,C,D,E>&  a_src,
                            BoxData<T,C,D,E>&  a_dst,
                            Bx&                a_box,
                            bool               a_initToZero,
                            T                  a_scale) 
{
  // Check input domain Boxes for proper nesting
  //TODO: Checks on inputs
  PROTO_ASSERT(!a_src.isAlias(a_dst),
               "Stencil::cudaApply(...) invalid. a_src and a_dst alias the same data.");
   
  //CH_TIMERS("Stencil::apply");
  // Initialize a_dst to zero within a_box if flag is set
  for (int ii = 0; ii < DIM; ii++)
  {
    PROTO_ASSERT((m_destShift[ii] >= 0) && (m_destShift[ii] < m_destRefratio[ii]), 
                 "m_destShift[ii] must be non-negative and less than m_destRefratio[ii] for all ii in [0,DIM)");
  }

  //make a copy of m_coefs and multiply by a_scale
  std::vector<T> coefs = m_coefs;
  if (a_scale != 1)
  {
    for (int ii = 0; ii < this->size(); ii++)
    {
      coefs[ii] *= a_scale;
    }
  }

  //m_coeff, m_offsets-lives on host
  //but m_offsets is in points and we need to make that 
  std::vector<T> & h_coeff  = coefs;
  std::vector<int>   h_offsetSrc, h_offsetDst;
  int boxpts;

  //first compute integer offsets
  getIntegerOffsets(h_offsetSrc, h_offsetDst, boxpts, a_box, a_src, a_dst);
  using std::cout;
  using std::endl;
  int srcSize = h_offsetSrc.size();
  int dstSize = h_offsetDst.size();
  int coeSize = h_coeff.size();

//  cout << "srcSize = " << srcSize << ", dstSize = " << dstSize << ", coeSize = " << coeSize << endl;
//  cout << "srcoffset  = ";
//  for(int isrc = 0; isrc < h_offsetSrc.size(); isrc++)
//  {
//    cout << h_offsetSrc[isrc] << " ";
//  }
//  cout << endl;
//  cout << "dstoffset  = ";
//  for(int idst = 0; idst < h_offsetDst.size(); idst++)
//  {
//    cout << h_offsetDst[idst] << " ";
//  }
//  cout << endl;

  size_t srcOffsetSize  = srcSize*sizeof(int);
  size_t dstOffsetSize  = dstSize*sizeof(int);
  size_t coeffSize      = coeSize*sizeof(T);
  //lives on device
  T*   d_coeff     = nullptr;
  int* d_offsetDst = nullptr;
  int* d_offsetSrc = nullptr;

  //get offsets and coefficients onto the device
  cudaMalloc(&d_coeff    ,  coeffSize);
  cudaMalloc(&d_offsetDst, dstOffsetSize);
  cudaMalloc(&d_offsetSrc, srcOffsetSize);

  cudaMemcpy(d_coeff    ,  h_coeff.data()    , coeffSize    , cudaMemcpyHostToDevice);
  cudaMemcpy(d_offsetSrc,  h_offsetSrc.data(), srcOffsetSize, cudaMemcpyHostToDevice);
  cudaMemcpy(d_offsetDst,  h_offsetDst.data(), dstOffsetSize, cudaMemcpyHostToDevice);

    

  // apply the stencil
  for (int ee = 0; ee < E; ee++)
  {    
    for (int dd = 0; dd < D; dd++)
    {    
      for (int cc = 0; cc < C; cc++)
      {    
        Point pt = a_box.low();
        int varoffsetSrc = a_src.index(pt, cc, dd, ee) - a_src.index(pt, 0, 0, 0);
        int varoffsetDst = a_dst.index(pt, cc, dd, ee) - a_dst.index(pt, 0, 0, 0);
//        cout << "cc dd ee varoffsetsrc varoffsetdst " << cc << " " << dd << " " << ee << " " << varoffsetSrc << " " << varoffsetDst << endl;
        int begin = 0;
        int end= boxpts;
        int n = coeSize;
        T* d_srcPtr = a_src.getRawPtr();
        T* d_dstPtr = a_dst.getRawPtr();
        int stride = a_box.size(0);
        int blocks = a_box.flatten(0).size();

        stencilIndexer<<<stride, blocks>>>(begin, end, n, varoffsetSrc, varoffsetDst, 
                                           a_initToZero, 
                                           d_srcPtr, d_dstPtr, d_coeff,
                                           d_offsetSrc, d_offsetDst);
        cudaError err = cudaGetLastError();
        if (err != cudaSuccess)
        {
          fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                  __FILE__, __LINE__, cudaGetErrorString(err));
        }
      }
    }
  }
  cudaFree(d_coeff    );
  cudaFree(d_offsetDst);
  cudaFree(d_offsetSrc);
  cudaDeviceSynchronize();
}
/**/
#endif
