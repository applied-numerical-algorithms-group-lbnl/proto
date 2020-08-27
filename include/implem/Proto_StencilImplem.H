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
template <unsigned int C, unsigned char D, unsigned char E>
void Stencil<T>::apply(const BoxData<T,C,D,E>&  a_src,
                       BoxData<T,C,D,E>&        a_dest,
                       const Box&               a_box,
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
//for debugging purposes only
//    hostApplySlow(a_src, a_dest, a_box, a_initToZero, a_scale);
#ifdef PROTO_BRICK
  hostApply(a_src, a_dest, a_box, a_initToZero, a_scale);
#else
  hostApply(a_src, a_dest, a_box, a_initToZero, a_scale);
#endif
#endif
    

    PR_FLOPS(this->numFlops(a_box));
}

template <typename T>
template<unsigned int C, unsigned char D, unsigned char E>
LazyStencil<T,C,D,E>
Stencil<T>::operator()(const BoxData<T,C,D,E>& a_src, T a_scale) const
{
    Box b = this->indexRange(a_src.box());
    return LazyStencil<T,C,D,E>(this, &a_src, b, a_scale);
}
template <typename T>
template<unsigned int C, unsigned char D, unsigned char E>
LazyStencil <T,C,D,E>
Stencil<T>::operator()(const BoxData<T,C,D,E>&  a_src, Box a_box, T a_scale) const
{
#ifndef NDEBUG
    Box b = this->indexRange(a_src.box());
    PROTO_ASSERT(b.contains(a_box),
            "Stencil::operator() invalid. User provided a Box larger than the valid range of the input data.");
#endif
    return LazyStencil<T, C, D, E>(this, &a_src, a_box, a_scale);
}

template <typename T, unsigned int C, unsigned char D, unsigned char E>
BoxData<T,C,D,E>&
operator|=(BoxData<T,C,D,E>& a_dest, LazyStencil<T,C,D,E>&& a_op)
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

template <class T, unsigned int C, unsigned char D, unsigned char E>
BoxData<T,C,D,E>&
operator+=(BoxData<T,C,D,E>& a_dest, LazyStencil<T,C,D,E>&& a_op)
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

#ifdef PROTO_BRICK

#endif

#ifndef PROTO_CUDA

template <typename T>
template <unsigned int C, unsigned char D, unsigned char E>
void Stencil<T>::hostApply(const BoxData<T,C,D,E>&  a_src,
                           BoxData<T,C,D,E>&        a_dest,
                           const Box&               a_box,
                           bool                     a_initToZero,
                           T                        a_scale) const
{
  if (m_coefs.size() == 0) return;

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
#ifndef PROTO_BRICK
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
#else
  // Source dimensions must be the same as dest
  auto strides = point2brickStrides(a_src.padded_box().high() - a_src.padded_box().low() + Point::Ones());
  unsigned* grid = a_src.getGrid();
  Point padded_box = a_src.padded_box().low();
  /*
  * Brick is limited by the code generator process:
  *   * no input/output shall have strides
  *   * input & output must have the same adjacency list
  */
  if (a_src.padded_box() == a_dest.padded_box() && m_srcRefratio == Point::Ones() && m_destRefratio == Point::Ones() &&
      m_destShift == Point::Zeros() && host_optimized != nullptr) {
    Point lb = a_box.low() - padded_box;
    Point rb = a_box.high() - padded_box + Point::Ones();

    auto base_point_f = [&strides, &grid, this, &coefs, &a_initToZero](Proto_brick &src, Proto_brick &dest,
                                                                       const Point &p) {
#if DIM == 3
      int tk = p[2] / BSIZE; int rk = p[2] % BSIZE;
#else
      int tk = 0;
      int rk = 0;
#endif
      int tj = p[1] / BSIZE;
      int rj = p[1] % BSIZE;
      int ti = p[0] / BSIZE;
      int ri = p[0] % BSIZE;
      auto b = DIM == 2 ? grid[ti + tj * strides[0]] : grid[ti + (tj + tk * strides[1]) * strides[0]];
      T res = 0.0;

      for (int ii = 0; ii < this->size(); ++ii) {
        auto br = src[b];
#if DIM == 3
        auto vk = br[rk + m_offsets[ii][2]];
#else
        auto &vk = br;
#endif
        auto vi = vk[rj + m_offsets[ii][1]][ri + m_offsets[ii][0]];
        res += coefs[ii] * vi;
      }
#if DIM == 3
      dest[b][rk][rj][ri] = (a_initToZero?0:dest[b][rk][rj][ri]) + res;
#else
      dest[b][rj][ri] = (a_initToZero ? 0 : dest[b][rj][ri]) + res;
#endif
    };

    /*
     * Call optimized version for Brick implementation
     * Note that optimized brick code is associated with only one pair of src/dst
     */
    typedef void optimized_t(Proto_brick &a_src,
                             Proto_brick &a_dest,
                             int brick_idx,
                             bool initToZero,
                             std::vector <T> &coefs);
    for (int ee = 0; ee < E; ee++)
      for (int dd = 0; dd < D; dd++)
        for (int cc = 0; cc < C; cc++) {
          auto src = a_src.getBrick(cc, dd, ee);
          auto dest = a_dest.getBrick(cc, dd, ee);

          auto brick_f = [&src, &dest, &strides, &grid, this, &a_initToZero, &coefs](const Point &p) {
#if DIM == 3
            int tk = p[2] / BSIZE; int rk = p[2] % BSIZE;
#else
            int tk = 0;
            int rk = 0;
#endif
            int tj = p[1] / BSIZE;
            int rj = p[1] % BSIZE;
            int ti = p[0] / BSIZE;
            int ri = p[0] % BSIZE;
            auto b = DIM == 2 ? grid[ti + tj * strides[0]] : grid[ti + (tj + tk * strides[1]) * strides[0]];
            ((optimized_t *) host_optimized)(src, dest, b, a_initToZero, coefs);
          };

          auto point_f = [&src, &dest, &base_point_f](const Point &p) {
            base_point_f(src, dest, p);
          };

          iterate_brick(lb, rb, point_f, brick_f);
        }
    return;
  }

  auto dest_strides = point2brickStrides(a_dest.padded_box().high() - a_dest.padded_box().low() + Point::Ones());
  unsigned* dest_grid = a_dest.getGrid();
  Point dest_padded_box = a_dest.padded_box().low();
  auto base_point_f = [&strides, &grid, this, &coefs, &a_initToZero, &padded_box,
                       &dest_strides, &dest_grid, &dest_padded_box]
      (Proto_brick &src, Proto_brick &dest, const Point &p) {
    T res = 0.0;
    {
      // Resolve source point
      Point spt = p * m_srcRefratio - padded_box;
  #if DIM == 3
      int tk = spt[2] / BSIZE; int rk = spt[2] % BSIZE;
  #else
      int tk = 0;
      int rk = 0;
  #endif
      int tj = spt[1] / BSIZE;
      int rj = spt[1] % BSIZE;
      int ti = spt[0] / BSIZE;
      int ri = spt[0] % BSIZE;
      auto b = DIM == 2 ? grid[ti + tj * strides[0]] : grid[ti + (tj + tk * strides[1]) * strides[0]];
      auto br = src[b];

      for (int ii = 0; ii < this->size(); ++ii) {
  #if DIM == 3
        auto vk = br[rk + m_offsets[ii][2]];
  #else
        auto &vk = br;
  #endif
        auto vi = vk[rj + m_offsets[ii][1]][ri + m_offsets[ii][0]];
        res += coefs[ii] * vi;
      }
    }
    {
      // resolve destination point
      Point dpt = p * m_destRefratio + m_destShift - dest_padded_box;
#if DIM == 3
      int tk = dpt[2] / BSIZE; int rk = dpt[2] % BSIZE;
#else
      int tk = 0;
      int rk = 0;
#endif
      int tj = dpt[1] / BSIZE;
      int rj = dpt[1] % BSIZE;
      int ti = dpt[0] / BSIZE;
      int ri = dpt[0] % BSIZE;
      auto b = DIM == 2 ? dest_grid[ti + tj * dest_strides[0]] :
               dest_grid[ti + (tj + tk * dest_strides[1]) * dest_strides[0]];

#if DIM == 3
      dest[b][rk][rj][ri] = (a_initToZero ? 0 : dest[b][rk][rj][ri]) + res;
#else
      dest[b][rj][ri] = (a_initToZero ? 0 : dest[b][rj][ri]) + res;
#endif
    }
  };

  /*
   * Regular stencil implementation will apply each offset
   */
  for (int ee = 0; ee < E; ee++)
    for (int dd = 0; dd < D; dd++)
      for (int cc = 0; cc < C; cc++) {
        auto src = a_src.getBrick(cc, dd, ee);
        auto dest = a_dest.getBrick(cc, dd, ee);
        Point p = Point::Zeros();
#if DIM == 3
        for (p[2] = a_box.low()[2]; p[2] <= a_box.high()[2]; ++p[2])
#endif
        for (p[1] = a_box.low()[1]; p[1] <= a_box.high()[1]; ++p[1])
          for (p[0] = a_box.low()[0]; p[0] <= a_box.high()[0]; ++p[0])
            base_point_f(src, dest, p);
      }

#endif
}

template <typename T>
template <unsigned int C, unsigned char D, unsigned char E>
void Stencil<T>::hostApplySlow(const BoxData<T,C,D,E>&  a_src,
                               BoxData<T,C,D,E>&        a_dest,
                               const Box&               a_box,
                               bool                     a_initToZero,
                               T                        a_scale) const                      
{
/**/

  if (m_coefs.size() == 0){return;}
   
  //CH_TIMERS("Stencil::apply");
  // Initialize a_dest to zero within a_box if flag is set
  for (int ii = 0; ii < DIM; ii++)
  {
    PROTO_ASSERT((m_destShift[ii] >= 0) && (m_destShift[ii] < m_destRefratio[ii]), 
                 "m_destShift[ii] must be non-negative and less than m_destRefratio[ii] for all ii in [0,DIM)");
  }

  if (a_initToZero)
  {
    for (int ee = 0; ee < E; ee++)
      for (int dd = 0; dd < D; dd++)
        for (int cc = 0; cc < C; cc++)
          for (auto iter = a_box.begin(); iter != a_box.end(); ++iter)
          {
            Point pt = (*iter)*m_destRefratio + m_destShift;
            T* val = a_dest.data(pt,cc,dd,ee);
            *val = 0;
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
  T destvalsave = 0;
  T destvalincr = 0;
  for (int ee = 0; ee < E; ee++)
    for (int dd = 0; dd < D; dd++)
      for (int cc = 0; cc < C; cc++)
      {    
        for (auto iter = a_box.begin(); iter != a_box.end(); ++iter)
        {
          for (int jj = 0; jj < this->size(); jj++)
          {
            Point dpt = (*iter)*m_destRefratio + m_destShift;
//begin debug
#if DIM==2
            Point ivdeb(55,17);
            int ideb = 1;
            if(dpt == ivdeb)
            {
              ideb = 1;
            }
#endif
//end debug
            if (!a_dest.box().contains(dpt)){continue;} //can happen when destShift is non-trivial
            Point spt = (*iter)*m_srcRefratio;
            T* destVal = a_dest.data(dpt,cc,dd,ee);
            const T* srcVal = a_src.data(spt,cc,dd,ee) + offsets[jj];
            T incr = coefs[jj]*(*srcVal);
            *destVal += incr;
//begin debug
#if DIM==2
            if(ideb == 1)
            {
              destvalsave  = *destVal;
              destvalincr  = incr;
            }
#endif
          }
        }
      }
/**/
}
#endif
//=======================================================================================
// LazyStencil ||
//=============++ 

template <typename T, unsigned int C, unsigned char D, unsigned char E>
LazyStencil<T,C,D,E>::LazyStencil(const Stencil<T>* a_stencil,
                                  const BoxData<T,C,D,E>* a_src,
                                  Box a_box,
                                  T a_scale)
{
    m_range = a_stencil->range(a_src->box()) & a_box.refine(a_stencil->destRatio());
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
        ret.m_stencil.push_back(&(m_s[ii]));
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
        ret.m_stencil.push_back(&(m_s[ii]));
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
#define PR_MAX_COEFFS 27
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
//  __shared__ T      s_coeff[PR_MAX_COEFFS];
//  __shared__ Point s_offset[PR_MAX_COEFFS];
//
//
//  if((threadIdx.x < a_ncoef)&&(threadIdx.x < PR_MAX_COEFFS) )
//    {
//      s_coeff[threadIdx.x] =a_coeff[threadIdx.x];
//      s_offset[threadIdx.x]=a_offset[threadIdx.x];
//    }
//  __syncthreads();

  if((idx< a_end))
  {
    Point boxPt = compute_box_point(a_boxLo, a_boxHi, idx);
    Point dstPt = boxPt;
    dstPt *= a_dstRefRat;
    dstPt += a_dstShift;
    unsigned int dstOff = Box::sindex(dstPt, a_dstLo, a_dstHi);
    Point srcSt = boxPt*a_srcRefRat;
    T& d = a_dst[dstOff];
    if(a_initToZero)
    {
      d = 0.;
    }
    for(int ivec = 0; ivec < a_ncoef; ivec++)
    {
      const     T & coef   = a_coeff[ivec];
      const Point & offset = a_offset[ivec];
      Point srcPt = srcSt + offset;
      unsigned int srcOff = Box::sindex(srcPt, a_srcLo, a_srcHi);
      
      d += a_scale*coef*a_src[srcOff];
    }
    dstPt[0]+=dstOff + idx;
  }
}
/**/
// Apply Stencil
template <typename T>
template <unsigned int C, unsigned char D, unsigned char E>
void Stencil<T>::cudaApply( BoxData<T,C,D,E> &  a_src,
                            BoxData<T,C,D,E> &  a_dst,
                            Box              &  a_box,
                            bool                a_initToZero,
                            T                   a_scale) 
{
  // Check input domain Boxes for proper nesting
   
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
        //if this fails, just make pr_max_coeffs bigger
        PROTO_ASSERT(ncoef <= PR_MAX_COEFFS, "you set max coeffs too small"); 

        const thrust::device_ptr<T>       coefptr = d_coeff.data();
        const thrust::device_ptr<Point> offsetptr = d_offset.data();
        size_t shared_mem = 0;
        cudaStream_t curstream = DisjointBoxLayout::getCurrentStream();
        stencilIndexer<<<blocks, stride, shared_mem, curstream>>>
          (begin, end,  ncoef, coefptr, offsetptr,
           a_initToZero, d_srcPtr, d_dstPtr, 
           boxLo, boxHi, srcLo, srcHi, dstLo, dstHi, 
           m_destShift, m_srcRefratio, m_destRefratio, a_scale);
#ifdef THRUST_DEBUG
        cudaError err = cudaGetLastError();
        if (err != cudaSuccess)
        {
          fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                  __FILE__, __LINE__, cudaGetErrorString(err));
        }
#endif
      }
    }
  }
}

template<typename T>
struct SVAL
{
  int offset;
  T   value;
};


template <typename T, int C>
__global__
void stencilIndexer2(int a_end, T* dstPtr, const T* srcPtr, Point srcStride, Point dstStride,
                     Point a_srcLo, Point a_srcHi, Point a_dstLo, Point a_dstHi,
                     Point a_dst, bool a_initToZero,
                     int ssize, thrust::device_ptr<T> coeffs,
                     thrust::device_ptr<Point> offsets, int ksize, T a_scale)
{
  int idx =  threadIdx.x+blockIdx.x*blockDim.x;
  int idy =  threadIdx.y+blockIdx.y*blockDim.y;
  extern __shared__ SVAL<T> srcIndices[];
  SVAL<T>* localIndices = srcIndices;
  if(idx < a_end)
    {
 
      if(threadIdx.x < ssize && threadIdx.y==0)
        {
          srcIndices[threadIdx.x].offset = Box::soffset(offsets[threadIdx.x], a_srcLo, a_srcHi);
          srcIndices[threadIdx.x].value = coeffs[threadIdx.x]*a_scale;
        }
       __syncthreads();

      a_dst[0]+=idx;
      a_dst[1]+=idy;
      unsigned int s_off = Box::sindex(a_dst, a_srcLo, a_srcHi);
 
      unsigned int d_off = Box::sindex(a_dst, a_dstLo, a_dstHi);
      T* dest = dstPtr + d_off;
      const T* src = srcPtr + s_off;
      for(int c=0; c<C; c++)
        {
          // C component loop
 
            for(int k=0; k<ksize; k++,src+=srcStride[DIM-1], dest+=dstStride[DIM-1])
              {
                if(a_initToZero) *dest=0;
                for(int s = 0; s<ssize; s++)
                  {
                    *dest+=srcIndices[s].value*src[srcIndices[s].offset];
                  }
                 s_off+=d_off*a_scale+localIndices[0].offset+localIndices[1].value + *dest;
              }
        }
    }
}

template <typename T>
template <unsigned int C, unsigned char D, unsigned char E>
void Stencil<T>::cudaApply2( const BoxData<T,C,D,E> &  a_src,
                            BoxData<T,C,D,E> &  a_dst,
                            const Box              &  a_box,
                            bool                a_initToZero,
                            T                   a_scale) 
{
  PR_TIME("cudaApply2");
  
  dim3 tblock(32, 4, 1);
  if(!m_isClosed)
    {
      d_coeff  = m_coefs;
      d_offset = m_offsets;
      m_isClosed = true;
    }
  // Let's ignore shift and ref ratio for the start here.
  int nx = a_box.high()[0]-a_box.low()[0]+1;
  int ny = a_box.high()[1]-a_box.low()[1]+1;
  dim3 grid((nx+32-1)/32,(ny+4-1)/4, 1);
  Point srcStride;
  Point dstStride;
  srcStride[0]=a_src.box().size(0);
  dstStride[0]=a_dst.box().size(0);
  for(int i=1; i<DIM; i++)
    {
      //these values change for refinement ratio changes
      srcStride[i]=a_src.box().size(i)*srcStride[i-1];
      dstStride[i]=a_dst.box().size(i)*dstStride[i-1];
    }
  int k=1;
#if DIM==3
  k= a_box.size(2);
#endif
  for (int ee = 0; ee < E; ee++)
    {    
      for (int dd = 0; dd < D; dd++)
        {
          const T* srcPtr = a_src.dataPtr(0, dd, ee);
          T* dstPtr = a_dst.dataPtr(0, dd, ee);
          cudaStream_t curstream = DisjointBoxLayout::getCurrentStream();
          stencilIndexer2<T, C><<<grid,
            tblock,
            m_coefs.size()*sizeof(SVAL<T>),
            curstream>>>(a_box.high()[0]-a_box.low()[0]+1, dstPtr, srcPtr, srcStride, dstStride,a_src.box().low(),
                         a_src.box().high(), a_dst.box().low(), a_dst.box().high(),
                         a_box.low(), a_initToZero,
                         m_coefs.size(), d_coeff.data(), d_offset.data(),k, a_scale );
        }
    }
  //PR_FLOPS(this->numFlops(a_box));
}


template <typename T>
template <unsigned int C, unsigned char D, unsigned char E>
void Stencil<T>::cudaApplyStream( BoxData<T,C,D,E> &  a_src,
                                  BoxData<T,C,D,E> &  a_dst,
                                  Box              &  a_box,
                                  bool                a_initToZero,
                                  T                   a_scale,
                                  cudaStream_t     &  a_stream,
                                  unsigned long long int & a_flops)
{
  PR_TIME("cudaApplyStream");
  // Check input domain Boxes for proper nesting
   
  //CH_TIMERS("Stencil::apply");
  // Initialize a_dst to zero within a_box if flag is set
  for (int ii = 0; ii < DIM; ii++)
  {
    PROTO_ASSERT((m_destShift[ii] >= 0) && (m_destShift[ii] < m_destRefratio[ii]), 
                 "m_destShift[ii] must be non-negative and less than m_destRefratio[ii] for all ii in [0,DIM)");
  }
//  static bool printed = false;
//  if(!printed)
//  {
//    printf("got to applystream \n");
//    printed = true;
//  }
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
        //if this fails, just make pr_max_coeffs bigger
        PROTO_ASSERT(ncoef <= PR_MAX_COEFFS, "you set max coeffs too small"); 

        const thrust::device_ptr<T>       coefptr = d_coeff.data();
        const thrust::device_ptr<Point> offsetptr = d_offset.data();
        size_t shared_mem = 0;
        stencilIndexer<<<blocks, stride, shared_mem, a_stream>>>
          (begin, end,  ncoef, coefptr, offsetptr,
           a_initToZero, d_srcPtr, d_dstPtr, 
           boxLo, boxHi, srcLo, srcHi, dstLo, dstHi, 
           m_destShift, m_srcRefratio, m_destRefratio, a_scale);
#ifdef THRUST_DEBUG
        cudaError err = cudaGetLastError();
        if (err != cudaSuccess)
        {
          fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
                  __FILE__, __LINE__, cudaGetErrorString(err));
        }
#endif
      }
    }
  }

  unsigned long long int numpts = a_box.size();
  unsigned long long int numflops = 3*size(); //multiplies by scale and  by weight and add

  a_flops = numpts*numflops;
}

/**/
#endif
