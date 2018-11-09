#include "PAssert.H"

//=================================================================================================
// Constructors ||
//==============++ 

Bx::Bx()
{
    m_low = Point::Zeros();
    m_high= Point::Ones()*(-1);
    m_size = 0;
}

Bx 
Bx::getFaceBox(int idir) const
{
  Point low  = m_low;
  Point high = m_high + Point::Basis(idir);
  return Bx(low, high);
}

Bx::Bx(const Bx& a_box)
{
    m_low = a_box.m_low;
    m_high = a_box.m_high;
    m_size = a_box.m_size;
}

Bx::Bx(const Point& a_lowCorner,const Point& a_highCorner)
{
    m_low = a_lowCorner;
    m_high = a_highCorner;
    recomputeSize();
}

Bx::Bx(const Point& a_pt)
{
    m_low = Point::Zeros();
    m_high = a_pt;
    recomputeSize();
}

void
Bx::recomputeSize()
{
    m_size = 1;
    for(int idir = 0; idir < DIM; idir++)
    {
        m_size *= size(idir);
    }
}

//=================================================================================================
// Static Methods ||
//================++ 

Bx Bx::Cube(int a_size)
{
    return Bx(Point::Zeros(),Point::Ones()*(a_size-1));
}

//=================================================================================================
// Accessors and Queries ||
//=======================++ 

std::size_t Bx::size(unsigned char a_dim) const
{
  PROTO_ASSERT((a_dim < DIM),
                "Bx::size(int a_dim) with a_dim = %i is invalid. a_dim should be in [0,DIM=%i).\n",
                 a_dim, DIM);
    if (m_high[a_dim] < m_low[a_dim]){return 0;}
    else {
        return (m_high[a_dim] - m_low[a_dim] + 1);
    }
}

bool Bx::contains(const Point& a_pt) const
{
    if (this->empty()){return false;}
    bool retval = true;
    for(int idir = 0; idir < DIM; idir++)
    {
        if(a_pt[idir] < m_low[idir])
        {
            retval = false;
        }
        if(a_pt[idir] > m_high[idir])
        {
            retval = false;
        }
    }
    return retval;
}
  
bool Bx::contains(const Bx& a_rhs) const
{
   return (this->contains(a_rhs.m_low) && this->contains(a_rhs.m_high));
}

bool Bx::onBoundary(const Point& a_p) const
{
  if (!(this->contains(a_p))){return false;}
  for (int ii = 0; ii < DIM; ii++)
  {
    if ((a_p[ii] == m_low[ii]) || (a_p[ii] == m_high[ii]))
    {
      return true;
    }
  }
  return false;
}

unsigned int Bx::index(const Point& a_pt) const
{
    PROTO_ASSERT(contains(a_pt),
    "Bx::index(Point a_pt) argument invalid. a_pt should be contained in *this Bx.\n");
    unsigned int factor = 1;
    unsigned int linIndex = a_pt[0] - m_low[0];
    for (unsigned char i = 1;i < DIM;i++)
    { 
        factor = factor*(m_high[i-1] - m_low[i-1]+1);
        linIndex = linIndex + (a_pt[i] - m_low[i])*factor;
    }
    return linIndex;
}

Point Bx::operator[](unsigned int a_index) const
{
  PROTO_ASSERT(((a_index < m_size)),
    "Bx::operator[](uint a_index) argument %i invalid. a_index should be in [0,%i).\n",
    a_index, m_size);
    
    int tuple[DIM];
    for (unsigned char i = 0;i < DIM; i++)
    {      
        int factor = (m_high[i] - m_low[i] + 1);
        int kred = a_index % factor;
        tuple[i] = kred + m_low[i];
        a_index = (a_index - kred)/factor;
    }
    return Point(tuple);
}

//=================================================================================================
// Operators ||
//===========++ 

Bx Bx::operator&(const Bx& a_rightBx) const
{
    int newLo[DIM];
    int newHi[DIM];
    for (int i = 0 ; i < DIM; i++)
    {
        newLo[i] = m_low[i];
        if (m_low[i] < a_rightBx.m_low[i])
        {
            newLo[i] = a_rightBx.m_low[i];
        }
        newHi[i] = m_high[i];
        if (m_high[i] > a_rightBx.m_high[i])
        {
            newHi[i] = a_rightBx.m_high[i];
        }
    }
    Point newLowCorner(newLo);
    Point newHighCorner(newHi); 
    for (int dir = 0; dir < DIM; dir++)
    {
        if (newHi[dir] < newLo[dir])
        {
            //if any of the dimensions are inside-out, return an empty (default) Bx
            Bx ret0;
            return ret0;
        }
    }
    Bx ret(newLowCorner,newHighCorner);
    return ret;
}

Bx Bx::operator&(Point a_pt) const
{
    if (contains(a_pt)){return *this;}
    else {
        Point L = m_low;
        Point H = m_high;
        for (int ii = 0; ii < DIM; ii++)
        {
            L[ii] = std::min(a_pt[ii], L[ii]);
            H[ii] = std::max(a_pt[ii], H[ii]);
        }
        return Bx(L,H);
    }
}

/*
void Bx::operator&=(const Bx& a_rightBx) 
{
    Bx retval = (*this) & a_rightBx;
    *this = retval;
}
*/

bool Bx::operator==(const Bx& a_rhsBx) const
{
    return ((m_low == a_rhsBx.m_low) 
       && (m_high == a_rhsBx.m_high));
}

Point Bx::operator%(const Point& a_pt) const
{
    return mod(a_pt);
}
//=================================================================================================
// Transformations ||
//=================++ 

Bx Bx::shift(int a_direction, int a_offset) const
{
  PROTO_ASSERT(((a_direction >= 0) && (a_direction < DIM)),
    "Bx::shift(int a_dir, int a_dist) with a_dir = %i is invalid. a_dir should be in [0,DIM=%i).\n",
    a_direction, DIM);
    Bx returnBx = Bx(*this);
    returnBx.m_low += Point::Basis(a_direction)*a_offset;
    returnBx.m_high += Point::Basis(a_direction)*a_offset;
    return returnBx;
}

Bx Bx::shift(const Point& a_pt) const
{
    Bx returnBx = Bx(*this);
    returnBx.m_low += a_pt;
    returnBx.m_high += a_pt;
    return returnBx;
}

Bx Bx::grow(const Point& a_offset) const
{
    Point lo = m_low;
    Point hi = m_high;
    lo -= a_offset;
    hi += a_offset;
    Bx returnBx(lo, hi);

    return returnBx;
}

Bx Bx::grow(int a_offset) const
{
    return grow(Point::Ones()*a_offset);
}

Bx Bx::coarsen(const Point& a_pt) const
{
    #ifndef NDEBUG
    for (int ii = 0; ii < DIM; ii++)
    {
        PROTO_ASSERT(a_pt[ii] > 0,
    "Bx::coarsen(Point a_pt) with a_pt[%i]=%i is invalid. All components of a_pt must be positive",
        ii,a_pt[ii]);
    }
    #endif
    Point lo = m_low;
    Point hi = m_high;
    lo /= a_pt;
    hi /= a_pt;
    Bx returnBx(lo, hi);

    return returnBx;
}

Bx Bx::coarsen(unsigned int a_nref) const
{
    PROTO_ASSERT(a_nref > 0,
    "Bx::coarsen(uint a_nref) with a_nref = %i is invalid. a_nref should be positive.", a_nref);

    return coarsen(Point::Ones()*a_nref);
}

Bx Bx::refine(const Point& a_pt) const
{
    #ifndef NDEBUG
    for (int ii = 0; ii < DIM; ii++)
    {
        PROTO_ASSERT(a_pt[ii] > 0,
    "Bx::refine(Point a_pt) with a_pt[%i]=%i is invalid. All components of a_pt must be positive.",
        ii,a_pt[ii]);
    }
    #endif
    Point lo = m_low;
    Point hi = m_high;
    lo *= a_pt;
    hi += Point::Ones();
    hi *= a_pt;
    hi -= Point::Ones();
    Bx returnBx(lo, hi);

    return returnBx;
}

Bx Bx::refine(unsigned int a_nref) const
{
    PROTO_ASSERT(a_nref > 0,
    "Bx::refine(uint a_nref) with a_nref = %i is invalid. a_nref should be positive.", a_nref);

    return refine(Point::Ones()*a_nref);
}

Bx Bx::edge(const Point& a_dir) const
{
  Point H = m_high; Point L = m_low;
  for (int ii = 0; ii < DIM; ii++)
  {
    if (a_dir[ii] > 0){L[ii] = H[ii];}
    if (a_dir[ii] < 0){H[ii] = L[ii];}
  }
  return Bx(L,H);
}

Bx Bx::flatten(const int a_dir, bool a_upper) const
{
  PROTO_ASSERT(((a_dir >= 0) && (a_dir < DIM)),
    "Bx::flatten(int a_dir) with a_dir = %i is invalid. a_dir should be in [0,DIM=%i).",
    a_dir, DIM);
    Point hi = m_high;
    Point lo = m_low;
    if (a_upper)
    {
        lo[a_dir] = hi[a_dir];
    } else {
        hi[a_dir] = lo[a_dir];
    }
    return Bx(lo,hi);
}

Bx Bx::extrude(int a_dir, int a_dist, bool a_upper) const
{
  PROTO_ASSERT(((a_dir >= 0) && (a_dir < DIM)),
    "Bx::flatten(int a_dir) with a_dir = %i is invalid. a_dir should be in [0,DIM=%i).",
    a_dir, DIM);
    Point hi = m_high;
    Point lo = m_low;
    if (a_upper)
    {
        hi[a_dir] += a_dist;
    } else {
        lo[a_dir] -= a_dist;
    }
    return Bx(lo,hi);
}

//=================================================================================================
// Utility ||
//=========++ 

BxIterator Bx::begin() const
{
    return BxIterator(*this,0);
}

BxIterator Bx::end() const
{
    return BxIterator(*this, this->size());
}

BxIterator Bx::rbegin() const
{
    return BxIterator(*this,this->size()-1);
}

BxIterator Bx::rend() const
{
    return BxIterator(*this, -1);
}

Point Bx::mod(const Point& a_pt) const
{
    PROTO_ASSERT(m_size > 0,
    "Bx::mod(Point a_pt) invalid for empty Bx.");
    Point p = a_pt - m_low;
    Point q = p % (m_high - m_low + Point::Ones());
    for (int ii = 0; ii < DIM; ii++)
    {
        if (q[ii] < 0)
        {
            q[ii] += this->size(ii);
        }
    }
    return q + m_low;
}

void Bx::print() const 
{
    std::cout<<*this<<std::endl;
}

