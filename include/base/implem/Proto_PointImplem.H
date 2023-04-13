
ACCEL_DECORATION
Point::Point()
{
    for (int j = 0; j < DIM; j++)
    {
        m_tuple[j]=0;
    }
}
ACCEL_DECORATION
Point::Point(const int (&a_tuple)[DIM])
{
    for (int j = 0; j < DIM; j++)
    {
        m_tuple[j]=a_tuple[j];
    }
}

ACCEL_DECORATION
Point::Point(const Point& a_pt)
{
    *this = a_pt;
}

ACCEL_DECORATION
Point Point::Ones(int a_scale)
{
    int ones[DIM];
    for (int k = 0 ; k < DIM; k++)
    {
        ones[k] = a_scale;
    }
    Point p(ones);
    return p;
}

ACCEL_DECORATION
Point Point::Zeros()
{
    int zeros[DIM];
    for (int k = 0 ; k < DIM; k++)
    {
        zeros[k] = 0;
    }
    return Point(zeros);
}

ACCEL_DECORATION
Point Point::Basis(int idir, int a_scale)
{
    Point retval = Point::Zeros();
    retval[idir] = a_scale;
    return retval;
}

ACCEL_DECORATION
Point Point::Basis(int idir, const Side::LoHiSide a_side)
{
    if (a_side == Side::Lo)
    {
        return Point::Basis(idir, -1);
    } else if (a_side == Side::Hi)
    {
        return Point::Basis(idir, 1);
    } else {
        return Point::Zeros();
    }
}

ACCEL_DECORATION
const int& Point::operator[](unsigned char a_index) const
{
    PROTO_ASSERT((a_index < DIM),
            "Point::operator[](uchar a_index) invalid for a_index = %i.\
            a_index must be in [0,DIM = %i)",
            a_index, DIM);
    return m_tuple[a_index];
}

ACCEL_DECORATION
int& Point::operator[](unsigned char a_index)
{
    PROTO_ASSERT(((a_index < DIM)),
            "Point::operator[](uchar a_index) invalid for a_index = %i.\
            a_index must be in [0,DIM = %i)",
            a_index, DIM);
    return m_tuple[a_index];
}



ACCEL_DECORATION
Point& Point::operator=(const Point& a_rhs)
{
    for (int ii = 0; ii < DIM; ii++)
    {
        m_tuple[ii] = a_rhs[ii];
    }
    return *this;
}

ACCEL_DECORATION
bool Point::operator<(const Point& a_rhs) const
{
    for (int ii = DIM-1; ii >= 0; ii--)
    {
        if (m_tuple[ii] > a_rhs[ii]){return false;}
        else if (m_tuple[ii] < a_rhs[ii]){return true;}
        //else they are equal, check next slowest index
    }
    //Points are equal; return false
    return false;
}

ACCEL_DECORATION
bool Point::operator==(const Point& a_rhsPoint) const 
{
    bool ret = true;
    for (int j = 0; j < DIM; j++)
    {
        if (m_tuple[j] != a_rhsPoint[j]) ret=false;
    }
    return ret;
}

ACCEL_DECORATION
bool Point::operator!=(const Point& a_rhsPoint) const 
{
    return !(*this == a_rhsPoint);
}

ACCEL_DECORATION
Point Point::operator+(const Point& a_rhsPoint) const
{
    int tuple[DIM];
    for (int j = 0; j < DIM; j++)
    {
        tuple[j] = m_tuple[j]+a_rhsPoint.m_tuple[j]; //a_tuple
    }
    Point ret(tuple);
    return ret;
}

ACCEL_DECORATION
Point Point::operator-(const Point& a_rhsPoint) const
{
    int tuple[DIM];
    for (int j = 0; j < DIM; j++)
    {
        tuple[j] = m_tuple[j]-a_rhsPoint.m_tuple[j]; 
    }
    Point ret(tuple);
    return ret;
}

ACCEL_DECORATION
Point Point::operator*(const Point& a_rhsPoint) const
{
    int tuple[DIM];
    for (int j = 0; j < DIM; j++)
    {
        tuple[j] = m_tuple[j]*a_rhsPoint.m_tuple[j]; //a_tuple
    }
    Point ret(tuple);
    return ret;
} 
ACCEL_DECORATION
Point Point::operator/(const Point& a_rhsPoint) const
{
    #ifndef NDEBUG
    for (int ii = 0; ii < DIM; ii++)
    {
        PROTO_ASSERT(a_rhsPoint[ii] != 0, 
        "Point::operator/(Point a_pt) invalid: division by zero (a_pt[%i] = 0)",ii);
    }
    // Division rounds towards -infinity, consistent with the geometric meaning coarsening.
    #endif
    int tuple[DIM];
    for (int j = 0; j < DIM; j++)
    {
        if (m_tuple[j]*a_rhsPoint[j] >= 0)
        {
            tuple[j] = m_tuple[j]/a_rhsPoint[j];
        }
        else if (-m_tuple[j]%a_rhsPoint[j] == 0)
          {
            tuple[j] = -(-m_tuple[j]/a_rhsPoint[j]);
          }
        else
          {
            tuple[j] = -(-m_tuple[j]/a_rhsPoint[j] + 1);
          }
    }
    Point ret(tuple);
    return ret;
}

ACCEL_DECORATION
Point Point::operator%(const Point& a_pt) const
{
    int arr[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
      if (m_tuple[ii] < 0)
        {  
          arr[ii] = (m_tuple[ii]+(-m_tuple[ii]/a_pt[ii]+1)*a_pt[ii]) % a_pt[ii];
        }
      else
        {
          arr[ii] = m_tuple[ii]%a_pt[ii];
        }
    }
    return Point(arr);
}

Point Point::operator+(int a_nref) const
{
    return (*this)+(Point::Ones(a_nref));
}    

Point Point::operator-(int a_nref) const
{
    return (*this)-(Point::Ones(a_nref));
}    

Point Point::operator*(int a_nref) const
{
    int tuple[DIM];
    for (int j = 0; j < DIM; j++)
    {
        tuple[j] = m_tuple[j]*a_nref;
    }
    Point ret(tuple);
    return ret;
}
    
Point Point::operator/(int a_nref) const
{
    PROTO_ASSERT(a_nref != 0,
    "Point::operator/(int a_nref) invalid for a_nref = %i. Division by zero.",
    a_nref);
    return (*this)/(Point::Ones()*a_nref);
}

Point Point::operator%(int a_nref) const
{
    PROTO_ASSERT(a_nref != 0,
    "Point::operator%%(int a_nref) invalid for a_nref = %i. Division by zero.",
    a_nref);
    return (*this)%(Point::Ones()*a_nref);
}

ACCEL_DECORATION
void Point::operator+=(const Point& a_rhsPoint)
{
    (*this) = (*this)+a_rhsPoint;
}

ACCEL_DECORATION
void Point::operator-=(const Point& a_rhsPoint)
{
    (*this) = (*this)-a_rhsPoint;
}

ACCEL_DECORATION
void Point::operator*=(const Point& a_rhsPoint)
{
    (*this) = (*this)*a_rhsPoint;
}

ACCEL_DECORATION
void Point::operator/=(const Point& a_rhsPoint)
{
    (*this) = (*this)/a_rhsPoint;
}

ACCEL_DECORATION
void Point::operator%=(const Point& a_pt)
{
    (*this) = (*this)%a_pt;
}

void Point::operator+=(int a_nref)
{
    (*this) = (*this)+a_nref;
}

void Point::operator-=(int a_nref)
{
    (*this) = (*this)-a_nref;
}
void Point::operator*=(int a_nref)
{
    (*this) = (*this)*a_nref;
}

void Point::operator/=(int a_nref)
{
    (*this) = (*this)/a_nref;
}

void Point::operator%=(int a_n)
{
    (*this) = (*this)%a_n;
}

int Point::dot(const Point& a_rhs) const
{
    return ((*this)*a_rhs).sum();
}

Point Point::coarsen(unsigned int a_refRatio) const
{
    PROTO_ASSERT(a_refRatio != 0, 
        "Error in Proto::Point::coarsen | refinement ratio must be non-zero.");
    Point retval = *this;
    for(int idir = 0; idir < DIM; idir++)
    {
        retval[idir] /= a_refRatio; 
    }
    return retval;
}

ACCEL_DECORATION
inline int Point::sum() const
{
    int s = 0;
    for (int ii = 0; ii < DIM; ii++) {
        s += m_tuple[ii];
    }
    return s;
}

int Point::min() const
{
    int val = m_tuple[0];
    for (int dir = 1; dir < DIM; dir++)
    {
        val = std::min(val, m_tuple[dir]);
    }
    return val;
}

int Point::max() const
{
    int val = m_tuple[0];
    for (int dir = 1; dir < DIM; dir++)
    {
        val = std::max(val, m_tuple[dir]);
    }
    return val;
}

Point Point::abs() const
{
    Point p;
    for (int d = 0; d < DIM; d++)
    {
        p[d] = std::abs((*this)[d]);
    }
    return p;
}

int Point::sumAbs() const
{
    int s = 0;
    for (int ii = 0; ii < DIM; ii++) {
        s += std::abs(m_tuple[ii]);
    }
    return s;
}

bool Point::isotropic() const
{
    for (int dir = 1; dir < DIM; dir++)
    {
        if (m_tuple[dir] != m_tuple[0]) {return false;}
    }
    return true;
}

int Point::codim() const
{
    int numZeros = 0;
    for (int d = 0; d < DIM; d++)
    {
        if ((*this)[d] != 0) { numZeros++; }
    }
    return numZeros;
}

std::array<Point, DIM> Point::parallelUnit() const
{
    std::array<Point, DIM> basis;
    for (int d = 0; d < DIM; d++)
    {
        if ((*this)[d] != 0)
        {
            basis[d] = Point::Basis(d);
        } else {
            basis[d] = Point::Zeros();
        }
    }
    return basis;
}

std::array<Point, DIM> Point::perpUnit() const
{
    std::array<Point, DIM> basis;
    for (int d = 0; d < DIM; d++)
    {
        if ((*this)[d] == 0)
        {
            basis[d] = Point::Basis(d);
        } else {
            basis[d] = Point::Zeros();
        }
    }
    return basis;
}

std::vector<Point> Point::split() const
{
    std::vector<Point> points;
    for (int ii = 0; ii < DIM; ii++)
    {
        int s = (*this)[ii];
        if (s != 0)
        {
            points.push_back(Point::Basis(ii, s));
        }
    }
    return points;
}

Point Point::dir() const
{
    Point d;
    for (int ii = 0; ii < DIM; ii++)
    {
        int di = (*this)[ii];
        d[ii] = di / std::abs(di);
    }
    return d;
}

size_t Point::linearSize() const
{
    size_t retval = DIM*sizeof(int);
    return retval;
}

void Point::linearIn(const char* a_buf)
{
    const int* intbuf = (int *) a_buf;
    for(int idir = 0; idir < DIM; idir++)
    {
        m_tuple[idir] = intbuf[idir];
    }
}


void Point::linearOut(char* a_buf) const
{
    int* intbuf = (int *) a_buf;
    for(int idir = 0; idir < DIM; idir++)
    {
        intbuf[idir] = m_tuple[idir];
    }
}

void Point::print() const
{
    using std::cout;
    using std::endl;
    cout << "(";
    for (int dir = 0; dir < DIM ; dir++)
    {
        cout << m_tuple[dir];
        if(dir < DIM-1) cout << ",";
    }
    cout <<   ")"  << endl;
}

std::ostream& operator <<(std::ostream& a_os, const Point& a_pt)
{
    a_os << "(" ;
    for (int dir = 0; dir < DIM; dir++) 
    {
        a_os << a_pt[dir];
        if (dir < DIM-1) 
        {
            a_os << ",";
        }
    }
    a_os << ")";
    return a_os;
}

std::istream& operator >>(std::istream& a_is,  Point& a_pt)
{
  char dum;
  a_is >> dum;
  for (int dir = 0; dir < DIM; dir++) 
  {
    a_is >> a_pt[dir];
    if (dir < DIM-1) 
    {
      a_is >> dum;
    }
  }
  a_is >> dum;
  return a_is;
}
