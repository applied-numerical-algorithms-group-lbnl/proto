
Point::Point()
{
    for (int j = 0; j < DIM; j++)
    {
        m_tuple[j]=0;
    }
}
Point::Point(const int (&a_tuple)[DIM])
{
    for (int j = 0; j < DIM; j++)
    {
        m_tuple[j]=a_tuple[j];
    }
}
Point::Point(const Point& a_pt)
{
    *this = a_pt;
}

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

Point Point::Zeros()
{
    int zeros[DIM];
    for (int k = 0 ; k < DIM; k++)
    {
        zeros[k] = 0;
    }
    return Point(zeros);
}

Point Point::Basis(int idir, int a_scale)
{
    Point retval = Point::Zeros();
    retval[idir] = a_scale;
    return retval;
}

const int& Point::operator[](unsigned char a_index) const
{
    PROTO_ASSERT((a_index < DIM),
            "Point::operator[](uchar a_index) invalid for a_index = %i.\
            a_index must be in [0,DIM = %i)",
            a_index, DIM);
    return m_tuple[a_index];
}

int& Point::operator[](unsigned char a_index)
{
    PROTO_ASSERT(((a_index < DIM)),
            "Point::operator[](uchar a_index) invalid for a_index = %i.\
            a_index must be in [0,DIM = %i)",
            a_index, DIM);
    return m_tuple[a_index];
}



Point& Point::operator=(const Point& a_rhs)
{
    for (int ii = 0; ii < DIM; ii++)
    {
        m_tuple[ii] = a_rhs[ii];
    }
    return *this;
}

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

bool Point::operator==(const Point& a_rhsPoint) const 
{
    bool ret = true;
    for (int j = 0; j < DIM; j++)
    {
        if (m_tuple[j] != a_rhsPoint[j]) ret=false;
    }
    return ret;
}

bool Point::operator!=(const Point& a_rhsPoint) const 
{
    return !(*this == a_rhsPoint);
}

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
Point Point::operator/(const Point& a_rhsPoint) const
{
    #ifndef NDEBUG
    for (int ii = 0; ii < DIM; ii++)
    {
        PROTO_ASSERT(a_rhsPoint[ii] != 0, 
        "Point::operator/(Point a_pt) invalid: division by zero (a_pt[%i] = 0)",ii);
    }
    #endif
    int tuple[DIM];
    for (int j = 0; j < DIM; j++)
    {
        if (m_tuple[j]*a_rhsPoint[j] >= 0)
        {
            tuple[j] = m_tuple[j]/a_rhsPoint[j];
        }
        else
        {
            tuple[j] = -(-m_tuple[j]/a_rhsPoint[j]);
        }
    }
    Point ret(tuple);
    return ret;
}

Point Point::operator%(const Point& a_pt) const
{
    int arr[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        arr[ii] = m_tuple[ii] % a_pt[ii];
    }
    return Point(arr);
}

Point Point::operator+(int a_nref) const
{
    return (*this)+(Point::Ones()*a_nref);
}    

Point Point::operator-(int a_nref) const
{
    return (*this)-(Point::Ones()*a_nref);
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

void Point::operator+=(const Point& a_rhsPoint)
{
    (*this) = (*this)+a_rhsPoint;
}

void Point::operator-=(const Point& a_rhsPoint)
{
    (*this) = (*this)-a_rhsPoint;
}

void Point::operator*=(const Point& a_rhsPoint)
{
    (*this) = (*this)*a_rhsPoint;
}

void Point::operator/=(const Point& a_rhsPoint)
{
    (*this) = (*this)/a_rhsPoint;
}

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

int Point::sum() const
{
    int s = 0;
    for (int ii = 0; ii < DIM; ii++) {
        s += m_tuple[ii];
    }
    return s;
}

int Point::absSum() const
{
    int s = 0;
    for (int ii = 0; ii < DIM; ii++) {
        s += std::abs(m_tuple[ii]);
    }
    return s;
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
