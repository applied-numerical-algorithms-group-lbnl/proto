
ProblemDomain::ProblemDomain(const Box& a_box, const Array<bool,DIM> a_isPeriodic)
{
    define(a_box,a_isPeriodic);
}

ProblemDomain::ProblemDomain(const Box& a_box, bool a_isPeriodic)
{
    Array<bool, DIM> periodic;
    for (int ii = 0; ii < DIM; ii++){periodic[ii] = a_isPeriodic;}
    define(a_box,periodic);
}

void ProblemDomain::define(const Box& a_box,const Array<bool,DIM> a_isPeriodic)
{
    m_isDefined = true;
    m_isPeriodic = a_isPeriodic;
    m_box = a_box;
    PROTO_ASSERT(a_box.low() == Point::Zeros(),
        "ProblemDomain constructor | Error: Box must have low corner at (0,...,0).");
}

Box ProblemDomain::operator&(Box a_box) const
{
    if (m_box.contains(a_box)) { return a_box; }
    else
    {
        Point lowCorner = a_box.low();
        Point highCorner = a_box.high();
        for (int dir = 0; dir < DIM; dir++)
        {
            if (!m_isPeriodic[dir])
            {
                lowCorner[dir] = std::max(lowCorner[dir],m_box.low()[dir]);
                highCorner[dir] =  std::min(highCorner[dir],m_box.high()[dir]);
            }
        }
        return Box(lowCorner,highCorner);
    }
}

bool ProblemDomain::contains(const Point& a_pt) const
{
    if (m_box.contains(a_pt)) { return true; }
    else {
        for (int dir = 0; dir < DIM; dir++)
        {
            if (!m_isPeriodic[dir])
            {
                if ((a_pt[dir] < m_box.low()[dir]) ||
                    (a_pt[dir] > m_box.high()[dir]))
                {
                    return false;
                }
            }
        }
    }
    return true;
}

Point ProblemDomain::image(const Point& a_pt) const
{
    PROTO_ASSERT(this->contains(a_pt),
        "ProblemDomain::image | Error: Point has no periodic image in this ProblemDomain.");
    return a_pt % (m_box.high() + Point::Ones());
}

Box ProblemDomain::image(const Box& a_box) const
{
    if (m_box.contains(a_box)) { return a_box; }
    else {
        Box truncated = (*this) & a_box;
        PROTO_ASSERT((m_box & truncated).empty(),
            "ProblemDomain::image | Error: Box input contains a periodic boundary.");
        Point lo = image(truncated.low());
        Point hi = image(truncated.high());
        return Box(lo, hi);
    }
}

bool ProblemDomain::coarsenable(Point a_boxSize) const
{
    return m_box.coarsenable(a_boxSize);
}

bool ProblemDomain::operator==(const ProblemDomain& a_domain) const
{
    return (m_box== a_domain.m_box) &&
           (m_isPeriodic==a_domain.m_isPeriodic) &&
           (m_isDefined==a_domain.m_isDefined);
}

bool ProblemDomain::operator!=(const ProblemDomain& a_domain) const
{
    return !(*this == a_domain);
}

ProblemDomain ProblemDomain::coarsen(Point a_refRatio) const
{
    PR_assert(this->coarsenable(a_refRatio));
    ProblemDomain retval;
    retval.m_box = m_box.coarsen(a_refRatio);
    retval.m_isPeriodic = m_isPeriodic;
    retval.m_isDefined = m_isDefined;
    PR_assert(m_isDefined);
    return retval;
}

ProblemDomain ProblemDomain::coarsen(int a_refRatio) const
{
    return coarsen(Point::Ones(a_refRatio));
}

ProblemDomain ProblemDomain::refine(Point a_refRatio) const
{
    ProblemDomain retval;
    retval.m_box = m_box.refine(a_refRatio);
    retval.m_isPeriodic = m_isPeriodic;
    retval.m_isDefined = m_isDefined;
    return retval;
}

ProblemDomain ProblemDomain::refine(int a_refRatio) const
{
    return refine(Point::Ones(a_refRatio));
}

std::ostream& operator<< (std::ostream& os, const ProblemDomain& a_pd)
{
    os << "ProblemDomain: Box = " << a_pd.box() <<
        " , periodicity = [" ;
    for (int d = 0; d < DIM; d++)
    {
        os<< a_pd.periodicity()[d] ;
        if (d < DIM-1) os << " , ";
    }
    os << "]";
    return os;
}
