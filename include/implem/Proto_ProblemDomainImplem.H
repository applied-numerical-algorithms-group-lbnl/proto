
ProblemDomain::ProblemDomain(const Box& a_bx, const array<bool,DIM> a_isPeriodic)
{
    PROTO_ASSERT(a_bx.low() == Point::Zeros(),
            "ProblemDomain constructor | Error: Box must have low corner at (0,...,0).");
    define(a_bx,a_isPeriodic);
}

ProblemDomain::ProblemDomain(const Box& a_bx, bool a_isPeriodic)
{
    std::array<bool, DIM> periodic;
    for (int ii = 0; ii < DIM; ii++){periodic[ii] = a_isPeriodic;}
    define(a_bx,periodic);
}

void ProblemDomain::define(const Box& a_bx,const array<bool,DIM> a_isPeriodic)
{
    m_isDefined = true;
    m_isPeriodic = a_isPeriodic;
    m_bx = a_bx;
    PROTO_ASSERT(m_bx.low() == Point::Zeros(),"Low corner of the problem domain not zero");
}

Box ProblemDomain::operator&(Box a_bx) const
{
    Point lowCorner = a_bx.low();
    Point highCorner = a_bx.high();
    for (int dir = 0; dir < DIM; dir++)
    {
        if (!m_isPeriodic[dir])
        {
            lowCorner[dir] = std::max(lowCorner[dir],m_bx.low()[dir]);
            highCorner[dir] =  std::min(highCorner[dir],m_bx.high()[dir]);
        }
    }
    return Box(lowCorner,highCorner);
}

bool ProblemDomain::image(Point& a_pt) const
{
    if (m_bx.empty()) { return false; }
    if (m_bx.contains(a_pt)) { return true; }

    Box domain = m_bx.shift(-m_bx.low());
    Point point = a_pt - m_bx.low();
    for (int ii = 0; ii < DIM; ii++)
    {
        if (point[ii] < domain.low()[ii] || point[ii] > domain.high()[ii])
        {
            if (!m_isPeriodic[ii]) { return false; }
            else {
                int a = point[ii]; int b = domain.size(ii);
                int pi = a % b;
                if (pi < 0) { pi += m_bx.size(ii); }
                point[ii] = pi;
            }
        }
    }
    a_pt = point + m_bx.low();
    return true;
}

Point ProblemDomain::shifted(const Point& a_pt) const
{
    PR_assert(  !( (*this) & Box(a_pt,a_pt)  ).empty() );
    return a_pt % (Point::Ones()*(m_bx.high() + 1));
}

ProblemDomain ProblemDomain::coarsen(Point a_refRatio) const
{
    PR_assert(this->coarsenable(a_refRatio));
    ProblemDomain retval;
    retval.m_bx = m_bx.coarsen(a_refRatio);
    retval.m_isPeriodic = m_isPeriodic;
    retval.m_isDefined = m_isDefined;
    PR_assert(m_isDefined);
    return retval;
}

ProblemDomain ProblemDomain::refine(Point a_refRatio) const
{
    ProblemDomain retval;
    retval.m_bx = m_bx.refine(a_refRatio);
    retval.m_isPeriodic = m_isPeriodic;
    retval.m_isDefined = m_isDefined;
    return retval;
}
