
NeighborIterator::NeighborIterator(const DisjointBoxLayout& a_dbl, const Box& a_bx)
{
    PROTO_ASSERT(!a_bx.empty(),
        "NeighborIterator::Constructor | Error: Input Box cannot be empty.");
    m_dbl = a_dbl;
    m_bx = a_bx;
    m_current = DataIndex(a_dbl.m_internals,0);
    m_end = DataIndex();
    m_dProbC = m_dbl.problemDomain().coarsen(m_dbl.boxSize());    

    // Set up BoxIterator.
    Box bxC = a_bx.coarsen(m_dbl.boxSize());
    Box bxPoints = m_dProbC & bxC;
    m_bxit = BoxIterator(bxPoints);
}

NeighborIterator&
NeighborIterator::begin()
{
    if (m_bxit.done()) { m_current = m_end; }
    m_bxit.begin();
    Point ptCoarseShifted = m_dProbC.shifted(*m_bxit);
    while ((m_dbl.find(ptCoarseShifted) == m_dbl.end()))
    {
        ++m_bxit;
        if (m_bxit.end()){ break; }
        ptCoarseShifted =m_dProbC.shifted(*m_bxit);
    }
    if (!m_bxit.end())
    {
        m_shift = (m_dProbC.shifted(*m_bxit) - *m_bxit)*m_dbl.boxSize();
        m_current = m_dbl.find(ptCoarseShifted);
    }
    else
    {
        m_current = m_end;
    }
    return *this;
}

NeighborIterator&
NeighborIterator::operator++()
{
    ++m_bxit;
    if (m_bxit.end())
    {
        m_current = m_end;
        return *this;
    }
    Point ptCoarseShifted = m_dProbC.shifted(*m_bxit);
    while ((m_dbl.find(ptCoarseShifted) == m_dbl.end()))
    {
        ++m_bxit;
        if (m_bxit.end()){break;}
        ptCoarseShifted = m_dProbC.shifted(*m_bxit);
    }
    if (!m_bxit.end())
    {
        m_shift = (m_dProbC.shifted(*m_bxit) - *m_bxit)*m_dbl.boxSize();
        m_current = m_dbl.find(ptCoarseShifted);
    }
    else
    {
        m_current = m_end;
    }
    return *this;
}

Point 
NeighborIterator::point() const
{
    return m_dProbC.shifted(*m_bxit);
}

Box
NeighborIterator::srcBox() const
{
    Box bx(*m_bxit,*m_bxit);
    return (bx.refine(m_dbl.boxSize()) & m_bx).shift(m_shift);
}

Box
NeighborIterator::destBox() const
{
    Box bx(*m_bxit,*m_bxit);
    return (bx.refine(m_dbl.boxSize()) & m_bx);
}

Point
NeighborIterator::shift()
{
    return m_shift;
}

DataIterator::DataIterator(const DisjointBoxLayout& a_dbl)
{
    //FIXME: Why is the iterator initialized this way? Shouldn't it return the same thing as begin()?
    m_dbl = a_dbl;
    m_current = DataIndex(m_dbl.m_internals,m_dbl.size());
    m_end = m_current;
}

DataIterator&
DataIterator::begin()
{
    m_current = DataIndex(m_dbl.m_internals, m_dbl.m_internals->startProc[procID()]);
    return *this;
}

bool
DataIterator::ok() const
{
    return (m_current != m_end);
}

DataIterator&
DataIterator::operator++()
{
    int currentInt = m_current.intIndex();
    if (currentInt <m_dbl.m_internals->startProc[procID()+1]-1)
    {
        m_current = DataIndex(m_dbl.m_internals,currentInt+1);
    }
    else
    {
        m_current = m_end;
    }
    return *this;
}

const DataIndex&
DataIterator::operator*() const
{
    return m_current;
}


Box 
DataIterator::box() const
{
    return m_dbl.box(m_current);
    /*
    Point pt = m_dbl[m_current];
    Point boxsize = m_dbl.boxSize();
    return Box(pt*boxsize,(pt + Point::Ones())*boxsize - Point::Ones());
    */
}

Point 
DataIterator::point() const
{
    return m_dbl.point(m_current);
    /*
    Point pt = m_dbl[m_current];
    return pt;
    */
}

DataIterator
DisjointBoxLayout::begin() const
{
    DataIterator retval(*this);
    retval.begin();
    return retval;
}

#if 0
FineFRIterator::FineFRIterator(
        const DataIndex& a_fineDit,
        const Point& a_refRatio)
{
}
    FRFineIterator&
FRFineIterator::begin()
{
}
    DataIndex&
FRFineIterator::end()
{
}
    FRFineIterator&
FRFineIterator::operator++()
{
}
    DataIndex&
FRFineIterator::operator()
{
}
    Box
FRFineIterator::box()
{
}
    Point
FRFineIterator::srcShift()
{
}
CoarseFRIterator::CoarseFRIterator(const DataIndex& a_coarseDit,
        const Point& a_refRatio)
{
}
    CoarseFRIterator&
CoarseFRIterator::begin()
{
}
    CoarseFRIterator&
CoarseFRIterator::end()
{
}
    CoarseFRIterator&
CoarseFRIterator::operator++()
{
}
    DataIndex
CoarseFRIterator::operator()
{
}
    Box
CoarseFRIterator::box()
{
}
    Point
CoarseFRIterator::srcShift()
{
}
    DataIndex
CoarseFRIterator::fromIndex()
{
}
#endif
