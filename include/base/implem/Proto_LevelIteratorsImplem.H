
NeighborIterator::NeighborIterator(const DisjointBoxLayout& a_dbl, const Box& a_bx)
{
    PROTO_ASSERT(!a_bx.empty(),
        "NeighborIterator::Constructor | Error: Input Box cannot be empty.");
    m_layout = a_dbl;
    m_bx = a_bx;
    m_current = DataIndex(a_dbl.m_internals,0);
    m_end = DataIndex();
    m_dProbC = m_layout.domain().coarsen(m_layout.boxSize());    

    // Set up BoxIterator.
    Box bxC = a_bx.coarsen(m_layout.boxSize());
    Box bxPoints = m_dProbC & bxC;
    m_bxit = BoxIterator(bxPoints);
}

NeighborIterator&
NeighborIterator::begin()
{
    if (!m_bxit.ok()) { m_current = m_end; }
    m_bxit.begin();
    Point ptCoarseShifted = m_dProbC.image(*m_bxit);
    while ((m_layout.find(ptCoarseShifted) == m_layout.end()))
    {
        ++m_bxit;
        if (!m_bxit.ok()){ break; }
        ptCoarseShifted =m_dProbC.image(*m_bxit);
    }
    if (m_bxit.ok())
    {
        m_shift = (m_dProbC.image(*m_bxit) - *m_bxit)*m_layout.boxSize();
        m_current = m_layout.find(ptCoarseShifted);
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
    if (!m_bxit.ok())
    {
        m_current = m_end;
        return *this;
    }
    Point ptCoarseShifted = m_dProbC.image(*m_bxit);
    while ((m_layout.find(ptCoarseShifted) == m_layout.end()))
    {
        ++m_bxit;
        if (!m_bxit.ok()){break;}
        ptCoarseShifted = m_dProbC.image(*m_bxit);
    }
    if (m_bxit.ok())
    {
        m_shift = (m_dProbC.image(*m_bxit) - *m_bxit)*m_layout.boxSize();
        m_current = m_layout.find(ptCoarseShifted);
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
    return m_dProbC.image(*m_bxit);
}

Box
NeighborIterator::srcBox() const
{
    Box bx(*m_bxit,*m_bxit);
    return (bx.refine(m_layout.boxSize()) & m_bx).shift(m_shift);
}

Box
NeighborIterator::destBox() const
{
    Box bx(*m_bxit,*m_bxit);
    return (bx.refine(m_layout.boxSize()) & m_bx);
}

Point
NeighborIterator::shift()
{
    return m_shift;
}

DataIterator::DataIterator(const DisjointBoxLayout& a_dbl)
{
    //FIXME: Why is the iterator initialized this way? Shouldn't it return the same thing as begin()?
    m_layout = a_dbl;
    m_current = DataIndex(m_layout.m_internals,m_layout.size());
    m_end = m_current;
}

DataIterator&
DataIterator::begin()
{
    m_current = DataIndex(m_layout.m_internals, m_layout.m_internals->startProc[procID()]);
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
    int currentInt = m_current.global();
    if (currentInt <m_layout.m_internals->startProc[procID()+1]-1)
    {
        m_current = DataIndex(m_layout.m_internals,currentInt+1);
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
    return m_layout.box(m_current);
}

Point 
DataIterator::point() const
{
    return m_layout.point(m_current);
}

const DisjointBoxLayout&
DataIterator::layout() const
{
    return m_layout;
}

bool
DataIterator::compatible(const DataIndex& a_index) const
{
    return m_layout.compatible(a_index);
}

bool
DataIterator::compatible(const DataIterator& a_iter) const
{
    return m_layout.compatible(a_iter);
}

DataIterator
DisjointBoxLayout::begin() const
{
    DataIterator retval(*this);
    retval.begin();
    return retval;
}

// This lives here because DataIterator is an incomplete
// type at the time when DisjointBoxLayout is defined.
// One of the drawbacks of header-only -CLG
bool 
DisjointBoxLayout::compatible(const DataIterator& a_iter) const
{
    return compatible(a_iter.layout());
}

