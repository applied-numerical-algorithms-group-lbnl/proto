DataIterator::DataIterator(const std::shared_ptr<DBLInternals> a_internals)
{
    m_internals = a_internals;
    m_current = DataIndex(a_internals, a_internals->startProc[procID()]);
    m_end = DataIndex(a_internals, a_internals->allBoxes.size());
}

DataIterator&
DataIterator::begin()
{
    m_current = DataIndex(m_internals, m_internals->startProc[procID()]);
    return *this;
}

DataIterator&
DataIterator::end()
{
    m_current = m_end;
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
    if (currentInt < m_internals->startProc[procID()+1]-1)
    {
        m_current = DataIndex(m_internals,currentInt+1);
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

/*
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
*/
bool
DataIterator::compatible(const DBLInternals& a_internals) const
{
    return m_internals->compatible(a_internals);
}

bool
DataIterator::compatible(const DataIndex& a_index) const
{
    return m_internals->compatible(*a_index.m_internals);
}

bool
DataIterator::compatible(const DataIterator& a_iter) const
{
    return m_internals->compatible(*a_iter.m_internals);
}

bool DataIndex::operator==(const DataIterator& a_iter) const
{
    return (a_iter == *this);
}

bool DataIndex::operator!=(const DataIterator& a_iter) const
{
    return (a_iter != *this);
}
