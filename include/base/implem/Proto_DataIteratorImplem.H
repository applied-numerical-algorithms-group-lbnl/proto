template<typename P>
DataIterator<P>::DataIterator(const std::shared_ptr<P> a_partition)
{
    m_partition = a_partition;
    //m_current = DataIndex<P>(a_partition, a_partition->procStartIndex(procID()));
    m_current = DataIndex<P>(a_partition, a_partition->globalIndex(0,procID()));
    m_end = DataIndex<P>(a_partition, a_partition->numBoxes());
}

template<typename P>
DataIterator<P>&
DataIterator<P>::begin()
{
    //m_current = DataIndex<P>(m_partition, m_partition->procStartIndex(procID()));
    m_current = DataIndex<P>(m_partition, m_partition->globalIndex(0,procID()));
    return *this;
}

template<typename P>
DataIterator<P>&
DataIterator<P>::end()
{
    m_current = m_end;
    return *this;
}

template<typename P>
bool
DataIterator<P>::ok() const
{
    return (m_current != m_end);
}

template<typename P>
unsigned int
DataIterator<P>::localSize() const
{
    return m_partition->numBoxes(procID());
    //return (m_partition->procEndIndex(procID()) - m_partition->procStartIndex(procID()));
}

template<typename P>
DataIterator<P>&
DataIterator<P>::operator++()
{
    int localIndex = m_current.local();
    if (localIndex < localSize()-1)
    {
        m_current = DataIndex<P>(m_partition, m_partition->globalIndex(localIndex+1,procID()));
    } else {
        m_current = m_end;
    }
    /*
    int currentInt = m_current.global();
    if (currentInt < m_partition->procEndIndex(procID())-1)
    {
        m_current = DataIndex<P>(m_partition, currentInt+1);
    }
    else
    {
        m_current = m_end;
    }
    */
    return *this;
}

template<typename P>
const DataIndex<P>
DataIterator<P>::operator[](unsigned int a_localIndex) const
{
    PROTO_ASSERT(a_localIndex < localSize(),
            "DataIterator::operator[] | Error: index out of bounds");
    //return DataIndex<P>(m_partition, m_partition->procStartIndex(procID()) + a_localIndex);
    return DataIndex<P>(m_partition, m_partition->globalIndex(a_localIndex, procID()));
}

template<typename P>
const DataIndex<P>&
DataIterator<P>::operator*() const
{
    return m_current;
}

template<typename P>
bool
DataIterator<P>::compatible(const P& a_partition) const
{
    return m_partition->compatible(a_partition);
}

template<typename P>
bool
DataIterator<P>::compatible(const DataIndex<P>& a_index) const
{
    return m_partition->compatible(*a_index.m_partition);
}

template<typename P>
bool
DataIterator<P>::compatible(const DataIterator<P>& a_iter) const
{
    return m_partition->compatible(*a_iter.m_partition);
}

template<typename P>
bool DataIndex<P>::operator==(const DataIterator<P>& a_iter) const
{
    return (a_iter == *this);
}

template<typename P>
bool DataIndex<P>::operator!=(const DataIterator<P>& a_iter) const
{
    return (a_iter != *this);
}
