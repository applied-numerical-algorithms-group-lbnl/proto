DataIterator::DataIterator(const std::shared_ptr<BoxPartition> a_partition)
{
    m_partition = a_partition;
    m_current = DataIndex(a_partition, a_partition->procStartIndex(procID()));
    m_end = DataIndex(a_partition, a_partition->numBoxes());
}

DataIterator&
DataIterator::begin()
{
    m_current = DataIndex(m_partition, m_partition->procStartIndex(procID()));
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
    if (currentInt < m_partition->procEndIndex(procID())-1)
    {
        m_current = DataIndex(m_partition, currentInt+1);
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

bool
DataIterator::compatible(const BoxPartition& a_partition) const
{
    return m_partition->compatible(a_partition);
}

bool
DataIterator::compatible(const DataIndex& a_index) const
{
    return m_partition->compatible(*a_index.m_partition);
}

bool
DataIterator::compatible(const DataIterator& a_iter) const
{
    return m_partition->compatible(*a_iter.m_partition);
}

bool DataIndex::operator==(const DataIterator& a_iter) const
{
    return (a_iter == *this);
}

bool DataIndex::operator!=(const DataIterator& a_iter) const
{
    return (a_iter != *this);
}
