MBDataPoint::MBDataPoint(MBIndex a_index, Point a_point, const MBDisjointBoxLayout& a_layout)
{
    PROTO_ASSERT(a_layout.defined(),
            "MBDataPoint::constructor | Error: Undefined layout");
    PROTO_ASSERT(a_index.compatible(a_layout.partition()),
            "MBDataPoint::constructor | Error: Invalid index / layout pair");
    m_index = a_index.global();
    m_point = a_point;
    m_block = -1;  // this means the block is the same as implied by a_index
    m_boundaryDir = Point::Zeros();
    m_inBoundary = false;
    m_partition = a_layout.partitionPtr();
}
MBDataPoint::MBDataPoint(MBIndex a_index, Point a_point, const MBDisjointBoxLayout& a_layout,
        Point a_boundaryDir, BlockIndex a_block)
{
    PROTO_ASSERT(a_block < a_layout.numBlocks(),
            "MBDataPoint::constructor | Error: Invalid block index");
    PROTO_ASSERT(a_layout.defined(),
            "MBDataPoint::constructor | Error: Undefined layout");
    PROTO_ASSERT(a_index.compatible(a_layout.partition()),
            "MBDataPoint::constructor | Error: Invalid index / layout pair");
    PROTO_ASSERT(a_boundaryDir.dot(a_boundaryDir) <= DIM,
            "MBDataPoint::constructor | Error: Invalid boundary direction");
    m_index = a_index.global();
    m_block = a_block;
    m_boundaryDir = a_boundaryDir;
    m_point = a_point;
    m_inBoundary = true;
    m_partition = a_layout.partitionPtr();
}

MBDataPoint::MBDataPoint(const MBDataPoint& a_rhs)
{
    m_index = a_rhs.m_index;
    m_boundaryDir = a_rhs.m_boundaryDir;
    m_point = a_rhs.m_point;
    m_block = a_rhs.m_block;
    m_inBoundary = a_rhs.m_inBoundary;
    m_partition = a_rhs.m_partition;
    PROTO_ASSERT(m_block < m_partition->numBlocks(),
            "MBDataPoint::constructor | Error: Invalid block index");
    PROTO_ASSERT(m_boundaryDir.dot(m_boundaryDir) <= DIM,
            "MBDataPoint::constructor | Error: Invalid boundary direction");
}

MBDataPoint MBDataPoint::shift(Point a_shift) const
{
    MBDataPoint ret(*this);
    ret.m_point = m_point + a_shift;
    return ret;
}

bool MBDataPoint::operator<(const MBDataPoint &a_rhs) const
{
    if (m_index != a_rhs.m_index)
    {
        return m_index < a_rhs.m_index;
    }
    else if (m_point != a_rhs.m_point)
    {
        return m_point < a_rhs.m_point;
    }
    if (m_inBoundary)
    {
        if (m_block != a_rhs.m_block)
        {
            return m_block < a_rhs.m_block;
        }
        else if (m_boundaryDir != a_rhs.m_boundaryDir)
        {
            return m_boundaryDir < a_rhs.m_boundaryDir;
        } else {
            return false;
        }
    }
    return false;
}

bool MBDataPoint::operator==(const MBDataPoint &a_rhs) const
{
    if (m_index != a_rhs.m_index)
    {
        return false;
    }
    if (m_point != a_rhs.m_point)
    {
        return false;
    }
    if (m_inBoundary)
    {
        if (m_block != a_rhs.m_block)
        {
            return false;
        }
        if (m_boundaryDir != a_rhs.m_boundaryDir)
        {
            return false;
        }
    }
    return true;
}
BlockIndex MBDataPoint::srcBlock() const
{
    return m_block >= 0 ? m_block : m_partition->getBlock(m_index); 
}

BlockIndex MBDataPoint::dstBlock() const
{
    return m_partition->getBlock(m_index); 
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
BoxData<T, C, MEM>& MBDataPoint::patch(MBLevelBoxData<T, C, MEM, CTR>& a_data) const
{
    if (m_inBoundary)
    {
        auto& boundData = a_data.bound(index(), m_boundaryDir, m_block);
        return *boundData.localData;
    } else {
        return a_data[index()];
    }
}
template<typename T, unsigned int C, MemType MEM, Centering CTR>
BoxData<T, C, MEM>& MBDataPoint::srcPatch(MBLevelBoxData<T, C, MEM, CTR>& a_data) const
{
    if (m_inBoundary)
    {
        auto& boundData = a_data.bound(index(), m_boundaryDir, m_block);
        return *boundData.adjData;
    } else {
        return a_data[index()];
    }
}
