MBDataPoint::MBDataPoint(MBIndex a_dstIndex, Point a_dstPoint, const MBDisjointBoxLayout& a_layout)
{
    PROTO_ASSERT(a_layout.defined(),
            "MBDataPoint::constructor | Error: Undefined layout");
    PROTO_ASSERT(a_dstIndex.compatible(a_layout.partition()),
            "MBDataPoint::constructor | Error: Invalid index / layout pair");

    m_dstIndex = a_dstIndex;
    m_dstPoint = a_dstPoint;
    m_dstBlock = a_layout.block(a_dstIndex);
    m_boundaryDir = Point::Zeros();

    m_srcIndex = m_dstIndex;
    m_srcPoint = m_dstPoint;
    m_srcBlock = m_dstBlock;
}
MBDataPoint::MBDataPoint(MBIndex a_dstIndex, Point a_dstPoint,
    MBIndex a_srcIndex, Point a_srcPoint,
    const MBDisjointBoxLayout& a_layout)
{
    PROTO_ASSERT(a_layout.defined(),
            "MBDataPoint::constructor | Error: Undefined layout");
    PROTO_ASSERT(a_dstIndex.compatible(a_layout.partition()),
            "MBDataPoint::constructor | Error: Invalid index / layout pair");
    PROTO_ASSERT(a_srcIndex.compatible(a_layout.partition()),
            "MBDataPoint::constructor | Error: Invalid index / layout pair");
    
    m_dstIndex = a_dstIndex;
    m_dstPoint = a_dstPoint;
    m_dstBlock = a_layout.block(a_dstIndex);

    m_srcIndex = a_srcIndex;
    m_srcPoint = a_srcPoint;
    m_srcBlock = a_layout.block(a_srcIndex);

    m_boundaryDir = a_layout.domain().graph().connectivity(m_dstBlock, m_srcBlock);
}

MBDataPoint::MBDataPoint(const MBDataPoint& a_rhs)
{
    m_dstIndex = a_rhs.m_dstIndex;
    m_dstPoint = a_rhs.m_dstPoint;
    m_dstBlock = a_rhs.m_dstBlock;

    m_srcIndex = a_rhs.m_srcIndex;
    m_srcPoint = a_rhs.m_srcPoint;
    m_srcBlock = a_rhs.m_srcBlock;

    m_boundaryDir = a_rhs.m_boundaryDir;
}

MBDataPoint MBDataPoint::shift(Point a_dstShift, const MBDisjointBoxLayout& a_layout) const
{
    MBDataPoint ret(*this);
    ret.m_dstPoint = m_dstPoint + a_dstShift;
    if (inBoundary())
    {
        auto R = a_layout.domain().graph().rotation(m_dstBlock, m_srcBlock);
        ret.m_srcPoint = m_srcPoint + R(a_dstShift);
    } else {
        ret.m_srcPoint = m_srcPoint + a_dstShift;
    }
    
    return ret;
}

bool MBDataPoint::operator<(const MBDataPoint &a_rhs) const
{
    if (m_dstIndex != a_rhs.m_dstIndex)
    {
        return m_dstIndex.global() < a_rhs.m_dstIndex.global();
    }
    if (m_dstPoint != a_rhs.m_dstPoint)
    {
        return m_dstPoint < a_rhs.m_dstPoint;
    }
    if (inBoundary())
    {
        if (m_srcIndex != a_rhs.m_srcIndex)
        {
            return m_srcIndex.global() < a_rhs.m_srcIndex.global();
        }
        if (m_srcPoint != a_rhs.m_srcPoint)
        {
            return m_srcPoint < a_rhs.m_srcPoint;
        }
    }
    return false;
}

bool MBDataPoint::operator==(const MBDataPoint &a_rhs) const
{
    if (m_dstIndex != a_rhs.m_dstIndex)
    {
        return false;
    }
    if (m_dstPoint != a_rhs.m_dstPoint)
    {
        return false;
    }
    if (inBoundary())
    {
        if (m_srcIndex != a_rhs.m_srcIndex)
        {
            return false;
        }
        if (m_srcPoint != a_rhs.m_srcPoint)
        {
            return false;
        }
    }
    return true;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
BoxData<T, C, MEM>& MBDataPoint::dstPatch(MBLevelBoxData<T, C, MEM, CTR>& a_data) const
{
    return a_data[dstIndex()];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
BoxData<T, C, MEM>& MBDataPoint::srcPatch(MBLevelBoxData<T, C, MEM, CTR>& a_data) const
{
    if (inBoundary())
    {
        Point patchBoundaryDir = a_data.layout().connectivity(m_dstIndex, m_srcIndex);
        auto& boundData = a_data.bound(dstIndex(), patchBoundaryDir, m_srcBlock);
        return (*boundData.localData);
    } else {
        return dstPatch(a_data);
    }
}
