template<int NBLOCK>
MBProblemDomain XPoint<NBLOCK>::Domain(
    int a_domainSize,
    int a_thickness)
{
    PROTO_ASSERT(NBLOCK > 2, "XPoint::Domain | Error: Map is not valide for 2 or fewer blocks");
    MBProblemDomain domain(NBLOCK);
    auto CCW = CoordPermutation::ccw();
    for (int ii = 0; ii < NBLOCK; ii++)
    {
        domain.defineBoundary(ii, (ii+1) % NBLOCK, 0, Side::Hi, CCW);
    }
    Point blockSize = Point::Ones(a_domainSize);
    for (int dir = 2; dir < DIM; dir++)
    {
        blockSize[dir] = a_thickness;
    }
    for (int bi = 0; bi < NBLOCK; bi++)
    {
        domain.defineDomain(bi, blockSize);
    }
    domain.close();
    return domain;
}

template<int NBLOCK>
template<MemType MEM>
MBLevelMap<MBMap_XPointRigid<NBLOCK, MEM>, MEM>
XPoint<NBLOCK>::Map(
        const MBDisjointBoxLayout&  a_layout,
        Point   a_ghost)
{
    MBLevelMap<MBMap_XPointRigid<NBLOCK, MEM>, MEM> map;
    map.define(a_layout, a_ghost);
    return map;
}

template<int NBLOCK>
template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBLevelMap<MBMap_XPointRigid<NBLOCK, MEM>, MEM>
XPoint<NBLOCK>::Map(
        const MBLevelBoxData<T,C,MEM,CTR>& a_data)
{
    MBLevelMap<MBMap_XPointRigid<NBLOCK, MEM>, MEM> map;
    map.define(a_data.layout(), a_data.ghost());
    return map;
}
