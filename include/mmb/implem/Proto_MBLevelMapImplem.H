template<MemType MEM>
MBLevelMap<MEM>::MBLevelMap(
        MBDisjointBoxLayout& a_layout,
        Array<Point, DIM+1>& a_ghost,
        int a_order)
{
    define(a_layout, a_ghost, a_order);
}

template<MemType MEM>
void MBLevelMap<MEM>::define(
        MBDisjointBoxLayout& a_layout,
        Array<Point, DIM+1>& a_ghost,
        int a_order)
{
    m_order = a_order;
    m_ghost = a_ghost;
    m_X.define(a_layout, a_ghost);
    m_J.define(a_layout, a_ghost);
    m_dx.resize(a_layout.numBlocks());

    for (int bi = 0; bi < a_layout.numBlocks(); bi++)
    {
        Box domainBox = a_layout.domain().blockDomain(bi).box();
        for (int dir = 0; dir < DIM; dir++)
        {
            m_dx[bi][dir] = 1.0/domainBox.size(dir);
        }
        std::cout << "MBLevelMap | block: " << bi << " | dx: " << m_dx[bi] << std::endl;
    }

    // user defined initialization
    init();

    for (auto iter : a_layout)
    {
        auto block = a_layout.block(iter);
        Box b_i = a_layout[iter];
        auto& X_i = m_X[iter];
        auto& J_i = m_J[iter];
        apply(X_i, J_i, iter);
    }
    m_X.fillBoundaries();
    m_X.exchange();
    m_J.fillBoundaries();
    m_J.exchange();
}

template<MemType MEM>
const MBDisjointBoxLayout& MBLevelMap<MEM>::layout() const { return m_X.layout(); }
    
template<MemType MEM>
const Array<double, DIM>& MBLevelMap<MEM>::dx(int a_blockID) const 
{
    return m_dx[a_blockID];
}

template<MemType MEM>
BoxData<double, DIM, MEM> MBLevelMap<MEM>::X(
        const Box& a_box,
        const Array<double, DIM>& a_dx) const
{
    BoxData<double, DIM, MEM> XOut(a_box);
    forallInPlace_p(
        [] PROTO_LAMBDA (Point& a_pt, Var<double, DIM, MEM>& a_x, Array<double, DIM> a_dx)
        {
            for (int dir = 0; dir < DIM; dir++)
            {
                a_x(dir) = a_pt[dir]*a_dx[dir];
            }
        }, XOut, a_dx);
    return XOut;
}

template<MemType MEM>
MBLevelBoxData<double, DIM, MEM, PR_NODE>& MBLevelMap<MEM>::map()
{
    return m_X;
}
template<MemType MEM>
const MBLevelBoxData<double, DIM, MEM, PR_NODE>& MBLevelMap<MEM>::map() const
{
    return m_X;
}
template<MemType MEM>
MBLevelBoxData<double, 1, MEM, PR_CELL>& MBLevelMap<MEM>::jacobian()
{
    return m_J;
}
template<MemType MEM>
const MBLevelBoxData<double, 1, MEM, PR_CELL>& MBLevelMap<MEM>::jacobian() const
{
    return m_J;
}

