template<MemType MEM>
MBMap<MEM>::MBMap(
        MBDisjointBoxLayout& a_layout,
        Array<Point, DIM+1>& a_ghost,
        int a_order = 4)
{
    define(a_layout, a_ghost, a_order);
}

template<MemType MEM>
void MBMap<MEM>::define(
        MBDisjointBoxLayout& a_layout,
        Array<Point, DIM+1>& a_ghost,
        int a_order = 4)
{
    m_order = a_order;
    m_ghost = a_ghost;
    m_X.define(a_layout, a_ghost);
    m_J.define(a_layout, a_ghost);

    // user defined initialization
    this->init();

    for (auto iter : a_layout)
    {
        auto block = a_layout.block(iter);
        Box b_i = a_layout[iter];
        auto& X_i = m_X[iter];
        auto& J_i = m_J[iter];
        Box domainBox = a_layout.domain().blockDomain(block).box();
        Array<double, DIM> h;
        for (int dir = 0; dir < DIM; dir++)
        {
            h[dir] = 1.0/domainBox.size(dir);
        }
        this->compute(X_i, J_i, b_i, h, block);
    }
}



template<MemType MEM>
const MBDisjointBoxLayout& MBMap<MEM>::layout() { return m_X.layout(); }

template<MemType MEM>
MBLevelBoxData<double, DIM, MEM, PR_NODE>& MBMap<MEM>::map()
{
    return m_X;
}
template<MemType MEM>
const MBLevelBoxData<double, DIM, MEM, PR_NODE>& MBMap<MEM>::map() const
{
    return m_X;
}
template<MemType MEM>
MBLevelBoxData<double, 1, MEM, PR_CELL>& MBMap<MEM>::jacobian()
{
    return m_J;
}
template<MemType MEM>
const MBLevelBoxData<double, 1, MEM, PR_CELL>& MBMap<MEM>::jacobian() const
{
    return m_J;
}

