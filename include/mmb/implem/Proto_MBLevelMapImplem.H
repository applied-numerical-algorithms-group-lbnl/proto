/*
    template<MemType MEM>
MBLevelMap<MEM>::MBLevelMap(
        MBDisjointBoxLayout& a_layout,
        Array<Point, DIM+1>& a_ghost,
        int a_order)
{
    define(a_layout, a_ghost, a_order);
}
*/
template<MemType MEM>
void MBLevelMap<MEM>::define(
        MBDisjointBoxLayout& a_layout,
        Array<Point, DIM+1>& a_ghost,
        int a_order)
{
    m_order = a_order;
    m_ghost = a_ghost;
    m_c2c = Stencil<double>::CornersToCells(a_order);
    m_X.define(a_layout, a_ghost);
    m_J.define(a_layout, a_ghost);
    m_X.setVal(0);
    m_J.setVal(0);
    m_dx.resize(a_layout.numBlocks());

    for (int bi = 0; bi < a_layout.numBlocks(); bi++)
    {
        Box domainBox = a_layout.domain().blockDomain(bi).box();
        for (int dir = 0; dir < DIM; dir++)
        {
            m_dx[bi][dir] = 1.0/domainBox.size(dir);
        }
    }

    // user defined initialization
    init();

    for (auto iter : a_layout)
    {
        auto block = a_layout.block(iter);
        Box b_i = a_layout[iter];
        auto& X_i = m_X[iter];
        auto& J_i = m_J[iter];
        apply(X_i, J_i, block);
    }
    m_X.fillBoundaries();
    m_X.exchange();
    m_J.fillBoundaries();
    m_J.exchange();
}

template<MemType MEM>
void MBLevelMap<MEM>::apply(
        BoxData<double, DIM, MEM>& a_X,
        BoxData<double, 1, MEM>& a_J,
        unsigned int a_block)
{
    FluxBoxData<double, DIM, MEM> NT(a_J.box());
    apply(a_X, a_J, NT, a_block);
}

template<MemType MEM>
void MBLevelMap<MEM>::apply(
        BoxData<double, DIM, MEM>& a_X,
        BoxData<double, 1, MEM>& a_J,
        FluxBoxData<double, DIM, MEM>& a_NT,
        unsigned int a_block)
{
    apply(a_X, a_J, a_block);
}

template<MemType MEM>
void MBLevelMap<MEM>::doApply(
        BoxData<double, DIM, MEM>& a_X,
        BoxData<double, 1, MEM>& a_J,
        unsigned int a_computeBlock,
        unsigned int a_outBlock)
{
    if (a_computeBlock != a_outBlock)
    {
        const auto& domain = layout().domain();
        Box BX = domain.convert(a_X.box(), a_outBlock, a_computeBlock, PR_NODE);
        Box BJ = domain.convert(a_J.box(), a_outBlock, a_computeBlock, PR_CELL);
        BoxData<double, DIM, MEM> X_tmp(BX);
        BoxData<double, 1, MEM> J_tmp(BJ);
        apply(X_tmp, J_tmp, a_computeBlock);
        Point conn = domain.graph().connectivity(a_computeBlock, a_outBlock);
        auto R = domain.graph().rotation(a_computeBlock, conn, a_outBlock);
        X_tmp.copyTo(a_X, R);
        J_tmp.copyTo(a_J, R);
    } else {
        apply(a_X, a_J, a_computeBlock);
    }
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
/*
template<MemType MEM>
Array<double, DIM> MBLevelMap<MEM>::operator()(
        const Point& a_cell,
        unsigned int a_block,
        Centering a_ctr)
{
    Box domainBox(a_cell, a_cell);
    auto index = layout().find(a_cell, a_block);
    PROTO_ASSERT(index != layout().end(),
            "MBLevelMap::operator() | Error: Point not found in block %u" << a_block);
    switch (a_ctr)
    {
        case PR_CELL:
            domainBox = m_c2c.domain(domainBox).grow(1); // c2c + convolve
            BoxData<double, DIM, MEM> X0(domainBox);
            BoxData<double, 1, MEM> J0(domainBox); // unused
            this->apply(X0, J0, index);
            BoxData<double, DIM, MEM> Xavg = m_c2c(X0);
            auto Xctr = Operator::deconvolve(Xavg);
            PROTO_ASSERT(Xctr.box().contains(a_cell),
                    "MBLevelMap::operator() | Error: Interpolated coordinates don't contain the requested cell");
            return Xctr(a_cell);
        case PR_NODE:
            BoxData<double, DIM, MEM> X0(domainBox);
            BoxData<double, 1, MEM> J0(domainBox); // unused
            this->apply(X0, J0, index);
            PROTO_ASSERT(X0.box().contains(a_cell),
                    "MBLevelMap::operator() | Error: Interpolated coordinates don't contain the requested cell");
            return X0(a_cell);
        default:
            MayDay<void>::Abort("MBLevelMap::operator() | Error: Case not implemented.");
    }
}

template<MemType MEM>
BoxData<double, DIM,MEM> MBLevelMap<MEM>::operator()(
        const Box& a_domain,
        unsigned int a_computeBlock,
        unsigned int a_domainBlock,
        Centering a_ctr)
{
*/
