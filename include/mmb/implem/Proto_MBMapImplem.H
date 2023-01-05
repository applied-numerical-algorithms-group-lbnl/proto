template<typename Func, MemType MEM>
MBMap<Func, MEM>::MBMap(
        const Func& a_func,
        MBDisjointBoxLayout& a_layout,
        Array<Point, DIM+1>& a_ghost,
        Point                a_boundGhost,
        unsigned int         a_order)
{
    define(a_func, a_layout, a_ghost, a_boundGhost, a_order);
}

template<typename Func, MemType MEM>
void MBMap<Func, MEM>::define(
        const Func& a_func,
        MBDisjointBoxLayout& a_layout,
        Array<Point, DIM+1>& a_ghost,
        Point                a_boundGhost,
        unsigned int         a_order)
{
    PR_TIME("MBMap::define");
    for (int d = 0; d < DIM+1; d++)
    {
        PROTO_ASSERT(a_ghost[d].min() >= 0,
                "MBMap::define | Error: Invalid ghost range");
    }
    PROTO_ASSERT(a_boundGhost.min() >= 0,
            "MBMap::define | Error: Invalid boundary ghost range");
    m_func = &a_func;
    m_order = a_order;
    auto ghost = a_ghost;
    for (int ii = 0; ii < DIM+1; ii++) { ghost[ii] = ghost[ii] + a_boundGhost; }
    m_X.define(a_layout, ghost, a_boundGhost);
    m_J.define(a_layout, ghost, a_boundGhost);
    
    const auto& domain = a_layout.domain();
    const auto& graph = domain.graph();
    for (auto iter : a_layout)
    {
        auto block = a_layout.block(iter);
        auto& patch = m_X[iter];
        Point blockSize = domain.blockDomain(block).box().sizes();
        
        BoxData<double, DIM> X0(patch.box().grow(a_ghost[0].max()+3));
        Centering centering = PR_NODE;
        forallInPlace_p(*m_func, X0, block, blockSize, centering);
        Array<BoxData<double, DIM>, DIM> NT;
        for (int d = 0; d < DIM; d++)
        {
            NT[d] = Operator::cofactor(X0, d);
        }
        auto J = Operator::jacobian(X0, NT);

        X0.copyTo(patch);
        J.copyTo(m_J[iter]);
    
        // manually fill the boundary data
        for (auto dir : Box::Kernel(1))
        {
            auto XBounds = m_X.bounds(iter, dir);
            auto JBounds = m_J.bounds(iter, dir);
            
            for (int bi = 0; bi < XBounds.size(); bi++) {
                auto XBound = XBounds[bi];
                auto JBound = JBounds[bi];
                auto adjBlock = a_layout.block(XBound.adjIndex);
                auto adjBlockSize = domain.blockDomain(adjBlock).box().sizes();
                Box adjDomain = XBound.adjData->box();
                adjDomain = adjDomain.extrude(Point::Ones(), 1); //hack related to Node centering

                // 3 ghost cells are needed to compute J with 0 ghost cells
                X0.define(adjDomain.grow(a_boundGhost.max()+3));
                forallInPlace_p(*m_func, X0, adjBlock, adjBlockSize, centering);
                for (int d = 0; d < DIM; d++)
                {
                    NT[d] = Operator::cofactor(X0, d);
                }
                auto J = Operator::jacobian(X0, NT);

                XBound.fill(X0);
                JBound.fill(J);
            }
        }
    }
}

template<typename Func, MemType MEM>
MBLevelBoxData<double, DIM, MEM, PR_NODE>&
MBMap<Func, MEM>::map()
{
    return m_X;
}

template<typename Func, MemType MEM>
const MBLevelBoxData<double, DIM, MEM, PR_NODE>&
MBMap<Func, MEM>::map() const
{
    return m_X;
}

template<typename Func, MemType MEM>
MBLevelBoxData<double, 1, MEM, PR_CELL>&
MBMap<Func, MEM>::jacobian()
{
    return m_J;
}

template<typename Func, MemType MEM>
const MBLevelBoxData<double, 1, MEM, PR_CELL>&
MBMap<Func, MEM>::jacobian() const
{
    return m_J;
}

template<typename Func, MemType MEM>
Array<double, DIM>
MBMap<Func, MEM>::operator()(const Array<double, DIM> a_X, unsigned int a_block) const
{
    PR_TIME("MBMap::operator()(Array, Block)");
    return m_func->map(a_X, a_block);
}

template<typename Func, MemType MEM>
Array<double, DIM>
MBMap<Func, MEM>::operator()(const MBDataPoint& a_pt, const Point a_dir) const
{
    return (*this)(a_pt.point, a_dir, m_X.layout().block(a_pt.index));
}

template<typename Func, MemType MEM>
Array<double, DIM>
MBMap<Func, MEM>::operator()(
        const Point& a_cell,
        const Point& a_dir,
        unsigned int a_block) const
{
    PR_TIME("MBMap::operator()(Cell, Dir, Block)");
    Array<double, DIM> X(a_cell);
    Array<double, DIM> dX(a_dir);
    dX += 1.0;
    dX *= 0.5;
    Array<double, DIM> blockSize = m_X.layout().domain().blockSize(a_block);
    X += dX;
    X /= blockSize;
    return m_func->map(X, a_block);
}

template<typename Func, MemType MEM>
BoxData<double, DIM>
MBMap<Func, MEM>::operator()(
        const Box& a_domain,
        unsigned int a_block,
        Centering a_ctr) const
{
    return (*this)(a_domain, a_block, a_block, a_ctr);
}

template<typename Func, MemType MEM>
BoxData<double, DIM>
MBMap<Func, MEM>::operator()(
        const Box& a_domain,
        unsigned int a_block,
        unsigned int a_refBlock,
        Centering a_ctr) const
{
    PR_TIME("MBMap::operator()(Box, Block, RefBlock, Centering)");
    // Compute x in a_block
    auto domain = m_X.layout().domain();
    Point blockSize = domain.blockSize(a_block);
    Box domainBox;
    switch (a_ctr)
    {
        case PR_CELL: domainBox = a_domain; break;
        case PR_NODE: domainBox = a_domain.extrude(Point::Ones(), 1); break;
        default:
            domainBox = a_domain.extrude(Point::Basis((int)a_ctr, 1), 1); break;
    }
    BoxData<double, DIM> xLocal(domainBox);
    forallInPlace_p(*m_func, xLocal, a_block, blockSize, a_ctr);
    
    // If a_refBlock is non-trivial, rotate the output
    if (a_block != a_refBlock)
    {
        Point conn = domain.graph().connectivity(a_block, a_refBlock);
        auto R = domain.graph().rotation(a_block, conn, a_refBlock);
        Box refBox = domain.convert(xLocal.box(), a_block, a_refBlock);
        xLocal.rotate(refBox, R);
    }
    return xLocal;
}

template<typename Func, MemType MEM>
bool MBMap<Func, MEM>::contains(
        const MBIndex& a_patch,
        const Point& a_cell,
        const Array<double, DIM>& a_X) const
{
    std::cout << "MBMap::contains is not yet implemented!!" << std::endl;
    return true;
}
