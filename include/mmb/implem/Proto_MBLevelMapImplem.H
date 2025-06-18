template<typename MAP, MemType MEM>
MBLevelMap<MAP,MEM>::MBLevelMap(
        const MBDisjointBoxLayout&  a_layout,
        const Array<Point, DIM+1>&  a_ghost,
        unsigned int                a_focalBlock)
{
    define(a_layout, a_ghost, a_focalBlock);
}

template<typename MAP, MemType MEM>
MBLevelMap<MAP,MEM>::MBLevelMap(
        const MBDisjointBoxLayout&  a_layout,
        Point                       a_ghost,
        unsigned int                a_focalBlock)
{
    define(a_layout, a_ghost, a_focalBlock);
}

template<typename MAP, MemType MEM>
void MBLevelMap<MAP,MEM>::define(
        const MBDisjointBoxLayout& a_layout,
        Point                      a_ghost,
        unsigned int               a_focalBlock)
{
    Array<Point, DIM+1> ghost;
    ghost.fill(a_ghost + Point::Ones(2)); //assumes 4th-5th order interpolation
    ghost[0] = a_ghost;
    define(a_layout, ghost, a_focalBlock);
}

template<typename MAP, MemType MEM>
void MBLevelMap<MAP,MEM>::define(
        const MBDisjointBoxLayout& a_layout,
        const Array<Point, DIM+1>& a_ghost,
        unsigned int               a_focalBlock)
{
    PR_TIME("MBLevelMap::define");
    m_layout = a_layout;
    m_ghost = a_ghost;
    m_focalBlock = a_focalBlock;
    m_ops.resize(a_layout.numBlocks());
    m_initialized = false;
    for (int bi = 0; bi < a_layout.numBlocks(); bi++)
    {
        Box domainBox = a_layout.domain().getBlock(bi).box();
        m_ops[bi] = std::make_shared<MAP>(dx(bi), bi, a_focalBlock);
    }
    m_c2c = Stencil<double>::CornersToCells(4);
    m_X = std::make_shared<MBLevelBoxData<double, DIM, MEM, PR_NODE>>(m_layout, m_ghost);
    m_J = std::make_shared<MBLevelBoxData<double, 1, MEM, PR_CELL>>(m_layout, m_ghost);
}

template<typename MAP, MemType MEM>
void MBLevelMap<MAP,MEM>::initialize() const
{
    m_X->setVal(0);
    m_J->setVal(0);
    
    for (auto iter : m_layout)
    {
        auto block = m_layout.block(iter);
        Box b_i = m_layout[iter];
        auto& X_i = (*m_X)[iter];
        auto& J_i = (*m_J)[iter];
        BoxData<double, 1, MEM> J0(X_i.box());
        apply(X_i, J0, block);
        J_i |= m_c2c(J0);
    }
    m_X->exchange();
    m_J->exchange();
    m_initialized = true;
}

template<typename MAP, MemType MEM>
MAP& MBLevelMap<MAP,MEM>::operator[](BlockIndex block)
{
    return *m_ops[block];
}

template<typename MAP, MemType MEM>
void MBLevelMap<MAP,MEM>::apply(
        BoxData<double, DIM, MEM>& a_X,
        BoxData<double, 1, MEM>& a_J,
        unsigned int a_block) const
{
    PR_TIME("MBLevelMap::apply");
    m_ops[a_block]->_apply(a_X, a_J);
}

template<typename MAP, MemType MEM>
void MBLevelMap<MAP,MEM>::apply(
        BoxData<double, DIM, MEM>& a_X,
        BoxData<double, 1, MEM>& a_J,
        FluxBoxData<double, DIM, MEM>& a_NT,
        unsigned int a_block) const
{
    PR_TIME("MBLevelMap::apply");
    m_ops[a_block]->_apply(a_X, a_J, a_NT);
}

template<typename MAP, MemType MEM>
void MBLevelMap<MAP,MEM>::doApply(
        BoxData<double, DIM, MEM>& a_X,
        BoxData<double, 1, MEM>& a_J,
        unsigned int a_computeBlock,
        unsigned int a_outBlock) const
{
    if (a_computeBlock != a_outBlock)
    {
        PR_TIME("MBLevelMap::doApply");
        const auto& domain = layout().domain();
        Box B = domain.convertBox(a_X.box(), a_outBlock, a_computeBlock, PR_NODE);
        BoxData<double, DIM, MEM> X_tmp(B);
        BoxData<double, 1, MEM> J_tmp(B);
        apply(X_tmp, J_tmp, a_computeBlock);
        Point conn = domain.graph().connectivity(a_computeBlock, a_outBlock);
        auto R = domain.graph().rotation(a_computeBlock, conn, a_outBlock);
        X_tmp.copyTo(a_X, R);
        J_tmp.copyTo(a_J, R);
    } else {
        apply(a_X, a_J, a_computeBlock);
    }
}

template<typename MAP, MemType MEM>
BoxData<double, DIM, MEM> MBLevelMap<MAP,MEM>::cellAveraged(
        const Box&   a_box,
        unsigned int a_computeBlock,
        unsigned int a_outBlock) const
{
    PR_TIME("MBLevelMap::cellAveraged(Box)");
    Box B0 = m_c2c.domain(a_box);
    BoxData<double, DIM, MEM> X0(B0);
    BoxData<double, 1, MEM> J0(B0); //unused
    doApply(X0, J0, a_computeBlock, a_outBlock);
    return m_c2c(X0);
}

template<typename MAP, MemType MEM>
BoxData<double, DIM, MEM> MBLevelMap<MAP,MEM>::cellCentered(
        const Box&   a_box,
        unsigned int a_computeBlock,
        unsigned int a_outBlock) const
{
    PR_TIME("MBLevelMap::cellCentered(Box)");
    auto X0 = cellAveraged(a_box.grow(1), a_computeBlock, a_outBlock);
    return Operator::deconvolve(X0);
}

template<typename MAP, MemType MEM>
Array<double, DIM> MBLevelMap<MAP,MEM>::cellAveraged(const MBDataPoint& a_point) const
{
    PR_TIME("MBLevelMap::cellAveraged(MBDataPoint)");
    auto computeBlock = a_point.srcBlock();
    auto domainBlock = a_point.dstBlock();
    Box B0(a_point.dstPoint(), a_point.dstPoint());
    auto X = cellAveraged(B0, computeBlock, domainBlock);
    return X.array(a_point.dstPoint());
}

template<typename MAP, MemType MEM>
Array<double, DIM> MBLevelMap<MAP,MEM>::cellCentered(const MBDataPoint& a_point) const
{
    PR_TIME("MBLevelMap::cellCentered(MBDataPoint)");
    auto computeBlock = a_point.srcBlock();
    auto domainBlock = a_point.dstBlock();
    Box B0(a_point.dstPoint(), a_point.dstPoint());
    auto X = cellCentered(B0, computeBlock, domainBlock);
    return X.array(a_point.dstPoint());
}

template<typename MAP, MemType MEM>
Array<double, DIM> MBLevelMap<MAP,MEM>::mappedCoords(
        const Point&          a_point,
        const BlockIndex&     a_block) const
{
    Array<double, DIM> offset;
    offset.fill(0);
    return mappedCoords(a_point, a_block, offset);
}

template<typename MAP, MemType MEM>
Array<double, DIM> MBLevelMap<MAP,MEM>::mappedCoords(
        const Point&          point,
        const BlockIndex&     block,
        const Array<double, DIM>&   offset) const
{
    auto dx = this->dx(block);
    Array<double, DIM> X = point;
    X += 0.5; //default is cell center
    X += offset;
    X *= dx;
    return X;
}

template<typename MAP, MemType MEM>
Array<double, DIM> MBLevelMap<MAP,MEM>::mappedCoords(
        const MBDataPoint&          a_point,
        const Array<double, DIM>&   a_offset) const
{
    return mappedCoords(a_point.dstPoint(), a_point.block(), a_offset);
}

template<typename MAP, MemType MEM>
Array<double, DIM> MBLevelMap<MAP,MEM>::mappedCoords(
        const MBDataPoint& a_point) const
{
    Array<double, DIM> offset;
    offset.fill(0);
    return mappedCoords(a_point, offset);
}

template<typename MAP, MemType MEM>
Array<double, DIM> MBLevelMap<MAP,MEM>::apply(
        const MBDataPoint&          a_point,
        const Array<double, DIM>&   a_offset) const
{
    auto xi = mappedCoords(a_point, a_offset);
    auto block = a_point.srcBlock();
    return op(block).apply(xi); 
}

template<typename MAP, MemType MEM>
Array<double, DIM> MBLevelMap<MAP,MEM>::apply(
        const MBDataPoint& a_point) const
{
    Array<double, DIM> offset;
    offset.fill(0);
    return apply(a_point, offset);
}


template<typename MAP, MemType MEM>
Array<double, DIM> MBLevelMap<MAP,MEM>::analyticCellCentered(const MBDataPoint& a_point) const
{
    return cellCentered(a_point);
    // TODO: We need to do this, but this isn't the right way. Turning this off 
    // until we get MBBoxOp online.
    /**
    auto computeBlock = a_point.srcBlock();
    Array<double,DIM> X;
    Array<double,DIM> thisDx = this->dx(a_point.srcBlock());
    Point pt = a_point.point();
    f_cubeSphereMapPoint<MEM>(pt,X,thisDx,CUBED_SPHERE_SHELL_R0,
                         CUBED_SPHERE_SHELL_R1,a_point.srcBlock(),
                         CUBED_SPHERE_SHELL_RADIAL_COORD);
    return X;
  }
  */
}

template<typename MAP, MemType MEM>
const MBDisjointBoxLayout& MBLevelMap<MAP,MEM>::layout() const { return m_layout; }
    
template<typename MAP, MemType MEM>
const Array<double, DIM> MBLevelMap<MAP,MEM>::dx(BlockIndex a_blockID) const 
{
    return layout().dx(a_blockID);
    //return m_dx[a_blockID];
}

template<typename MAP, MemType MEM>
MAP& MBLevelMap<MAP,MEM>::op(BlockIndex a_blockID) const
{
    return *m_ops[a_blockID];
}

template<typename MAP, MemType MEM>
BoxData<double, DIM, MEM> MBLevelMap<MAP,MEM>::X(
        const Box& a_box,
        const Array<double, DIM>& a_dx) const
{
    PR_TIME("MBLevelMap::X");
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

template<typename MAP, MemType MEM>
MBLevelBoxData<double, DIM, MEM, PR_NODE>& MBLevelMap<MAP,MEM>::map()
{
    if (!m_initialized) { initialize(); }
    return *m_X;
}
template<typename MAP, MemType MEM>
const MBLevelBoxData<double, DIM, MEM, PR_NODE>& MBLevelMap<MAP,MEM>::map() const
{
    // PROTO_ASSERT(m_initialized,
    //     "MBLevelMap::map() | Error: Coordinates not initialized (try calling MBLevelMap::initialize())");
    if (!m_initialized) { initialize(); }
    return *m_X;
}
template<typename MAP, MemType MEM>
MBLevelBoxData<double, 1, MEM, PR_CELL>& MBLevelMap<MAP,MEM>::jacobian()
{
    if (!m_initialized) { initialize(); }
    return *m_J;
}
template<typename MAP, MemType MEM>
const MBLevelBoxData<double, 1, MEM, PR_CELL>& MBLevelMap<MAP,MEM>::jacobian() const
{
    // PROTO_ASSERT(m_initialized,
    // "MBLevelMap::map() | Error: Coordinates not initialized (try calling MBLevelMap::initialize())");
    if (!m_initialized) { initialize(); }
    return *m_J;
}

template<typename MAP, MemType MEM>
Point MBLevelMap<MAP,MEM>::convertPoint(Point srcPoint, BlockIndex srcBlock, BlockIndex dstBlock) const
{
    auto xSrc = mappedCoords(srcPoint, srcBlock);
    auto X = op(srcBlock).apply(xSrc);
    auto xDst = op(dstBlock).inverse(X);
    Point dstPoint;
    Box dstDomainBox = layout().getBlock(dstBlock).domain().box();
    Point dstDomainSizes = dstDomainBox.sizes();
    for (int ii = 0; ii < DIM; ii++) { dstPoint[ii] = std::floor(xDst[ii]*dstDomainSizes[ii]); }
    return dstPoint;
}
