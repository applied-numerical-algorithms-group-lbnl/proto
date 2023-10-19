template<template<MemType> typename MAP, MemType MEM>
MBLevelMap<MAP,MEM>::MBLevelMap(
        const MBDisjointBoxLayout&  a_layout,
        const Array<Point, DIM+1>&  a_ghost,
        unsigned int                a_focalBlock)
{
    define(a_layout, a_ghost, a_focalBlock);
}

template<template<MemType> typename MAP, MemType MEM>
void MBLevelMap<MAP,MEM>::define(
        const MBDisjointBoxLayout& a_layout,
        const Array<Point, DIM+1>& a_ghost,
        unsigned int               a_focalBlock)
{
    m_focalBlock = a_focalBlock;
    m_c2c = Stencil<double>::CornersToCells(4);
    m_X = std::make_shared<MBLevelBoxData<double, DIM, MEM, PR_NODE>>(a_layout, a_ghost);
    m_J = std::make_shared<MBLevelBoxData<double, 1, MEM, PR_CELL>>(a_layout, a_ghost);
    m_X->setVal(0);
    m_J->setVal(0);
    m_dx.resize(a_layout.numBlocks());
    m_ops.resize(a_layout.numBlocks());
    for (int bi = 0; bi < a_layout.numBlocks(); bi++)
    {
        Box domainBox = a_layout.domain().blockDomain(bi).box();
        for (int dir = 0; dir < DIM; dir++)
        {
            m_dx[bi][dir] = 1.0/domainBox.size(dir);
        }
        m_ops[bi].define(m_dx[bi], bi, a_focalBlock);
    }

    for (auto iter : a_layout)
    {
        auto block = a_layout.block(iter);
        Box b_i = a_layout[iter];
        auto& X_i = (*m_X)[iter];
        auto& J_i = (*m_J)[iter];
        apply(X_i, J_i, block);
    }
    m_X->exchange();
    m_J->exchange();
}

template<template<MemType> typename MAP, MemType MEM>
void MBLevelMap<MAP,MEM>::apply(
        BoxData<double, DIM, MEM>& a_X,
        BoxData<double, 1, MEM>& a_J,
        unsigned int a_block) const
{
    m_ops[a_block]._apply(a_X, a_J);
}

template<template<MemType> typename MAP, MemType MEM>
void MBLevelMap<MAP,MEM>::apply(
        BoxData<double, DIM, MEM>& a_X,
        BoxData<double, 1, MEM>& a_J,
        FluxBoxData<double, DIM, MEM>& a_NT,
        unsigned int a_block) const
{
    m_ops[a_block]._apply(a_X, a_J, a_NT);
}

template<template<MemType> typename MAP, MemType MEM>
void MBLevelMap<MAP,MEM>::doApply(
        BoxData<double, DIM, MEM>& a_X,
        BoxData<double, 1, MEM>& a_J,
        unsigned int a_computeBlock,
        unsigned int a_outBlock) const
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

template<template<MemType> typename MAP, MemType MEM>
BoxData<double, DIM, MEM> MBLevelMap<MAP,MEM>::cellAveraged(
        const Box&   a_box,
        unsigned int a_computeBlock,
        unsigned int a_outBlock) const
{
    Box B0 = m_c2c.domain(a_box);
    BoxData<double, DIM, MEM> X0(B0);
    BoxData<double, 1, MEM> J0(B0.extrude(Point::Ones(), -1)); //unused
    doApply(X0, J0, a_computeBlock, a_outBlock);
    return m_c2c(X0);
}

template<template<MemType> typename MAP, MemType MEM>
BoxData<double, DIM, MEM> MBLevelMap<MAP,MEM>::cellCentered(
        const Box&   a_box,
        unsigned int a_computeBlock,
        unsigned int a_outBlock) const
{
    auto X0 = cellAveraged(a_box.grow(1), a_computeBlock, a_outBlock);
    return Operator::deconvolve(X0);
}

template<template<MemType> typename MAP, MemType MEM>
Array<double, DIM> MBLevelMap<MAP,MEM>::cellAveraged(const MBDataPoint& a_point) const
{
    auto computeBlock = a_point.srcBlock();
    auto domainBlock = a_point.dstBlock();
    Box B0(a_point.point, a_point.point);
    auto X = cellAveraged(B0, computeBlock, domainBlock);
    return X.array(a_point.point);
}

template<template<MemType> typename MAP, MemType MEM>
Array<double, DIM> MBLevelMap<MAP,MEM>::cellCentered(const MBDataPoint& a_point) const
{
    auto computeBlock = a_point.srcBlock();
    auto domainBlock = a_point.dstBlock();
    Box B0(a_point.point, a_point.point);
    auto X = cellCentered(B0, computeBlock, domainBlock);
    return X.array(a_point.point);
}


template<template<MemType> typename MAP, MemType MEM>
const MBDisjointBoxLayout& MBLevelMap<MAP,MEM>::layout() const { return m_X->layout(); }
    
template<template<MemType> typename MAP, MemType MEM>
const Array<double, DIM>& MBLevelMap<MAP,MEM>::dx(int a_blockID) const 
{
    return m_dx[a_blockID];
}

template<template<MemType> typename MAP, MemType MEM>
BoxData<double, DIM, MEM> MBLevelMap<MAP,MEM>::X(
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

template<template<MemType> typename MAP, MemType MEM>
MBLevelBoxData<double, DIM, MEM, PR_NODE>& MBLevelMap<MAP,MEM>::map()
{
    return *m_X;
}
template<template<MemType> typename MAP, MemType MEM>
const MBLevelBoxData<double, DIM, MEM, PR_NODE>& MBLevelMap<MAP,MEM>::map() const
{
    return *m_X;
}
template<template<MemType> typename MAP, MemType MEM>
MBLevelBoxData<double, 1, MEM, PR_CELL>& MBLevelMap<MAP,MEM>::jacobian()
{
    return *m_J;
}
template<template<MemType> typename MAP, MemType MEM>
const MBLevelBoxData<double, 1, MEM, PR_CELL>& MBLevelMap<MAP,MEM>::jacobian() const
{
    return *m_J;
}
