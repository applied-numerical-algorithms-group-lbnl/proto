template<template<MemType> typename MAP, MemType MEM>
MBLevelMap<MAP,MEM>::MBLevelMap(
        const MBDisjointBoxLayout&  a_layout,
        const Array<Point, DIM+1>&  a_ghost,
        unsigned int                a_focalBlock)
{
    define(a_layout, a_ghost, a_focalBlock);
}

template<template<MemType> typename MAP, MemType MEM>
MBLevelMap<MAP,MEM>::MBLevelMap(
        const MBDisjointBoxLayout&  a_layout,
        Point                       a_ghost,
        unsigned int                a_focalBlock)
{
    define(a_layout, a_ghost, a_focalBlock);
}

template<template<MemType> typename MAP, MemType MEM>
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

template<template<MemType> typename MAP, MemType MEM>
void MBLevelMap<MAP,MEM>::define(
        const MBDisjointBoxLayout& a_layout,
        const Array<Point, DIM+1>& a_ghost,
        unsigned int               a_focalBlock)
{
    PR_TIME("MBLevelMap::define");
    m_ghost = a_ghost;
    m_focalBlock = a_focalBlock;
    m_c2c = Stencil<double>::CornersToCells(4);
    m_X = std::make_shared<MBLevelBoxData<double, DIM, MEM, PR_NODE>>(a_layout, m_ghost);
    m_J = std::make_shared<MBLevelBoxData<double, 1, MEM, PR_CELL>>(a_layout, m_ghost);
    m_X->setVal(0);
    m_J->setVal(0);
    //m_dx.resize(a_layout.numBlocks());
    m_ops.resize(a_layout.numBlocks());
    
    for (int bi = 0; bi < a_layout.numBlocks(); bi++)
    {
        Box domainBox = a_layout.domain().getBlock(bi).box();
        m_ops[bi].define(dx(bi), bi, a_focalBlock);
    }
    

    for (auto iter : a_layout)
    {
        auto block = a_layout.block(iter);
        Box b_i = a_layout[iter];
        auto& X_i = (*m_X)[iter];
        auto& J_i = (*m_J)[iter];
        BoxData<double, 1, MEM> J0(X_i.box());
        apply(X_i, J0, block);
        J_i |= m_c2c(J0);
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
    PR_TIME("MBLevelMap::apply");
    m_ops[a_block]._apply(a_X, a_J);
}

template<template<MemType> typename MAP, MemType MEM>
void MBLevelMap<MAP,MEM>::apply(
        BoxData<double, DIM, MEM>& a_X,
        BoxData<double, 1, MEM>& a_J,
        FluxBoxData<double, DIM, MEM>& a_NT,
        unsigned int a_block) const
{
    PR_TIME("MBLevelMap::apply");
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
        PR_TIME("MBLevelMap::doApply");
        const auto& domain = layout().domain();
        Box B = domain.convert(a_X.box(), a_outBlock, a_computeBlock, PR_NODE);
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

template<template<MemType> typename MAP, MemType MEM>
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

template<template<MemType> typename MAP, MemType MEM>
BoxData<double, DIM, MEM> MBLevelMap<MAP,MEM>::cellCentered(
        const Box&   a_box,
        unsigned int a_computeBlock,
        unsigned int a_outBlock) const
{
    PR_TIME("MBLevelMap::cellCentered(Box)");
    auto X0 = cellAveraged(a_box.grow(1), a_computeBlock, a_outBlock);
    return Operator::deconvolve(X0);
}

template<template<MemType> typename MAP, MemType MEM>
Array<double, DIM> MBLevelMap<MAP,MEM>::cellAveraged(const MBDataPoint& a_point) const
{
    PR_TIME("MBLevelMap::cellAveraged(MBDataPoint)");
    auto computeBlock = a_point.srcBlock();
    auto domainBlock = a_point.dstBlock();
    Box B0(a_point.point(), a_point.point());
    auto X = cellAveraged(B0, computeBlock, domainBlock);
    return X.array(a_point.point());
}

template<template<MemType> typename MAP, MemType MEM>
Array<double, DIM> MBLevelMap<MAP,MEM>::cellCentered(const MBDataPoint& a_point) const
{
    PR_TIME("MBLevelMap::cellCentered(MBDataPoint)");
    auto computeBlock = a_point.srcBlock();
    auto domainBlock = a_point.dstBlock();
    Box B0(a_point.point(), a_point.point());
    auto X = cellCentered(B0, computeBlock, domainBlock);
    return X.array(a_point.point());
}

template<template<MemType> typename MAP, MemType MEM>
Array<double, DIM> MBLevelMap<MAP,MEM>::mappedCoords(
        const MBDataPoint&          a_point,
        const Array<double, DIM>&   a_offset) const
{
    auto block = a_point.srcBlock();
    auto dx = this->dx(block);
    Array<double, DIM> X = a_point.point();
    X += 0.5; //default is cell center
    X += a_offset;
    X *= dx;
    return X;
}

template<template<MemType> typename MAP, MemType MEM>
Array<double, DIM> MBLevelMap<MAP,MEM>::mappedCoords(
        const MBDataPoint& a_point) const
{
    Array<double, DIM> offset;
    offset.fill(0);
    return mappedCoords(a_point, offset);
}

template<template<MemType> typename MAP, MemType MEM>
Array<double, DIM> MBLevelMap<MAP,MEM>::apply(
        const MBDataPoint&          a_point,
        const Array<double, DIM>&   a_offset) const
{
    auto xi = mappedCoords(a_point, a_offset);
    auto block = a_point.srcBlock();
    return op(block).apply(xi); 
}

template<template<MemType> typename MAP, MemType MEM>
Array<double, DIM> MBLevelMap<MAP,MEM>::apply(
        const MBDataPoint& a_point) const
{
    Array<double, DIM> offset;
    offset.fill(0);
    return apply(a_point, offset);
}


template<template<MemType> typename MAP, MemType MEM>
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

template<template<MemType> typename MAP, MemType MEM>
const MBDisjointBoxLayout& MBLevelMap<MAP,MEM>::layout() const { return m_X->layout(); }
    
template<template<MemType> typename MAP, MemType MEM>
const Array<double, DIM> MBLevelMap<MAP,MEM>::dx(BlockIndex a_blockID) const 
{
    return layout().dx(a_blockID);
    //return m_dx[a_blockID];
}

template<template<MemType> typename MAP, MemType MEM>
MAP<MEM>& MBLevelMap<MAP,MEM>::op(BlockIndex a_blockID) const
{
    return m_ops[a_blockID];
}

template<template<MemType> typename MAP, MemType MEM>
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
