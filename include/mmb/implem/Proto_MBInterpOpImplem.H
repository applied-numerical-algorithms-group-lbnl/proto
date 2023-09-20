
template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_binomialPower_tmp(Var<T, 1, MEM>& a_xp, const Var<T, DIM, MEM>& a_x, Point a_p)
{
    a_xp(0) = 1.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        if (a_p[dir] != 0)
        {
            a_xp(0) *= pow(a_x(dir), a_p[dir]);
        }
    }
}
PROTO_KERNEL_END(f_binomialPower_tmp, f_binomialPower)

template<template<MemType> typename MAP, MemType MEM>
MBPointInterpOp::MBPointInterpOp(
        MBDataPoint a_dst,
        Point a_ghost,
        MBLevelMap<MAP, MEM>& a_map,
        const std::vector<Point>& a_footprint,
        unsigned int a_order)
{
    PR_TIME("MBPointInterpOp::constructor");
    m_dst = a_dst;
    m_srcs = a_map.map().interpFootprint(a_dst.point, a_ghost, a_footprint, a_dst.index);
    
    // References
    auto& layout = a_map.map().layout();
    unsigned int dstBlock = a_dst.dstBlock();
    
    // Compute Exponents
    Box spanBox;
    for (auto si : m_srcs)
    {
        spanBox &= si.point;
    }
    Point expBoxSize = Point::Ones(a_order);
    // Cull exponents if the footprint is "flat" in a given direction
    // This allows for sub-DIM dimensional operators
    for (int ii = 0; ii < DIM; ii++)
    {
        if (spanBox.size(ii) == 1) { expBoxSize[ii] = 1; }
    }
    m_exponents.clear();
    for (auto bi : Box(expBoxSize))
    {
        if (bi.sum() < a_order)
        {
            m_exponents.push_back(bi);
        }
    }
    int N = m_exponents.size();
    int M = m_srcs.size();

    // Rg = 1/M*sum_j(sqrt(<xj - xg>^2))
    // Cjp = <(x-xg)^p/Rg>_j  |  p = nomial, j = source cell
    // Mp = <(x-xg)^p/Rg>_g   |  p = nomial
    // Sj = sum_p(Mp*Cinv_pj) = m_coefs[j]
    
    // Compute Rg
    // This is the approximate radius in physical space of the operator
    // Normalizing by Rg helps keep the matrices well conditioned
    double Rg = 0.0;
    auto xg = a_map.cellCentered(a_dst);
    for (int ii = 0; ii < M; ii++)
    {
        auto xj = a_map.cellCentered(m_srcs[ii]);
        auto dist = xj - xg;
        Rg += dist.norm();
    }
    Rg /= M;
    // Cache Moment Data
    // We cache < J * ((x-xg)/Rg)^p > for each exponent p
    auto C2C = Stencil<double>::CornersToCells(a_order);
    std::unordered_map<int, Box> srcBoxes;
    std::unordered_map<int, std::vector<BoxData<double>>> srcData; 
    Box B0(a_dst.point, a_dst.point);
    srcBoxes[a_dst.dstBlock()] = B0;
    for (auto si : m_srcs)
    {
        int bi = si.srcBlock();
        if (srcBoxes.find(bi) == srcBoxes.end())
        {
            srcBoxes[bi] = Box(si.point, si.point); 
        } else {
            srcBoxes[bi] &= si.point;
        }
    }
    for (auto& [block, box] : srcBoxes)
    {
        BoxData<double, DIM, MEM> Xi(C2C.domain(box).grow(1));
        BoxData<double, 1, MEM> Ji(box.grow(1));
        a_map.doApply(Xi, Ji, block, dstBlock); 
        Xi -= xg;
        Xi /= Rg;
        for (int jj = 0; jj < N; jj++)
        {
            auto Xi_p = forall<double, 1>(f_binomialPower, Xi, m_exponents[jj]);
            BoxData<double> Xi_p_avg = C2C(Xi_p);
            Ji.setVal(1.0);
            srcData[block].push_back(Operator::cellProduct(Ji, Xi_p_avg));
        }
    }
   
    // Compute Matrices
    m_C.define(M,N);
    for (int ii = 0; ii < M; ii++)
    {
        auto& src = m_srcs[ii];
        auto srcBlock = src.srcBlock();
        for (int jj = 0; jj < N; jj++)
        {
            m_C(ii,jj) = srcData[srcBlock][jj](src.point);
        }
    }
    m_D.define(1,N);
    for (int jj = 0; jj < N; jj++)
    {
        m_D(0,jj) = srcData[dstBlock][jj](a_dst.point);
    }
    
    // Compute "Stencil"
    m_Cinv = m_C.inverse();
    m_S = m_D*m_Cinv;
    m_coefs.clear();
    m_coefs.resize(M);
    for (int ii = 0; ii < M; ii++)
    {
        m_coefs[ii] = m_S(0,ii);
    }
}

template<typename T, unsigned int C, MemType MEM>
void MBPointInterpOp::apply(
        MBLevelBoxData<T, C, MEM>& a_dst,
        MBLevelBoxData<T, C, MEM>& a_src)
{
    PR_TIME("MBPointInterpOp::apply");
    auto dst = a_dst[m_dst]; //dst is a Var<T, C, MEM>
    for (int ii = 0; ii < size(); ii++)
    {
        auto src = a_src[m_srcs[ii]]; //src is a Var<T, C, MEM>
        for (int cc = 0; cc < C; cc++)
        {
            if (ii == 0) { dst(cc) = 0.0; }
            dst(cc) += m_coefs[ii]*src(cc);
        }
    }
}

template<typename T, unsigned int C, MemType MEM>
Matrix<double> MBPointInterpOp::coefs(
        MBLevelBoxData<T, C, MEM>& a_src,
        int a_srcComponent)
{
    PROTO_ASSERT(a_srcComponent < C, "Component index %i is out of bounds.", a_srcComponent);
    int M = m_Cinv.M();
    int N = m_Cinv.N();
    PROTO_ASSERT(N == m_srcs.size(),"");
    Matrix<double> Y0(N,1);
    for (int ii = 0; ii < N; ii++)
    {
        Y0(ii,0) = a_src[m_srcs[ii]](a_srcComponent);
    }
    return m_Cinv*Y0;
}

unsigned int MBPointInterpOp::size() const
{
    return m_coefs.size();
}

MBInterpOp::MBInterpOp(Point a_ghost, unsigned int a_order)
{
    PR_TIME("MBInterpOp::constructor");
    m_order = a_order;
    m_ghost = a_ghost;
    m_ops.clear();
}

template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp::define(
        MBLevelMap<MAP,MEM>&            a_map,
        const std::vector<Point>&   a_footprint,
        int                         a_block)
{
    const auto& layout = a_map.map().layout();
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        if ((a_block >= 0) && (block != a_block)) {continue; }
        Box blockDomainBox = layout.domain().blockDomain(block).box(); 
        Box patchBox = layout[iter];
        for (auto dir : Box::Kernel(1))
        {
            if (layout.domain().graph().isTriplePoint(block, dir))
            {
                std::set<Point> fset;
                auto dirs = dir.split();
                for (auto fi : a_footprint)
                {
                    fset.insert(fi);
                    for (auto di : dirs)
                    {
                        fset.insert(fi-di);
                    }
                }
                std::vector<Point> extFootprint(fset.begin(), fset.end());
                Box boundBox = patchBox.adjacent(m_ghost*dir);
                if (blockDomainBox.contains(boundBox))
                {
                    continue;
                }
                for (auto bi : boundBox)
                {
                    MBDataPoint dstDataPoint(iter, bi, layout);
                    m_ops.push_back(
                            MBPointInterpOp(dstDataPoint, m_ghost, a_map, extFootprint, 4));
                }
                continue;
            }
            auto bounds = a_map.map().bounds(iter, dir);
            for (auto bound : bounds)
            {
                Box boundBox = patchBox.adjacent(m_ghost*dir);
                if (blockDomainBox.contains(boundBox))
                {
                    continue;
                }
                for (auto bi : boundBox)
                {
                    MBDataPoint dstDataPoint(iter, bi, layout);
                    m_ops.push_back(
                            MBPointInterpOp(dstDataPoint, m_ghost, a_map, a_footprint, 4));
                }
            }
        }
    }
}

template<typename T, unsigned int C, MemType MEM>
void MBInterpOp::apply(
        MBLevelBoxData<T, C, MEM>& a_dst,
        MBLevelBoxData<T, C, MEM>& a_src)
{
    PR_TIME("MBInterpOp::apply");
    for (auto& op : m_ops)
    {
        auto block = a_dst.layout().block(op.target().index);
        op.apply(a_dst, a_src);
    }
}

template<typename T, unsigned int C, unsigned int P, MemType MEM>
void MBInterpOp::coefs(
        MBLevelBoxData<T, P, MEM>& a_coefs,
        MBLevelBoxData<T, C, MEM>& a_src,
        int a_srcComponent)
{
    for (auto& op : m_ops)
    {
        auto poly = op.coefs(a_src, a_srcComponent);
        for (int pi = 0; pi < P; pi++)
        {
            a_coefs[op.target()](pi) = poly(pi,0);
        }
    }
}

template<template<MemType> typename Map, typename T, unsigned int C, MemType MEM, Centering CTR>
void interpBoundaries(MBLevelBoxData<T, C, MEM, CTR>& a_data, unsigned int a_order)
{
    Map<MEM> map;
    map.define(a_data.layout(), a_data.ghost());
    interpBoundaries<Map, T, C, MEM, CTR>(a_data, map, a_order);
}

template<template<MemType> typename Map, typename T, unsigned int C, MemType MEM, Centering CTR>
void interpBoundaries(
        MBLevelBoxData<T, C, MEM, CTR>& a_data,
        Map<MEM>& a_map,
        unsigned int a_order)
{
    std::vector<Point> footprint;
    Box K = Box::Kernel(2);
    for (auto pi : K)
    {
        if (pi.abs().sum() <= 2)
        {
            footprint.push_back(pi);
        }
    }
    MBInterpOp op(a_data.ghost()[0], a_order);
    op.define(a_map, footprint);
    op.apply(a_data, a_data);
}






