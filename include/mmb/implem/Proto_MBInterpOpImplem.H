
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

MBPointInterpOp::MBPointInterpOp(
        MBDataPoint              a_dst,
        std::vector<MBDataPoint> a_srcs,
        std::vector<double>      a_coefs)
{
    PROTO_ASSERT(a_srcs.size() == a_coefs.size(),
            "MBPointInterpOp::constructor | Error: source and coefs should be the same size.");
    m_dst = a_dst;
    m_srcs = a_srcs;
    m_coefs = a_coefs;
}

template<template<MemType> typename MAP, MemType MEM>
MBPointInterpOp::MBPointInterpOp(
        MBDataPoint a_dst,
        Point a_ghost,
        const MBLevelMap<MAP, MEM>& a_map,
        std::vector<Point> a_footprint,
        unsigned int a_order)
{
    PR_TIME("MBPointInterpOp::constructor");
    m_dst = a_dst;
    m_srcs = a_map.map().interpFootprint(a_dst.point, a_ghost, a_footprint, a_dst.index);
     
    // References
    auto& layout = a_map.map().layout();
    unsigned int dstBlock = a_dst.dstBlock();
    Stencil<double> av2 = 1.0*Shift::Zeros();
    Stencil<double> HOCorners = (1.0*Shift::Zeros() + (-1./12.)*(Stencil<double>::Laplacian()));
    for (int dir = 0; dir < DIM; dir++)
    {
        av2=av2*(.5*Shift::Zeros() + .5*Shift::Basis(dir,1));
    }
    
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
    std::vector<Point> exponents;
    //m_exponents.clear();
    for (auto bi : Box(expBoxSize))
    {
        if (bi.sum() < a_order)
        {
            exponents.push_back(bi);
        }
    }

    int N = exponents.size();
    int M = m_srcs.size();
    
    if (N > M)
    {
        pr_out() << "Underdetermined system detected" << std::endl;
        pr_out() << "\tN = " << N << ", M = " << M << std::endl;
        pr_out() << "\tDst Point: " << a_dst << std::endl;
        pr_out() << "\tSrc Points: " << std::endl;
        for (auto si : m_srcs)
        {
            pr_out() << "\t\t" << si << std::endl;
        }
    }
    //std::cout << "N: " << N << " | M: " << M << std::endl;

    // Rg = 1/M*sum_j(sqrt(<xj - xg>^2))
    // Cjp = <(x-xg)^p/Rg>_j  |  p = nomial, j = source cell
    // Mp = <(x-xg)^p/Rg>_g   |  p = nomial
    // Sj = sum_p(Mp*Cinv_pj) = m_coefs[j]
    
    // Compute Rg
    // This is the approximate radius in physical space of the operator
    // Normalizing by Rg helps keep the matrices well conditioned
    double Rg = 0.0;
    //auto xg = a_map.cellCentered(a_dst);
    auto xg = a_map.apply(a_dst);
    
    for (int ii = 0; ii < M; ii++)
    {
        //auto xj = a_map.analyticCellCentered(m_srcs[ii]);
        auto xj = a_map.apply(m_srcs[ii]);
        auto dist = xj - xg;
        Rg += dist.norm();
    }
    Rg /= M;
    // Cache Moment Data
    // We cache < J * ((x-xg)/Rg)^p > for each exponent p
    //auto C2C = Stencil<double>::CornersToCells(a_order);
    auto C2C = Stencil<double>::CornersToCells(4);
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
        BoxData<double, 1, MEM> Ji(Xi.box());
        a_map.doApply(Xi, Ji, block, dstBlock); 
        Xi -= xg;
        Xi /= Rg;
        for (int jj = 0; jj < N; jj++)
        {
            auto Xi_p = forall<double, 1>(f_binomialPower, Xi, exponents[jj]);
            BoxData<double> Xi_p_avg = C2C(Xi_p);
            srcData[block].push_back(std::move(Xi_p_avg));
        }
    }  
    // Compute Matrices
    Matrix<double> C(M,N);
    //m_C.define(M,N);
    for (int ii = 0; ii < M; ii++)
    {
        auto& src = m_srcs[ii];
        auto srcBlock = src.srcBlock();
        for (int jj = 0; jj < N; jj++)
        {
            C(ii,jj) = srcData[srcBlock][jj](src.point);
        }
    }
    Matrix<double> D(1,N);
    //m_D.define(1,N);
    for (int jj = 0; jj < N; jj++)
    {
        D(0,jj) = srcData[dstBlock][jj](a_dst.point);
    }
    
    // Compute "Stencil"
    auto Cinv = C.inverse();
    //m_Cinv = m_C.inverse();
    auto S = D*Cinv;
    //m_S = m_D*m_Cinv;
    m_coefs.clear();
    m_coefs.resize(M);
    //std::cout << "\nComputed operator for dest " << a_dst << std::endl;
    for (int ii = 0; ii < M; ii++)
    {
        m_coefs[ii] = S(0,ii);
        //std::cout << "\t" << m_srcs[ii].point - a_dst.point << " : " << m_coefs[ii] << std::endl;
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
    /*
    PROTO_ASSERT(a_srcComponent < C, "Component index %i is out of bounds.", a_srcComponent);
    int M = m_Cinv.M();
    int N = m_Cinv.N();
    PROTO_ASSERT(N == m_srcs.size(),"");
    Matrix<double> Y0(N,1);
    for (int ii = 0; ii < N; ii++)
    {
        Y0(ii,0) = a_src[m_srcs[ii]](a_srcComponent);
    }
    */
    //return m_Cinv*Y0;
    PROTO_ASSERT(false,
            "MBPointInterp::coefs | Error: Deprecated (no longer implemented)");
    return Matrix<double>(1,1);
}

unsigned int MBPointInterpOp::size() const
{
    return m_coefs.size();
}

bool MBPointInterpOp::validate() const
{
    for (auto si : m_srcs)
    {
        if (!si.grid().defined()) { return false; }
    }
    return true;
}

void MBPointInterpOp::print() const
{
    pr_out() << "--------------------------------" << std::endl;
    pr_out() << m_dst << std::endl;
    for (int ii = 0; ii < size(); ii++) {
        pr_out() << "\tcoef: " << m_coefs[ii] << " | " << m_srcs[ii] << std::endl;
    }
}

MBInterpOp::MBInterpOp(Point a_ghost, unsigned int a_order)
{
    PR_TIME("MBInterpOp::constructor");
    m_order = a_order;
    m_ghost = a_ghost;
    m_ops.clear();
    m_points.clear();
}

template<template<MemType> typename MAP, MemType MEM>
MBInterpOp::MBInterpOp(
        const MBLevelMap<MAP,MEM>&    a_map,
        int                           a_order)
{
    for (int bi = 0; bi < a_map.layout().numBlocks(); bi++)
    {
        define(a_map, a_order, bi);
    }
}

template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp::define(
        const MBLevelMap<MAP,MEM>&    a_map,
        int                     a_order,
        int                     a_block)
{
    std::vector<Point> footprint;
    switch (a_order)
    {
        case 4:
            for (auto pi : Box::Kernel(2))
            {
                if (pi.abs().sum() <= 2)
                {
                    footprint.push_back(pi);
                }
            }
            break;
        case 5:
            for (auto pi : Box::Kernel(2))
            {
                if (pi.abs().sum() <= 3)
                {
                    footprint.push_back(pi);
                }
            }
            break;
        default:
            MayDay<void>::Error("MBInterpOp::define | Function undefined for order %i", a_order);
            break;
    }
    define(a_map, footprint, a_order, a_block);
}

template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp::define(
        const MBLevelMap<MAP,MEM>&  a_map,
        std::vector<Point>          a_footprint,
        int                         a_order,
        int                         a_block)
{
    MBInterpLayout interpLayout;
    for (auto fi : a_footprint)
    {
        interpLayout.addPoint(fi);
    }
    define(a_map, interpLayout, a_order, a_block);
}

template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp::define(
        const MBLevelMap<MAP,MEM>&  a_map,
        MBInterpLayout              a_layout,
        int                         a_order,
        int                         a_block)
{
    m_layout = a_layout;
    m_ghost = a_map.ghost()[0];
    const auto& layout = a_map.map().layout();
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        if ((a_block >= 0) && (block != a_block)) {continue; }
        Box blockDomainBox = layout.domain().getBlock(block).box(); 
        Box patchBox = layout[iter];
        for (auto dir : Box::Kernel(1))
        {
            Box boundBox = patchBox.adjacent(m_ghost*dir);
            if (blockDomainBox.contains(boundBox))
            {
                continue;
            }
            if (layout.domain().graph().isTriplePoint(block, dir))
            {
                for (auto bi : boundBox)
                {
                    // check if this point is already done
                    auto key = std::make_tuple(bi, iter.global(), block);
                    if (m_points.count(key)) {continue; }

                    Point corner = patchBox.edge(dir,1).low();
                    std::set<Point> fset;
                    for (int ii = 0; ii < DIM; ii++)
                    {
                        if (dir[ii] == 0) { continue; }
                        Point s = Point::Basis(ii,-1)*(bi-corner);
                        for (auto fi : a_layout.footprint())
                        {
                            fset.insert(fi + s);
                        }
                    }
                    std::vector<Point> extFootprint(fset.begin(), fset.end());
                    MBDataPoint dstDataPoint(iter, bi, layout);
                    m_ops.push_back(MBPointInterpOp( dstDataPoint, m_ghost, a_map, extFootprint, a_order));
                    auto& op = m_ops[m_ops.size()-1];
                    m_points.insert(key);
                    
                    //do copy
                    for (int dir = 0; dir < DIM; dir++)
                    {
                        if (a_layout.doCopy(dir))
                        {
                            int colHeight = boundBox.sizes()[dir];
                            Box boundColumn = Box(bi,bi).grow(Point::Basis(dir,colHeight));
                            boundColumn &= boundBox;
                            for (auto bj : boundColumn)
                            {
                                if (bi == bj) {continue; }
                                auto key_j = std::make_tuple(bj, iter.global(), block);
                                if (m_points.count(key_j)) {continue; }
                                Point shift = bj - bi;
                                MBDataPoint dst_j = dstDataPoint.shift(shift);
                                std::vector<MBDataPoint> srcs_j;
                                for (auto si : op.sources())
                                {
                                    srcs_j.push_back(si.shift(shift));
                                }
                                MBPointInterpOp op_j(dst_j, srcs_j, op.coefs());
                                m_ops.push_back(op_j);
                                m_points.insert(key_j);    
                            }
                        }
                    }
                }
                continue;
            }
            auto bounds = a_map.map().bounds(iter, dir);
            std::set<Point> fset;
            for (auto fi : a_layout.footprint())
            {
                fset.insert(fi);
            }
            for (auto bound : bounds)
            {
                //handle points adjacent to triple points
                for (auto bi : boundBox)
                {
                    auto key = std::make_tuple(bi, iter.global(), block);
                    if (m_points.count(key)) {continue; }
                    for (int ii = 0; ii < DIM; ii++)
                    {
                        Point d2 = Point::Basis(ii);
                        if (dir[ii] != 0) {continue; }
                        bool onLoEdge = boundBox.edge(-d2).contains(bi);
                        bool onHiEdge = boundBox.edge( d2).contains(bi);
                        if (onLoEdge && layout.domain().graph().isTriplePoint(block, dir - d2))
                        {
                            for (auto fi : a_layout.footprint())
                            {
                                fset.insert(fi + d2);
                            }
                        }
                        if (onHiEdge && layout.domain().graph().isTriplePoint(block, dir + d2))
                        {
                            for (auto fi : a_layout.footprint())
                            {
                                fset.insert(fi - d2);
                            }
                        }
                    }
                    std::vector<Point> extFootprint(fset.begin(), fset.end());
                    MBDataPoint dstDataPoint(iter, bi, layout);
                    m_ops.push_back(MBPointInterpOp(dstDataPoint, m_ghost, a_map, extFootprint, a_order));
                    auto& op = m_ops[m_ops.size()-1];
                    m_points.insert(key);
                    
                    //do copy
                    for (int dir = 0; dir < DIM; dir++)
                    {
                        if (a_layout.doCopy(dir))
                        {
                            int colHeight = boundBox.sizes()[dir];
                            Box boundColumn = Box(bi,bi).grow(Point::Basis(dir,colHeight));
                            boundColumn &= boundBox;
                            for (auto bj : boundColumn)
                            {
                                if (bi == bj) {continue; }
                                auto key_j = std::make_tuple(bj, iter.global(), block);
                                if (m_points.count(key_j)) {continue; }
                                Point shift = bj - bi;
                                MBDataPoint dst_j = dstDataPoint.shift(shift);
                                std::vector<MBDataPoint> srcs_j;
                                for (auto si : op.sources())
                                {
                                    srcs_j.push_back(si.shift(shift));
                                }
                                MBPointInterpOp op_j(dst_j, srcs_j, op.coefs());
                                m_ops.push_back(op_j);
                                m_points.insert(key_j);    
                            }
                        }
                    }
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

void MBInterpOp::print() const
{
    pr_out() << "Number of PointInterpOps: " << m_ops.size() << std::endl;
    for (auto& op : m_ops)
    {
        op.print(); 
    }
}

bool MBInterpOp::validate() const
{
    for (auto& op : m_ops)
    {
        if (!op.validate()) { return false; }
    }
    return true;
}

template<template<MemType> typename MAP, typename T, unsigned int C, MemType MEM, Centering CTR>
void interpBoundaries(MBLevelBoxData<T, C, MEM, CTR>& a_data, unsigned int a_order)
{
    MBLevelMap<MAP, MEM> map(a_data.layout(), a_data.ghost());
    interpBoundaries<MAP, T, C, MEM, CTR>(a_data, map, a_order);
}

template<template<MemType> typename MAP, typename T, unsigned int C, MemType MEM, Centering CTR>
void interpBoundaries(
        MBLevelBoxData<T, C, MEM, CTR>& a_data,
        MBLevelMap<MAP, MEM>& a_map,
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






