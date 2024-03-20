
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
    m_srcs = a_map.map().interpFootprint(a_dst.point(), a_ghost, a_footprint, a_dst.index());
    
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
        spanBox &= si.point();
    }
    Point expBoxSize = Point::Ones(a_order);
    // Cull exponents if the footprint is "flat" in a given direction
    // This allows for sub-DIM dimensional operators
    for (int ii = 0; ii < DIM; ii++)
    {
        if (spanBox.size(ii) == 1) { expBoxSize[ii] = 1; }
    }
    std::vector<Point> exponents;
    for (auto bi : Box(expBoxSize))
    {
        if (bi.sum() < a_order)
        {
            exponents.push_back(bi);
        }
    }

    int N = exponents.size();
    int M = m_srcs.size();
   
#if PR_VERBOSE > 1
    pr_out() << "\t\t\t\t\tAdding MBPointInterpOp | point: " << m_dst.point() << " | N: " << N << " | M: " << M << std::endl;
#endif

    if (N > M)
    {
        pr_out() << "====================================================================" << std::endl;
        pr_out() << "Underdetermined system detected" << std::endl;
        pr_out() << "\tN = " << N << ", M = " << M << std::endl;
        pr_out() << "\tDst Point: " << a_dst << std::endl;
        pr_out() << "\tGhost: " << a_ghost << std::endl;
        #if PR_VERBOSE > 1
        pr_out() << "\tFootprint: " << std::endl;
        for (auto fi : a_footprint)
        {
            pr_out() << "\t\t" << fi << std::endl;
        }
        pr_out() << "\tSrc Points: " << std::endl;
        for (auto si : m_srcs)
        {
            pr_out() << "\t\t" << si << std::endl;
        }
        pr_out() << "\tExponents: " << std::endl;
        for (auto ei : exponents)
        {
            pr_out() << "\t\t" << ei << std::endl;
        }
    
        #endif
        pr_out() << "====================================================================" << std::endl;
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
      PR_TIMERS("compute_Rg");
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
    Box B0(a_dst.point(), a_dst.point());
    srcBoxes[a_dst.dstBlock()] = B0;
    for (auto si : m_srcs)
    {
        PR_TIMERS("compute_moment_domains");
        int bi = si.srcBlock();
        if (srcBoxes.find(bi) == srcBoxes.end())
        {
            srcBoxes[bi] = Box(si.point(), si.point()); 
        } else {
            srcBoxes[bi] &= si.point();
        }
    }
    
    for (auto& [block, box] : srcBoxes)
      {
        PR_TIMERS("compute_moments");
        BoxData<double, DIM, MEM> Xi(C2C.domain(box).grow(1));
        BoxData<double, 1, MEM> Ji(Xi.box());
        a_map.doApply(Xi, Ji, block, dstBlock); 
        Xi -= xg;
        Xi /= Rg;
        for (int jj = 0; jj < N; jj++)
        {
            auto Xi_p = forall<double, 1>(f_binomialPower, Xi, exponents[jj]);
            {
                BoxData<double> temp = HOCorners(Xi_p);
                {
                    BoxData<double> Xi_p_avg = av2(temp);
                    srcData[block].push_back(std::move(Xi_p_avg));
                }
            }
        }
    }  
    // Compute Matrices
    Matrix<double> C(M,N);
    for (int ii = 0; ii < M; ii++)
    {
        PR_TIMERS("build_C_matrix");
        auto& src = m_srcs[ii];
        auto srcBlock = src.srcBlock();
        for (int jj = 0; jj < N; jj++)
        {
            C(ii,jj) = srcData[srcBlock][jj](src.point());
        }
    }
    Matrix<double> D(1,N);
    for (int jj = 0; jj < N; jj++)
    {
        PR_TIMERS("build_D_matrix");
        D(0,jj) = srcData[dstBlock][jj](a_dst.point());
    }
    
    // Compute "Stencil"
    auto Cinv = C.inverse();
    auto S = D*Cinv;
    m_coefs.resize(M);
    for (int ii = 0; ii < M; ii++)
    {
        m_coefs[ii] = S(0,ii);
    }
}

template<typename T, unsigned int C, MemType MEM>
void MBPointInterpOp::apply(
        MBLevelBoxData<T, C, MEM>& a_dst,
        MBLevelBoxData<T, C, MEM>& a_src)
{
  //PR_TIME("MBPointInterpOp::apply");
    auto dst = a_dst[m_dst]; //dst is a Var<T, C, MEM>
    for (int ii = 0; ii < size(); ii++)
    {
      //PR_TIMERS("size loop");
        auto src = a_src[m_srcs[ii]]; //src is a Var<T, C, MEM>
        for (int cc = 0; cc < C; cc++)
        {
          //PR_TIME("increment");
            if (ii == 0) { dst(cc) = 0.0; }
            dst(cc) += m_coefs[ii]*src(cc);
        }
    }
}

template<typename T, unsigned int C, MemType MEM>
void MBInterpOp::printErrorPoints(
        MBLevelBoxData<T, C, MEM>& a_err,
        T                          a_eps)
{
    auto& layout = a_err.layout();
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        auto& err_i = a_err[iter];
        for (auto pi : err_i.box())
        {
            for (int cc = 0; cc < C; cc++)
            {
                if (err_i(pi, cc) > a_eps)
                {   
                    for (auto& op : m_ops)
                    {
                        auto& dst = op->target();
                        if (dst.point() == pi && dst.dstBlock() == block)
                        {
                            std::cout << "\nError Point found | error: " << err_i(pi, cc) << " | target: " << pi << " | block: " << block << std::endl;
                            std::cout << "Source Points: " << std::endl;
                            for (auto& si : op->sources())
                            {
                                std::cout << "\t" << si << std::endl;
                            }
                            break;
                        }
                    }
                    break;
                }
                
            }
        }
    }
}


unsigned int MBPointInterpOp::size() const
{
    return m_coefs.size();
}

void MBPointInterpOp::print() const
{
    pr_out() << "--------------------------------" << std::endl;
    pr_out() << m_dst << std::endl;
    for (int ii = 0; ii < size(); ii++) {
        pr_out() << "\tcoef: " << m_coefs[ii] << " | " << m_srcs[ii] << std::endl;
    }
}

MBInterpOp::MBInterpOp()
{
    PR_TIME("MBInterpOp::constructor");
    m_copyAxes.fill(false);
}

MBInterpOp::MBInterpOp(Point a_ghost, unsigned int a_order)
{
    PR_TIME("MBInterpOp::constructor");
    m_order = a_order;
    m_ghost = a_ghost;
    m_copyAxes.fill(false);
}

template<template<MemType> typename MAP, MemType MEM>
MBInterpOp::MBInterpOp(
        const MBLevelMap<MAP,MEM>&    a_map,
        int                           a_order)
{
    m_copyAxes.fill(false);
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
    Box spanBox;
    for (auto fi : a_footprint)
    {
        spanBox &= fi;
    }
    Point footprintDims = Point::Zeros();
    for (int dd = 0; dd < DIM; dd++)
    {
        if (spanBox.sizes()[dd] > 1)
        {
            footprintDims[dd] = 1;
        }
    }
#if PR_VERBOSE > 1
    pr_out() << "MBInterpOp Constructor" << std::endl;
    pr_out() << "Order: " << a_order << " | Block: " << a_block << std::endl;
    pr_out() << "Base Footprint: " << std::endl;
    for (auto fi : a_footprint)
    {
        pr_out() << fi << ", ";
    }
    pr_out() << std::endl;
#endif
    PR_TIME("MBInterpOp::define");
    m_ghost = a_map.ghost()[0];
    
    const auto& layout = a_map.map().layout();
    const auto& patchDomain = layout.partition().domain();
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        Point patch_i = layout.point(iter);
        if ((a_block >= 0) && (block != a_block)) {continue; }
        Box blockDomainBox = layout.domain().getBlock(block).box(); 
        Box patchBox = layout[iter];
#if PR_VERBOSE > 1
        pr_out() << "\tPatch: " << patch_i << " | Block: " << block << std::endl;
        pr_out() << "\tPatch Box: " << patchBox << " | Block Box: " << blockDomainBox << std::endl;
#endif
        for (auto dir : Box::Kernel(1))
        {
            Point patch_j = patch_i + dir; 
#if PR_VERBOSE > 1
            pr_out() << "\t\tDir: " << dir << " | Neighbor: " << patch_j << std::endl;
            if (patchDomain.inInterior(patch_j, block))
            {
                pr_out() << "\t\t\tNeighbor is in the interior -> skipping" << std::endl;
            }
            if (patchDomain.inDomainBoundary(patch_j, block))
            {
                pr_out() << "\t\t\tNeighbor is in the domain boundary -> skipping" << std::endl;
            }
#endif
            if (patchDomain.inInterior(patch_j, block)) { continue; }
            if (patchDomain.inDomainBoundary(patch_j, block)) { continue; }
            Box boundBox = patchBox.adjacent(m_ghost*dir);
            if (patchDomain.inTriplePoint(patch_j, block))
            {
#if PR_VERBOSE > 1
                pr_out() << "\t\t\tNeighbor is inside a triple point region: " << boundBox << std::endl;
#endif
                for (auto bi : boundBox)
                {
#if PR_VERBOSE > 1
                    pr_out() << "\t\t\tChecking point: " << bi << std::endl;
#endif
                    PR_TIME("add_triple_point");
                    MBPointID key(bi, iter.global());
#if PR_VERBOSE > 1
                    if (m_points.count(key))
                    {
                        pr_out() << "\t\t\t\tOperator already exists -> skipping" << std::endl;
                    }
#endif
                    if (m_points.count(key)) { continue; }
                    Point corner = patchBox.edge(dir,1).low();
                    std::set<Point> fset;
                    for (int ii = 0; ii < DIM; ii++)
                    {
                        if (dir[ii] == 0) { continue; }
                        Point s = Point::Basis(ii,-1)*(bi-corner)*footprintDims;
                        for (auto fi : a_footprint)
                        {
                            fset.insert(fi + s);
                        }
                    }
                    std::vector<Point> extFootprint(fset.begin(), fset.end());
#if PR_VERBOSE > 1
                    pr_out() << "\t\t\t\tExtended Footprint: " << std::endl;
                    for (auto fi : extFootprint)
                    {
                        pr_out() << "\t\t\t\t" << fi;
                        if (std::find(a_footprint.begin(), a_footprint.end(), fi) == a_footprint.end())
                        {
                            pr_out() << " <- added";
                        }
                        pr_out() << std::endl;
                    }
#endif
                    MBDataPoint dstDataPoint(iter, bi, layout);
                    auto op = std::make_shared<MBPointInterpOp>(dstDataPoint, m_ghost, a_map, extFootprint, a_order);
                    m_ops.push_back(op);
                    m_points.insert(key);
                    
                    for (int dd = 0; dd < DIM; dd++)
                    {
                        PR_TIME("copy_along_axis");
                        if (m_copyAxes[dd])
                        {
                            int colHeight = boundBox.sizes()[dd];
                            Box boundColumn = Box(bi,bi).grow(Point::Basis(dd,colHeight-1));
                            boundColumn &= boundBox;
                            for (auto bj : boundColumn)
                            {
                                if (bi == bj) {continue; }
                                MBPointID key_j(bj, iter.global());
                                if (m_points.count(key_j)) {continue; }
                                Point shift = bj - bi;
                                MBDataPoint dst_j = dstDataPoint.shift(shift);
                                std::vector<MBDataPoint> srcs_j;
                                for (auto si : op->sources())
                                {
                                    srcs_j.push_back(si.shift(shift));
                                }
                                auto coefs = op->coefs();
                                auto op_j = std::make_shared<MBPointInterpOp>(dst_j, srcs_j, coefs);
                                m_ops.push_back(op_j);
                                m_points.insert(key_j);    
                            }
                        }
                    }
                }
                continue;
            }
            PROTO_ASSERT(patchDomain.inBlockBoundary(patch_j, block),
                    "MBInterpOp::define | Data corruption.");
            auto bounds = a_map.map().bounds(iter, dir);
#if PR_VERBOSE > 1
            pr_out() << "\t\t\tNeighbor is inside a normal boundary region: " << boundBox << " | number of boundaries: " << bounds.size() << std::endl;
#endif
            std::set<Point> fset;
            for (auto fi : a_footprint)
            {
                fset.insert(fi);
            }
            for (auto bound : bounds)
            {

                //handle points adjacent to triple points
                for (auto bi : boundBox)
                {
#if PR_VERBOSE > 1
                    pr_out() << "\t\t\tChecking point: " << bi << std::endl;
#endif
                    PR_TIME("add_non_triple_point");
                    MBPointID key(bi, iter.global());
#if PR_VERBOSE > 1
                    if (m_points.count(key))
                    {
                        pr_out() << "\t\t\t\tOperator already exists -> skipping" << std::endl;
                    }
#endif
                    if (m_points.count(key)) { continue; }
                    for (int ii = 0; ii < DIM; ii++)
                    {
                        PR_TIME("adjacent_triple_point_logic");
                        Point d2 = Point::Basis(ii);
                        bool onLoEdge = boundBox.edge(-d2).contains(bi);
                        bool onHiEdge = boundBox.edge( d2).contains(bi);
                        bool hasLoTriplePoint = patchDomain.inTriplePoint(patch_j - d2, block);
                        bool hasHiTriplePoint = patchDomain.inTriplePoint(patch_j + d2, block);
                        if (onLoEdge && hasLoTriplePoint)
                        {
#if PR_VERBOSE > 1
                            if (m_points.count(key))
                            {
                                pr_out() << "\t\t\t\tDetected adjacent triple point patch (lo): " << patch_j - d2 << std::endl;
                            }
#endif
                            for (auto fi : a_footprint)
                            {
                                fset.insert((fi + d2)*footprintDims);
                            }
                        }
                        if (onHiEdge && hasHiTriplePoint)
                        {
#if PR_VERBOSE > 1
                            if (m_points.count(key))
                            {
                                pr_out() << "\t\t\t\tDetected adjacent triple point patch (hi): " << patch_j + d2 << std::endl;
                            }
#endif
                            for (auto fi : a_footprint)
                            {
                                fset.insert((fi - d2)*footprintDims);
                            }
                        }
                    }
                    std::vector<Point> extFootprint(fset.begin(), fset.end());
#if PR_VERBOSE > 1
                    pr_out() << "\t\t\t\tExtended Footprint: " << std::endl;
                    for (auto fi : extFootprint)
                    {
                        pr_out() << "\t\t\t\t" << fi;
                        if (std::find(a_footprint.begin(), a_footprint.end(), fi) == a_footprint.end())
                        {
                            pr_out() << " <- added";
                        }
                        pr_out() << std::endl;
                    }
#endif
                    MBDataPoint dstDataPoint(iter, bi, layout);
                    auto op = std::make_shared<MBPointInterpOp>(dstDataPoint, m_ghost, a_map, extFootprint, a_order);
                    m_ops.push_back(op);
                    m_points.insert(key);
                    for (int dd = 0; dd < DIM; dd++)
                    {
                        PR_TIME("copy_along_axis");
                        if (m_copyAxes[dd])
                        {
                            int colHeight = boundBox.sizes()[dd];
                            Box boundColumn = Box(bi,bi).grow(Point::Basis(dd,colHeight-1));
                            boundColumn &= boundBox;
                            for (auto bj : boundColumn)
                            {
                                if (bi == bj) {continue; }
                                MBPointID key_j(bj, iter.global());
                                if (m_points.count(key_j)) {continue; }
                                Point shift = bj - bi;
                                MBDataPoint dst_j = dstDataPoint.shift(shift);
                                std::vector<MBDataPoint> srcs_j;
                                for (auto si : op->sources())
                                {
                                    srcs_j.push_back(si.shift(shift));
                                }
                                
                                auto op_j = std::make_shared<MBPointInterpOp>(dst_j, srcs_j, op->coefs());
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
        auto block = a_dst.layout().block(op->target().index());
        op->apply(a_dst, a_src);
    }
}

/*
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
*/

void MBInterpOp::print() const
{
    pr_out() << "Number of PointInterpOps: " << m_ops.size() << std::endl;
    for (auto& op : m_ops)
    {
        op->print(); 
    }
}

void MBInterpOp::copyAxis(unsigned int a_axis, bool a_flag)
{
    m_copyAxes[a_axis] = a_flag;
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






