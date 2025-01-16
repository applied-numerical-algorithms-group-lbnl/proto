
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
    m_order = a_order;
    m_dst = a_dst;
    m_srcs = a_map.map().interpFootprint(a_dst.point(), a_ghost, a_footprint, a_dst.index());
    
    computeExponents();
    validateMatrixDimensions(a_footprint);
    auto momentData = computeMoments(a_map);
    computeCoefficients(momentData);
}

void MBPointInterpOp::computeExponents()
{
    Box sourceSpanBox;
    for (auto si : m_srcs)
    {
        sourceSpanBox &= si.point();
    }
    Point maxExponents = Point::Ones(m_order);

    // If the span of the footprint only contains a single layer
    // of data along the ith axis, assume that the interpolating
    // polynomial is in a plane orthogonal to that axis
    for (int ii = 0; ii < DIM; ii++)
    {
        if (sourceSpanBox.size(ii) == 1) { maxExponents[ii] = 1; }
    }
    m_exponents.clear();
    for (auto bi : Box(maxExponents))
    {
        if (bi.sum() < m_order)
        {
            m_exponents.push_back(bi);
        }
    }
}

void MBPointInterpOp::validateMatrixDimensions(
    const std::vector<Point>& footprint) const
{
    int N = m_exponents.size();
    int M = m_srcs.size();
    if (N > M)
    {
        pr_out() << "====================================================================" << std::endl;
        pr_out() << "Underdetermined system detected" << std::endl;
        pr_out() << "\tN = " << N << ", M = " << M << std::endl;
        pr_out() << "\tDst Point: " << m_dst << std::endl;
        //pr_out() << "\tGhost: " << a_ghost << std::endl;
#if PR_VERBOSE > 1
        pr_out() << "\tFootprint: " << std::endl;
        for (auto fi : footprint)
        {
            pr_out() << "\t\t" << fi << std::endl;
        }
        pr_out() << "\tSrc Points: " << std::endl;
        for (auto si : m_srcs)
        {
            pr_out() << "\t\t" << si << std::endl;
        }
        pr_out() << "\tExponents: " << std::endl;
        for (auto ei : m_exponents)
        {
            pr_out() << "\t\t" << ei << std::endl;
        }
#endif
        pr_out() << "====================================================================" << std::endl;
    }
}

template<template<MemType> typename MAP, MemType MEM>
std::unordered_map<BlockIndex, std::vector<BoxData<double>>>
MBPointInterpOp::computeMoments(
    const MBLevelMap<MAP, MEM>& map) const
{
    // the radius in physical space of the neighborhood containing
    // the source data points
    auto xCenter = map.apply(m_dst);
    double sourceDataRadius = computeSourceDataRadius(map, xCenter);

    std::unordered_map<BlockIndex, Box> momentDomains;

    std::unordered_map<BlockIndex, std::vector<BoxData<double>>> momentData;

    // The domain of the block containing the destination point must
    // at minimum include that point
    Box B0(m_dst.point(), m_dst.point());
    momentDomains[m_dst.dstBlock()] = B0;

    for (auto si : m_srcs)
    {
        int bi = si.srcBlock();
        if (momentDomains.find(bi) == momentDomains.end())
        {
            momentDomains[bi] = Box(si.point(), si.point());
        } else {
            momentDomains[bi] &= si.point();
        }
    }

    auto C2C = Stencil<double>::CornersToCells(4);
    for (auto& [block, momentDomain] : momentDomains)
    {
        BoxData<double, DIM, MEM> Xi(C2C.domain(momentDomain).grow(1));
        BoxData<double, 1, MEM> Ji(Xi.box()); //Ji is unused
        map.doApply(Xi, Ji, block, m_dst.dstBlock()); 
        Xi -= xCenter;
        Xi /= sourceDataRadius;
        for (auto ei : m_exponents)
        {
            BoxData<double> moment = forall<double, 1>(f_binomialPower, Xi, ei);
            BoxData<double> momentAvg = C2C(moment);
            momentData[block].push_back(std::move(momentAvg));
        }
    }  

    return momentData;
}

template<template<MemType> typename MAP, MemType MEM>
double MBPointInterpOp::computeSourceDataRadius(
    const MBLevelMap<MAP, MEM>& map,
    Array<double, DIM>& xCenter) const
{
    double radius = 0.0;
    
    for (auto& src : m_srcs)
    {
        auto xSrc = map.apply(src);
        auto dist = xSrc - xCenter;
        radius += dist.norm();
    }
    radius /= (m_srcs.size());
    return radius;
}

void MBPointInterpOp::computeCoefficients(
    std::unordered_map<BlockIndex, std::vector<BoxData<double>>>& momentData)
{
    int N = m_exponents.size();
    int M = m_srcs.size();
    Matrix<double> C(M,N);
    for (int ii = 0; ii < M; ii++)
    {
        PR_TIMERS("build_C_matrix");
        auto& src = m_srcs[ii];
        auto srcBlock = src.srcBlock();
        for (int jj = 0; jj < N; jj++)
        {
            C(ii,jj) = momentData[srcBlock][jj](src.point());
        }
    }

    Matrix<double> D(1,N);
    auto dstBlock = m_dst.dstBlock();
    for (int jj = 0; jj < N; jj++)
    {
        PR_TIMERS("build_D_matrix");
        D(0,jj) = momentData[dstBlock][jj](m_dst.point());
    }
    
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

template<template<MemType> typename MAP, MemType MEM>
MBInterpOp<MAP,MEM>::MBInterpOp()
{
    PR_TIME("MBInterpOp<MAP,MEM>::constructor");
    m_copyAxes.fill(false);
}

template<template<MemType> typename MAP, MemType MEM>
MBInterpOp<MAP,MEM>::MBInterpOp(Point a_ghost, unsigned int a_order)
{
    PR_TIME("MBInterpOp<MAP,MEM>::constructor");
    m_order = a_order;
    m_ghost = a_ghost;
    m_copyAxes.fill(false);
}

template<template<MemType> typename MAP, MemType MEM>
MBInterpOp<MAP,MEM>::MBInterpOp(
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
void MBInterpOp<MAP,MEM>::define(
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
            MayDay<void>::Error("MBInterpOp<MAP,MEM>::define | Function undefined for order %i", a_order);
            break;
    }
    define(a_map, footprint, a_order, a_block);
}

template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp<MAP,MEM>::define(
        const MBLevelMap<MAP,MEM>&  a_map,
        std::vector<Point>          a_footprint,
        int                         a_order,
        int                         a_block)
{
    PR_TIME("MBInterpOp<MAP,MEM>::define");

    m_ghost = a_map.ghost()[0];
    m_map = a_map;
    m_baseFootprint = a_footprint;

    Box spanBox;
    for (auto fi : a_footprint)
    {
        spanBox &= fi;
    }
    m_footprintDims = Point::Zeros();
    for (int dd = 0; dd < DIM; dd++)
    {
        if (spanBox.sizes()[dd] > 1)
        {
            m_footprintDims[dd] = 1;
        }
    }

    const auto& layout = a_map.map().layout();
    const auto& patchDomain = layout.partition().domain();
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        Point patch_i = layout.point(iter);
        if ((a_block >= 0) && (block != a_block)) {continue; }
        Box blockDomainBox = layout.domain().getBlock(block).box(); 
        Box patchBox = layout[iter];

        for (auto dir : Box::Kernel(1))
        {
            Point patch_j = patch_i + dir; 

            if (patchDomain.inInterior(patch_j, block)) { continue; }
            if (patchDomain.inDomainBoundary(patch_j, block)) { continue; }
            Box boundBox = patchBox.adjacent(m_ghost*dir);
            if (patchDomain.inTriplePoint(patch_j, block))
            {
                for (auto bi : boundBox)
                {
                    PR_TIME("add_triple_point");
                    MBPointID key(bi, iter.global());

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
                    "MBInterpOp<MAP,MEM>::define | Data corruption.");
            auto bounds = a_map.map().bounds(iter, dir);

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

                    PR_TIME("add_non_triple_point");
                    MBPointID key(bi, iter.global());

                    if (m_points.count(key)) { continue; }
                    for (int ii = 0; ii < DIM; ii++)
                    {
                        PR_TIME("adjacent_triple_point_logic");
                        Point d2 = Point::Basis(ii);
                        bool onLoEdge = boundBox.edge(-d2).containsPoint(bi);
                        bool onHiEdge = boundBox.edge( d2).containsPoint(bi);
                        bool hasLoTriplePoint = patchDomain.inTriplePoint(patch_j - d2, block);
                        bool hasHiTriplePoint = patchDomain.inTriplePoint(patch_j + d2, block);
                        if (onLoEdge && hasLoTriplePoint)
                        {

                            for (auto fi : a_footprint)
                            {
                                fset.insert((fi + d2)*footprintDims);
                            }
                        }
                        if (onHiEdge && hasHiTriplePoint)
                        {

                            for (auto fi : a_footprint)
                            {
                                fset.insert((fi - d2)*footprintDims);
                            }
                        }
                    }
                    std::vector<Point> extFootprint(fset.begin(), fset.end());

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

template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp<MAP,MEM>::addPointOps(const MBIndex& index, Point dir)
{
    PatchID patch = m_map.layout().point(index);
    PatchID adjPatch = patch + dir;
    BlockIndex block = m_map.layout().block(index);
    auto& patchDomain = m_map.layout().patchDomain();

    if (patchDomain.inInterior(adjPatch, block)) { return; }
    else if (patchDomain.inDomainBoundary(adjPatch, block)) { return; }
    else if (patchDomain.inTriplePoint(adjPatch, block))
    {
        addTriplePointOps(index, dir);
    } else {
        PROTO_ASSERT(patchDomain.inBlockBoundary(adjPatch, block),
            "MBInterpOp<MAP,MEM>::addPointOps | Data corruption.");
        addBlockBoundaryOps(index, dir);
    }

    
}
template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp<MAP,MEM>::addTriplePointOps(const MBIndex& index, Point dir)
{
    Box patchBox = m_map.layout()[index];
    Box boundBox = patchBox.adjacent(m_ghost*dir);
    for (auto boundPoint : boundBox)
    {
        if (this->containsOp(boundPoint, index)) { continue; }
        auto footprint = getTriplePointFootprint(boundPoint, index, dir);
        addPointOp(boundPoint, index, footprint);
    }
}
template<template<MemType> typename MAP, MemType MEM>
std::vector<Point> MBInterpOp<MAP,MEM>::getTriplePointFootprint(
    Point centerPoint,
    const MBIndex& index,
    Point dir) const
{
    Box patchBox = m_map.layout()[index];
    Point corner = patchBox.high();
    Point positionFromCorner = (centerPoint - corner);
    int maxDistFromCorner = (positionFromCorner*m_footprintDims).abs().max();
    std::set<Point> footprint;
    for (int dd = 0; dd < DIM; dd++)
    {
        if (dir[dd] == 0) { continue; }
        Point sourcePositionFromCorner = 
            Point::Basis(dd, maxDistFromCorner*dir[dd]);
        Point sourcePositionFromCenter = 
            (corner + sourcePositionFromCorner - centerPoint);
        for (auto fi : m_baseFootprint)
        {
            footprint.insert((fi + sourcePositionFromCenter)*m_footprintDims);
        }
    }
    return std::vector<Point>(footprint.begin(), footprint.end());
}

template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp<MAP,MEM>::addBlockBoundaryOps(const MBIndex& index, Point dir)
{
    addTriplePointAdjacentOps(index, dir);
}


template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp<MAP,MEM>::addTriplePointAdjacentOps(const MBIndex& index, Point dir)
{
    PatchID patch = m_map.layout().point(index);
    PatchID boundPatch = patch + dir;
    auto& patchDomain = m_map.layout().patchDomain();
    Box patchBox = m_map.layout()[index];
    Box boundBox = patchBox.adjacent(m_ghost*dir);
    for (auto adjDir : Point::DirectionsOfCodim(1))
    {
        PatchID boundPatchAdj = boundPatch + adjDir;
        if (patchDomain.inTriplePoint(boundPatchAdj))
        {
            for (auto boundPoint : boundBox.edge(adjDir, 1))
            {

            }
        }
    }
}

template<template<MemType> typename MAP, MemType MEM>
bool MBInterpOp<MAP,MEM>::containsOp(Point dataPoint, const MBIndex& index) const
{
    MBPointID key(dataPoint, index.global());
    return (m_points.count(key) > 0);
}

template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp<MAP,MEM>::addPointOp(
    Point dataPoint,
    const MBIndex& index,
    const std::vector<Point>& footprint)
{
    PROTO_ASSERT(!containsOp(dataPoint, index),
        "MBInterpOp::adPointOp | Error: Trying to add duplicate op");
    
    MBDataPoint dstDataPoint(index, dataPoint, m_map.layout());
    auto op = std::make_shared<MBPointInterpOp>(
        dstDataPoint, m_ghost, m_map, footprint, m_order);
    
    m_ops.push_back(op);
    MBPointID key(dataPoint, index.global());
    m_points.insert(key);

    copyAlongSymmetryAxes(*op);
}

template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp<MAP,MEM>::copyAlongSymmetryAxes(const MBPointInterpOp& pointOp)
{
    Point dataPoint = pointOp.target().point();
    MBIndex index = pointOp.target().index();
    Box patchBox = m_map.layout()[index];
    for (int dd = 0; dd < DIM; dd++)
    {
        PR_TIME("MBInterpOp::copyAlongSymmetryAxis");

        if (!m_copyAxes[dd]) { continue; }

        Point low = dataPoint;
        Point high = dataPoint;
        low[dd] = patchBox.low()[dd];
        high[dd] = patchBox.high()[dd];
        Box axisBox(low, high);

        for (auto axisPoint : axisBox)
        {
            if (containsOp(axisPoint, index)) { continue; }
            Point shift = axisPoint - dataPoint;
            addShiftedOp(pointOp, shift);
        }
    }
}

template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp<MAP,MEM>::addShiftedOp(
    const MBPointInterpOp& baseOp,
    Point shift)
{
    PROTO_ASSERT(!containsOp(dataPoint, index),
        "MBInterpOp::adPointOp | Error: Trying to add duplicate op");

    MBDataPoint shiftDst = pointOp.target().shift(shift);
    std::vector<MBDataPoint> shiftSrc;
    for (auto si : pointOp.sources())
    {
        shiftSrc.push_back(si.shift(shift));
    }
    auto coefs = pointOp.coefs();
    auto shiftedOp = std::make_shared<MBPointInterpOp>(shiftDst, shiftSrc, coefs);
    
    m_ops.push_back(shiftedOp);
    MBPointID key(shiftDst.point(), index.global());
    m_points.insert(key);
}

template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp<MAP,MEM>::insert(std::shared_ptr<MBPointInterpOp> a_op)
{
    Point dst = a_op->target().point();
    int globalIndex = a_op->target().index();
    MBPointID key(dst, globalIndex);
    PROTO_ASSERT(m_points.count(key) == 0,
            "MBInterpOp<MAP,MEM>::insert | Tried to insert a duplicate MBPointInterpOp");

    m_ops.push_back(a_op);
    m_points.insert(key);
}
template<template<MemType> typename MAP, MemType MEM>
template<typename T, unsigned int C, MemType MEM>
void MBInterpOp<MAP,MEM>::apply(
        MBLevelBoxData<T, C, MEM>& a_dst,
        MBLevelBoxData<T, C, MEM>& a_src)
{
    PR_TIME("MBInterpOp<MAP,MEM>::apply");
    for (auto& op : m_ops)
    {
        auto block = a_dst.layout().block(op->target().index());
        op->apply(a_dst, a_src);
    }
}
template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp<MAP,MEM>::print() const
{
    pr_out() << "Number of PointInterpOps: " << m_ops.size() << std::endl;
    for (auto& op : m_ops)
    {
        op->print(); 
    }
}
template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp<MAP,MEM>::copyAxis(unsigned int a_axis, bool a_flag)
{
    m_copyAxes[a_axis] = a_flag;
}
template<template<MemType> typename MAP, MemType MEM>
template<typename T, unsigned int C, MemType MEM>
void MBInterpOp<MAP,MEM>::printErrorPoints(
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






