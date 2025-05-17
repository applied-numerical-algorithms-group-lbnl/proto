
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

template<typename MAP, MemType MEM>
MBPointInterpOp::MBPointInterpOp(
        MBDataPoint a_dst,
        Point a_ghost,
        const MBLevelMap<MAP, MEM>& a_map,
        std::set<Point> a_footprint,
        unsigned int a_order)
{
    PR_TIME("MBPointInterpOp::constructor");
    m_order = a_order;
    m_dst = a_dst;
    computeSourcePoints(a_footprint, a_map.layout(), a_ghost);
    computeExponents();
    validateMatrixDimensions(a_footprint);
    auto momentData = computeMoments(a_map);
    computeCoefficients(momentData);
}

void MBPointInterpOp::computeSourcePoints(
    const std::set<Point>& footprint,
    const MBDisjointBoxLayout& layout,
    Point ghostSizes)
{
    PR_TIME("MBPointInterpOp::computeSourcePoints");
    auto index = m_dst.index();
    auto center = m_dst.point();
    Box patchBox = layout[index];
    patchBox = patchBox.grow(ghostSizes);
    auto patch = layout.patch(index);
    Box domainBox = layout.domain().getBlock(patch.block).box();
    m_srcs.clear();

    for (auto fi : footprint)
    {
        Point pi = center + fi;
        if (domainBox.containsPoint(pi) && patchBox.containsPoint(pi))
        {
            MBDataPoint data(index, pi, layout);
            m_srcs.push_back(data);
        } else if (!domainBox.containsPoint(pi))
        {
            Point dir = layout[index].whichBoundaryContains(pi);
            auto adjIndices = layout.adjacentIndices(patch, dir);
            for (auto adjIndex : adjIndices)
            {
                BlockIndex adjBlock = layout.block(adjIndex);
                MBDataPoint data(index, pi, layout, dir, adjBlock);
                m_srcs.push_back(data);
            }
        }
    }
    // std::cout << "dst point: " << m_dst.point() << " | block: " << m_dst.dstBlock() << std::endl;
    // for (auto si : m_srcs)
    // {
    //     std::cout << "\tsrc point: " << si.point() << " | block: " << si.srcBlock() << std::endl;
    // }
    // std::cout << std::endl;
    PROTO_ASSERT(m_srcs.size() > 0,
        "MBPointInterpOp::computeSourcePoints | Error: Could not compute any source data points for interpolation to point %s in block %i",
        m_dst.point().str(), m_dst.dstBlock());
}

void MBPointInterpOp::computeExponents()
{
    PR_TIME("MBPointInterpOp::computeExponents");
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
    const std::set<Point>& footprint) const
{
    PR_TIME("MBPointInterpOp::validateMatrixDimensions");
    int N = m_exponents.size();
    int M = m_srcs.size();
    if (N == 0)
    {
        std::cout << "MBInterpOp | Invalid system detected: N = 0 (see pout for details)" << std::endl;
        pr_out() << "====================================================================" << std::endl;
        pr_out() << "Invalid system detected" << std::endl;
        pr_out() << "\tNo exponents were computed" << std::endl;
        pr_out() << "\tN = " << N << ", M = " << M << std::endl;
        pr_out() << "\tDst Point: " << m_dst << std::endl;
        pr_out() << "====================================================================" << std::endl;
    }
    if (M == 0)
    {
        std::cout << "MBInterpOp | Invalid system detected: M = 0 (see pout for details)" << std::endl;
        pr_out() << "====================================================================" << std::endl;
        pr_out() << "Invalid system detected" << std::endl;
        pr_out() << "\tNo source points were found" << std::endl;
        pr_out() << "\tN = " << N << ", M = " << M << std::endl;
        pr_out() << "\tDst Point: " << m_dst << std::endl;
        pr_out() << "====================================================================" << std::endl;
    }
    if (N > M)
    {
        std::cout << "MBInterpOp | Undetermined system detected: N > M (see pout for details)" << std::endl;
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

template<typename MAP, MemType MEM>
std::unordered_map<BlockIndex, std::vector<BoxData<double>>>
MBPointInterpOp::computeMoments(
    const MBLevelMap<MAP, MEM>& map) const
{
    PR_TIME("MBPointInterpOp::computeMoments");
    

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
    
    // the radius in physical space of the neighborhood containing
    // the source data points
    auto xCenter = map.apply(m_dst);
    double sourceDataRadius = computeSourceDataRadius(map, xCenter);
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

template<typename MAP, MemType MEM>
double MBPointInterpOp::computeSourceDataRadius(
    const MBLevelMap<MAP, MEM>& map,
    Array<double, DIM>& xCenter) const
{
    PR_TIME("MBPointInterpOp::computeSourceDataRadius");
    double radius = 0.0;
    for (auto& src : m_srcs)
    {
        auto xSrc = map.apply(src);
        auto dist = xSrc - xCenter;
        radius += dist.norm();
    }
    radius /= (m_srcs.size());
    PROTO_ASSERT(!(radius == 0 || std::isnan(radius) || std::isinf(radius)),
        "MBPointInterpOp::computeSourceDataRadius | Error: Invalid radius");

    return radius;
}

void MBPointInterpOp::computeCoefficients(
    std::unordered_map<BlockIndex, std::vector<BoxData<double>>>& momentData)
{
    PR_TIME("MBPointInterpOp::computeCoefficients");
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
    //#if PR_DEBUG
    double coefSum = std::accumulate(m_coefs.begin(), m_coefs.end(), 0.0);
    if (abs(coefSum - 1.0) > 1e-12)
    {
        std::cout << "debug handle" << std::endl;
    }
    PROTO_ASSERT(abs(coefSum - 1.0) < 1e-12,
        "MBPointInterpOp::computeCoeffcients | Error: Coefficients must sum to 1. Observed: %f for target point %s in block %i",
        coefSum, m_dst.point().str().c_str(), m_dst.dstBlock());
    //TODO: Check for coefficients with values greater than some value
    //#endif
}

template<typename T, unsigned int C, MemType MEM>
void MBPointInterpOp::apply(
        MBLevelBoxData<T, C, MEM>& a_dst,
        MBLevelBoxData<T, C, MEM>& a_src)
{
    PR_TIME("MBPointInterpOp::apply");
    auto dst = a_dst[m_dst]; //dst is a Var<T, C, MEM>
    PROTO_ASSERT(a_dst[m_dst.index()].containsAddress(&dst(0)), "Data Corruption");
    for (int ii = 0; ii < size(); ii++)
    {
        
        auto src = a_src[m_srcs[ii]]; //src is a Var<T, C, MEM>
        for (int cc = 0; cc < C; cc++)
        {
            if (ii == 0) { dst(cc) = 0.0; }
            auto& patch = m_srcs[ii].patch(a_src);
            PROTO_ASSERT(patch.containsAddress(&src(cc)), "Data Corruption");
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
void MBPointInterpOp::writeFootprint(std::string filename) const 
{
    Box B;
    B &= target().point();
    for (auto si : m_srcs)
    {
        B &= si.point();
    }
    BoxData<double,1> D(B);
    D.setVal(0);
    for (auto si : m_srcs)
    {
        D(si.point()) = 1 + si.srcBlock();
    }
    D(target().point()) = -1;
    HDF5Handler h5;
    h5.writePatch(D, filename);
}
MBInterpOp::MBInterpOp()
{
    PR_TIME("MBInterpOp::constructor");
}

MBInterpOp::MBInterpOp(Point a_ghost, unsigned int a_order)
{
    PR_TIME("MBInterpOp::constructor");
    m_order = a_order;
    m_ghost = a_ghost;
}

template<typename MAP, MemType MEM>
MBInterpOp::MBInterpOp(
        const MBLevelMap<MAP,MEM>&    a_map,
        int a_order)
{
    PR_TIME("MBInterpOp::constructor");
    define(a_map, a_order);
}

template<typename MAP, MemType MEM>
void MBInterpOp::define(
        const MBLevelMap<MAP,MEM>&    a_map,
        int a_order)
{
    MBInterpLayout interpLayout(a_map.layout(), a_order);
    define(a_map, interpLayout, a_order);
}

template<typename MAP, MemType MEM>
void MBInterpOp::define(
        const MBLevelMap<MAP,MEM>&  a_map,
        MBInterpLayout              a_interpLayout,
        int                         a_order)
{
    PR_TIME("MBInterpOp::define");
    PROTO_ASSERT(a_interpLayout.layout().compatible(a_map.layout()),
        "MBInterpOp::define | Error: Incompatible layouts among inputs");

    m_order = a_order;
    m_ghost = a_map.ghost()[0];
    m_interpLayout = a_interpLayout;
    for (auto iter : m_interpLayout.layout())
    {
        for (auto dir : Box::Kernel(1))
        {
           addPointOps(iter, dir, a_map);
        }
    }
}

template<typename MAP, MemType MEM>
void MBInterpOp::define(
        const std::vector<std::shared_ptr<MBLevelMap<MAP,MEM>>>&  a_maps,
        MBInterpLayout              a_interpLayout,
        int                         a_order)
{
    PR_TIME("MBInterpOp::define");
    PROTO_ASSERT(a_maps.size() == a_interpLayout.layout().numBlocks(),
        "MBInterpOp::define | Error: Incorrect number of maps. There must be one map per block");
    for (auto bi = 0; bi < a_interpLayout.layout().numBlocks(); bi++)
    {
        PROTO_ASSERT(a_interpLayout.layout().compatible(a_maps[bi]->layout()),
            "MBInterpOp::define | Error: Incompatible layouts among inputs");
        PROTO_ASSERT(a_maps[bi]->ghost() == a_maps[0]->ghost(),
            "MBInterpOp::define | Error: Input maps have different ghost regions");
    }

    m_order = a_order;
    m_ghost = a_maps[0]->ghost()[0];
    m_interpLayout = a_interpLayout;
    for (auto iter : m_interpLayout.layout())
    {
        auto block = m_interpLayout.layout().block(iter);
        auto map = a_maps[block];
        for (auto dir : Box::Kernel(1))
        {
           addPointOps(iter, dir, *map);
        }
    }
}

template<typename MAP, MemType MEM>
void MBInterpOp::addPointOps(const MBIndex& index, Point dir, const MBLevelMap<MAP,MEM>& map)
{
    PR_TIME("MBInterpOp::addPointOps");
    auto& layout = m_interpLayout.layout();
    auto patch = layout.patch(index);
    PatchID adjPatch = patch.point + dir;
    BlockIndex block = patch.block;
    auto& patchDomain = layout.patchDomain();

    if (patchDomain.isPointInInterior(adjPatch, block)) { return; }
    else if (patchDomain.isPointInDomainBoundary(adjPatch, block)) { return; }
    else if (layout.isRefinementBoundary(patch, dir)) { return; }
    else {
        Box patchBox = layout[index];
        Box boundBox = patchBox.adjacent(m_ghost * dir);
        for (auto boundPoint : boundBox)
        {
            if (this->containsOp(boundPoint, index))
            {
                continue;
            }
            auto footprint = m_interpLayout.footprint(boundPoint, dir, index);
            addPointOp(boundPoint, index, footprint, map);
        }
    }
}

bool MBInterpOp::containsOp(Point dataPoint, const MBIndex& index) const
{
    PR_TIME("MBInterpOp::containsOp");
    MBPointID key(dataPoint, index.global());
    return (m_points.count(key) > 0);
}

template<typename MAP, MemType MEM>
void MBInterpOp::addPointOp(
    Point dataPoint,
    const MBIndex& index,
    const std::set<Point>& footprint,
    const MBLevelMap<MAP,MEM>& map)
{
    PR_TIME("MBInterpOp::addPointOp");
    PROTO_ASSERT(!containsOp(dataPoint, index),
        "MBInterpOp::adPointOp | Error: Trying to add duplicate op");
    
    MBDataPoint dstDataPoint(index, dataPoint, m_interpLayout.layout());
    auto op = std::make_shared<MBPointInterpOp>(
        dstDataPoint, m_ghost, map, footprint, m_order);
    
    m_ops.push_back(op);
    MBPointID key(dataPoint, index.global());
    m_points.insert(key);

    copyAlongSymmetryAxes(*op);
}

void MBInterpOp::copyAlongSymmetryAxes(const MBPointInterpOp& pointOp)
{
    PR_TIME("MBInterpOp::copyAlongSymmetryAxis");
    Point dataPoint = pointOp.target().point();
    MBIndex index = pointOp.target().index();
    Box patchBox = m_interpLayout.layout()[index];
    Point boundaryDir = patchBox.whichBoundaryContains(dataPoint);
    Box patchBoundBox = patchBox.adjacent(boundaryDir * m_ghost);
    for (int dd = 0; dd < DIM; dd++)
    {
        PR_TIME("MBInterpOp::copyAlongSymmetryAxis");

        if (!m_interpLayout.isCopied(dd)) { continue; }

        Point low = dataPoint;
        Point high = dataPoint;
        low[dd] = patchBoundBox.low()[dd];
        high[dd] = patchBoundBox.high()[dd];
        Box axisBox(low, high);
        PROTO_ASSERT(axisBox.containsPoint(dataPoint),
            "MBInterpOp::copyAlongSymmetryAxes | Error: Data Corruption");
        for (auto axisPoint : axisBox)
        {
            if (containsOp(axisPoint, index)) { continue; }
            Point shift = axisPoint - dataPoint;
            addShiftedOp(pointOp, shift);
        }
    }
}

void MBInterpOp::addShiftedOp(
    const MBPointInterpOp& baseOp,
    Point shift)
{
    PR_TIME("MBInterpOp::addShiftedOp");
    PROTO_ASSERT(containsOp(baseOp.target().point(), baseOp.target().index()),
        "MBInterpOp::adPointOp | Error: Trying to add duplicate op");
    if (shift == Point::Zeros()) { return; }
    MBDataPoint shiftDst = baseOp.target().shift(shift);
    std::vector<MBDataPoint> shiftSrc;
    for (auto si : baseOp.sources())
    {
        shiftSrc.push_back(si.shift(shift));
    }
    auto coefs = baseOp.coefs();
    auto shiftedOp = std::make_shared<MBPointInterpOp>(shiftDst, shiftSrc, coefs);
    
    m_ops.push_back(shiftedOp);
    MBPointID key(shiftDst.point(), baseOp.target().index().global());
    m_points.insert(key);
}

void MBInterpOp::insert(std::shared_ptr<MBPointInterpOp> a_op)
{
    PR_TIME("MBInterpOp::insert");
    Point dst = a_op->target().point();
    int globalIndex = a_op->target().index();
    MBPointID key(dst, globalIndex);
    PROTO_ASSERT(m_points.count(key) == 0,
            "MBInterpOp::insert | Tried to insert a duplicate MBPointInterpOp");

    m_ops.push_back(a_op);
    m_points.insert(key);
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
void MBInterpOp::print() const
{
    pr_out() << "Number of PointInterpOps: " << m_ops.size() << std::endl;
    for (auto& op : m_ops)
    {
        op->print(); 
    }
}
void MBInterpOp::writeFootprint(std::string filename) const
{
    std::map<int, int> counters;
    for (auto& op : m_ops)
    {
        int block = op->target().dstBlock();
        int ii = counters[block]++;
        std::string suffix = ("B" + std::to_string(block) + "_I" + std::to_string(ii));
        op->writeFootprint(filename+suffix);
    }
}
// void MBInterpOp::copyAxis(unsigned int a_axis, bool a_flag)
// {
//     m_copyAxes[a_axis] = a_flag;
// }
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


template<typename MAP, typename T, unsigned int C, MemType MEM, Centering CTR>
void interpBoundaries(MBLevelBoxData<T, C, MEM, CTR>& a_data, unsigned int a_order)
{
    MBLevelMap<MAP, MEM> map(a_data.layout(), a_data.ghost());
    interpBoundaries<MAP, T, C, MEM, CTR>(a_data, map, a_order);
}

template<typename MAP, typename T, unsigned int C, MemType MEM, Centering CTR>
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






