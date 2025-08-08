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
        const MBInterpLayout& a_interpLayout,
        unsigned int a_order)
{
    PR_TIME("MBPointInterpOp::constructor");
    m_order = a_order;
    m_dst = a_dst;
    Box patchBox = a_map.layout()[m_dst.dstIndex()];
    auto footprint = a_interpLayout.footprint(m_dst.dstPoint(), m_dst.dstIndex());
    m_srcs = a_interpLayout.computeSourcePoints(footprint, a_map, m_dst.dstPoint(), m_dst.dstIndex());
    computeExponents();
    validateMatrixDimensions(footprint);
    auto momentData = computeMoments(a_map);
    try {
        computeCoefficients(momentData);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        writeLevelFootprint(a_map, "ERROR_FOOTPRINT");
        std::abort();
    }
    
}

// void MBPointInterpOp::computeSourcePoints(
//     const std::set<Point>& footprint,
//     const MBDisjointBoxLayout& layout,
//     Point ghostSizes)
// {
//     PR_TIME("MBPointInterpOp::computeSourcePoints");
//     auto index = m_dst.dstIndex();
//     auto center = m_dst.dstPoint();
//     Box patchBox = layout[index];
//     patchBox = patchBox.grow(ghostSizes);
//     auto patch = layout.patch(index);
//     Box domainBox = layout.domain().getBlock(patch.block).box();
//     m_srcs.clear();

//     for (auto fi : footprint)
//     {
//         Point pi = center + fi;
//         if (domainBox.containsPoint(pi) && patchBox.containsPoint(pi))
//         {
//             MBDataPoint data(index, pi, layout);
//             m_srcs.push_back(data);
//         } else if (!domainBox.containsPoint(pi))
//         {
//             Point dir = layout[index].whichBoundaryContains(pi);
//             auto adjIndices = layout.adjacentIndices(patch, dir);
//             for (auto adjIndex : adjIndices)
//             {
//                 BlockIndex adjBlock = layout.block(adjIndex);
//                 MBDataPoint data(index, pi, layout, dir, adjBlock);
//                 m_srcs.push_back(data);
//             }
//         }
//     }
//     // std::cout << "dst point: " << m_dst.point() <<" | patch: " << layout.point(m_dst.index()) << " | block: " << m_dst.dstBlock() << std::endl;
//     // for (auto si : m_srcs)
//     // {
//     //     std::cout << "\tsrc point: " << si.point() << " | block: " << si.srcBlock() << std::endl;
//     // }
//     // std::cout << std::endl;
//     PROTO_ASSERT(m_srcs.size() > 0,
//         "MBPointInterpOp::computeSourcePoints | Error: Could not compute any source data points for interpolation to point %s in block %i",
//         m_dst.dstPoint().str(), m_dst.dstBlock());
// }

void MBPointInterpOp::computeExponents()
{
    PR_TIME("MBPointInterpOp::computeExponents");
    Box sourceSpanBox;
    for (auto si : m_srcs)
    {
        sourceSpanBox &= si.dstPoint();
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
    Box B0(m_dst.dstPoint(), m_dst.dstPoint());
    momentDomains[m_dst.dstBlock()] = B0;

    for (auto si : m_srcs)
    {
        int bi = si.srcBlock();
        if (momentDomains.find(bi) == momentDomains.end())
        {
            momentDomains[bi] = Box(si.dstPoint(), si.dstPoint());
        } else {
            momentDomains[bi] &= si.dstPoint();
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
            C(ii,jj) = momentData[srcBlock][jj](src.dstPoint());
        }
    }

    Matrix<double> D(1,N);
    auto dstBlock = m_dst.dstBlock();
    for (int jj = 0; jj < N; jj++)
    {
        PR_TIMERS("build_D_matrix");
        D(0,jj) = momentData[dstBlock][jj](m_dst.dstPoint());
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
    double error = abs(coefSum - 1);
    //TODO: This check should be much more strict
    if (error > 1e-1)
    {
        #if PR_VERBOSE > 1
        pr_out() << "Built operator for point " << m_dst << std::endl;
        for (int ii = 0; ii < m_srcs.size(); ii++)
        {
            pr_out() << "\tcoef: " << m_coefs[ii] << " | src: " << m_srcs[ii] << std::endl;
        }
        pr_out() << "\tCoefficients sum to " << coefSum << std::endl;
        #endif
        char errorMsg[100];
        sprintf(errorMsg, "Coefficients must sum to 1. Error: %.2e for target point %s in block %i.", error, m_dst.dstPoint().str().c_str(), m_dst.dstBlock());
        throw std::runtime_error(std::string(errorMsg));
    }
    PROTO_ASSERT(error < 1e-1,
        "MBPointInterpOp::computeCoeffcients | Error: Coefficients must sum to 1. Error: %.2e for target point %s in block %i. See pout for details.",
        error, m_dst.dstPoint().str().c_str(), m_dst.dstBlock());
    
    
}

template<typename T, unsigned int C, MemType MEM>
void MBPointInterpOp::apply(
        MBLevelBoxData<T, C, MEM>& a_dst,
        MBLevelBoxData<T, C, MEM>& a_src)
{
    PR_TIME("MBPointInterpOp::apply");
    auto dst = a_dst[m_dst]; //dst is a Var<T, C, MEM>
    PROTO_ASSERT(a_dst[m_dst.dstIndex()].containsAddress(&dst(0)), "Data Corruption");
    for (int ii = 0; ii < size(); ii++)
    {
        auto src = a_src[m_srcs[ii]]; //src is a Var<T, C, MEM>
        for (int cc = 0; cc < C; cc++)
        {
            if (ii == 0) { dst(cc) = 0.0; }
            auto& patch = m_srcs[ii].srcPatch(a_src);
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
    B &= target().dstPoint();
    for (auto si : m_srcs)
    {
        B &= si.dstPoint();
    }
    BoxData<double,1> D(B);
    D.setVal(0);
    for (auto si : m_srcs)
    {
        D(si.dstPoint()) = 5 + si.srcBlock();
    }
    D(target().dstPoint()) = -1;
    HDF5Handler h5;
    h5.writePatch(D, filename);
}
template<typename MAP, MemType MEM>
void MBPointInterpOp::writeLevelFootprint(const MBLevelMap<MAP, MEM>& map, std::string filename) const 
{
    MBLevelMap<MAP, MEM> _map;
    _map.define(map.layout(), map.ghost()[1]);
    MBLevelBoxData<double, 1, HOST> data(map.layout(), map.ghost()[1]);
    data.setVal(0);
    data.setBoundary(0);
    for (auto si : m_srcs)
    {
        if (si.inBoundary())
        {
            data[si.srcIndex()](si.srcPoint()) = 1.0;
        } 
        data[si.dstIndex()](si.dstPoint()) = 1.0;
    }
    auto var = data[m_dst];
    var(0) = -1;
    HDF5Handler h5;
    h5.writeMBLevel(_map, data, filename);
}