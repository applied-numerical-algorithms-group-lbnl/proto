
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

template<MemType MEM>
MBPointInterpOp::MBPointInterpOp(
        MBDataPoint a_dst,
        Point a_ghost,
        MBLevelMap<MEM>& a_map,
        const std::vector<Point>& a_footprint,
        unsigned int a_order)
{
    PR_TIME("MBPointInterpOp::constructor");
#if PR_VERBOSE > 1
    pout() << "DST_POINT: " << a_dst.point << " | BLOCK: " << a_dst.dstBlock() << std::endl;
#endif
    m_dst = a_dst;
    m_srcs = a_map.map().interpFootprint(a_dst.point, a_ghost, a_footprint, a_dst.index);
#if PR_VERBOSE > 1
    pout() << "SRC_POINTS: " << std::endl;
    for (auto si : m_srcs)
    {
        pout() << "\tpoint: " << si.point << " | block: " << si.srcBlock() << std::endl;
    }
    pout() << std::endl;
#endif

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
#if PR_VERBOSE > 1
    pout() << "SPAN_BOX: " << spanBox << std::endl;
    pout() << "EXPONENTS: ";
    for (auto ei : m_exponents)
    {
        pout() << ei << ", ";
    }
    pout() << std::endl;
#endif
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
#if PR_VERBOSE > 1
    pout() << "Xg: " << xg << " | point: " << a_dst.point << std::endl;
#endif
    for (int ii = 0; ii < M; ii++)
    {
        auto xj = a_map.cellCentered(m_srcs[ii]);
        auto dist = xj - xg;
#if PR_VERBOSE > 1
        pout() << "\tXj: " << xj << " | point: " << m_srcs[ii].point << " | dist: " << dist << std::endl;
#endif
        Rg += dist.norm();
    }
    Rg /= M;
    Rg = 1.0; // debugging
#if PR_VERBOSE > 1
    pout() << "Rg: " << Rg << std::endl;
#endif
   
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
#if PR_VERBOSE > 1
    pout() << "Cached Data Regions: " << std::endl;
#endif
    for (auto& [block, box] : srcBoxes)
    {
#if PR_VERBOSE > 1
        pout() << "block: " << block << " | box: " << box << std::endl;
#endif
        BoxData<double, DIM, MEM> Xi(C2C.domain(box).grow(1));
        BoxData<double, 1, MEM> Ji(box.grow(1));
        a_map.doApply(Xi, Ji, block, dstBlock); 
#if PR_VERBOSE > 1
        pout() << "X: " << std::endl;
        Xi.printData();
        pout() << "J: " << std::endl;
        Ji.printData();
#endif
        Xi -= xg;
        Xi /= Rg;
#if PR_VERBOSE > 1
        pout() << "X (normalized): " << std::endl;
        Xi.printData();
#endif
        for (int jj = 0; jj < N; jj++)
        {
            auto Xi_p = forall<double, 1>(f_binomialPower, Xi, m_exponents[jj]);
            BoxData<double> Xi_p_avg = C2C(Xi_p);
            Ji.setVal(1);
            srcData[block].push_back(Operator::cellProduct(Ji, Xi_p_avg));
        
#if PR_VERBOSE > 1
            pout() << "moment " << m_exponents[jj] << std::endl;
            Xi_p_avg.printData();
            pout() << "moment " << m_exponents[jj] << " * J" << std::endl;
            srcData[block][srcData[block].size()-1].printData();
#endif
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
/*
template<MemType MEM>
MBInterpOp::MBInterpOp(
        MBLevelMap<MEM>& a_map,
        const std::vector<Point>& a_footprint,
        Point a_ghost,
        unsigned int a_order)
{
    PR_TIME("MBInterpOp::constructor");
    m_footprint = a_footprint;
    m_exponents.clear();
    for (auto bi : Box::Cube(a_order))
    {
        if (bi.sum() < a_order - 1)
        {
            m_exponents.push_back(bi);
        }
    }
    m_order = a_order;
    m_ghost = a_ghost;
    m_ops.clear();
    const auto& layout = a_map.map().layout();
    const auto& domain = layout.domain();
    const auto& graph = layout.domain().graph();
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        Box domainBox = layout.domain().blockDomain(block).box(); 
        for (auto dir : Box::Kernel(1))
        {
           if (dir == Point::Zeros()) {continue;}
           if (a_map.map().bounds(iter, dir).size() > 0)
           {
               Box patchBox = layout[iter];
               Box ghostBox = patchBox.adjacent(dir*a_ghost);
               for (auto gi : ghostBox)
               {
                   MBDataPoint p(iter, gi, layout);
                   m_ops.push_back(MBPointInterpOp(p, a_ghost, a_map, a_footprint, a_order));
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
*/






