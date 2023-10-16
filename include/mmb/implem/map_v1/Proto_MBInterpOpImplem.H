
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

template<typename Func>
MBPointInterpOp::MBPointInterpOp(
        MBDataPoint a_dst,
        Point a_ghost,
        MBMap<Func>& a_map,
        const std::vector<Point>& a_footprint,
        unsigned int a_order)
{
    PR_TIME("MBPointInterpOp::constructor");
#if PR_VERBOSE > 1
    pr_out() << "DST_POINT: " << a_dst.point << " | BLOCK: " << a_dst.dstBlock() << std::endl;
#endif
    m_dst = a_dst;
    m_srcs = a_map.map().interpFootprint(a_dst.point, a_ghost, a_footprint, a_dst.index);


#if PR_VERBOSE > 1
    pr_out() << "SRC_POINTS: ";
    for (auto si : m_srcs)
    {
        pr_out() << "\tpoint: " << si.point << " | block: " << si.srcBlock() << std::endl;
    }
    pr_out() << std::endl;
#endif
    // References
    auto& layout = a_map.map().layout();
    unsigned int dstBlock = a_dst.dstBlock();
    
    // Compute footprint span
#if 0
    Box spanBox;
    for (auto si : m_srcs)
    {
        spanBox &= si.point;
    }
    Point ignoreDims = Point::Zeros();
    for (int ii = 0; ii < DIM; ii++)
    {
        if (spanBox.size(ii) == 1) {ignoreDims[ii] = 1; }
    }
#endif
    // Compute exponents
    m_exponents.clear();
    //for (auto bi : Box((a_order-1)*ignoreDims))
    for (auto bi : Box::Cube(a_order))
    {
        if (bi.sum() < a_order)
        {
            m_exponents.push_back(bi);
        }
    }
    int N = m_exponents.size();
    int M = m_srcs.size();

    // Rg = 1/M*sum_j(sqrt(<xj - xg>^2))
    // Cjp = <(x-xg)^p/Rg>_j  |  p = binomial, j = source cell
    // Mp = <(x-xg)^p/Rg>_g   |  p = binomial
    // Sj = sum_p(Mp*Cinv_pj = m_coefs[j]
    
    // Compute Rg
    double Rg = 0.0;
    auto xg = a_map(a_dst);
    for (int ii = 0; ii < M; ii++)
    {
        auto xj = a_map(m_srcs[ii]);
        auto dist = xj - xg;
        Rg += dist.norm();
    }
    Rg /= M;
   
    // Cache Moment Data
    auto C2C = Stencil<double>::CornersToCells(4);
    std::unordered_map<int, Box> srcBoxes;
    std::unordered_map<int, std::vector<BoxData<double>>> srcData; 
    std::unordered_map<int, BoxData<double>*> srcJ;
    Box B0 = C2C.domain(Box::Kernel(1).shift(a_dst.point));
    srcBoxes[a_dst.dstBlock()] = B0;
    srcJ[a_dst.dstBlock()] = &a_dst.patch(a_map.jacobian());
    for (auto si : m_srcs)
    {
        int bi = si.srcBlock();
        auto B = Box::Kernel(1).shift(si.point);
        B = C2C.domain(B);
        if (srcBoxes.find(bi) == srcBoxes.end())
        {
            srcBoxes[bi] = B; 
            srcJ[bi] = &si.patch(a_map.jacobian());
        } else {
            srcBoxes[bi] &= B.low();
            srcBoxes[bi] &= B.high();
        }
    }
    int ind = 0;
    for (auto& [block, box] : srcBoxes)
    {
        auto Xi = a_map(box, block, dstBlock, PR_NODE);
        Xi -= xg;
        Xi /= Rg;
#if PR_VERBOSE==3
        if (dstBlock == 0)
        {
            pr_out() << "Cached Jacobian for data from block " << block << " to block " << dstBlock << std::endl;
            srcJ[block]->printData(4);
        }
#endif
        for (int jj = 0; jj < N; jj++)
        {
            auto Xi_p = forall<double, 1>(f_binomialPower, Xi, m_exponents[jj]);
            BoxData<double> Xi_p_avg = C2C(Xi_p);
            //BoxData<double> ones(srcJ[block]->box());
            //ones.setVal(1);
            srcData[block].push_back(Operator::cellProduct(*srcJ[block], Xi_p_avg));
            //srcData[block].push_back(Operator::cellProduct(ones, Xi_p_avg));
        }
        ind++;
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

template<typename Func>
MBInterpOp::MBInterpOp(
        MBMap<Func>& a_map,
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







