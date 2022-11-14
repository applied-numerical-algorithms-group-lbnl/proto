/*
int choose(int n, int k)
{
    if (k == 0){return 1;}
    return (n * choose(n-1,k-1)/k);
}
*/

template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_binomialPower_tmp(Var<T, 1, MEM>& a_xp, const Var<T, DIM, MEM>& a_x, Point a_p)
{
    a_xp(0) = 1.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        a_xp(0) *= pow(a_x(dir), a_p[dir]);
    }
}
PROTO_KERNEL_END(f_binomialPower_tmp, f_binomialPower)

template<typename Func>
MBPointInterpOp::MBPointInterpOp(
        MBDataPoint a_dst,
        MBMap<Func>& a_map,
        const std::vector<Point>& a_footprint,
        unsigned int a_order)
{
    PR_TIME("MBPointInterpOp::constructor");
    m_dst = a_dst;
    m_srcs = a_map.map().interpFootprint(a_dst.point, a_footprint, a_dst.index);

    // References
    auto& layout = a_map.map().layout();
    unsigned int dstBlock = layout.block(a_dst.index);

    // Compute exponents
    std::vector<Point> exponents;
    for (auto bi : Box::Cube(a_order))
    {
        if (bi.sum() < a_order - 1)
        {
            exponents.push_back(bi);
        }
    }
    int N = exponents.size();
    int M = m_srcs.size();

    // ** Rg = 1/M*sum_j(sqrt(<xj - xg>^2))
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
    
    // Compute Matrices
    Matrix<double> C(M,N);
    std::array<BoxData<double, DIM>, DIM> NT;
    BoxData<double> J(Box::Kernel(1));
    for (int ii = 0; ii < M; ii++)
    {
        // Get data at (and near) the source point
        auto& src = m_srcs[ii];
        unsigned int srcBlock = layout.block(src.index);
        // 3 ghost cells are needed in all directions to compute J
        // 4 ghost cells are needed to compute the quadrature with J
        Box srcBox = Box::Kernel(4).shift(src.point);
        auto xi = a_map(srcBox, srcBlock);
        auto& Ji = src.patch(a_map.jacobian());
        // subtract xg
        xi -= xg;
        // divide by Rg
        xi /= Rg;
        
        Point shift = src.point - Point::Ones() - J.box().low();
        J.shift(shift);
        Ji.copyTo(J);
        for (int jj = 0; jj < N; jj++)
        {
            // compute binomial power
            BoxData<double> xi_p(xi.box());
            forallInPlace(f_binomialPower, xi_p, xi, exponents[jj]);
        
            BoxData<double> xi_p_avg = Stencil<double>::cornersToCells(4)(xi_p);
            auto Jx = Operator::cellProduct(J, xi_p_avg);
            C(ii,jj) = Jx(src.point);
        }
    }
    Matrix<double> D(1,N);

    for (int jj = 0; jj < N; jj++)
    {
        unsigned int dstBlock = layout.block(a_dst.index);
        Box dstBox = Box::Kernel(4).shift(a_dst.point);
        auto xi = a_map(dstBox, dstBlock);
        // compute metric terms
        for (int d = 0; d < DIM; d++)
        {
            NT[d] = Operator::cofactor(xi, d);
        }
        J = Operator::jacobian(xi, NT);
        xi -= xg;
        xi /= Rg;
        BoxData<double> xi_p(xi.box());
        forallInPlace(f_binomialPower, xi_p, xi, exponents[jj]);
        BoxData<double> xi_p_avg = Stencil<double>::cornersToCells(4)(xi_p);
        auto Jx = Operator::cellProduct(J, xi_p_avg);
        D(0,jj) = Jx(a_dst.point);
    }

    // Compute "Stencil"
    auto Cinv = C.inverse();
    auto S = D*Cinv;
    m_coefs.clear();
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

template<typename Func>
MBInterpOp::MBInterpOp(
        MBMap<Func>& a_map,
        const std::vector<Point>& a_footprint,
        Point a_maxGhostSize,
        unsigned int a_order)
{
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
               Box ghostBox = patchBox.adjacent(dir*a_maxGhostSize);
               for (auto gi : ghostBox)
               {
                   MBDataPoint p(iter, gi);
                   MBPointInterpOp op(p, a_map, a_footprint, a_order);
                   m_ops.push_back(op);
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
    for (auto op : m_ops)
    {
        auto block = a_dst.layout().block(op.target().index);
        // check that dst has the data destination point
        op.apply(a_dst, a_src);
    }
}





