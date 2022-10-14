
int choose(int n, int k)
{
    if (k == 0){return 1;}
    return (n * choose(n-1,k-1)/k);
}

template<typename T>
MBInterpOp<T>::MBInterpOp(
        MBDataPoint a_dst,
        MBIndex& a_patch,
        MBMap<T>& a_map,
        const std::vector<Point>& a_footprint,
        unsigned int order)
{
    PROTO_ASSERT(a_srcs.size() == a_coefs.size(),
            "MBInterpOp::constructor | Error: \
            The source and coefficient arrays must be the same size.");
    m_dst = a_dst;
    // TODO: InterpFootprint should probably be a member of MBInterpOp
    m_srcs = a_map.map().interpFootprint(a_dst.point, a_footprint, a_patch);
    
    // Compute LS Matrix
    int N = choose(DIM + order - 1, DIM);
    int M = m_srcs.size();
    Matrix<double> A(M,N);
    Matrix<double> b(M,1);
    std::set<Point> binomials;
    for (auto bi : Box::Cube(order))
    {
        if (bi.sum() < order - 1)
        {
            binomials.insert(bi);
        }
    }
    PROTO_ASSERT(binomials.size() == N,
        "MBInterpOp::constructor | Error: incorrect number of binomial terms \
        (There should be %i coefficients for DIM = %i and order = %u)",N, DIM, order);
    auto& dstData = m_dst.patch(a_map.map());
    auto dstDataAvg = Stencil<double>::cornersToFaces(order)(dstData);
    for (int ii = 0; ii < M; ii++)
    {
        auto& srcData = m_srcs[ii].patch(a_map.map());
        auto srcDataAvg = Stencil<double>::cornersToFaces(order)(srcData);
        double R = 0;
        double Aij = 0;
        for (int dir = 0; dir < DIM; dir++)
        {
            double x_dst = dstDataAvg(m_dst.point, dir);
            double x_src = srcDataAvg(m_src[ii].point, dir);
            R += pow(x_dst - x_src, 2);
        }
        R = sqrt(R);
        for (int jj = 0; jj < N; jj++)
        {
            for (int dir = 0; dir < DIM; dir++)
            {
                double x0 = dstDataAvg(m_dst.point, dir);
                double xL = srcData(m_srcs[ii].point, dir);
                double xH = srcData(m_srcs[ii].point + Point::Basis(dir), dir);
                double p = binomials[jj][dir];
                Aij += pow(xH - x0, p+1);
                Aij -= pow(xL - x0, p+1);
                Aij /= (p+1);
            }
            A(ii,jj) = Aij;
        }
    }
}

template<typename T>
template<unsigned int C, MemType MEM>
void MBInterpOp<T>::execute(
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
            dst(cc) += m_coefs[ii]*src(cc)
        }
    }
}

template<typename T>
unsigned int MBInterpOp<T>::size() const
{
    return m_coefs.size();
}
