/*
int choose(int n, int k)
{
    if (k == 0){return 1;}
    return (n * choose(n-1,k-1)/k);
}
*/
template<typename T>
MBInterpOp<T>::MBInterpOp(
        MBDataPoint a_dst,
        MBMap<T>& a_map,
        const std::vector<Point>& a_footprint,
        unsigned int a_order)
{
    m_dst = a_dst;
    m_srcs = a_map.map().interpFootprint(a_dst.point, a_footprint, a_dst.index);
    
    // Compute LS Matrix
    std::vector<Point> binomials;
    for (auto bi : Box::Cube(a_order))
    {
        if (bi.sum() < a_order - 1)
        {
            binomials.push_back(bi);
        }
    }
    int N = binomials.size();
    int M = m_srcs.size();
    Matrix<double> A(M,N);
    Matrix<double> b(M,1);
    auto& dstData = m_dst.patch(a_map.map());
    BoxData<double, 3> dstDataAvg = Stencil<double>::cornersToCells(a_order)(dstData);
    for (int ii = 0; ii < M; ii++)
    {
        auto& srcData = m_srcs[ii].patch(a_map.map());
        BoxData<double, 3> srcDataAvg = Stencil<double>::cornersToCells(a_order)(srcData);
        double R = 0;
        for (int dir = 0; dir < DIM; dir++)
        {
            double x_dst = dstDataAvg(m_dst.point, dir);
            double x_src = srcDataAvg(m_srcs[ii].point, dir);
            R += pow(x_dst - x_src, 2);
        }
        R = sqrt(R);
        for (int jj = 0; jj < N; jj++)
        {
            double Aij = 1;
            if (ii != jj) {Aij = 0;}
            for (int dir = 0; dir < DIM; dir++)
            {
                // Compute Aij
                /*
                double x0 = dstDataAvg(m_dst.point, dir);
                double xL = srcData(m_srcs[ii].point, dir);
                double xH = srcData(m_srcs[ii].point + Point::Basis(dir), dir);
                double p = binomials[jj][dir];
                Aij += pow(xH - x0, p+1);
                Aij -= pow(xL - x0, p+1);
                Aij /= (p+1);
                */
            }
            A(ii,jj) = Aij;
        }
        b(ii,0) = 1;
    }
    solveLS(A,b);

    m_coefs.clear();
    m_coefs.resize(M);
    for (int ii = 0; ii < M; ii++)
    {
        m_coefs[ii] = b(ii,0);
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
            dst(cc) += m_coefs[ii]*src(cc);
        }
    }
}

template<typename T>
unsigned int MBInterpOp<T>::size() const
{
    return m_coefs.size();
}
