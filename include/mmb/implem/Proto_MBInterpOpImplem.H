
template<typename T>
MBInterpOp<T>::MBInterpOp(
        MBDataPoint a_dst,
        std::vector<MBDataPoint> a_srcs,
        std::vector<T> a_coefs)
{
    PROTO_ASSERT(a_srcs.size() == a_coefs.size(),
            "MBInterpOp::constructor | Error: \
            The source and coefficient arrays must be the same size.");
    m_dst = a_dst;
    m_srcs = a_srcs;
    m_coefs = a_coefs;
}

template<typename T>
template<unsigned int C, MemType MEM>
void execute(
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

inline unsigned int size() const
{
    return m_coefs.size();
}
