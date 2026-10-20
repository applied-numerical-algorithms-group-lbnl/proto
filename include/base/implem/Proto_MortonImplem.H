Morton& Morton::instance()
{
    static Morton m;
    return m;
}

uint64_t Morton::index(const Point& a_pt)
{
    uint64_t retval = 0;
#if DIM > 1
    bool cached=true;
    for(int d=0;d<DIM;++d)
        if(a_pt[d]<0 || a_pt[d]>=MORTONSIZE) cached=false;
    if(cached)
    {
        for(int d=0;d<DIM;++d) retval+=instance().m_morton1D[d][a_pt[d]];
    }
    else
    {
        // Sparse AMR can reach large logical patch coordinates without
        // allocating a large mesh. Keep the small-coordinate lookup fast.
        constexpr unsigned bits=64/DIM;
        for(unsigned d=0;d<DIM;++d)
        {
            PROTO_ASSERT(a_pt[d]>=0 && uint64_t(a_pt[d])<(uint64_t(1)<<bits),
                "Morton::index | Patch coordinate exceeds the 64-bit Morton range");
            const uint64_t coordinate=a_pt[d];
            for(unsigned bit=0;bit<bits;++bit)
                retval|=((coordinate>>bit)&1)<<(DIM*bit+d);
        }
    }
#else
    retval = a_pt[0];
#endif
    return retval;
}

void Morton::sort(vector<Point>& a_pts)
{
    vector<pair<Point,uint64_t> > sorter;
    for (int k = 0; k < a_pts.size();++k)
    {
        sorter.push_back(pair<Point,uint64_t>(a_pts[k], Morton::index(a_pts[k])));
    }

    std::sort(sorter.begin(),sorter.end(),compareSecond);

    for (int k = 0; k < sorter.size();++k)
    {
        a_pts[k] = sorter[k].first;
    }
}

Morton::Morton()
{
    PR_TIMERS("Morton define");
    PROTO_ASSERT(DIM < 7,
            "Morton::Constructor | Constructor is not defined for DIM >= 7");
    for (int d = 0; d < DIM; d++)
    {
        m_morton1D[d]=vector<uint64_t>(MORTONSIZE);
    }
#if DIM > 1
    uint64_t mask0 = 1;
    for (uint64_t i = 0; i < MORTONSIZE; i++)
    {
        for (uint64_t d = 0; d < DIM;d++)
        {
            m_morton1D[d][i] = 0;
        }
        for (uint64_t logi = 0; logi < LOGSIZE; logi++)
        {
            for (uint64_t d = 0;d < DIM; d++) 
            {
                m_morton1D[d][i] += 
                    ((i >> logi) & mask0) << (DIM*logi + d);
            }
        }
    }
#endif
}
