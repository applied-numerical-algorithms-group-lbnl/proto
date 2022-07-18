
MBProblemDomain::MBProblemDomain()
{
    m_closed = false;
    m_domains.clear();
    m_domainsDefined.clear();
}

MBProblemDomain::MBProblemDomain(unsigned int a_numBlocks)
{
    m_closed= false;
    m_domains.clear();
    m_domains.resize(a_numBlocks);
    m_graph = std::make_shared<MBGraph>(a_numBlocks);
    m_domainsDefined.clear();
    m_domainsDefined.resize(a_numBlocks, false);
}

void MBProblemDomain::defineDomain(unsigned int a_blockID, Point a_domainSize)
{
    PROTO_ASSERT(!m_closed,
            "MBProblemDomain::defineDomain | Error: MBProblemDomain is closed.");
    std::array<bool, DIM> periodicity;
    periodicity.fill(false);
    Box domainBox(a_domainSize);
    m_domains[a_blockID].define(domainBox, periodicity);
    m_domainsDefined[a_blockID] = true;
}

void MBProblemDomain::defineBoundary(
        unsigned int        a_srcBlock,
        unsigned int        a_dstBlock,
        unsigned int        a_dir,
        Side::LoHiSide      a_side,
        CoordPermutation&   a_rotation)
{
    PROTO_ASSERT(!m_closed,
            "MBProblemDomain::defineDomain | Error: MBProblemDomain is closed.");
    m_graph->defineBoundary(a_srcBlock, a_dstBlock, a_dir, a_side, a_rotation);
}

void MBProblemDomain::defineBoundary(
        unsigned int        a_srcBlock,
        unsigned int        a_dstBlock,
        Point               a_dir,
        CoordPermutation&   a_rotation)
{
    PROTO_ASSERT(!m_closed,
            "MBProblemDomain::defineDomain | Error: MBProblemDomain is closed.");
    m_graph->defineBoundary(a_srcBlock, a_dstBlock, a_dir, a_rotation);
}

const ProblemDomain& MBProblemDomain::blockDomain(unsigned int a_blockID) const
{
    return m_domains[a_blockID];
}

MBProblemDomain MBProblemDomain::coarsen(Point a_refRatio) const
{
    PROTO_ASSERT(coarsenable(a_refRatio),
            "MBProblemDomain::coarsen | Error: Problem domain is not coarsenable.");
    MBProblemDomain coarsened(size());
    for (int ii = 0; ii < size(); ii++)
    {
        coarsened.m_domains[ii] = m_domains[ii].coarsen(a_refRatio);
        coarsened.m_domainsDefined[ii] = true;
    }
    coarsened.m_graph = m_graph;
    return coarsened;
}

MBProblemDomain MBProblemDomain::coarsen(std::vector<Point> a_refRatios) const
{
    MBProblemDomain coarsened(size());
    for (int ii = 0; ii < size(); ii++)
    {
        PROTO_ASSERT(coarsenable(ii, a_refRatios[ii]),
                "MBProblemDomain::coarsen | Error: Problem domain is not coarsenable by inputs.");
        coarsened.m_domains[ii] = m_domains[ii].coarsen(a_refRatios[ii]);
        coarsened.m_domainsDefined[ii] = true;
    }
    coarsened.m_graph = m_graph;
    return coarsened;
}

MBProblemDomain MBProblemDomain::refine(Point a_refRatio) const
{
    MBProblemDomain refined(size());
    for (int ii = 0; ii < size(); ii++)
    {
        refined.m_domains[ii] = m_domains[ii].refine(a_refRatio);
        refined.m_domainsDefined[ii] = true;
    }
    refined.m_graph = m_graph;
    return refined;
}

bool MBProblemDomain::coarsenable(Point a_refRatio) const
{
    for (auto d : m_domains)
    {
        if (!d.coarsenable(a_refRatio)) { return false; }
    }
    return true;
}

bool MBProblemDomain::coarsenable(unsigned int a_block, Point a_refRatio) const
{
    return m_domains[a_block].coarsenable(a_refRatio);
}

bool MBProblemDomain::defined()
{
    for (auto def : m_domainsDefined)
    {
        if (!def) { return false; }
    }
    return true;
}

void MBProblemDomain::close() const
{
    // TODO: validate conformity of blocks
    m_closed = true;

}
