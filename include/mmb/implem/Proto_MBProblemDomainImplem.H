
MBProblemDomain::MBProblemDomain()
{
    m_domains.clear();
    m_domainsDefined.clear();
    m_closed = false;
}

MBProblemDomain::MBProblemDomain(unsigned int a_numBlocks)
{
    m_closed = false;
    m_domains.clear();
    m_domains.resize(a_numBlocks);
    m_graph = std::make_shared<MBGraph>(a_numBlocks);
    m_domainsDefined.clear();
    m_domainsDefined.resize(a_numBlocks, false);
}

MBProblemDomain::MBProblemDomain(std::shared_ptr<MBGraph> a_graph)
{
    m_closed = false;
    m_domains.clear();
    m_domains.resize(a_graph->size());
    m_graph = a_graph;
    m_domainsDefined.clear();
    m_domainsDefined.resize(a_graph->size(), false);
}

void MBProblemDomain::defineDomain(BlockIndex a_blockID, Point a_domainSize)
{
    PROTO_ASSERT(!m_closed,
        "MBProblemDomain::defineDomain | Error: ProblemDomain is closed.");
    Array<bool, DIM> periodicity;
    periodicity.fill(false);
    Box domainBox(a_domainSize);
    m_domains[a_blockID].define(domainBox, periodicity);
    m_domainsDefined[a_blockID] = true;
}

void MBProblemDomain::defineBoundary(
        BlockIndex        a_srcBlock,
        BlockIndex        a_dstBlock,
        unsigned int        a_dir,
        Side::LoHiSide      a_side,
        CoordPermutation&   a_rotation)
{
    PROTO_ASSERT(!m_closed,
        "MBProblemDomain::defineBoundary | Error: ProblemDomain is closed.");
    m_graph->defineBoundary(a_srcBlock, a_dstBlock, a_dir, a_side, a_rotation);
}

void MBProblemDomain::defineBoundary(
        BlockIndex        a_srcBlock,
        BlockIndex        a_dstBlock,
        Point               a_dir,
        CoordPermutation&   a_rotation)
{
    PROTO_ASSERT(!m_closed,
        "MBProblemDomain::defineBoundary | Error: ProblemDomain is closed.");
    m_graph->defineBoundary(a_srcBlock, a_dstBlock, a_dir, a_rotation);
}

const ProblemDomain& MBProblemDomain::getBlock(BlockIndex a_blockID) const
{
    return m_domains[a_blockID];
}

Point MBProblemDomain::blockSize(BlockIndex a_blockID) const
{
    return getBlock(a_blockID).box().sizes();
}

double MBProblemDomain::dv(BlockIndex a_blockID) const
{
    double ret = 1.0;
    Point sizes = blockSize(a_blockID);
    for (int dir = 0; dir < DIM; dir++)
    {
        ret /= ((double)sizes[dir]);
    }
    return ret;
}

std::vector<MBGraphArc> MBProblemDomain::boundaries(
        BlockIndex a_srcBlock,
        Point        a_dir) const
{
    return m_graph->boundaries(a_srcBlock, a_dir);
}

Point MBProblemDomain::convert(
        Point a_point,
        BlockIndex a_srcBlock,
        BlockIndex a_dstBlock) const
{
    PROTO_ASSERT(a_srcBlock < size(),
            "MBProblemDomain::convert | Error: srcBlock is out of bounds.");
    PROTO_ASSERT(a_dstBlock < size(),
            "MBProblemDomain::convert | Error: dstBlock is out of bounds.");
    Point n01 = m_graph->connectivity(a_srcBlock, a_dstBlock);
    PROTO_ASSERT(n01.codim() <= 2,
            "MBProblemDomain::convert | Error: src and dst blocks are do not \
share a codimension 1 boundary.");
    
    Box srcBox = m_domains[a_srcBlock].box();
    Box dstBox = m_domains[a_dstBlock].box();
    auto R01 = m_graph->rotation(a_srcBlock, n01, a_dstBlock);
    
    Point d01, d10, v01, v10, s10;
    
    // convert the normal n01 to a codimension-DIM boundary in dst
    // which is contained in the codimension-1 boundary n01
    for (int dir = 0; dir < DIM; dir++)
    {
        if (n01[dir] == 0)
        {
            d01[dir] = -1;
        } else {
            d01[dir] = n01[dir];
        }
    }
    // find the corresponding codimension-DIM boundary in src
    d10 = m_graph->reverseArc(a_srcBlock, a_dstBlock, d01);
    
    // construct position vectors to the common corner defined
    // by d01 and d10
    for (int dir = 0; dir < DIM; dir++)
    {
        v01[dir] = d01[dir] > 0 ? srcBox.size(dir) : 0;
        v10[dir] = d10[dir] > 0 ? dstBox.size(dir) : 0;
    }

    // compute the displacement between origins in dst coordinates
    s10 = v10 - R01(v01);

    // convert the input point to dst coordinates and add to the origin
    // displacement
    return s10 + R01(a_point);
}

Box MBProblemDomain::convert(
        Box a_box,
        BlockIndex a_srcBlock,
        BlockIndex a_dstBlock,
        Centering a_ctr) const
{
    Point L, H;
    switch (a_ctr)
    {
        case PR_CELL:
        {
            L = convert(a_box.low(), a_srcBlock, a_dstBlock);
            H = convert(a_box.high() + Point::Ones(), a_srcBlock, a_dstBlock);
            Box B(L, L);
            B = B & H;
            return Box(B.low(), B.high()-Point::Ones());
        }
        case PR_NODE:
        {
            L = convert(a_box.low(), a_srcBlock, a_dstBlock);
            H = convert(a_box.high(), a_srcBlock, a_dstBlock);
            Box B(L, L);
            B = B & H;
            return B;
        }
        default:
        {
            MayDay<void>::Abort("Method Not Implemented");
        }
    }
    return Box(L,H); //suppresses a compiler warning
}

MBProblemDomain
MBProblemDomain::coarsen(Point a_refRatio) const
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

MBProblemDomain
MBProblemDomain::coarsen(const std::vector<Point>& a_refRatios) const
{
    PROTO_ASSERT(size() == a_refRatios.size(),
            "MBProblemDomain::coarsen | Error: Wrong number of refinement ratios");
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

MBProblemDomain
MBProblemDomain::refine(Point a_refRatio) const
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

MBProblemDomain
MBProblemDomain::refine(const std::vector<Point>& a_refRatios) const
{
    PROTO_ASSERT(size() == a_refRatios.size(),
            "MBProblemDomain::refine | Error: Wrong number of refinement ratios");
    MBProblemDomain refined(size());
    for (int ii = 0; ii < size(); ii++)
    {
        refined.m_domains[ii] = m_domains[ii].refine(a_refRatios[ii]);
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

bool MBProblemDomain::coarsenable(BlockIndex a_block, Point a_refRatio) const
{
    return m_domains[a_block].coarsenable(a_refRatio);
}

bool MBProblemDomain::defined() const
{
    for (auto def : m_domainsDefined)
    {
        if (!def) { return false; }
    }
    return true;
}

bool MBProblemDomain::operator==(const MBProblemDomain& a_rhs) const
{
    bool equal = true;
    equal &= ((*m_graph) == (*a_rhs.m_graph));
    for (int bi = 0; bi < numBlocks(); bi++)
    {
        equal &= (m_domains[bi] == a_rhs.m_domains[bi]);
    }
    return equal;
}

void MBProblemDomain::close() const
{
    if (m_closed) { return; }
    m_graph->close();
    PROTO_ASSERT(defined(),
        "MBProblemDomain::close | Error: Some blocks have undefined domains.");
    for (BlockIndex bi = 0; bi < size(); bi++)
    {
        Box box_i = getBlock(bi).box();
        for (auto dir_ij : Box::Kernel(1))
        {
            if (dir_ij == Point::Zeros()) { continue; }
            if (dir_ij.codim() > 2) {continue; } //TODO: make this work for higher codim
            for (auto bound : boundaries(bi, dir_ij))
            {
                BlockIndex bj = bound.dstBlock;
                auto dir_ji = bound.dstToSrc;
                Box box_j = getBlock(bj).box();

                Box box_ij = box_i.edge(dir_ij,1);
                Box box_ji = box_j.edge(dir_ji,1);
                PROTO_ASSERT(bound.R(box_ij.sizes()).abs() == box_ji.sizes().abs(),
                    "MBProblemDomain::close | Error: Non-conformal boundary between blocks %u and %u.",
                    bi, bj);
            }
        }
    }
    m_closed = true;
}

std::vector<MBPatchID_t> MBProblemDomain::patches(int a_boxSize) const
{
    return patches(Point::Ones(a_boxSize));
}
std::vector<MBPatchID_t> MBProblemDomain::patches(Point a_boxSize) const
{
    std::vector<Point> boxSizes(numBlocks(), a_boxSize);
    return patches(boxSizes);
}

std::vector<MBPatchID_t> MBProblemDomain::patches(const std::vector<Point>& a_boxSizes) const
{
    PROTO_ASSERT(a_boxSizes.size() == numBlocks(),
            "MBProblemDomain::patches | Error: Wrong number of box sizes specified");
    std::vector<MBPatchID_t> patchVect;

    for (int bi = 0; bi < numBlocks(); bi++)
    {
        PROTO_ASSERT(coarsenable(bi, a_boxSizes[bi]),
                "MBProblemDomain::patches | Error: Domain is not coarsenable by box size");
        Box patchDomain = m_domains[bi].box().coarsen(a_boxSizes[bi]);
        for (auto pi : patchDomain)
        {
            patchVect.push_back(MBPatchID_t(pi, bi));
        }
    }
    return patchVect;
}

Point MBProblemDomain::getDir(Point a_point, BlockIndex a_block) const
{
    Box patchDomain = getBlock(a_block).box();
    Point hiDiff = a_point - patchDomain.high();
    Point loDiff = patchDomain.low() - a_point;
    Point dir = Point::Zeros();
    for (int di = 0; di < DIM; di++)
    {
        if (hiDiff[di] > 0)
        {
            dir[di] = 1;
        } else if (loDiff[di] > 0)
        {
            dir[di] = -1;
        } else {
            dir[di] = 0;
        }
    }
    return dir;
}

bool MBProblemDomain::inTriplePoint(
        Point       a_point,
        BlockIndex  a_block) const
{
    return graph().isTriplePoint(a_block, getDir(a_point, a_block));
}

bool MBProblemDomain::inBlockBoundary(
        Point       a_point,
        BlockIndex  a_block) const
{
    return graph().isBlockBoundary(a_block, getDir(a_point, a_block));
}

bool MBProblemDomain::inDomainBoundary(
        Point       a_point,
        BlockIndex  a_block) const
{
    return graph().isDomainBoundary(a_block, getDir(a_point, a_block));
}

bool MBProblemDomain::inInterior(
        Point       a_point,
        BlockIndex  a_block) const
{
    return getBlock(a_block).box().contains(a_point);
}
