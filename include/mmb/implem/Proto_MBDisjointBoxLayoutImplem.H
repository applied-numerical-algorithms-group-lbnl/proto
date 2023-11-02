MBDisjointBoxLayout::MBDisjointBoxLayout(
        MBProblemDomain& a_domain,
        Point a_boxSize)
{
    m_partition = std::make_shared<MBBoxPartition>();
    std::vector<Point> boxSizes(a_domain.numBlocks(), a_boxSize);
    define(a_domain, boxSizes);
}

MBDisjointBoxLayout::MBDisjointBoxLayout(
        MBProblemDomain& a_domain,
        const std::vector<Point>& a_boxSizes)
{
    m_partition = std::make_shared<MBBoxPartition>();
    define(a_domain, a_boxSizes);
}

MBDisjointBoxLayout::MBDisjointBoxLayout(
        MBProblemDomain& a_domain,
        const std::vector<MBPatchID_t>& a_patches,
        const std::vector<Point>& a_boxSizes)
{
    m_partition = std::make_shared<MBBoxPartition>();
    define(a_domain, a_patches, a_boxSizes);
}

MBDisjointBoxLayout::MBDisjointBoxLayout(
        std::shared_ptr<MBBoxPartition> a_partition,
        const std::vector<Point>& a_boxSizes)
{
    define(a_partition, a_boxSizes);
}

void MBDisjointBoxLayout::define(
        MBProblemDomain& a_domain, 
        const std::vector<Point>& a_boxSizes)
{
    PROTO_ASSERT(a_boxSizes.size() == a_domain.size(),
            "MBDisjointBoxLayout::define | Error: incorrect number of boxSize inputs.");
    
    std::vector<MBPatchID_t> patches;
    for (unsigned int bi = 0; bi < a_domain.numBlocks(); bi++)
    {
        Box patchBox = a_domain.getBlock(bi).box().coarsen(a_boxSizes[bi]);
        for (auto pi : patchBox)
        {
            patches.push_back(MBPatchID_t(pi, bi));
        }
    }
    define(a_domain, patches, a_boxSizes);
}

void MBDisjointBoxLayout::define(
        MBProblemDomain& a_domain,
        const std::vector<MBPatchID_t>& a_patches,
        const std::vector<Point>& a_boxSizes)
{
    a_domain.close();
    PROTO_ASSERT(a_boxSizes.size() == a_domain.size(),
            "MBDisjointBoxLayout::define | Error: incorrect number of boxSize inputs.");

    m_domain = a_domain;
    m_partition->define(a_domain.coarsen(a_boxSizes), a_patches);
    m_layouts.clear();
    m_layouts.resize(a_domain.numBlocks());
    
    for (unsigned int bi = 0; bi < a_domain.numBlocks(); bi++)
    {
        m_layouts[bi].define(m_partition->blockPartition(bi), a_boxSizes[bi]);
    }
}

void MBDisjointBoxLayout::define(
        std::shared_ptr<MBBoxPartition> a_partition,
        const std::vector<Point>& a_boxSizes)
{
    m_partition = a_partition;
    m_domain = a_partition->domain().refine(a_boxSizes);
    m_layouts.clear();
    m_layouts.resize(m_domain.size());
    for (int bi = 0; bi < m_domain.size(); bi++)
    {
        m_layouts[bi].define(m_partition->blockPartition(bi), a_boxSizes[bi]);
    }
}

const MBProblemDomain& MBDisjointBoxLayout::domain() const
{
    return m_domain;
}

unsigned int MBDisjointBoxLayout::procID(const MBIndex& a_index) const
{
    return m_partition->partition()[a_index.global()].second;
}

unsigned int MBDisjointBoxLayout::numBlocks() const
{
    return m_partition->numBlocks();
}

unsigned int MBDisjointBoxLayout::numBoxes() const 
{
    return m_partition->numBoxes();
}

unsigned int MBDisjointBoxLayout::numBoxes(unsigned int a_proc) const 
{
    return m_partition->numBoxes(a_proc);
}

unsigned int MBDisjointBoxLayout::numProcs() const
{
    return m_partition->numProcs();
}

DataIndex<BoxPartition>
MBDisjointBoxLayout::blockIndex(const MBIndex& a_index) const
{
    Point p = point(a_index);
    unsigned int b = block(a_index);
    return m_layouts[b].find(p); 
}

Box MBDisjointBoxLayout::operator[](const MBIndex& a_index) const
{
    return box(a_index);
}

MBIterator MBDisjointBoxLayout::begin() const
{
    MBIterator iter(m_partition);
    return iter.begin();
}

MBIterator MBDisjointBoxLayout::end() const
{
    MBIterator iter(m_partition);
    return iter.end();
}

Point MBDisjointBoxLayout::point(const MBIndex& a_index) const
{
    auto data = m_partition->partition()[a_index.global()];
    auto block = data.first.second;
    auto patch = data.first.first;
    return patch;
}

Box MBDisjointBoxLayout::box(const MBIndex& a_index) const
{
    auto data = m_partition->partition()[a_index.global()];
    auto block = data.first.second;
    auto patch = data.first.first;
    auto boxSize = m_layouts[block].boxSize();
    return Box(patch, patch).refine(boxSize);
}

unsigned int MBDisjointBoxLayout::block(const MBIndex& a_index) const
{
    auto data = m_partition->partition()[a_index.global()];
    auto block = data.first.second;
    return block;
}

DisjointBoxLayout& MBDisjointBoxLayout::getBlock(unsigned int a_block)
{
    PROTO_ASSERT(a_block < numBlocks(),
            "MBDisjointBoxLayout::layout | Error: Block index is out of bounds \
            (numBlocks() = %u)", numBlocks());
    return m_layouts[a_block];
}

const DisjointBoxLayout& MBDisjointBoxLayout::getBlock(unsigned int a_block) const
{
    PROTO_ASSERT(a_block < numBlocks(),
            "MBDisjointBoxLayout::layout | Error: Block index is out of bounds \
            (numBlocks() = %u)", numBlocks());
    return m_layouts[a_block];
}

DisjointBoxLayout MBDisjointBoxLayout::levelLayout() const
{
    return m_layouts[0];
}

MBDisjointBoxLayout::MBIndex MBDisjointBoxLayout::find(Point a_point, unsigned int a_block) const
{
    auto patchDomain = m_layouts[a_block].patchDomain();
    if (!patchDomain.contains(a_point))
    {
        return *end();
    }
    auto globalIndex = m_partition->find(a_point, a_block);
    if (globalIndex == numBoxes())
    {
        return *end();
    } else {
        return MBIndex(m_partition, globalIndex);
    }
}

size_t MBDisjointBoxLayout::offset(unsigned int a_proc) const
{
    return m_partition->procStartIndex(a_proc);
}

bool MBDisjointBoxLayout::compatible(const MBDisjointBoxLayout& a_rhs) const
{
    return m_partition->compatible(*a_rhs.m_partition);
}

bool MBDisjointBoxLayout::isBlockBoundary(
        const MBIndex&  a_index,
        Point           a_dir,
        unsigned int    a_block) const
{
    unsigned int thisBlock = block(a_index);
    Point thisPatch = point(a_index);
    
    Point conn = domain().graph().connectivity(thisBlock, a_block);
    if (conn == Point::Zeros()) {return false;}
    Point invConn = domain().graph().reverseDir(thisBlock, a_block, conn);
    auto R = domain().graph().rotation(thisBlock, conn, a_block);

    Box localBoundPatches = getBlock(thisBlock).patchDomain().box().adjacent(conn);
    Box adjBoundPatches = getBlock(a_block).patchDomain().box().edge(invConn);
    Point neighbor = thisPatch + a_dir;

    Point adjPatch = R.rotateCell(
            neighbor, localBoundPatches, adjBoundPatches);

    auto adjIndex = find(adjPatch, a_block);
    if (adjIndex == end()){return false;}
    return true;
}

std::vector<MBIndex> MBDisjointBoxLayout::bounds(
        const MBIndex& a_index) const
{
    std::vector<MBIndex> ret;
    auto b0 = block(a_index);
    for (unsigned int bi = 0; bi < numBlocks(); bi++)
    {
        auto conn = domain().graph().connectivity(b0, bi);
        if (bi != b0 && conn == Point::Zeros()) { continue; }
        Box patchDomain = getBlock(bi).patchDomain().box();
        for (auto pi : patchDomain)
        {
            auto index = find(pi, bi);
            if (index != *end() && (connectivity(a_index, index) != Point::Zeros()))
            {
                ret.push_back(index);
            }
        }
    }
    return ret;
}

Point MBDisjointBoxLayout::connectivity(
        const MBIndex& a_i1,
        const MBIndex& a_i2) const
{
    unsigned int b1 = block(a_i1);
    unsigned int b2 = block(a_i2);
    Point p1 = point(a_i1);
    Point p2 = point(a_i2);

    // patches on same block
    if (b1 == b2)
    {
        if (Box::Kernel(1).shift(p1).contains(p2))
        {
            return (p2 - p1);
        }
    }
    
    // patches on different blocks

    Point conn12 = domain().graph().connectivity(b1, b2);
    if (conn12 == Point::Zeros()) { return Point::Zeros(); }
    Point conn21 = domain().graph().connectivity(b2, b1);
    PROTO_ASSERT(conn21 != Point::Zeros(),
            "MBDisjointBoxLayout::connectivity | Error: Data corruption.");
    Box adjPatch1 = getBlock(b1).patchDomain().box().adjacent(conn12);
    Box adjPatch2 = getBlock(b2).patchDomain().box().edge(conn21);
    auto R = domain().graph().rotation(b1,conn12,b2);

    Box dirs = Box::Kernel(1).grow(-conn12.abs()).shift(conn12);
    //std::cout << "checking " << p1 << " -> " << p2 << " | conn12: " << conn12 << " | conn21: " << conn21;
    //std::cout << " | dirs: " << dirs << " | adjPatch1: " << adjPatch1 << " | adjPatch2: " << adjPatch2 << std::endl;
    for (auto dir : dirs)
    {
        Point n1 = p1 + dir;
        if (!adjPatch1.contains(n1)){ continue; }
        Point n2 = R.rotateCell(
                n1, adjPatch1, adjPatch2);
        if (n2 == p2) {return dir; }
    }
    return Point::Zeros();
}

bool MBDisjointBoxLayout::isBlockBoundary(
        const MBIndex& a_i1,
        const MBIndex& a_i2) const
{
    if (block(a_i1) == block(a_i2)) {return false;}
    if (connectivity(a_i1, a_i2) == Point::Zeros()) {return false;}
    return true;
}

unsigned int MBDisjointBoxLayout::size() const
{
    return m_partition->numBoxes();
}

unsigned int MBDisjointBoxLayout::localSize() const
{
    return m_partition->numBoxes(Proto::procID());
}

std::vector<Point> MBDisjointBoxLayout::boxSizes() const
{
    std::vector<Point> ret;
    for (int bi = 0; bi < numBlocks(); bi++)
    {
        ret.push_back(m_layouts[bi].boxSize());
    }
    return ret;
}

MBDisjointBoxLayout MBDisjointBoxLayout::coarsen(const Point& a_refRatio) const
{
    std::vector<Point> refRatios(numBlocks(), a_refRatio);
    return this->coarsen(refRatios);
}

MBDisjointBoxLayout MBDisjointBoxLayout::coarsen(const std::vector<Point>& a_refRatios) const
{
    PROTO_ASSERT(a_refRatios.size() == numBlocks(),
            "MBDisjointBoxLayout::coarsen | Error: Wrong number of ref ratios");
    std::vector<Point> crseBoxSize;
    std::vector<Point> fineBoxSize = this->boxSizes();
    for (int bi = 0; bi < numBlocks(); bi++)
    {
        PROTO_ASSERT(fineBoxSize[bi] % a_refRatios[bi] == Point::Zeros(),
            "MBDisjointBoxLayout::coarsen | Error: Layout is not coarsenable by ref ratio.");
        crseBoxSize.push_back(fineBoxSize[bi] / a_refRatios[bi]);
    }
    auto crseDomain = m_domain.coarsen(a_refRatios);
    return MBDisjointBoxLayout(crseDomain, crseBoxSize);
}

std::ostream& operator<<(std::ostream& a_os, const MBIndex& a_di)
{
    auto& partition = a_di.partition().partition();;
    a_os << "block: " << partition[a_di.global()].first.second;
    a_os << " | patch: " << partition[a_di.global()].first.first;
    a_os << " | localIndex: " << a_di.local();
    a_os << " | globalIndex: " << a_di.global() << std::endl;
    return a_os;
}
