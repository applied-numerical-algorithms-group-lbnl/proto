MBDisjointBoxLayout::MBDisjointBoxLayout(
        MBProblemDomain&    a_domain,
        Point                               a_boxSize)
{
    std::vector<Point> boxSizes(a_domain.numBlocks(), a_boxSize);
    define(a_domain, boxSizes);
}

MBDisjointBoxLayout::MBDisjointBoxLayout(
        MBProblemDomain&    a_domain,
        std::vector<Point>                  a_boxSizes)
{
    define(a_domain, a_boxSizes);
}

MBDisjointBoxLayout::MBDisjointBoxLayout(
        MBProblemDomain&  a_domain,
        std::vector<MBPatchID_t>          a_patches,
        std::vector<Point>                a_boxSizes)
{
    define(a_domain, a_patches, a_boxSizes);
}

MBDisjointBoxLayout::MBDisjointBoxLayout(
        std::shared_ptr<MBBoxPartition> a_partition,
        std::vector<Point>              a_boxSizes)
{
    define(a_partition, a_boxSizes);
}

void MBDisjointBoxLayout::define(
        MBProblemDomain& a_domain, 
        std::vector<Point>               a_boxSizes)
{
    PROTO_ASSERT(a_boxSizes.size() == a_domain.size(),
            "MBDisjointBoxLayout::define | Error: incorrect number of boxSize inputs.");
    
    std::vector<MBPatchID_t> patches;
    for (BlockIndex bi = 0; bi < a_domain.numBlocks(); bi++)
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
        MBProblemDomain&            a_domain,
        std::vector<MBPatchID_t>    a_patches,
        std::vector<Point>          a_boxSizes)
{
    a_domain.close();
    PROTO_ASSERT(a_boxSizes.size() == a_domain.size(),
            "MBDisjointBoxLayout::define | Error: incorrect number of boxSize inputs.");

    m_domain = a_domain;
    m_partition = std::make_shared<MBBoxPartition>(a_domain.coarsen(a_boxSizes), a_patches);
    m_layouts.clear();
    m_layouts.resize(a_domain.numBlocks());
    for (BlockIndex bi = 0; bi < a_domain.numBlocks(); bi++)
    {
        m_layouts[bi].define(m_partition->blockPartition(bi), a_boxSizes[bi]);
    }
    m_defined = true;
}

void MBDisjointBoxLayout::define(
        std::shared_ptr<MBBoxPartition> a_partition,
        std::vector<Point>              a_boxSizes)
{
    a_partition->domain().close();
    m_partition = a_partition;
    m_domain = a_partition->domain().refine(a_boxSizes);
    m_layouts.clear();
    m_layouts.resize(m_domain.numBlocks());
    for (int bi = 0; bi < m_domain.size(); bi++)
    {
        m_layouts[bi].define(m_partition->blockPartition(bi), a_boxSizes[bi]);
    }
    m_defined = true;
}

const MBProblemDomain& MBDisjointBoxLayout::domain() const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::domain | Error: Layout is not defined");
    return m_domain;
}

unsigned int MBDisjointBoxLayout::procID(const MBIndex& a_index) const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::procID | Error: Layout is not defined");
    return m_partition->partition()[a_index.global()].second;
}

unsigned int MBDisjointBoxLayout::numBlocks() const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::numBlocks | Error: Layout is not defined");
    return m_partition->numBlocks();
}

unsigned int MBDisjointBoxLayout::numBoxes() const 
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::numBoxes | Error: Layout is not defined");
    return m_partition->numBoxes();
}

unsigned int MBDisjointBoxLayout::numBoxes(unsigned int a_proc) const 
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::numBoxes | Error: Layout is not defined");
    return m_partition->numBoxes(a_proc);
}

unsigned int MBDisjointBoxLayout::numProcs() const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::numProcs | Error: Layout is not defined");
    return m_partition->numProcs();
}

std::vector<MBPatchID_t> MBDisjointBoxLayout::patches() const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::patches | Error: Layout is not defined");
    std::vector<MBPatchID_t> patchList;
    for (auto item : m_partition->partition())
    {
        patchList.push_back(item.first);
    }
    return patchList;
}

DataIndex<BoxPartition>
MBDisjointBoxLayout::blockIndex(const MBIndex& a_index) const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::blockIndex | Error: Layout is not defined");
    Point p = point(a_index);
    BlockIndex b = block(a_index);
    return m_layouts[b].find(p); 
}

Box MBDisjointBoxLayout::operator[](const MBIndex& a_index) const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::operator[] | Error: Layout is not defined");
    return box(a_index);
}

MBIterator MBDisjointBoxLayout::begin() const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::begin | Error: Layout is not defined");
    MBIterator iter(m_partition);
    return iter.begin();
}

MBIterator MBDisjointBoxLayout::end() const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::end | Error: Layout is not defined");
    MBIterator iter(m_partition);
    return iter.end();
}

Point MBDisjointBoxLayout::point(const MBIndex& a_index) const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::point | Error: Layout is not defined");
    auto data = m_partition->partition()[a_index.global()];
    auto block = data.first.second;
    auto patch = data.first.first;
    return patch;
}

Box MBDisjointBoxLayout::box(const MBIndex& a_index) const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::box | Error: Layout is not defined");
    auto& data = m_partition->partition()[a_index.global()];
    auto block = data.first.second;
    auto patch = data.first.first;
    auto boxSize = m_layouts[block].boxSize();
    return Box(patch, patch).refine(boxSize);
}

unsigned int MBDisjointBoxLayout::block(const MBIndex& a_index) const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::block | Error: Layout is not defined");
    auto& data = m_partition->partition()[a_index.global()];
    auto block = data.first.second;
    return block;
}

DisjointBoxLayout& MBDisjointBoxLayout::getBlock(BlockIndex a_block)
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::getBlock | Error: Layout is not defined");
    PROTO_ASSERT(a_block < numBlocks(),
            "MBDisjointBoxLayout::layout | Error: Block index is out of bounds \
            (numBlocks() = %u)", numBlocks());
    return m_layouts[a_block];
}

const DisjointBoxLayout&
MBDisjointBoxLayout::getBlock(BlockIndex a_block) const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::getBlock | Error: Layout is not defined");
    PROTO_ASSERT(a_block < numBlocks(),
            "MBDisjointBoxLayout::layout | Error: Block index is out of bounds \
            (numBlocks() = %u)", numBlocks());
    return m_layouts[a_block];
}

MBDisjointBoxLayout::MBIndex MBDisjointBoxLayout::find(Point a_point, BlockIndex a_block) const
{
    auto globalIndex = m_partition->find(a_point, a_block);
    if (globalIndex == numBoxes()) { return *end(); }
    else { return MBIndex(m_partition, globalIndex); }
    /*
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::find | Error: Layout is not defined");
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
    */
}

size_t MBDisjointBoxLayout::offset(unsigned int a_proc) const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::offset | Error: Layout is not defined");
    return m_partition->procStartIndex(a_proc);
}

bool MBDisjointBoxLayout::compatible(const MBDisjointBoxLayout& a_rhs) const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::compatible | Error: Layout is not defined");
    return m_partition->compatible(*a_rhs.m_partition);
}

bool MBDisjointBoxLayout::isBlockBoundary(
        const MBIndex&  a_index,
        Point           a_dir,
        unsigned int    a_block) const
{
  //PR_TIME("MBLevelBoxData::isBlockBoundary");
    PROTO_ASSERT(m_defined,
            "MBDisjointBoxLayout::isBlockBoundary | Error: Layout is not defined");
    return m_partition->isBlockBoundary(a_index.global(), a_dir, a_block);
}

bool MBDisjointBoxLayout::isBlockBoundary(
        const MBIndex&  a_index,
        Point           a_dir) const
{
    PROTO_ASSERT(m_defined,
            "MBDisjointBoxLayout::isBlockBoundary | Error: Layout is not defined");
    return m_partition->isBlockBoundary(a_index.global(), a_dir);
}

bool MBDisjointBoxLayout::isDomainBoundary(
        const MBIndex&  a_index,
        Point           a_dir) const
{
    PROTO_ASSERT(m_defined,
            "MBDisjointBoxLayout::isDomainBoundary | Error: Layout is not defined");
    return m_partition->isDomainBoundary(a_index, a_dir);
}

bool MBDisjointBoxLayout::isInteriorBoundary(
        const MBIndex&  a_index,
        Point           a_dir) const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::isInteriorBoundary | Error: Layout is not defined");
    return m_partition->isInteriorBoundary(a_index.global(), a_dir);
}
std::vector<MBIndex> MBDisjointBoxLayout::bounds(
        const MBIndex& a_index) const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::bounds | Error: Layout is not defined");
    std::vector<MBIndex> ret;
    auto b0 = block(a_index);
    for (BlockIndex bi = 0; bi < numBlocks(); bi++)
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
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::connectivity | Error: Layout is not defined");
    BlockIndex b1 = block(a_i1);
    BlockIndex b2 = block(a_i2);
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
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::isBlockBoundary | Error: Layout is not defined");
    if (block(a_i1) == block(a_i2)) {return false;}
    if (connectivity(a_i1, a_i2) == Point::Zeros()) {return false;}
    return true;
}

unsigned int MBDisjointBoxLayout::size() const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::size | Error: Layout is not defined");
    return m_partition->numBoxes();
}

unsigned int MBDisjointBoxLayout::localSize() const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::localSize | Error: Layout is not defined");
    return m_partition->numBoxes(Proto::procID());
}

std::vector<Point> MBDisjointBoxLayout::boxSizes() const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::boxSizes | Error: Layout is not defined");
    std::vector<Point> ret;
    for (int bi = 0; bi < numBlocks(); bi++)
    {
        ret.push_back(m_layouts[bi].boxSize());
    }
    return ret;
}

MBDisjointBoxLayout
MBDisjointBoxLayout::coarsen(const Point& a_refRatio) const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::coarsen | Error: Layout is not defined");
    std::vector<Point> refRatios(numBlocks(), a_refRatio);
    return this->coarsen(refRatios);
}

MBDisjointBoxLayout
MBDisjointBoxLayout::coarsen(const std::vector<Point>& a_refRatios) const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::coarsen | Error: Layout is not defined");
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
    MBDisjointBoxLayout crseLayout(crseDomain, patches(), crseBoxSize);
    return crseLayout;
}

void MBDisjointBoxLayout::print() const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::print | Error: Layout is not defined");
    for (auto iter : *this)
    {
        auto bi = block(iter);
        auto pi = point(iter);
        pr_out() << "block: " << bi << " | patch: " << pi << " | box: " << (*this)[iter] << std::endl;
    }
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
