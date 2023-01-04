MBDisjointBoxLayout::MBDisjointBoxLayout(
        MBProblemDomain& a_domain,
        Point a_boxSize)
{
    std::vector<Point> boxSizes(a_domain.numBlocks(), a_boxSize);
    define(a_domain, boxSizes);
}

MBDisjointBoxLayout::MBDisjointBoxLayout(
        MBProblemDomain& a_domain,
        const std::vector<Point>& a_boxSizes)
{
    define(a_domain, a_boxSizes);
}

MBDisjointBoxLayout::MBDisjointBoxLayout(
        MBProblemDomain& a_domain,
        const std::vector<MBPatchID_t>& a_patches,
        const std::vector<Point>& a_boxSizes)
{
    define(a_domain, a_patches, a_boxSizes);
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
        Box patchBox = a_domain.blockDomain(bi).box().coarsen(a_boxSizes[bi]);
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
    m_partition = std::make_shared<MBBoxPartition>(a_domain.coarsen(a_boxSizes), a_patches);
    m_layouts.clear();
    m_layouts.resize(a_domain.numBlocks());
    
    for (unsigned int bi = 0; bi < a_domain.numBlocks(); bi++)
    {
        m_layouts[bi].define(m_partition->blockPartition(bi), a_boxSizes[bi]);
    }
}

MBProblemDomain MBDisjointBoxLayout::domain() const
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

DisjointBoxLayout& MBDisjointBoxLayout::blockLayout(unsigned int a_block)
{
    PROTO_ASSERT(a_block < numBlocks(),
            "MBDisjointBoxLayout::layout | Error: Block index is out of bounds \
            (numBlocks() = %u)", numBlocks());
    return m_layouts[a_block];
}

const DisjointBoxLayout& MBDisjointBoxLayout::blockLayout(unsigned int a_block) const
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

std::ostream& operator<<(std::ostream& a_os, const MBIndex& a_di)
{
    auto& partition = a_di.partition().partition();;
    a_os << "block: " << partition[a_di.global()].first.second;
    a_os << " | patch: " << partition[a_di.global()].first.first;
    a_os << " | localIndex: " << a_di.local();
    a_os << " | globalIndex: " << a_di.global() << std::endl;
    return a_os;
}
