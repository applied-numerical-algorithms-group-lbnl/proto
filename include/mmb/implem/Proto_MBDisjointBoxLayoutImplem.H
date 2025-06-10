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
        std::vector<MBPoint>          a_patches,
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
    
    std::vector<MBPoint> patches;
    for (BlockIndex bi = 0; bi < a_domain.numBlocks(); bi++)
    {
        Box patchBox = a_domain.getBlock(bi).box().coarsen(a_boxSizes[bi]);
        for (auto pi : patchBox)
        {
            patches.push_back(MBPoint(pi, bi));
        }
    }
    define(a_domain, patches, a_boxSizes);
}

void MBDisjointBoxLayout::define(
        MBProblemDomain&            a_domain,
        std::vector<MBPoint>    a_patches,
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

std::vector<MBPoint> MBDisjointBoxLayout::patches() const
{
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::patches | Error: Layout is not defined");
    std::vector<MBPoint> patchList;
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
    validateIndex(a_index, "point");
    return m_partition->getPoint(a_index.global());
}

Box MBDisjointBoxLayout::box(const MBIndex& a_index) const
{
    validateIndex(a_index, "box");
    auto patch = point(a_index);
    auto boxSize = m_layouts[block(a_index)].boxSize();
    return Box(patch, patch).refine(boxSize);
}

Array<double, DIM> MBDisjointBoxLayout::dx(BlockIndex a_block) const
{
    Point sizes = blockDomainBox(a_block).sizes();
    Array<double, DIM> h;
    h.fill(1.0);
    for (int dd = 0; dd < DIM; dd++) { h[dd] /= sizes[dd]; }
    return h;
}

Array<double, DIM> MBDisjointBoxLayout::dx(MBIndex a_index) const
{
    return dx(block(a_index));
}

unsigned int MBDisjointBoxLayout::block(const MBIndex& a_index) const
{
    validateIndex(a_index, "block");
    return m_partition->getBlock(a_index.global());
}

MBPoint MBDisjointBoxLayout::patch(const MBIndex& index) const
{
    validateIndex(index, "patch");
    return m_partition->getPatch(index.global());
}

DisjointBoxLayout& MBDisjointBoxLayout::getBlock(BlockIndex a_block)
{
    #ifndef NDEBUG
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::getBlock | Error: Layout is not defined");
    PROTO_ASSERT(a_block < numBlocks(),
            "MBDisjointBoxLayout::layout | Error: Block index is out of bounds \
            (numBlocks() = %u)", numBlocks());
    #endif
    return m_layouts[a_block];
}

const DisjointBoxLayout&
MBDisjointBoxLayout::getBlock(BlockIndex a_block) const
{
    validateSelf("getBlock");
    validateBlock(a_block, "getBlock");
    return m_layouts[a_block];
}

MBDisjointBoxLayout::MBIndex MBDisjointBoxLayout::find(Point a_point, BlockIndex a_block) const
{
    auto globalIndex = m_partition->find(a_point, a_block);
    if (globalIndex == numBoxes()) { return *end(); }
    else { return MBIndex(m_partition, globalIndex); }
}

size_t MBDisjointBoxLayout::offset(unsigned int a_proc) const
{
    validateSelf("offset");
    return m_partition->procStartIndex(a_proc);
}

bool MBDisjointBoxLayout::compatible(const MBDisjointBoxLayout& a_rhs) const
{
    validateSelf("compatible");
    return m_partition->compatible(*a_rhs.m_partition);
}

bool MBDisjointBoxLayout::compatible(const MBIndex& a_index) const
{
    validateSelf("compatible");
    return m_partition->compatible(a_index.partition());
}

Point MBDisjointBoxLayout::connectivity(
    const MBPoint &patch_1,
    const MBPoint &patch_2) const
{
    validateSelf("connectivity(patch, patch)");
    return m_partition->connectivity(patch_1, patch_2);
}

Point MBDisjointBoxLayout::connectivity(
    const MBIndex &i1,
    const MBIndex &i2) const
{
    auto p1 = patch(i1);
    auto p2 = patch(i2);
    return connectivity(p1, p2);
}

bool MBDisjointBoxLayout::isBlockBoundary(
    const MBPoint &patch_1,
    const MBPoint &patch_2) const
{
    validateSelf("isBlockBoundary(patch, patch)");
    if (patch_1.block == patch_2.block)
    {
        return false;
    }
    if (connectivity(patch_1, patch_2) == Point::Zeros())
    {
        return false;
    }
    return true;
}

bool MBDisjointBoxLayout::isBlockBoundary(
    const MBIndex &patch_1,
    const MBIndex &patch_2) const
{
    return isBlockBoundary(patch(patch_1), patch(patch_2));
}

bool MBDisjointBoxLayout::isBlockBoundary(
    const MBPoint& patch,
    Point          dir) const
{
    validateSelf("isBlockBoundary");
    return m_partition->isBlockBoundary(patch, dir);
}

bool MBDisjointBoxLayout::isBlockBoundary(
    const MBIndex& patch,
    Point          dir) const
{
    return isBlockBoundary(this->patch(patch), dir);
}

bool MBDisjointBoxLayout::isBlockBoundary(
        const MBPoint&  patch,
        Point           dir,
        BlockIndex      adjBlock) const
{
    validateSelf("isBlockBoundary");
    return m_partition->isBlockBoundary(patch, dir, adjBlock);
}



bool MBDisjointBoxLayout::isPatchOnBlockBoundary(
    const MBPoint& patch) const 
{
    validateSelf("isPatchOnBlockBoundary");
    for (auto dir : Point::Directions())
    {
        if (isPatchOnBlockBoundary(patch, dir)) { return true; }
    }
    return false;
}

bool MBDisjointBoxLayout::isPatchOnBlockBoundary(
    const MBPoint& patch,
    Point a_dir) const 
{
    validateSelf("isPatchOnBlockBoundary");
    return patchDomain().isPointOnBlockBoundary(patch.point, patch.block, a_dir);
}

bool MBDisjointBoxLayout::isPatchInBlockBoundary(
    const MBPoint& patch) const
{
    validateSelf("isPatchInBlockBoundary");
    return patchDomain().isPointInBlockBoundary(patch.point, patch.block);
}

bool MBDisjointBoxLayout::isPatchInBlockBoundary(
    const MBPoint& patch,
    Point dir) const
{
    validateSelf("isPatchInBlockBoundary");
    return patchDomain().isPointInBlockBoundary(patch.point, patch.block, dir);
}

bool MBDisjointBoxLayout::isDomainBoundary(
        const MBPoint& patch,
        Point          dir) const
{
    validateSelf("isDomainBoundary");
    return m_partition->isDomainBoundary(patch, dir);
}

bool MBDisjointBoxLayout::isDomainBoundary(
    const MBIndex& patch,
    Point          dir) const
{
    return isDomainBoundary(this->patch(patch), dir);
}

bool MBDisjointBoxLayout::isPatchOnDomainBoundary(
    const MBPoint& patch) const 
{
    validateSelf("isPatchOnDomainBoundary");
    for (auto dir : Point::Directions())
    {
        if (isPatchOnDomainBoundary(patch, dir)) { return true; }
    }
    return false;
}

bool MBDisjointBoxLayout::isPatchOnDomainBoundary(
    const MBPoint& patch,
    Point          dir) const 
{
    validateSelf("isPatchOnDomainBoundary");
    return patchDomain().isPointOnDomainBoundary(patch.point, patch.block, dir);
}

bool MBDisjointBoxLayout::isPatchInDomainBoundary(
    const MBPoint& patch) const
{
    validateSelf("isPatchInDomainBoundary");
    return patchDomain().isPointInDomainBoundary(patch.point, patch.block);
}
bool MBDisjointBoxLayout::isPatchInDomainBoundary(
    const MBPoint& patch,
    Point dir) const
{
    validateSelf("isPatchInDomainBoundary");
    return patchDomain().isPointInDomainBoundary(patch.point, patch.block, dir);
}

bool MBDisjointBoxLayout::isRefinementBoundary(
    const MBPoint& patch,
    Point          dir) const
{
    validateSelf("isRefinementBoundary");
    return m_partition->isRefinementBoundary(patch, dir);
}
bool MBDisjointBoxLayout::isPatchOnRefinementBoundary(
    const MBPoint& patch) const 
{
    validateSelf("isPatchOnRefinementBoundary");
    for (auto dir : Point::Directions())
    {
        if (isRefinementBoundary(patch, dir)) { return true; }
    }
    return false;
}

bool MBDisjointBoxLayout::isInteriorBoundary(
        const MBPoint&  patch,
        Point           dir) const
{
    validateSelf("isInteriorBoundary");
    return m_partition->isInteriorBoundary(patch, dir);
}

bool MBDisjointBoxLayout::isPatchInTriplePointRegion(
    const MBPoint& patch) const
{
    validateSelf("isPatchInTriplePointRegion");
    return patchDomain().isPointInTriplePoint(patch.point, patch.block);
}
bool MBDisjointBoxLayout::isPatchInTriplePointRegion(
    const MBPoint& patch,
    Point dir) const
{
    validateSelf("isPatchInTriplePointRegion");
    return patchDomain().isPointInTriplePoint(patch.point, patch.block, dir);
}
std::vector<MBIndex> MBDisjointBoxLayout::adjacentIndices( const MBPoint& patch) const
{
    validateSelf("adjacentIndices");
    std::vector<MBIndex> indices;
    auto adjacentPatches = m_partition->adjacentPatches(patch);
    for (auto& adjPatch : adjacentPatches)
    {
        auto index = find(adjPatch.point, adjPatch.block);
        if (index != end())
        {
            indices.push_back(index);
        }
    }
    return indices;
}

std::vector<MBIndex> MBDisjointBoxLayout::adjacentIndices(
        const MBPoint& patch,
        Point dir) const
{
    validateSelf("adjacentIndices");
    std::vector<MBIndex> indices;
    auto adjacentPatches = m_partition->adjacentPatches(patch);
    for (auto& adjPatch : adjacentPatches)
    {
        if (connectivity(patch, adjPatch) != dir) { continue; }
        auto index = find(adjPatch.point, adjPatch.block);
        if (index != end())
        {
            indices.push_back(index);
        }
    }
    return indices;
}

#if 0
inline std::vector<MBIndex> MBDisjointBoxLayout::findAdjacentIndices(
    MBIndex index,
    Point dir) const
{
    PR_TIME("MBDisjointBoxLayout::findAdjacentIndices");
    auto allBounds = bounds(index);
    std::vector<MBIndex> dirBounds;
    for (auto bound : allBounds)
    {
        if (connectivity(index, bound) == dir)
        {
            dirBounds.push_back(bound);
        }
    }
    return dirBounds;
}

MBIndex MBDisjointBoxLayout::findAdjacentIndex(
    PatchID a_patch,
    BlockIndex a_block,
    Point a_dir) const
{
    PROTO_ASSERT(a_dir.codim() == 1,
        "MBDisjointBoxLayout::findAdjacentIndex | Error: This function is only well defined for codimension-1 directions");
    PatchID adjacentPatch = a_patch + a_dir;
    auto& patchDomain = partition().domain();
    PROTO_ASSERT(patchDomain.getBlock(a_block).box().containsPoint(a_patch),
        "MBDisjointBoxLayout::findAdjacentIndex | Error: Base patch is not inside the patch domain");
    MBIndex adjacentIndex = *this->end();

    if (patchDomain.isPointInInterior(adjacentPatch, a_block))
    {
        // INTERIOR case (normal case)
        adjacentIndex = find(adjacentPatch, a_block);
    } else if (patchDomain.isPointInBlockBoundary(adjacentPatch, a_block))
    {
        // BLOCK BOUNDARY case
        BlockIndex adjacentBlock = patchDomain.graph().adjacent(a_block, a_dir);
        PatchID adjacentPatchLocalCoords = 
            patchDomain.convertPoint(adjacentPatch, a_block, adjacentBlock, PR_CELL);
        adjacentIndex = find(adjacentPatchLocalCoords, adjacentBlock);
    }
    else if (patchDomain.isPointInDomainBoundary(adjacentPatch, a_block))
    {
        // EXTERIOR case (no CF boundary)
    } else {
        // INVALID case
        PROTO_ASSERT(false, "MBDisjointBoxLayout::findAdjacentIndex | Error: Data corruption.");
    }  
    return adjacentIndex;
}

std::vector<MBIndex> MBDisjointBoxLayout::bounds(
        const MBIndex& a_index) const
{
    PR_TIME("MBDisjointBoxLayout::bounds");
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::bounds | Error: Layout is not defined");
    std::vector<MBIndex> ret;
    auto adjacent = m_partition->adjacentPatches(MBPoint(point(a_index), block(a_index)));
    for (auto& adj : adjacent)
    {
        auto index = find(adj.point, adj.block);
        if (index != end())
        {
            ret.push_back(index);
        }
    }
    return ret;
}
#endif



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

Point MBDisjointBoxLayout::boxSize(BlockIndex a_block) const
{
    PROTO_ASSERT(a_block >= 0 && a_block < numBlocks(),
            "MBDisjointBoxLayout::boxSize | Error: invalid block index.");
    return boxSizes()[a_block];
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

void MBDisjointBoxLayout::validateSelf(std::string methodName) const
{
    #ifndef NDEBUG
    PROTO_ASSERT(m_defined, "MBDisjointBoxLayout::%s | Error: Incompatible index", methodName);
    #endif
}

void MBDisjointBoxLayout::validateIndex(const MBIndex& index, std::string methodName) const
{
    #ifndef NDEBUG
    PROTO_ASSERT(compatible(index), "MBDisjointBoxLayout::%s | Error: Incompatible index", methodName.c_str());
    #endif
}

void MBDisjointBoxLayout::validateBlock(const BlockIndex& block, std::string methodName) const
{
    #ifndef NDEBUG
    PROTO_ASSERT(block < numBlocks(),
        "MBDisjointBoxLayout::%s | Error: Block index %i is out of bounds (numBlocks = %i)",
        methodName.c_str(), block, numBlocks());
    #endif
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
    a_os << "block: " << partition[a_di.global()].first.block;
    a_os << " | patch: " << partition[a_di.global()].first.point;
    a_os << " | localIndex: " << a_di.local();
    a_os << " | globalIndex: " << a_di.global() << std::endl;
    return a_os;
}
