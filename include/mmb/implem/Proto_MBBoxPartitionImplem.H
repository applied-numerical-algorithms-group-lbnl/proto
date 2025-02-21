MBBoxPartition::MBBoxPartition(
        const MBProblemDomain& a_patchDomain)
{
    std::vector<MBPoint> patches;
    // creates the BoxPartition with no patches in it
    define(a_patchDomain, patches);
}

MBBoxPartition::MBBoxPartition(
        const MBProblemDomain& a_patchDomain,
        const std::vector<MBPoint>& a_patches)
{
    define(a_patchDomain, a_patches);
}

void MBBoxPartition::define(
        const MBProblemDomain& a_patchDomain,
        const std::vector<MBPoint>& a_patches)
{
    a_patchDomain.close();
    m_patchDomain = a_patchDomain;
    if (m_blockPartitions.size() != a_patchDomain.size())
    {
        for (unsigned int bi = 0; bi < a_patchDomain.size(); bi++)
        {
            const auto& domain = a_patchDomain.getBlock(bi);
            m_blockPartitions.push_back(std::make_shared<BoxPartition>(domain));
        }
    }
    define(a_patches);
}

void MBBoxPartition::define(
        const std::vector<MBPoint>& a_patches)
{
    std::vector<std::vector<Point>> blockPatches;
    blockPatches.resize(m_patchDomain.numBlocks());
    for (auto pi : a_patches)
    {
        blockPatches[pi.block].push_back(pi.point);
    }
    
    for (unsigned int bi = 0; bi < m_patchDomain.size(); bi++)
    {
        m_blockPartitions[bi]->define(
                m_patchDomain.getBlock(bi), blockPatches[bi], 0, numProc());
    }
    loadBalance();
    findBoundaries();
    findAdjacencies();
}

void MBBoxPartition::buildLocalMaps()
{
    uint64_t N = numBoxes();
    uint64_t n = N / Proto::numProc();
    uint64_t r = N % Proto::numProc();
    std::vector<uint64_t> unassignedBoxes(Proto::numProc(), n);
    for (int ii = 0; ii < r; ii++)
    {
        unassignedBoxes[ii]++;
    }
    unsigned int proc = 0;
    for (BlockIndex bi = 0; bi < numBlocks(); bi++)
    {
        uint64_t numAssigned = 0;
        uint64_t blockSize = m_blockPartitions[bi]->numBoxes();
        std::vector<std::pair<int, unsigned int>> blockAssignment;
        while (numAssigned < blockSize)
        {
            uint64_t ni = min(unassignedBoxes[proc], blockSize - numAssigned);
            blockAssignment.push_back(std::pair<int, unsigned int>(proc,ni));
            unassignedBoxes[proc] -= ni;
            numAssigned += ni;
            if (unassignedBoxes[proc] == 0) { proc++; }
        }
        m_blockPartitions[bi]->loadAssign(blockAssignment);
    }
}

void MBBoxPartition::buildGlobalMaps()
{
    m_partition.clear();
    m_indexMap.clear();
    m_procMap.clear();
    uint64_t globalIndex = 0;
    unsigned int currentProc = 0;
    uint64_t procRangeStart = 0;
    for (BlockIndex bi = 0; bi < numBlocks(); bi++)
    {
        auto patchData = m_blockPartitions[bi]->partition();
        for (auto patchData_i : patchData)
        {
            Point patch = patchData_i.first;
            unsigned int  proc  = patchData_i.second;
            if (proc != currentProc)
            {
                std::pair<uint64_t, uint64_t> procRange(procRangeStart, globalIndex);
                m_procMap[currentProc] = procRange;
                currentProc = proc;
                procRangeStart = globalIndex;
            }
            MBPoint mbPatch(patch, bi);
            m_partition.push_back(std::pair<MBPoint, unsigned int>(mbPatch, proc));
            m_indexMap[mbPatch] = globalIndex;
            globalIndex++;
        }
        std::pair<uint64_t, uint64_t> procRange(procRangeStart, numBoxes());
        m_procMap[currentProc] = procRange;
    }
}

void MBBoxPartition::loadBalance()
{
    buildLocalMaps();
    buildGlobalMaps();
}

void MBBoxPartition::findBoundaries()
{
    typedef std::pair<uint64_t, Point> PatchDir;
    m_blockBoundaries.clear();
    m_domainBoundaries.clear();
    auto& graph = m_patchDomain.graph();
    for (uint64_t ind_i = 0; ind_i < numBoxes(); ind_i++)
    {
        auto& partitionData = m_partition[ind_i];
        int   block_i = partitionData.first.block;
        auto  patch_i = partitionData.first.point;
        for (auto dir : Point::Directions())
        {
            PatchBound key; key.index = ind_i; key.dir = dir;
            bool foundBlockBoundary = false; 
            auto patch_ji = patch_i + dir;
            if (m_patchDomain.getBlock(block_i).box().containsPoint(patch_ji))
            {
                // patch_j is interior to the block
                continue; 
            } else if (m_patchDomain.isPointInTriplePoint(patch_i + dir, block_i))
            {
                std::vector<uint64_t> emptyVector;
                m_blockBoundaries[key] = emptyVector;
                foundBlockBoundary = true;
            }
            else
            {
                for (BlockIndex block_j = 0; block_j < numBlocks(); block_j++)
                {
                    if (block_j == block_i)
                    {
                        continue;
                    }
                    Point conn_ij = graph.connectivity(block_i, block_j);
                    if (conn_ij.dot(dir) <= 0)
                    {
                        continue;
                    }
                    Point conn_ji = graph.reverseDir(block_i, block_j, conn_ij);
                    auto R = graph.rotation(block_i, conn_ij, block_j);

                    Box localBoundPatches = m_patchDomain.getBlock(block_i).box().adjacent(conn_ij);
                    Box adjBoundPatches = m_patchDomain.getBlock(block_j).box().edge(conn_ji);

                    Point patch_jj = R.rotateCell(patch_ji, localBoundPatches, adjBoundPatches);

                    if (!localBoundPatches.containsPoint(patch_ji))
                    {
                        continue;
                    }

                    PROTO_ASSERT(adjBoundPatches.containsPoint(patch_jj),
                                 "MBBoxPartition::findBoundaries | Error: Data Corruption.");

                    auto ind_j = find(patch_jj, block_j);
                    m_blockBoundaries[key].push_back(ind_j);
                    foundBlockBoundary = true;
                    // Note that if patch_jj does not exist in the refinement, the block boundary is
                    // still added but the associated patch index will be equal to numBoxes(). This
                    // is the case where a block boundary is also a refinement boundary.
                }
            }

            if (!foundBlockBoundary)
            {
                m_domainBoundaries.insert(key);
            }
        }
    }
    findRefinementBoundaries();
}

void MBBoxPartition::findRefinementBoundaries()
{
    m_refinementBoundaries.clear();
    for (uint64_t index = 0; index < numBoxes(); index++)
    {
        auto& partitionData = m_partition[index];
        int   block = partitionData.first.block;
        auto  patch = partitionData.first.point;
        for (auto dir : Point::Directions())
        {
            PatchBound key; key.index = index; key.dir = dir;
            if (isDomainBoundary(index, dir)) { continue; }

            auto adjPatch = patch + dir;
            Box patchDomainBox = m_patchDomain.getBlock(block).box();
            if (patchDomainBox.containsPoint(adjPatch))
            {
                if (find(adjPatch, block) == numBoxes())
                {
                    m_refinementBoundaries.insert(key);
                }
            } else {
                PROTO_ASSERT(m_blockBoundaries.count(key),
                    "MBBoxPartition::findRefinementBoundaries | Error: Data Corruption");
                auto adjIndices = m_blockBoundaries[key];
                for (BlockIndex adjIndex : adjIndices)
                {
                    if (adjIndex == numBoxes()) 
                    {
                        m_refinementBoundaries.insert(key);
                    }
                    // For codimension 2+ boundaries, it's possible that some adjacent blocks
                    // will constitute a refinement boundary and others will not. This function
                    // takes the convention that if there is at least one refinement boundary, then
                    // this (patch, dir) pair constitutes a refinement boundary. Due to the nature
                    // of regridding and coarse / fine boundary computations, I don't believe that
                    // this distinction will matter for now. - 2/3/2025
                }
            }
        }
    }
}

void MBBoxPartition::findAdjacencies() {
    
    m_patchAdjacencies.clear();
    for (BlockIndex block = 0; block < numBlocks(); block++)
    {
        for (PatchID pi : m_patchDomain.getBlock(block).box())
        {
            MBPoint patch(pi, block);
            bool exists = contains(pi, block);
            MBPatchInfo patchInfo(patch, exists);

            getBoundaryInfo(patchInfo);
            
            m_patchAdjacencies[patch] = patchInfo;
        }
    }
    updateRefinementBoundaries();
}

void MBBoxPartition::getBoundaryInfo(MBPatchInfo& patchInfo)
{
    BlockIndex block = patchInfo.patch.block;
    Point pi = patchInfo.patch.point;
    Box blockDomainBox = m_patchDomain.getBlock(block).box();
    auto& graph = m_patchDomain.graph();
    for (auto dir : Point::Directions())
    {
        MBPatchBoundInfo boundInfo;
        PatchID pjLocal = pi + dir;
        if (blockDomainBox.containsPoint(pjLocal))
        {
            MBPoint adjPatch(pjLocal, block);
            boundInfo.adjacentPatches.insert(adjPatch);
            boundInfo.boundaryType = PR_INTERIOR_BOUNDARY;
        } else {
            Point boundaryDir = blockDomainBox.whichBoundaryContains(pjLocal);
            if (graph.isBlockBoundary(block, boundaryDir))
            {
                for (BlockIndex adjBlock : graph.adjacentBlocks(block, boundaryDir))
                {
                    PatchID pjAdjacent = m_patchDomain.convertPoint(pjLocal, block, adjBlock, PR_CELL);
                    MBPoint adjPatch(pjAdjacent, adjBlock);
                    boundInfo.adjacentPatches.insert(adjPatch);
                }
                boundInfo.boundaryType = PR_BLOCK_BOUNDARY;
            } else if (graph.isDomainBoundary(block, boundaryDir))
            {
                boundInfo.boundaryType = PR_DOMAIN_BOUNDARY;
            }
        }
        boundInfo.isRefinementBoundary = false; //this information is updated in updateRefinementBoundaries
        patchInfo.boundaries[dir] = boundInfo;
    }
    
}

void MBBoxPartition::updateRefinementBoundaries()
{
    for (auto& [patch, patchInfo] : m_patchAdjacencies)
    {
        bool patchExists = contains(patch.point, patch.block);
        for (auto& [dir, boundInfo] : patchInfo.boundaries)
        {
            for (auto adjPatch : boundInfo.adjacentPatches)
            {
                bool adjPatchExists = contains(adjPatch.point, adjPatch.block);
                boundInfo.isRefinementBoundary = (patchExists ^ adjPatchExists);

            }
        }
        patchInfo.exists = patchExists;
    }
}

BlockIndex MBBoxPartition::block(uint64_t a_index) const
{
    auto& partitionData = m_partition[a_index];
    return partitionData.first.block;
}

bool MBBoxPartition::compatible(const MBBoxPartition& a_rhs)
{
    if (this == &a_rhs) { return true; }
    if (numBlocks() != a_rhs.numBlocks()) { return false; }
    for (BlockIndex bi = 0; bi < numBlocks(); bi++)
    {
        if (!m_blockPartitions[bi]->compatible(*a_rhs.m_blockPartitions[bi]))
        {
            return false;
        }
    }
    return true;
}
unsigned int MBBoxPartition::numProcs() const
{
    return m_procMap.size();
}
unsigned int MBBoxPartition::numBoxes() const
{
    int N = 0;
    for (auto p : m_blockPartitions)
    {
        N += p->numBoxes();
    }
    return N;
}
unsigned int MBBoxPartition::numBoxes(unsigned int a_proc)
{
    return m_procMap[a_proc].second - m_procMap[a_proc].first;
}

BlockIndex MBBoxPartition::numBlocks() const
{
    return m_blockPartitions.size();
}

uint64_t MBBoxPartition::procStartIndex(unsigned int a_proc) const
{
    auto data = m_procMap.find(a_proc);
    if (data == m_procMap.end()) { return numBoxes(); }
    return (*data).second.first;
}
uint64_t MBBoxPartition::procEndIndex(unsigned int a_proc) const
{
    auto data = m_procMap.find(a_proc);
    if (data == m_procMap.end()) { return numBoxes(); }
    return (*data).second.second;
}
uint64_t MBBoxPartition::localIndex(uint64_t a_globalIndex, unsigned int a_proc) const
{
    return a_globalIndex - procStartIndex(a_proc);
}
uint64_t MBBoxPartition::globalIndex(uint64_t a_localIndex, unsigned int a_proc) const
{
    return a_localIndex + procStartIndex(a_proc);
}


inline std::shared_ptr<BoxPartition> MBBoxPartition::blockPartition(BlockIndex a_block) const
{
    return m_blockPartitions[a_block];
}

uint64_t MBBoxPartition::find(Point a_patch, BlockIndex a_block) const
{
    auto blockDomain = m_patchDomain.getBlock(a_block);
    Point image = blockDomain.image(a_patch);
    // if the image isn't in the domain, return failure
    if (!blockDomain.contains(image)) { return numBoxes(); }
    MBPoint patchID(image, a_block);
    auto data = m_indexMap.find(patchID);
    // if the box isn't in the partition, return failure
    if (data == m_indexMap.end()) { return numBoxes(); }
    return m_indexMap[patchID];
}

bool MBBoxPartition::contains(Point a_patch, BlockIndex a_block) const
{
    return (find(a_patch, a_block) < numBoxes());
}

void MBBoxPartition::validatePatch(const MBPoint& patch, std::string methodName) const
{
#ifndef NDEBUG
    PROTO_ASSERT(m_patchDomain.getBlock(patch.block).box().containsPoint(patch.point),
                 "%s | Error: %s is not a valid patch in block %i",
                 methodName, patch.point.str(), patch.block);
#endif
}

Point MBBoxPartition::connectivity(const MBPoint& src, const MBPoint& dst) const
{
    validatePatch(src, "MBBoxPartition::connectivity");
    validatePatch(dst, "MBBoxPartition::connectivity");
    for (auto& [dir, boundInfo] : m_patchAdjacencies[src].boundaries)
    {
        if (boundInfo.adjacentPatches.count(dst) > 0)
        {
            return dir;
        }
    }
    return Point::Zeros();
}

std::set<MBPoint> MBBoxPartition::adjacentPatches(const MBPoint& patch) const {
    validatePatch(patch, "MBBoxPartition::adjacentPatches");
    std::set<MBPoint> adjacent;
    for (auto& [dir, boundInfo] : m_patchAdjacencies[patch].boundaries)
    {
        for (auto adjPatch : boundInfo.adjacentPatches)
        adjacent.insert(adjPatch);
    }
    return adjacent;
}

bool MBBoxPartition::isBlockBoundary(
        uint64_t  a_index,
        Point     a_dir) const
{
    PatchBound key;
    key.index = a_index;
    key.dir = a_dir;
    return m_blockBoundaries.find(key) != m_blockBoundaries.end();
}

bool MBBoxPartition::isBlockBoundary(
        uint64_t        a_index,
        Point           a_dir,
        BlockIndex      a_adjBlock) const
{
    if (a_dir == Point::Zeros()) { return false; }
    PatchBound key;
    key.index = a_index;
    key.dir = a_dir;
    
    if (m_blockBoundaries.find(key) != m_blockBoundaries.end())
    {
        for (auto& adjInd : m_blockBoundaries[key])
        {
            auto& partitionData = m_partition[adjInd];
            BlockIndex block = partitionData.first.block;
            if (block == a_adjBlock) { return true; }
        }
    }
    return false;
}
bool MBBoxPartition::isDomainBoundary(
        uint64_t  a_index,
        Point           a_dir) const
{
    PatchBound key;
    key.index = a_index;
    key.dir = a_dir;
    return m_domainBoundaries.find(key) != m_domainBoundaries.end();
}

bool MBBoxPartition::isInteriorBoundary(
        uint64_t  a_index,
        Point     a_dir) const
{
    if (isDomainBoundary(a_index, a_dir)) { return false; }
    if (isBlockBoundary(a_index, a_dir)) { return false; }
    return true;
}

bool MBBoxPartition::isRefinementBoundary(
    uint64_t index,
    Point dir) const
{
    PatchBound key;
    key.index = index;
    key.dir = dir;
    return m_refinementBoundaries.find(key) != m_domainBoundaries.end();
}

void MBBoxPartition::print()
{
    using Proto::pr_out;
    pr_out() << "MBBoxPartition" << std::endl;
    for (BlockIndex bi = 0; bi < numBlocks(); bi++)
    {
        pr_out() << "Block " << bi << std::endl;
        m_blockPartitions[bi]->print();
    }
    pr_out() << "Global Partition Data" << std::endl;
    pr_out() << "\tNum Boxes:  " << numBoxes() << std::endl;
    pr_out() << "\tNum Blocks: " << numBlocks() << std::endl; 
    pr_out() << "\tNum Procs:  " << numProcs() << std::endl;
    pr_out() << "\tData By Processor: " << std::endl;
    for (auto procData : m_procMap)
    {
        unsigned int proc = procData.first;
        pr_out() << "\t\tProc: " << proc << " | Num Boxes: " << numBoxes(proc);
        pr_out() << " | Start Index: " << procStartIndex(proc);
        pr_out() << " | End Index: " << procEndIndex(proc) << std::endl;
    }
    pr_out() << "\tData Partition: " << std::endl;
    for (auto patchData : m_partition)
    {
        auto patchID = patchData.first;
        auto proc = patchData.second;
        pr_out() << "\t\tBlock: " << patchID.block << " | Patch: " << patchID.point;
        pr_out() << " | Proc: " << proc << std::endl;
    }
}

