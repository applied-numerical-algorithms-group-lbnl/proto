MBBoxPartition::MBBoxPartition(
        const MBProblemDomain& a_patchDomain)
{
    std::vector<MBPatchID_t> patches;
    // creates the BoxPartition with no patches in it
    define(a_patchDomain, patches);
}

MBBoxPartition::MBBoxPartition(
        const MBProblemDomain& a_patchDomain,
        const std::vector<MBPatchID_t>& a_patches)
{
    define(a_patchDomain, a_patches);
}

void MBBoxPartition::define(
        const MBProblemDomain& a_patchDomain,
        const std::vector<MBPatchID_t>& a_patches)
{
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
        const std::vector<MBPatchID_t>& a_patches)
{
    std::vector<std::vector<Point>> blockPatches;
    blockPatches.resize(m_patchDomain.numBlocks());
    for (auto pi : a_patches)
    {
        blockPatches[pi.second].push_back(pi.first);
    }
    
    for (unsigned int bi = 0; bi < m_patchDomain.size(); bi++)
    {
        m_blockPartitions[bi]->define(
                m_patchDomain.getBlock(bi), blockPatches[bi], 0, numProc());
    }
    loadBalance();
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
            MBPatchID_t mbPatch(patch, bi);
            m_partition.push_back(std::pair<MBPatchID_t, unsigned int>(mbPatch, proc));
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
    findBoundaries();
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
        int   block_i = partitionData.first.second;
        auto  patch_i = partitionData.first.first;
        for (auto dir : Box::Kernel(1))
        {
            if (dir == Point::Zeros()) { continue; }
            auto patch_ji = patch_i + dir;
            if (m_patchDomain.getBlock(block_i).box().contains(patch_ji))
            {
                // patch_j is interior to the block
                continue; 
            }
            PatchBound key; key.index = ind_i; key.dir = dir;
            bool foundBlockBoundary = false; 
            for (BlockIndex block_j = 0; block_j < numBlocks(); block_j++)
            {
                if (block_j == block_i) { continue; }
                Point conn_ij = graph.connectivity(block_i, block_j);
                if (conn_ij.dot(dir) <= 0) { continue; }
                Point conn_ji = graph.reverseDir(block_i, block_j, conn_ij);
                auto R = graph.rotation(block_i, conn_ij, block_j);

                Box localBoundPatches = m_patchDomain.getBlock(block_i).box().adjacent(conn_ij);
                Box adjBoundPatches = m_patchDomain.getBlock(block_j).box().edge(conn_ji);

                Point patch_jj = R.rotateCell(patch_ji, localBoundPatches, adjBoundPatches);

                if (!localBoundPatches.contains(patch_ji)) { continue; }

                PROTO_ASSERT(adjBoundPatches.contains(patch_jj),
                        "MBBoxPartition::findBoundaries | Error: Data Corruption.");
                
                auto ind_j = find(patch_jj, block_j);
                if (ind_j < numBoxes())
                {
                    m_blockBoundaries[key].push_back(ind_j);
                    foundBlockBoundary = true;
                }
            }
            if (!foundBlockBoundary)
            {
                m_domainBoundaries.insert(key);
            }
        }
    }
}

BlockIndex MBBoxPartition::block(uint64_t a_index) const
{
    auto& partitionData = m_partition[a_index];
    return partitionData.first.second;
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
    PROTO_ASSERT(blockDomain.contains(a_patch),
            "MBBoxPartition::find | Error: \
            Patch is not contained within the problem domain.");
    Point image = blockDomain.image(a_patch);
    MBPatchID_t patchID(image, a_block);
    auto data = m_indexMap.find(patchID);
    if (data == m_indexMap.end()) { return numBoxes(); }
    return m_indexMap[patchID];
}

bool MBBoxPartition::isBlockBoundary(
        uint64_t  a_index,
        Point           a_dir) const
{
    if (a_dir == Point::Zeros()) { return false; }
    for (BlockIndex bi = 0; bi < numBlocks(); bi++)
    {
        if (isBlockBoundary(a_index, a_dir, bi)) { return true; }
    }
    return false;
}

bool MBBoxPartition::isBlockBoundary(
        uint64_t  a_index,
        Point           a_dir,
        BlockIndex      a_block) const
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
            BlockIndex block = partitionData.first.second;
            if (block == a_block) { return true; }
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
        Point           a_dir) const
{
    if (isDomainBoundary(a_index, a_dir)) { return false; }
    if (isBlockBoundary(a_index, a_dir)) { return false; }
    return true;
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
        pr_out() << "\t\tBlock: " << patchID.second << " | Patch: " << patchID.first;
        pr_out() << " | Proc: " << proc << std::endl;
    }
}

