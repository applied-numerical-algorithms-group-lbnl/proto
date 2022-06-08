
MBBoxPartition::MBBoxPartition(
        const MBProblemDomain& a_domain,
        const std::vector<std::shared_ptr<BoxPartition>> a_blockPartitions)
{
    define(a_domain, a_blockPartitions);
}

void MBBoxPartition::define(
        const MBProblemDomain& a_domain,
        const std::vector<std::shared_ptr<BoxPartition>> a_blockPartitions)
{
    PROTO_ASSERT(a_domain.size() == a_blockPartitions.size(),
            "MBBoxPartition::Define | Error: Exactly one BoxPartition must be \
            specified for each block.");
    m_domain = a_domain;
    m_blockPartitions.clear();
    for (auto partition : a_blockPartitions)
    {
        m_blockPartitions.push_back(partition);
    }
    loadBalance();
}

void MBBoxPartition::buildLocalMaps()
{
    globalIndex_t N = numBoxes();
    globalIndex_t n = N / numProcs();
    globalIndex_t r = N % numProcs();
    std::vector<globalIndex_t> unassignedBoxes(numProcs(), n);
    for (int ii = 0; ii < r; ii++)
    {
        unassignedBoxes[ii]++;
    }
    procID_t proc = 0;
    for (blockID_t bi = 0; bi < numBlocks(); bi++)
    {
        globalIndex_t numAssigned = 0;
        globalIndex_t blockSize = m_blockPartitions[bi]->numBoxes();
        std::vector<procID_t, globalIndex_t> blockAssignment;
        while (numAssigned < blockSize)
        {
            globalIndex_t ni = min(unassignedBoxes[proc], blockSize - numAssigned);
            procAssign.push_back(std::pair<procID_t, globalIndex_t>(proc,ni));
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
    globalIndex_t globalIndex = 0;
    procID_t currentProc = 0;
    globalIndex_t procRangeStart = 0;
    for (blockID_t bi = 0; bi < numBlocks(); bi++)
    {
        auto patchData = m_blockPartitions[block]->partition();
        for (auto patchData_i : patchData)
        {
            patchID_t patch = patchData_i.first;
            procID_t  proc  = patchData_i.second;
            if (proc != currentProc)
            {
                std::pair<globalIndex_t, globalIndex_t> procRange(procRangeStart, globalIndex);
                m_procMap[currentProc] = procRange;
                currentProc = proc;
                procRangeStart = globalIndex;
            }
            MBPatchID_t mbPatch(patch, block);
            m_partition.push_back(std::pair<MBPatchID_t, procID_t>(mbPatch, proc));
            m_indexMap[mbPatch] = globalIndex;
            globalIndex++;
        }
    }
}

void MBBoxPartition::loadBalance()
{
    buildLocalMaps();
    buildGlobalMaps();
}

bool MBBoxPartition::compatible(const MBBoxPartition& a_rhs)
{
    if (this == &a_rhs) { return true; }
    if (numBlocks() != a_rhs.numBlocks()) { return false; }
    for (blockID_t bi = 0; bi < numBlocks(); bi++)
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
unsigned int MBBoxPartition::numBoxes(procID_t a_proc)
{
    return m_procMap[a_proc].second - m_procMap[a_proc].first;
}

unsigned int MBBoxPartition::numBlocks() const
{
    return m_blockPartitions.size();
}

globalIndex_t MBBoxPartition::procStartIndex(procID_t a_proc) const
{
    auto data = m_procMap.find(a_proc);
    if (data == m_procMap.end()) { return numBoxes(); }
    return (*data).first;
}
globalIndex_t MBBoxPartition::procEndIndex(procID_t a_proc) const
{
    auto data = m_procMap.find(a_proc);
    if (data == m_procMap.end()) { return numBoxes(); }
    return (*data).second;
}

globalIndex_t MBBoxPartition::find(Point a_patch, blockID_t a_block)
{
    auto blockDomain = m_patchDomains.domain(a_block);
    PROTO_ASSERT(blockDomain.contains(a_patch),
            "MBBoxPartition::find | Error: \
            Patch is not contained within the problem domain.");
    Point image = blockDomain.image(a_patch);
    MBPatchID patchID(image, a_block);
    auto data = m_indexMap.find(patchID);
    if (data == m_indexMap.end()) { return numBoxes(); }
    return m_indexMap[patchID];
}


void MBBoxPartition::print() const
{
    pout() << "MBBoxPartition" << std::endl;
    for (blockID_t bi = 0; bi < numBlocks(); bi++)
    {
        pout() << "Block " << bi << std::endl;
        m_blockPartitions[bi]->print();
    }
    pout() << "Global Partition Data" << std::endl;
    pout() << "\tNum Boxes:  " << numBoxes() << std::endl;
    pout() << "\tNum Blocks: " << numBlocks() << std::endl; 
    pout() << "\tNum Procs:  " << numProcs() << std::endl;
    pout() << "\tData By Processor: " << std::endl;
    for (auto procData : m_procMap)
    {
        procID_t proc = procData.first;
        pout() << "\t\tProc: " << proc << " | Num Boxes: " << numBoxes(proc);
        pout() << " | Start Index: " << procStartIndex(proc);
        pout() << " | End Index: " << procEndIndex(proc) << std::endl;
    }
    pout() << "\tData Partition: " << std::endl;
    for (auto patchData : m_partition)
    {
        auto patchID = patchData.first;
        auto proc = patchData.second;
        pout() << "\t\tBlock: " << patchID.second << " | Patch: " << patchID.first;
        pout() << " | Proc: " << proc << std::endl;
    }
}

