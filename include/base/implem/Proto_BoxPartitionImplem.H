BoxPartition::BoxPartition(
        const ProblemDomain& a_patchDomain,
        const std::vector<Point>& a_patches)
{
    define(a_patchDomain, a_patches, 0, numProc());
}

BoxPartition::BoxPartition(
        const ProblemDomain& a_patchDomain,
        const std::vector<Point>& a_patches,
        unsigned int a_startProc,
        unsigned int a_endProc)
{
    define(a_patchDomain, a_patches, a_startProc, a_endProc);
}

void BoxPartition::define(
        const ProblemDomain& a_patchDomain,
        const std::vector<Point>& a_patches,
        unsigned int a_startProc,
        unsigned int a_endProc)
{
    for (auto patch : a_patches)
    {
        PROTO_ASSERT(a_patchDomain.box().contains(patch),
                "BoxPartion::Constructor | Error: \
                All patches must be contained in the patchDomain.");
    }

    m_patchDomain = a_patchDomain;
    std::vector<Point> sortedPatches = a_patches;
    Morton::sort(sortedPatches);
    loadBalance(sortedPatches, a_startProc, a_endProc);
}

void BoxPartition::loadBalance(
        std::vector<Point>& a_patches,
        unsigned int a_startProc,
        unsigned int a_endProc)
{
    PROTO_ASSERT(a_startProc < numProc(),
            "BoxPartition::loadBalance | Error: \
            start proc %u exceeds the number of processes %u.",
            a_startProc, numProc());
    PROTO_ASSERT(a_startProc < a_endProc,
            "BoxPartition::loadBalance | Error: \
            start proc %u is larger than the end proc %u.",
            a_startProc, a_endProc);
    
    // convert vector to a set to preclude the possibility of duplicates
    int nsegs  = a_endProc - a_startProc;
    int length = a_patches.size() / nsegs;
    int remainder = a_patches.size() % nsegs;
    int globalIndex = 0;
    m_procIndices.clear();
    m_indexMap.clear();
    m_partition.clear();
    for (int procNum = 0; procNum < nsegs; procNum++)
    {
        m_procIndices.push_back(globalIndex);
        int segLength = length;
        if (procNum < remainder)
        {
            segLength = length+1;
        }
        for (int ii = 0; ii < segLength; ii++)
        {
            Point patch = a_patches[globalIndex];
            uint64_t mortonIndex = Morton::index(patch);
            m_indexMap.insert(std::pair<uint64_t, int>(mortonIndex, globalIndex));
            auto procAssign = pair<Point, unsigned int>(patch, procNum);
            m_partition.push_back(procAssign);
            globalIndex++;
        }
    }
    m_procIndices.push_back(a_patches.size());
}

bool BoxPartition::compatible(const BoxPartition& a_rhs)
{
    if (this == &a_rhs) { return true; }
    auto lhsPartition = m_partition;
    auto rhsPartition = a_rhs.partition();
    if (lhsPartition.size() != rhsPartition.size())
    {
        return false;
    } 
    for (int ii = 0; ii < lhsPartition.size(); ii++)
    {
        auto& lhs_i = lhsPartition[ii];
        auto& rhs_i = rhsPartition[ii];
        if (lhs_i.first != rhs_i.first) 
        {
            return false;
        }
        if (lhs_i.second != rhs_i.second) 
        {
            return false; 
        }
    }
    return true;
}
