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

void BoxPartition::assign(
        std::vector<Point>& a_patches,
        unsigned int a_globalIndex,
        unsigned int a_proc,
        unsigned int a_num)
{
    int global = a_globalIndex;
    for (int ii = 0; ii < a_num; ii++, global++)
    {
        Point patch = a_patches[global];
        uint64_t mortonIndex = Morton::index(patch);
        m_indexMap.insert(std::pair<uint64_t, int>(mortonIndex, global));
        auto procAssign = pair<Point, unsigned int>(patch, a_proc);
        m_partition.push_back(procAssign);
    }
}

template<typename... vals>
void BoxPartition::loadBalance(
        std::vector<Point>& a_patches,
        vals... a_args)
{

}


void BoxPartition::loadBalance(
        std::vector<Point>& a_patches,
        unsigned int a_startProc,
        unsigned int a_endProc)
{
    // convert vector to a set to preclude the possibility of duplicates
    m_procIndices.clear();
    m_indexMap.clear();
    m_partition.clear();
    m_startProc = a_startProc;
    if (a_startProc >= numProc() || a_startProc >= a_endProc) { return; } //nothing else to do
    
    int nsegs  = a_endProc - a_startProc;
    int length = a_patches.size() / nsegs;
    int remainder = a_patches.size() % nsegs;
    int globalIndex = 0;
    for (int procNum = a_startProc; procNum < a_startProc + nsegs; procNum++)
    {
        m_procIndices.push_back(globalIndex);
        int segLength = length;
        if (procNum < remainder)
        {
            segLength = length+1;
        }
        assign(a_patches, globalIndex, procNum, segLength);
        globalIndex += segLength;
        /*
        for (int ii = 0; ii < segLength; ii++)
        {
            Point patch = a_patches[globalIndex];
            uint64_t mortonIndex = Morton::index(patch);
            m_indexMap.insert(std::pair<uint64_t, int>(mortonIndex, globalIndex));
            auto procAssign = pair<Point, unsigned int>(patch, procNum);
            m_partition.push_back(procAssign);
            globalIndex++;
        }
        */
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

unsigned int BoxPartition::numProcs() const
{
    if (m_procIndices.size() == 0) { return 0; }
    return m_procIndices.size()-1;
}

unsigned int BoxPartition::numBoxes() const
{
    return m_partition.size();
}

unsigned int BoxPartition::numBoxes(unsigned int a_proc) const
{
    return procEndIndex(a_proc) - procStartIndex(a_proc);
}

unsigned int BoxPartition::procStartIndex(unsigned int a_proc) const
{
    if (a_proc < m_startProc || a_proc >= m_startProc + numProcs())
    {
        return numBoxes();
    }
    unsigned int localProc = a_proc - m_startProc;
    return m_procIndices[localProc];
}

unsigned int BoxPartition::procEndIndex(unsigned int a_proc) const
{
    if (a_proc < m_startProc || a_proc >= m_startProc + numProcs())
    {
        return numBoxes();
    }
    unsigned int localProc = a_proc - m_startProc;
    return m_procIndices[localProc + 1];
}

unsigned int BoxPartition::patchIndex(Point a_pt)
{
    PROTO_ASSERT(m_patchDomain.contains(a_pt),
            "BoxPartition::patchIndex | Error: \
            The input patch is not contained in the domain or any of its periodic images.");
    Point image = m_patchDomain.image(a_pt);
    return m_indexMap[Morton::index(image)];
}

unsigned int BoxPartition::find(Point a_pt)
{
    PROTO_ASSERT(m_patchDomain.contains(a_pt),
            "BoxPartition::patchIndex | Error: \
            The input patch is not contained in the domain or any of its periodic images.");
    Point image = m_patchDomain.image(a_pt);
    auto mortonIndex = Morton::index(image);
    auto iter = m_indexMap.find(mortonIndex);
    if (iter == m_indexMap.end())
    {
        return numBoxes();
    }
    return m_indexMap[mortonIndex];
}

void BoxPartition::print() const
{
    pout() << "Box Partition: " << std::endl;
    pout() << "\tPatch Domain: " << m_patchDomain.box() << std::endl;
    pout() << "\tNum Boxes: " << numBoxes() << std::endl;
    pout() << "\tStart Proc: " << startProc() << std::endl;
    pout() << "\tEnd Proc: " << endProc() << std::endl;
    pout() << "\tNum Procs: " << numProcs() << std::endl;
    pout() << "\tData By Processor: " << std::endl;
    for (int ii = startProc(); ii < endProc(); ii++)
    {
        pout() << "\t\tProc: " << ii << " | NumBoxes: " << numBoxes(ii);
        pout() << " | Start Index: " << procStartIndex(ii);
        pout() << " | End Index: " << procEndIndex(ii) << std::endl;
    }
    pout() << "\tData Partition: " << std::endl;
    for (auto item : m_partition)
    {
        pout() << "\t\tPatch: " << item.first << " | Proc: " << item.second << std::endl; 
    }
}

template<typename... vals>
inline void BoxPartition::unpack(
        std::vector<Point>& a_patches,
        unsigned int a_globalIndex,
        unsigned int a_proc,
        unsigned int a_num,
        vals... a_args)
{
    assign(a_patches, a_globalIndex, a_proc, a_num); 
    unpack(a_patches, a_globalIndex + a_num, a_args...);
}

template<>
inline void BoxPartition::unpack(
        std::vector<Point>& a_patches,
        unsigned int a_globalIndex,
        unsigned int a_proc,
        unsigned int a_num)
{
    assign(a_patches, a_globalIndex, a_proc, a_num); 
}
