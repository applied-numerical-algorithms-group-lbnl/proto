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

BoxPartition::BoxPartition(
        const ProblemDomain& a_patchDomain,
        const std::vector<std::pair<Point, unsigned int>>& a_patches)
{
    define(a_patchDomain, a_patches);
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
    m_doLoadBalance = true; 
    m_procSorted = true;
    m_patchDomain = a_patchDomain;
    std::vector<Point> sortedPatches = a_patches;
    Morton::sort(sortedPatches);
    loadBalance(sortedPatches, a_startProc, a_endProc);
}

bool BoxPartition::checkProcSorting(const std::vector<std::pair<Point, unsigned int>> a_data)
{
    if (a_data.size() == 0) { return true; }
    unsigned int proc = a_data[0].second;
    for (auto data : a_data)
    {
        if (data.second < proc) { return false; }
        if (data.second > proc) { proc = data.second; }
    }
    return true;
}

void BoxPartition::define(
        const ProblemDomain& a_patchDomain,
        const std::vector<std::pair<Point, unsigned int>>& a_patches)
{
    m_procSorted = checkProcSorting(a_patches);
    m_doLoadBalance = false;
    m_patchDomain = a_patchDomain;
    m_indexMap.clear();
    m_partition.clear();
    m_partition = a_patches;
    m_procMap.clear();
    m_procMapUnsorted.clear();
    if (a_patches.size() == 0) { return; } //nothing else to do
    
    unsigned int procIndexStart = 0;
    unsigned int currentProc = a_patches[0].second;
    for (int global = 0; global < a_patches.size(); global++)
    {
        auto patch = a_patches[global].first;
        m_indexMap.insert(std::pair<uint64_t, int>(Morton::index(patch), global));
        
        if (m_procSorted) {
            if (a_patches[global].second > currentProc)
            {
                m_procMap[currentProc] =
                    std::pair<unsigned int, unsigned int>(procIndexStart, global);
                procIndexStart = global;
                currentProc = a_patches[global].second;
            }
        } else {
            m_procMapUnsorted[a_patches[0].second].push_back(global);
        }
    }
}



void BoxPartition::assign(
        std::vector<Point>& a_patches,
        unsigned int a_globalIndex,
        int          a_proc,
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
void BoxPartition::loadAssign(
        std::vector<Point>& a_patches,
        vals... a_args)
{
    if (!m_doLoadBalance) {return; }
    m_procMap.clear();
    m_indexMap.clear();
    m_partition.clear();
    unpack(a_patches, 0, a_args...);
}

void BoxPartition::loadAssign(
        std::vector<std::pair<int, unsigned int>>& a_assignment)
{
    if (!m_doLoadBalance) {return;}
    std::vector<Point> patches;
    for (auto data : m_partition)
    {
        patches.push_back(data.first);
    }
    loadAssign(a_assignment, patches);
}
void BoxPartition::loadAssign(
        std::vector<std::pair<int, unsigned int>>& a_assignment,
        std::vector<Point>& a_patches)
{
    if (!m_doLoadBalance) {return;}
    m_procMap.clear();
    m_indexMap.clear();
    m_partition.clear();
    int globalIndex = 0;
    for (auto item : a_assignment)
    {
        int proc = item.first;
        unsigned int numBoxes = item.second;
        m_procMap[proc] = std::pair<int, int>(globalIndex, globalIndex + numBoxes);
        assign(a_patches, globalIndex, proc, numBoxes);
        globalIndex += numBoxes;
    }
    if (globalIndex < a_patches.size())
    {
        MayDay<void>::Warning("BoxPartition::loadAssign | Warning: \
                Incomplete load assignment.");
    }
}

void BoxPartition::loadBalance(
        std::vector<Point>& a_patches,
        unsigned int a_startProc,
        unsigned int a_endProc)
{
    if (!m_doLoadBalance) {return;}
    m_procMap.clear();
    m_indexMap.clear();
    m_partition.clear();
    if (a_startProc >= numProc() || a_startProc >= a_endProc) { return; } //nothing else to do
    
    int nsegs  = a_endProc - a_startProc;
    int length = a_patches.size() / nsegs;
    int remainder = a_patches.size() % nsegs;
    int globalIndex = 0;
    for (int procNum = a_startProc; procNum < a_startProc + nsegs; procNum++)
    {
        int segLength = length;
        if (procNum < remainder)
        {
            segLength = length+1;
        }
        m_procMap[procNum] = std::pair<int, int>(globalIndex, globalIndex + segLength);
        assign(a_patches, globalIndex, procNum, segLength);
        globalIndex += segLength;
    }
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
    if (m_procSorted)
    {
        return m_procMap.size();
    } else {
        return m_procMapUnsorted.size();
    }
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
    PROTO_ASSERT(m_procSorted,
            "BoxPartition::procStartIndex() | Error: This function is invalid for unsorted objects.");
    auto data = m_procMap.find(a_proc);
    if (data == m_procMap.end()) { return numBoxes(); }
    return (*data).second.first;
}

unsigned int BoxPartition::procEndIndex(unsigned int a_proc) const
{
    PROTO_ASSERT(m_procSorted,
            "BoxPartition::procEndIndex() | Error: This function is invalid for unsorted objects.");
    auto data = m_procMap.find(a_proc);
    if (data == m_procMap.end()) { return numBoxes(); }
    return (*data).second.second;
}

unsigned int BoxPartition::globalIndex(unsigned int a_localIndex, unsigned int a_proc)
{
    if (m_procSorted)
    {
        return procStartIndex(a_proc) + a_localIndex;
    } else {
        return m_procMapUnsorted[a_proc][a_localIndex];
    }
}

unsigned int BoxPartition::localIndex(unsigned int a_globalIndex, unsigned int a_proc)
{
    if (m_procSorted)
    {
        return a_globalIndex - procStartIndex(a_proc);
    } else {
        auto indices = m_procMapUnsorted[a_proc];
        for (int ii = 0; ii < indices.size(); ii++)
        {
            if (indices[ii] == a_globalIndex) {return ii; }
        }
        MayDay<void>::Error("BoxPartition::localIndex | Error: global index is not present on this proc");
    }
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
    pr_out() << "Box Partition: " << std::endl;
    pr_out() << "\tPatch Domain: " << m_patchDomain.box() << std::endl;
    pr_out() << "\tNum Boxes: " << numBoxes() << std::endl;
    pr_out() << "\tNum Procs: " << numProcs() << std::endl;
    pr_out() << "\tData By Processor: " << std::endl;
    if (m_procSorted)
    {
        for (auto procData : m_procMap)
        {
            unsigned int proc = procData.first;
            pr_out() << "\t\tProc: " << proc << " | NumBoxes: " << numBoxes(proc);
            pr_out() << " | Start Index: " << procStartIndex(proc);
            pr_out() << " | End Index: " << procEndIndex(proc) << std::endl;
        }
    }
    pr_out() << "\tData Partition: " << std::endl;
    for (auto item : m_partition)
    {
        pr_out() << "\t\tPatch: " << item.first << " | Proc: " << item.second << std::endl; 
    }
}

template<typename... vals>
inline void BoxPartition::unpack(
        std::vector<Point>& a_patches,
        unsigned int a_globalIndex,
        int          a_proc,
        unsigned int a_num,
        vals... a_args)
{
    m_procMap[a_proc] = std::pair<int, int>(a_globalIndex, a_globalIndex + a_num);
    assign(a_patches, a_globalIndex, a_proc, a_num); 
    unpack(a_patches, a_globalIndex + a_num, a_args...);
}

template<>
inline void BoxPartition::unpack(
        std::vector<Point>& a_patches,
        unsigned int a_globalIndex,
        int          a_proc,
        unsigned int a_num)
{
    m_procMap[a_proc] = std::pair<int, int>(a_globalIndex, a_globalIndex + a_num);
    assign(a_patches, a_globalIndex, a_proc, a_num); 
}

