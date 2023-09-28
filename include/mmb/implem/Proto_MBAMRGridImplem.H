MBAMRGrid::MBAMRGrid(
        const MBProblemDomain&          a_coarseDomain,
        const std::vector<Point>&       a_coarseBoxSizes,
        const std::vector<Point>&       a_refRatios)
    : m_partition(a_coarseDomain.coarsen(a_coarseBoxSizes), a_refRatios)
{
    define(a_coarseDomain, a_coarseBoxSizes, a_refRatios);
}

MBAMRGrid::MBAMRGrid(
        const MBProblemDomain&          a_coarseDomain,
        const std::vector<MBPatchID_t>& a_coarsePatches,
        const std::vector<Point>&       a_coarseBoxSizes,
        const std::vector<Point>&       a_refRatios)
    : m_partition(a_coarseDomain.coarsen(a_coarseBoxSizes), a_refRatios)
{
    define(a_coarseDomain, a_coarsePatches, a_coarseBoxSizes, a_refRatios);
}

void MBAMRGrid::define(
        const MBProblemDomain&          a_coarseDomain,
        const std::vector<Point>&       a_coarseBoxSizes,
        const std::vector<Point>&       a_refRatios)
{
    std::vector<MBPatchID_t> patches;
    for (unsigned int bi = 0; bi < a_coarseDomain.numBlocks(); bi++)
    {
        Box patchBox = a_coarseDomain.blockDomain(bi).box().coarsen(a_coarseBoxSizes[bi]);
        for (auto pi : patchBox)
        {
            patches.push_back(MBPatchID_t(pi, bi));
        }
    }
    define(a_coarseDomain, patches, a_coarseBoxSizes, a_refRatios);
}

void MBAMRGrid::define(
        const MBProblemDomain&          a_coarseDomain,
        const std::vector<MBPatchID_t>& a_coarsePatches,
        const std::vector<Point>&       a_coarseBoxSizes,
        const std::vector<Point>&       a_refRatios)
{
    m_refRatios = a_refRatios;
    m_layouts.clear();
    std::vector<Point> boxSizes = a_coarseBoxSizes;
    std::vector<std::vector<DisjointBoxLayout>> gridLayouts;
    gridLayouts.resize(a_coarseDomain.numBlocks());
    for (int li = 0; li <= a_refRatios.size(); li++)
    {
        MBDisjointBoxLayout layout(m_partition.getLevel(li), boxSizes);
        m_layouts.push_back(layout);
        for (int bi = 0; bi < a_coarseDomain.numBlocks(); bi++)
        {
            gridLayouts[bi].push_back(layout.getBlock(bi));
        }
        if (li == a_refRatios.size()) { continue; }
    }
    m_grids.clear();
    for (int bi = 0; bi < a_coarseDomain.numBlocks(); bi++)
    {
        AMRGrid grid(gridLayouts[bi], a_refRatios, a_refRatios.size()+1);
        m_grids.push_back(grid);
    }
    setPatches(a_coarsePatches);
    m_partition.loadBalance();

}

void MBAMRGrid::setPatches(
        const std::vector<MBPatchID_t>& a_patches)
{
    std::vector<MBPatchID_t> patches = a_patches;
    std::vector<MBPatchID_t> finePatches;
    for (int li = 0; li < numLevels(); li++)
    {
        setPatches(li, patches);
        if (li < m_refRatios.size())
        {
            finePatches.clear();
            for (auto pi : patches)
            {
                Box b(pi.first, pi.first);
                b = b.refine(m_refRatios[li]);
                for (auto pj : b)
                {
                    finePatches.push_back(MBPatchID_t(pj, pi.second));
                }
            }
            patches = finePatches;
        }
    }
}

void MBAMRGrid::setPatches(
        unsigned int a_level,
        const std::vector<MBPatchID_t>& a_patches)
{
    auto part = m_partition.getLevel(a_level);
    part->define(a_patches);
    part->loadBalance();
}


const AMRGrid& MBAMRGrid::getBlock(unsigned int a_block) const
{
    PROTO_ASSERT(a_block < numBlocks(),
        "MBAMRGrid::getBlock | Error: block index out of bounds");
    return m_grids[a_block];
}
AMRGrid& MBAMRGrid::getBlock(unsigned int a_block)
{
    PROTO_ASSERT(a_block < numBlocks(),
        "MBAMRGrid::getBlock | Error: block index out of bounds");
    return m_grids[a_block];
}
const MBDisjointBoxLayout& MBAMRGrid::getLevel(unsigned int a_level) const
{
    PROTO_ASSERT(a_level < numLevels(),
        "MBAMRGrid::getLevel | Error: level index out of bounds");
    return m_layouts[a_level];
}

MBDisjointBoxLayout& MBAMRGrid::getLevel(unsigned int a_level)
{
    PROTO_ASSERT(a_level < numLevels(),
        "MBAMRGrid::getLevel | Error: level index out of bounds");
    return m_layouts[a_level];
}
