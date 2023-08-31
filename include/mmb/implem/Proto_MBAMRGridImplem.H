MBAMRGrid::MBAMRGrid(
        const MBProblemDomain&          a_coarseDomain,
        const std::vector<Point>&       a_coarseBoxSizes,
        const std::vector<Point>&       a_refRatios)
    : m_partition(a_coarseDomain, a_refRatios)
{
    define(a_coarseDomain, a_coarseBoxSizes, a_refRatios);
}

MBAMRGrid::MBAMRGrid(
        const MBProblemDomain&          a_coarseDomain,
        const std::vector<MBPatchID_t>& a_coarsePatches,
        const std::vector<Point>&       a_coarseBoxSizes,
        const std::vector<Point>&       a_refRatios)
    : m_partition(a_coarseDomain, a_refRatios)
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
    
}

const AMRGrid& MBAMRGrid::getBlock(unsigned int a_block) const
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
