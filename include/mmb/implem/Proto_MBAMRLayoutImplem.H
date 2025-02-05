MBAMRLayout::MBAMRLayout(
        const MBProblemDomain&          a_coarseDomain,
        const std::vector<Point>&       a_coarseBoxSizes,
        const std::vector<Point>&       a_refRatios)
{
    define(a_coarseDomain, a_coarseBoxSizes, a_refRatios);
}

MBAMRLayout::MBAMRLayout(
        const MBProblemDomain&          a_coarseDomain,
        const std::vector<MBPoint>& a_coarsePatches,
        const std::vector<Point>&       a_coarseBoxSizes,
        const std::vector<Point>&       a_refRatios)
{
    define(a_coarseDomain, a_coarsePatches, a_coarseBoxSizes, a_refRatios);
}

MBAMRLayout::MBAMRLayout(
        const MBProblemDomain&      a_coarseDomain,
        Point                       a_coarseBoxSize,
        Point                       a_refRatio,
        int                         a_numLevels)
{
    std::vector<Point> coarseBoxSizes;
    std::vector<Point> refRatios;
    for (BlockIndex bi = 0; bi < a_coarseDomain.numBlocks(); bi++)
    {
        coarseBoxSizes.push_back(a_coarseBoxSize);
    }
    for (int li = 0; li < a_numLevels-1; li++)
    {
        refRatios.push_back(a_refRatio);
    }
    define(a_coarseDomain, coarseBoxSizes, refRatios);
}
void MBAMRLayout::define(
        const MBProblemDomain&          a_coarseDomain,
        const std::vector<Point>&       a_coarseBoxSizes,
        const std::vector<Point>&       a_refRatios)
{
    std::vector<MBPoint> patches;
    for (BlockIndex bi = 0; bi < a_coarseDomain.numBlocks(); bi++)
    {
        Box patchBox = a_coarseDomain.getBlock(bi).box().coarsen(a_coarseBoxSizes[bi]);
        for (auto pi : patchBox)
        {
            patches.push_back(MBPoint(pi, bi));
        }
    }
    define(a_coarseDomain, patches, a_coarseBoxSizes, a_refRatios);
}

void MBAMRLayout::define(
        const MBProblemDomain&          a_coarseDomain,
        const std::vector<MBPoint>& a_coarsePatches,
        const std::vector<Point>&       a_coarseBoxSizes,
        const std::vector<Point>&       a_refRatios)
{
    auto crsPatchDomain = a_coarseDomain.coarsen(a_coarseBoxSizes);
    m_partition = std::make_shared<MBAMRBoxPartition>(crsPatchDomain, a_refRatios);
    m_refRatios = a_refRatios;
    m_layouts.clear();
    std::vector<Point> boxSizes = a_coarseBoxSizes;
    std::vector<std::vector<DisjointBoxLayout>> gridLayouts;
    gridLayouts.resize(a_coarseDomain.numBlocks());
    for (int li = 0; li <= a_refRatios.size(); li++)
    {
        MBDisjointBoxLayout layout(m_partition->getLevel(li), boxSizes);
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
    m_partition->loadBalance();

}

void MBAMRLayout::setPatches(
        const std::vector<MBPoint>& a_patches)
{
    std::vector<MBPoint> patches = a_patches;
    std::vector<MBPoint> finePatches;
    for (int li = 0; li < numLevels(); li++)
    {
        setPatches(li, patches);
        if (li < m_refRatios.size())
        {
            finePatches.clear();
            for (auto pi : patches)
            {
                Box b(pi.point, pi.point);
                b = b.refine(m_refRatios[li]);
                for (auto pj : b)
                {
                    finePatches.push_back(MBPoint(pj, pi.block));
                }
            }
            patches = finePatches;
        }
    }
}

void MBAMRLayout::setPatches(
        int a_level,
        const std::vector<MBPoint>& a_patches)
{
    auto part = m_partition->getLevel(a_level);
    part->define(a_patches);
    part->loadBalance();
}


const AMRGrid& MBAMRLayout::getBlock(BlockIndex a_block) const
{
    PROTO_ASSERT(a_block < numBlocks(),
        "MBAMRLayout::block | Error: block index out of bounds");
    return m_grids[a_block];
}
AMRGrid& MBAMRLayout::getBlock(BlockIndex a_block)
{
    PROTO_ASSERT(a_block < numBlocks(),
        "MBAMRLayout::block | Error: block index out of bounds");
    return m_grids[a_block];
}
const MBDisjointBoxLayout& MBAMRLayout::getLevel(int a_level) const
{
    PROTO_ASSERT(a_level < numLevels(),
        "MBAMRLayout::level | Error: level index out of bounds");
    return m_layouts[a_level];
}

MBDisjointBoxLayout& MBAMRLayout::getLevel(int a_level)
{
    PROTO_ASSERT(a_level < numLevels(),
        "MBAMRLayout::level | Error: level index out of bounds");
    return m_layouts[a_level];
}

MBIterator MBAMRLayout::begin(int a_level) const
{
    PROTO_ASSERT(a_level < numLevels(),
        "MBAMRLayout::level | Error: level index out of bounds");
    return m_layouts[a_level].begin();    

}

MBIterator MBAMRLayout::end(int a_level) const
{
    PROTO_ASSERT(a_level < numLevels(),
        "MBAMRLayout::level | Error: level index out of bounds");
    return m_layouts[a_level].end();    

}

MBIndex MBAMRLayout::find(Point a_point, BlockIndex a_block, int a_level) const
{
    PROTO_ASSERT(a_level < numLevels(),
        "MBAMRLayout::level | Error: level index out of bounds");

    return m_layouts[a_level].find(a_point, a_block);
}


