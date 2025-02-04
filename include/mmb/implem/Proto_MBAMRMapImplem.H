
template<template<MemType> typename MAP, MemType MEM>
MBAMRMap<MAP,MEM>::MBAMRMap(
        const MBAMRLayout&            a_grid,
        const Array<Point, DIM+1>&  a_ghost,
        unsigned int                a_focalBlock)
{
    define(a_grid, a_ghost, a_focalBlock);
}

template<template<MemType> typename MAP, MemType MEM>
MBAMRMap<MAP,MEM>::MBAMRMap(
        const MBAMRLayout&            a_grid,
        Point                       a_ghost,
        unsigned int                a_focalBlock)
{
    define(a_grid, a_ghost, a_focalBlock);
}

template<template<MemType> typename MAP, MemType MEM>
void MBAMRMap<MAP,MEM>::define(
        const MBAMRLayout&            a_grid,
        Point                       a_ghost,
        unsigned int                a_focalBlock)
{
    Array<Point, DIM+1> ghost;
    ghost.fill(a_ghost + Point::Ones(2)); //assumes 4th - 5th order interpolation
    ghost[0] = a_ghost;
    define(a_grid, ghost, a_focalBlock);
}

template<template<MemType> typename MAP, MemType MEM>
void MBAMRMap<MAP,MEM>::define(
        const MBAMRLayout&            a_grid,
        const Array<Point, DIM+1>&  a_ghost,
        unsigned int                a_focalBlock)
{
    m_grid = &a_grid;
    m_maps.clear();
    std::vector<std::shared_ptr<MBLevelBoxData<double, DIM, MEM, PR_NODE>>> xLevels;
    std::vector<std::shared_ptr<MBLevelBoxData<double, 1, MEM, PR_CELL>>> jLevels;

    for (int li = 0; li < m_grid->numLevels(); li++)
    {
        auto levelMap = std::make_shared<MBLevelMap<MAP,MEM>>(
                m_grid->getLevel(li), a_ghost, a_focalBlock);
        m_maps.push_back(levelMap);
        xLevels.push_back(levelMap->mapData());
        jLevels.push_back(levelMap->jacobianData());
    }
    m_X.define(a_grid, xLevels);
    m_J.define(a_grid, jLevels);
}

template<template<MemType> typename MAP, MemType MEM>
void MBAMRMap<MAP,MEM>::initialize()
{
    for (int li = 0; li < m_grid->numLevels(); li++)
    {
        m_maps[li]->initialize();
    }
}

template<template<MemType> typename MAP, MemType MEM>
const MBLevelMap<MAP,MEM>& MBAMRMap<MAP,MEM>::getLevel(unsigned int a_level) const
{
    PROTO_ASSERT(a_level < m_grid->numLevels(),
            "MBLevelMap::getLevel | Error: Level %u is out of bounds", a_level);
    return *m_maps[a_level];
}
