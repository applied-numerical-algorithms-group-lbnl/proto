template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBAMRData<T,C,MEM,CTR>::MBAMRData(){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBAMRData<T,C,MEM,CTR>::MBAMRData(
        const MBAMRGrid& a_grid,
        Array<Point, DIM+1> a_ghost,
        Point a_boundGhost)
{
    define(a_grid, a_ghost, a_boundGhost);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::define(
        const MBAMRGrid& a_grid,
        Array<Point, DIM+1> a_ghost,
        Point a_boundGhost)
{
    m_ghost = a_ghost;
    m_grid = &a_grid;
    m_levelData.clear();
    m_amrData.clear();
    for (int li = 0; li < a_grid.numLevels(); li++)
    {
        auto& layout = a_grid.getLevel(li);
        m_levelData.push_back(
                std::make_shared<MBLevelBoxData<T,C,MEM,CTR>>(layout, a_ghost, a_boundGhost));
    }
    //TODO: Create AMRData
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::setVal(T a_value)
{
    for (int li = 0; li < numLevels(); li++)
    {
        getLevel(li).setVal(a_value);
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::setRandom(T a_low, T a_high)
{
    for (int li = 0; li < numLevels(); li++)
    {
        getLevel(li).setRandom(a_low, a_high);
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::exchange()
{
    for (int li = 0; li < numLevels(); li++)
    {
        getLevel(li).exchange();
    }
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBLevelBoxData<T,C,MEM,CTR>&
MBAMRData<T,C,MEM,CTR>::getLevel(unsigned int a_level)
{
    PROTO_ASSERT(a_level < numLevels(),
            "MBAMRData::getLevel | Error: level %u is out of bounds", a_level);
    return *(m_levelData[a_level]);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const MBLevelBoxData<T,C,MEM,CTR>&
MBAMRData<T,C,MEM,CTR>::getLevel(unsigned int a_level) const
{
    PROTO_ASSERT(a_level < numLevels(),
            "MBAMRData::getLevel | Error: level %u is out of bounds", a_level);
    return *(m_levelData[a_level]);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
AMRData<T,C,MEM,CTR>&
MBAMRData<T,C,MEM,CTR>::getBlock(unsigned int a_block)
{
    PROTO_ASSERT(a_block < numBlocks(),
            "MBAMRData::getBLock | Error: block %u is out of bounds", a_block);
    return *(m_amrData[a_block]);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const AMRData<T,C,MEM,CTR>&
MBAMRData<T,C,MEM,CTR>::getBlock(unsigned int a_block) const
{
    PROTO_ASSERT(a_block < numBlocks(),
            "MBAMRData::getBLock | Error: block %u is out of bounds", a_block);
    return *(m_amrData[a_block]);
}

// Operators
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::averageDown(){}        

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::increment(
        const MBAMRData<T,C,MEM,CTR>& a_data, T a_scale){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::operator*=(T a_value){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::operator+=(T a_value){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
T MBAMRData<T,C,MEM,CTR>::integrate(Array<T,DIM> a_cdx, unsigned int a_c){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
T MBAMRData<T,C,MEM,CTR>::integrate(T a_cdx, unsigned int a_c){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
T MBAMRData<T,C,MEM,CTR>::integrateAbs(Array<T,DIM> a_cdx, unsigned int a_c){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
T MBAMRData<T,C,MEM,CTR>::integrateAbs(T a_cdx, unsigned int a_c){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
T MBAMRData<T,C,MEM,CTR>::absMax(unsigned int a_c){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::copyTo(MBAMRData& a_rhs){}
