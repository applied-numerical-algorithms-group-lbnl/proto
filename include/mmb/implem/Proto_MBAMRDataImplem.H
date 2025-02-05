template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBAMRData<T,C,MEM,CTR>::MBAMRData(){}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBAMRData<T,C,MEM,CTR>::MBAMRData(
        const MBAMRLayout& a_grid,
        Array<Point, DIM+1> a_ghost,
        Point a_boundGhost)
{
    define(a_grid, a_ghost, a_boundGhost);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBAMRData<T,C,MEM,CTR>::MBAMRData(
        const MBAMRLayout& a_grid,
        Point a_ghost,
        Point a_boundGhost)
{
    Array<Point, DIM+1> ghost;
    ghost.fill(a_ghost + Point::Ones(2)); //assumes 4th order interpolation
    ghost[0] = a_ghost;
    define(a_grid, ghost, a_boundGhost);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBAMRData<T,C,MEM,CTR>::MBAMRData(
        const MBAMRLayout& a_grid,
        std::vector<std::shared_ptr<MBLevelBoxData<T,C,MEM,CTR>>> a_levelData)
{
    define(a_grid, a_levelData);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::define(
        const MBAMRLayout& a_grid,
        Array<Point, DIM+1> a_ghost,
        Point a_boundGhost)
{
    m_ghost = a_ghost;
    m_grid =  a_grid;
    m_levelData.clear();
    for (int li = 0; li < a_grid.numLevels(); li++)
    {
        auto& layout = a_grid.getLevel(li);
        auto levelData = std::make_shared<MBLevelBoxData<T,C,MEM,CTR>>(
                layout, a_ghost, a_boundGhost);
        m_levelData.push_back(levelData);
    }
    defineAMRData();
    PROTO_ASSERT(validate(),
        "MBAMRData::define | Error: Validation failed");
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
bool MBAMRData<T,C,MEM,CTR>::validate() const
{
    bool valid = true;
    valid &= (numLevels() == m_levelData.size()); 
    valid &= (numBlocks() == m_amrData.size());
    if (!valid) { std::cout << "failed on initial sizing: " << std::endl;} 
    for (BlockIndex bi = 0; bi < numBlocks(); bi++)
    {
        for (int lvl = 0; lvl < numLevels(); lvl++)
        {
            auto p1 = &(*m_amrData[bi])[lvl];
            auto p2 = &m_levelData[lvl]->getBlock(bi);
            valid &= (p1 == p2);
            if (!valid) { std::cout << "failed on lvl: " << lvl << " | block: " << bi << std::endl;}
        }
    }
    return valid; 
}


template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::define(
        const MBAMRLayout& a_grid,
        std::vector<std::shared_ptr<MBLevelBoxData<T,C,MEM,CTR>>> a_levelData)
{
    PROTO_ASSERT(a_grid.numLevels() > 0,
            "MBAMRData::define (alias) | Error: Grid has zero levels");
    PROTO_ASSERT(a_levelData.size() == a_grid.numLevels(),
            "MBAMRData::define (alias) | Error: incorrect number of levels");
    m_ghost = a_levelData[0]->ghost();
    m_grid = a_grid;
    m_levelData = a_levelData;
    defineAMRData();
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::defineAMRData()
{
    std::vector<std::vector<std::shared_ptr<LevelBoxData<T,C,MEM,CTR>>>> amrData;
    amrData.resize(numBlocks());
    for (int li = 0; li < numLevels(); li++)
    {
        for (int bi = 0; bi < numBlocks(); bi++)
        {
            amrData[bi].push_back(m_levelData[li]->m_data[bi]);
        }
    }
    m_amrData.clear();
    for (int bi = 0; bi < numBlocks(); bi++)
    {
        m_amrData.push_back(
                std::make_shared<AMRData<T,C,MEM,CTR>>(amrData[bi]));
    }
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
            "MBAMRData::level | Error: level %u is out of bounds", a_level);
    return *(m_levelData[a_level]);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const MBLevelBoxData<T,C,MEM,CTR>&
MBAMRData<T,C,MEM,CTR>::getLevel(unsigned int a_level) const
{
    PROTO_ASSERT(a_level < numLevels(),
            "MBAMRData::level | Error: level %u is out of bounds", a_level);
    return *(m_levelData[a_level]);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
AMRData<T,C,MEM,CTR>&
MBAMRData<T,C,MEM,CTR>::getBlock(unsigned int a_block)
{
    PROTO_ASSERT(a_block < numBlocks(),
            "MBAMRData::block | Error: block %u is out of bounds", a_block);
    return *(m_amrData[a_block]);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const AMRData<T,C,MEM,CTR>&
MBAMRData<T,C,MEM,CTR>::getBlock(unsigned int a_block) const
{
    PROTO_ASSERT(a_block < numBlocks(),
            "MBAMRData::block | Error: block %u is out of bounds", a_block);
    return *(m_amrData[a_block]);
}

// Operators
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::averageDown(){}        

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::interpolate(unsigned int a_level)
{
    if (a_level == 0) { return; } //nothing to interpolate from
    for (BlockIndex bi = 0; bi < grid().numBlocks(); bi++)
    {
        auto refRatio = m_amrData[bi]->grid().refRatio(a_level-1);
        auto I = InterpStencil<T>::Constant(refRatio);
        auto& fineLevel = (*m_amrData[bi])[a_level];
        auto& crseLevel = (*m_amrData[bi])[a_level-1];

        DisjointBoxLayout cfLayout = fineLevel.layout().coarsen(refRatio);
        LevelBoxData<T,C,MEM,CTR> tmp(cfLayout, ghost()[0]);
        crseLevel.copyTo(tmp);
        tmp.exchange();
        for (auto iter : tmp.layout())
        {
            auto& crse_i = tmp[iter];
            auto& fine_i = fineLevel[iter];
            pr_out() << "-----------------------------------------" << std::endl; 
            pr_out() << "stencil span: " << I.span() << " | domain: " << I.domain(fine_i.box()) << " | range: " << I.range(crse_i.box()) << std::endl;
            pr_out() << "crse (input): " << std::endl;
            crse_i.printData();
            
            fine_i |= I(crse_i);
            pr_out() << "fine (output): " << std::endl;
            fine_i.printData();
        }
    }
}        

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void MBAMRData<T,C,MEM,CTR>::interpolate()
{
    for (unsigned int li = 1; li < numLevels(); li++)
    {
        interpolate(li);
    }
}

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
