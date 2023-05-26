
template<typename T, unsigned int C, MemType MEM, Centering CTR>
LevelBC<T,C,MEM,CTR>::LevelBC(DisjointBoxLayout& a_layout)
{
    define(a_layout);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void LevelBC<T,C,MEM,CTR>::define(DisjointBoxLayout& a_layout)
{
    m_layout = a_layout;
    this->init();
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void LevelBC<T,C,MEM,CTR>::setNumConstants(int a_numData) const
{
    m_constData.clear();
    m_constData.resize(a_numData);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void LevelBC<T,C,MEM,CTR>::setConstant(int a_index, double a_data) const
{
    PROTO_ASSERT(a_index < m_constData.size(),
    "LevelBC::setConstant | Error: stage %i is out of bounds. Call setNumConstants first.");

    m_constData[a_index] = a_data;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
double LevelBC<T,C,MEM,CTR>::getConstant(int a_index) const
{
    PROTO_ASSERT(a_index < m_constData.size(),
    "LevelBC::getConstant | Error: stage %i is out of bounds. Call setNumConstants first.");

    return m_constData[a_index];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void LevelBC<T,C,MEM,CTR>::setNumLevelData(int a_numStages) const
{
    m_fixedData.clear();
    m_fixedData.resize(a_numStages);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void LevelBC<T,C,MEM,CTR>::setLevelData(
        int a_index, const LevelBoxData<T,C,MEM,CTR>& a_data) const
{
    PROTO_ASSERT(a_index < m_fixedData.size(),
    "LevelBC::setFixedData | Error: stage %i is out of bounds. Call setNumFixedData first.");

    m_fixedData[a_index] = &a_data;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const LevelBoxData<T,C,MEM,CTR>& LevelBC<T,C,MEM,CTR>::getLevelData(int a_index) const
{
    PROTO_ASSERT(a_index < m_fixedData.size(),
    "LevelBC::getFixedData | Error: stage %i is out of bounds. Call setNumFixedData first.");

    return *m_fixedData[a_index];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void LevelBC<T,C,MEM,CTR>::setNumRKStages(int a_numStages) const
{
    m_rkStageData.clear();
    m_rkStageData.resize(a_numStages);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void LevelBC<T,C,MEM,CTR>::setRKStageData(
        int a_stage, const LevelBoxData<T,C,MEM,CTR>& a_data) const
{
    PROTO_ASSERT(a_stage < m_rkStageData.size(),
    "LevelBC::setRKStageData | Error: stage %i is out of bounds. Call setNumRKStages first.");

    m_rkStageData[a_stage] = &a_data;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const LevelBoxData<T,C,MEM,CTR>& LevelBC<T,C,MEM,CTR>::getRKStageData(int a_stage) const
{
    PROTO_ASSERT(a_stage < m_rkStageData.size(),
    "LevelBC::setRKStageData | Error: stage %i is out of bounds. Call setNumRKStages first.");

    return *m_rkStageData[a_stage];
}
