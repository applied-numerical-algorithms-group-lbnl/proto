
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void LevelBC<T,C,MEM,CTR>::setNumFixedData(int a_numStages)
{
    m_fixedData.clear();
    m_fixedData.resize(a_numStages);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void LevelBC<T,C,MEM,CTR>::setFixedData(int a_index, const LevelBoxData<T,C,MEM,CTR>& a_data)
{
    PROTO_ASSERT(a_index < m_fixedData.size(),
    "LevelBC::setFixedData | Error: stage %i is out of bounds. Call setNumFixedData first.");

    m_fixedData[a_index] = &a_data;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const LevelBoxData<T,C,MEM,CTR>& LevelBC<T,C,MEM,CTR>::getFixedData(int a_index)
{
    PROTO_ASSERT(a_index < m_fixedData.size(),
    "LevelBC::getFixedData | Error: stage %i is out of bounds. Call setNumFixedData first.");

    return *m_fixedData[a_index];
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void LevelBC<T,C,MEM,CTR>::setNumRKStages(int a_numStages)
{
    m_rkStageData.clear();
    m_rkStageData.resize(a_numStages);
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void LevelBC<T,C,MEM,CTR>::setRKStageData(int a_stage, const LevelBoxData<T,C,MEM,CTR>& a_data)
{
    PROTO_ASSERT(a_stage < m_rkStageData.size(),
    "LevelBC::setRKStageData | Error: stage %i is out of bounds. Call setNumRKStages first.");

    m_rkStageData[a_stage] = &a_data;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
const LevelBoxData<T,C,MEM,CTR>& LevelBC<T,C,MEM,CTR>::getRKStageData(int a_stage)
{
    PROTO_ASSERT(a_stage < m_rkStageData.size(),
    "LevelBC::setRKStageData | Error: stage %i is out of bounds. Call setNumRKStages first.");

    return *m_rkStageData[a_stage];
}
