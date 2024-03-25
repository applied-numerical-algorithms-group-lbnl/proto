
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
void LevelBC<T,C,MEM,CTR>::setTime(T a_time) const
{
    m_time = a_time;
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
