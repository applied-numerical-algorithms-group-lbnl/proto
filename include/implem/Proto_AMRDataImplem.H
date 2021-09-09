
template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void AMRData<T, C, MEM, CTR>::define(AMRGrid& a_grid, Point a_ghost)
{
    m_ghost = a_ghost;
    m_grid = a_grid;
    m_data.clear();
    for (int ii = 0; ii < a_grid.numLevels(); ii++)
    {
        auto level = std::make_shared<LevelBoxData<T, C, MEM, CTR>>(a_grid[ii], a_ghost);
        m_data.push_back(level);
    }
    m_defined = true;
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
LevelBoxData<T, C, MEM, CTR>&
AMRData<T, C, MEM, CTR>::operator[](unsigned int a_level)
{
    PROTO_ASSERT(a_level < m_data.size(),
        "AMRData::operator[] | Error: level %u is out of bounds.", a_level);
    return *m_data[a_level];
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
const LevelBoxData<T, C, MEM, CTR>&
AMRData<T, C, MEM, CTR>::operator[](unsigned int a_level) const
{
    PROTO_ASSERT(a_level < m_data.size(),
        "AMRData::operator[] | Error: level %u is out of bounds.", a_level);
    return *m_data[a_level];
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void
AMRData<T, C, MEM, CTR>::setToZero()
{
    for (int ii = 0; ii < m_data.size(); ii++)
    {
        m_data[ii]->setToZero();
    }
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
template<typename Func, typename... Srcs>
void
AMRData<T, C, MEM, CTR>::initialize(double a_dx0, Func& a_func, Srcs... a_srcs)
{
    double dx = a_dx0;
    for (int ii = 0; ii < m_data.size(); ii++)
    {
        m_data[ii]->initialize(a_func, dx, a_srcs...);
        dx /= PR_AMR_REFRATIO;
    }
}

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void
AMRData<T, C, MEM, CTR>::regrid(unsigned int a_level)
{
}
