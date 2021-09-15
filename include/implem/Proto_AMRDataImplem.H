
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

template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
void
AMRData<T, C, MEM, CTR>::computeTags(
    LevelTagData& a_tags,
    unsigned int a_level,
    Point a_bufferSize,
    T a_threshold) const
{
    PROTO_ASSERT(a_level != maxLevels()-1,
        "AMRData::computeTags | Error: Cannot generate tags on the finest level");
    PROTO_ASSERT(a_level < maxLevels()-1,
        "AMRData::computeTags | Error: Level %u is out of bounds", a_level);
    PROTO_ASSERT(a_level < m_data.size(),
        "AMRData::computeTags | Error: No data defined on level %u.", a_level);

    a_tags.define(m_grid[a_level], a_bufferSize);

    std::vector<Stencil<T>> D;
    for (int dir = 0; dir < DIM; dir++)
    {
        Stencil<T> Si = ((T)1)*Shift::Basis(dir, 1) - ((T)1)*Shift::Basis(dir, -1);
        D.push_back(Si);
    }

    for (auto iter = m_grid[a_level].begin(); iter.ok(); ++iter)
    {
        const auto& data = operator[](a_level)[*iter];
        auto& tags = a_tags[*iter];
        BoxData<T, DIM> diffs(iter.box()); 
        for (int dir = 0; dir < DIM; dir++)
        {
            BoxData<T, 1> diff = slice(diffs, dir);
            diff |= D[dir](data);
        }
        forallInPlace(
        [] PROTO_LAMBDA (Var<char, 1>& v_tags, Var<T, DIM>& v_diffs, T v_threshold)
        {
            T diffSq = 0;
            for (int dir = 0; dir < DIM; dir++)
            {
                diffSq += v_diffs(dir)*v_diffs(dir);
            }
            T diff = sqrt(diffSq);
            if (diff > v_threshold)
            {
                v_tags(0) = 1;
            } else {
                v_tags(0) = 0;
            }
        }, tags, diffs, a_threshold);
    }
    
    AMRGrid::buffer(a_tags, a_bufferSize);
}

