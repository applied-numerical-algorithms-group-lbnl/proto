
template<class OpType, typename T,unsigned int C,MemType MEM>
AMRRK4::AMRRK4(
        AMRData<T,C,MEM>&   a_data,
        std::array<T, DIM>  a_dx,
        int      a_timeRefRatio)
{
    define(a_data, a_dx, a_timeRefRatio);
}

template<class OpType, typename T,unsigned int C,MemType MEM>
void AMRRK4::define(
        AMRData<T,C,MEM>&   a_data,
        std::array<T, DIM>  a_dx,
        int      a_timeRefRatio)
{
    define(a_data, a_dx, a_timeRefRatio, 0);
}

template<class OpType, typename T,unsigned int C,MemType MEM>
void AMRRK4::define(
        AMRData<T,C,MEM>&   a_data,
        std::array<T, DIM>  a_dx,
        int                 a_timeRefRatio,
        unsigned int        a_level)
{
    m_data = &a_data;
    m_level = a_level;
    m_time = 0.0;
    m_dx = a_dx;

    auto& layout = m_data->grid()[m_level];
    // if there are finer levels than m_level
    if (m_level < m_data->numLevels()-1)
    {
        m_spaceRefRatio = m_data->grid().refRatio(m_level);
        m_timeRefRatio = m_spaceRefRatio.max(); //CFL restriction
        std::array<T, DIM> dxFine;
        for (int dir = 0; dir < DIM; dir++)
        {
            dxFine[dir] = a_dx[dir] / m_spaceRefRatio[dir];
        }
        m_fine = std::make_shared<AMRRK4<OpType, T, C, MEM>>(a_data, dxFine, m_timeRefRatio);
        m_UTaylor.resize(RKORDER);
        for (int ii = 0; ii < RKORDER; ii++)
        {
            m_UTaylor[ii] = std::make_shared<LevelBoxData<T, C, MEM>>(layout, Point::Zeros()); 
        }
    }
    if (m_level > 0)
    {
        m_spaceRefRatioCrse = m_data->grid().refRatio(m_level-1);
        m_timeRefRatioCrse = m_spaceRefRatioCrse.max(); //CFL restriction
        m_cfInterp.define(OpType::order() + 1, m_spaceRefRatioCrse);

        DisjointBoxLayout crseLayout = layout.coarsen(m_spaceRefRatioCrse);
        m_UTaylorCrse.resize(RKORDER);
        Point interpGhost = OpType::ghost()/m_spaceRefRatioCrse + m_cfInterp.ghost();
        for (int ii = 0; ii < RKORDER; ii++)
        {
            m_UTaylorCrse[ii] = std::make_shared<LevelBoxData<T, C, MEM>>(crseLayout, interpGhost);
        }
    }
}

template<class OpType, typename T,unsigned int C,MemType MEM>
void AMRRK4::advance(const T& a_dt)
{
    auto& layout = m_data->grid()[m_level];
    LevelBoxData<T, C, MEM> UStage(layout, OpType::ghost());
    (*m_data)[m_level].copyTo(UStage);
    LevelBoxData<T, C, MEM> UTotal(layout, Point::Zeros());
    UTotal.setToZero();
    
    for (int stage = 0; stage < RKSTAGES; stage++)
    {
        // Interpolate Ghost Cells
    }




}
