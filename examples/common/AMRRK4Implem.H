
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
    m_op.define(a_dx);

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
        m_UTaylor.define(layout, Point::Zeros());
    }
    if (m_level > 0)
    {
        m_spaceRefRatioCrse = m_data->grid().refRatio(m_level-1);
        m_timeRefRatioCrse = m_spaceRefRatioCrse.max(); //CFL restriction
        m_cfInterp.define(OpType::order() + 1, m_spaceRefRatioCrse);

        DisjointBoxLayout crseLayout = layout.coarsen(m_spaceRefRatioCrse);
        Point interpGhost = OpType::ghost()/m_spaceRefRatioCrse + m_cfInterp.ghost();
        m_UTaylorCrse.define(crseLayout, interpGhost);
    }
}

template<class OpType, typename T,unsigned int C,MemType MEM>
void AMRRK4::advance(const T& a_dt, const T& a_chi)
{
    auto& layout = m_data->grid()[m_level];
    LevelBoxData<T, C, MEM> UStage(layout, OpType::ghost());
    (*m_data)[m_level].copyTo(UStage);
    LevelBoxData<T, C, MEM> RHSTotal(layout, Point::Zeros());
    RHSTotal.setToZero();
    auto& U0 = (*m_data)[m_level];
    Stencil<T> I = 1.0*Shift::Zeros();
    for (int stage = 0; stage < RKSTAGES; stage++)
    {
        // Interpolate Ghost Cells
        if (m_level > 0)
        {
            T chiStage  = a_chi + m_stageTimes[stage] / m_timeRefRatioCrse;
            T timeStage = m_time + a_dt*m_stageTimes[stage];
            for (auto iter = layout.begin(); iter.ok(); ++iter)
            {
                // Only interpolate non-interior patches
                if (layout.onLevelBoundary(iter.point()))
                {
                    // store the interior data
                    BoxData<T, C, MEM> UStageTemp(layout[*iter]);
                    auto& UStage_i = UStage[*iter];
                    UStageTemp |= I(UStage_i); //stand-in for copyTo

                    // interpolate in time
                    auto UTimeInterp = forall<T, C>
                        (RKTimeInterp, m_UTayolorCrse[*iter], chiStage);
                    
                    // interpolate boundary conditions
                    m_cfInterp.apply(UStage_i, UTimeInterp);

                    // replace interior data
                    UStage_i |= I(UStageTemp); //stand-in for copyTo
                }
            }
        }
        UStage.exchange(); //fixes up ghost cells after interpolation
        
        // Do the stage update
        T dtStageNext = a_dt*m_stageTimes[stage+1];
        T dtStage = a_dt*m_stageTimes[stage];
        for (auto iter = layout.begin(); iter.ok(); ++iter)
        {
            auto& UStage_i = UStage[*iter];
            auto& RHS_i = RHSTotal[*iter];
            auto& UTaylor_i = m_UTaylor[*iter];
            auto& U0_i = U0[*iter];

            std::array<BoxData<T, C, MEM>, DIM> fluxes;
            for (int dir = 0; dir < DIM; dir++)
            {
                fluxes[dir].define(layout[*iter].grow(dir, Side::Hi, 1), 0);
            }
            BoxData<T, C, MEM> kStage(layout[*iter], 0);
            m_op(kStage, fluxes, UStage_i);
            forallInPlace
                (stageUpdate, UStage_i, RHS_i, U0_i, kStage, dtStageNext, m_stageWgt[stage]); 
            if (m_level < (m_data->grid().numLevels() - 1))
            {
                kStage *= a_dt;
                forallInPlace(RKTaylorCoefs, UTaylor_i, kStage, U0_i, stage);
            }
            if (m_level > 0)
            {
            }
        }

    }




}
