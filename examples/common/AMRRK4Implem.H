template <class OpType, typename T,unsigned int C,MemType MEM>
void AMRRK4<OpType,T,C,MEM>::define(
        AMRData<T, C, MEM>& a_data,
        T                   a_dt,
        int                 a_level)
{
    int numLevels = a_data.grid().numLevels();
    
    m_data = &a_data;
    
    m_timeRefRatios.resize(numLevels - 1);
    m_spaceRefRatios.resize(numLevels - 1);
    m_ops.resize(numLevels);
    m_cfInterp.resize(numLevels - 1);
    m_UStage.resize(numLevels);
    m_RHSTotal.resize(numLevels);
    m_UTaylor.resize(numLevels);
    m_UTaylorC.resize(numLevels - 1);
    m_fluxRegisters.resize(numLevels - 1);
}

template <class OpType, typename T,unsigned int C,MemType MEM>
void AMRRK4<OpType,T,C,MEM>::advance(int a_level, T a_dt, T a_chi)
{
    auto& layout    = m_data->grid[a_level];
    auto& U0 = (*m_data)[a_level];
    auto& UStage = m_UStage[a_level];
    auto& RHSTotal = m_RHSTotal[a_level];
    
    T timeRefRatio = m_timeRefRatios[a_level];

    U0.copyTo(UStage);
    RHSTotal.setToZero();

    for (int stage = 0; stage < RKSTAGES; stage++)
    {
        if (a_level > 0)
        {
            interpBoundaries(a_level, stage, a_chi);
        }
        stageUpdate(a_level, stage, a_dt);
    } // end stage loop

    // Update U0
    for (auto iter = layout.begin(); iter.ok(); ++iter)
    {
        RHSTotal[*iter] *= a_dt;
        U0[*iter] += RHStotal[*iter];
    }

    // Recursive Advance
    if (a_level < m_data->numLevels() - 1)
    {
        m_UTaylor[a_level].copyTo(m_UTaylorC[a_level]);
        for (int nstep = 0; nstep < timeRefRatio; nstep++)
        {
            T chi = nstep*1.0/timeRefRatio;
            T dtFine = a_dt/timeRefRatio;
            advance(a_level+1, dtFine, chi);
        }
        Proto::averageDown((*m_data)[a_level], (*m_data)[a_level+1], m_spaceRefRatios[a_level]);
        // 1/dx factor is done at the increment stage to account for anisotropy
        m_fluxRegisters[a_level].reflux(U0, 1.0); 
    }
    
    // Do Regridding
    if (a_level == m_regridBaseLevel)
    {
        regrid(a_level);
    }
}

template <class OpType, typename T,unsigned int C,MemType MEM>
void AMRRK4<OpType,T,C,MEM>::interpBoundaries(int a_level, int a_stage, T a_chi)
{
    PROTO_ASSERT(a_level > 0,
        "AMRRK4::interpBoundaries | Error: Interpolation is only done for level > 0.");
    
    auto& layout    = m_data->grid[a_level];
    auto& UStage    = m_UStage[a_level];
    auto& UTaylorC  = m_UTaylorC[a_level-1];
    auto& cfInterp  = m_cfInterp[a_level-1];

    T chi = a_chi + m_stageTimes[a_stage]/m_timeRefRatios[a_level-1];
    Stencil<T> I = 1.0*Shift::Zeros();

    for (auto iter = layout.begin(); iter.ok(); ++iter)
    {
        if (layout.onLevelBoundary(iter.point()))
        {
            // get data
            auto& UStage_i   = UStage[*iter];
            auto& UTaylorC_i = UTaylorC[*iter];
            BoxData<T, C, MEM> UStageTemp(layout[*iter]);
            
            // save copy of interior
            UStageTemp |= I(Ustage_i);

            // interpolate in time
            auto UTimeInterp = forall<T, C>(f_RKTimeInterp, UTaylorC_i, chi);
        
            // interpolate boundary conditions
            cfInterp.apply(UStage_i, UTimeInterp);
        
            // copy back the interior
            UStage_i |= I(UStageTemp);
        }
    }
    UStage.exchange();
}

template <class OpType, typename T,unsigned int C,MemType MEM>
void AMRRK4<OpType,T,C,MEM>::stageUpdate(int a_level, int a_stage, T a_dt)
{
    auto& layout    = m_data->grid[a_level];
    auto& UStage    = m_UStage[a_level];
    auto& U0        = (*m_data)[a_level];
    auto& RHSTotal  = m_RHSTotal[a_level];
    auto& UTaylor   = m_UTaylor[a_level];
    auto& op        = m_ops[a_level];
    
    T dtStageNext = a_dt*m_stageTimes[a_stage+1];
    T stageWeight = m_stageWeights[a_stage];
    T dtWeight    = a_dt*m_stageWeights[a_stage];

    for (auto iter = layout.begin(); iter.ok(); ++iter)
    {
        // get data
        auto& UStage_i   = UStage[*iter];
        auto& U0_i       = U0[*iter];
        auto& RHSTotal_i = RHSTotal[*iter];
        auto& UTaylor_i  = UTaylor[*iter];

        std::arry<BoxData<T, C, MEM>, DIM> fluxes;
        for (int dir = 0; dir < DIM; dir++)
        {
            fluxes[dir].define(layout[*iter].grow(dir, Side::Hi, 1), 0.0);
        }
        BoxData<T,C,MEM> kStage(layout[*iter], 0.0);
        
        // apply spatial operator to compute k
        op(kStage, fluxes, UStage_i);
        
        // do stage update
        forallInPlace(f_stageUpdate, UStage_i, RHSTotal_i, U0_i, kStage, dtStageNext, stageWeight);
    
        if (a_level < m_data->numLevels() - 1)
        {
            // store stage information
            kStage *= a_dt;
            forallInPlace(f_RKTaylorCoefs, UTaylor_i, kStage, U0_i, a_stage);

            // increment flux register
            for (int dir = 0; dir < DIM; dir++)
            {
                T cdx = op.dx()[dir];
                m_fluxRegisters[a_level].incrementCoarse(fluxes[dir], *iter, dtWeight/cdx, dir);
            }
        }
        if (a_level > 0)
        {
            for (int dir = 0; dir < DIM; dir++)
            {
                T cdx = m_ops[a_level-1].dx()[dir];
                m_fluxRegisters[a_level-1].incrementFine(fluxes[dir], *iter, dtWeight/cdx, dir);
            }
        }
    } // end patch loop
}

template <class OpType, typename T,unsigned int C,MemType MEM>
void AMRRK4<OpType,T,C,MEM>::regrid(int a_baseLevel)
{
    AMRGrid newGrids = m_data->grid();
    for (int level = m_data->numLevels()-2; level >= a_baseLevel; level--)
    {
        auto& layout = newGrids[level];

        // Compute tags on this level
        LevelTagData tags(layout, Point::Ones(m_bufferSize));
        tags.setToZero();
        for (auto iter = layout.begin(); iter.ok(); ++iter)
        {
            auto& tags_i = tags[*iter];
            auto& U0_i   = (*m_data)[level][*iter];
            m_ops[level].generageTags(tags_i, U0_i);
        }

        // If there are 2+ finer levels
        if (level+2 < m_data->numLevels())
        {
            newGrids.addFinerGrids(tags, level);
        }
        AMRGrid::buffer(tags, Point::Ones(m_bufferSize));
        newGrids.regrid(tags, level);
    }
    for (int level = max(2, a_baseLevel); level < newGrids.numLevels(); level++)
    {
        newGrids.enforceNesting2(level);
    }
    m_data->regrid(newGrids, a_baseLevel, OpType::order()+1);
    
    define(*m_data, m_ops[a_baseLevel].dx(), a_baseLevel);
}














