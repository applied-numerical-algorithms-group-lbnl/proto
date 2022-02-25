
PROTO_KERNEL_START
void
f_RKTaylorCoefsF(Var<double,RKORDER*NUMCOMPS>& a_UTaylor,
              Var<double,NUMCOMPS>& a_kstage,
              Var<double,NUMCOMPS>& a_U,
              int a_stage
              )
{
  for (int n = 0; n < NUMCOMPS;n++)
    {
      if (a_stage == 0) 
        {
          a_UTaylor(n*RKORDER) = a_U(n);
          a_UTaylor(1 + n*RKORDER) = a_kstage(n);
          a_UTaylor(2 + n*RKORDER) = -(3.0/2)*a_kstage(n);
          a_UTaylor(3 + n*RKORDER) = (2.0/3)*a_kstage(n);
        }
      if (a_stage == 1)
        {
          a_UTaylor(2 + n*RKORDER) += a_kstage(n);
          a_UTaylor(3 + n*RKORDER) += (-2.0/3)*a_kstage(n);
        }
      if (a_stage == 2)
        {
          a_UTaylor(2 + n*RKORDER) += a_kstage(n);
          a_UTaylor(3 + n*RKORDER) += (-2.0/3)*a_kstage(n);
          
        }
      if (a_stage == 3)
        {
          a_UTaylor(2 + n*RKORDER) += (-1.0/2)*a_kstage(n);
          a_UTaylor(3 + n*RKORDER) += (2.0/3)*a_kstage(n);
        }
    }
}
PROTO_KERNEL_END(f_RKTaylorCoefsF, f_RKTaylorCoefs)

PROTO_KERNEL_START
void
f_RKTimeInterpF(Var<double,NUMCOMPS>& a_U,
              Var<double,RKORDER*NUMCOMPS>& a_UTaylor,
              double& a_chi
              )
{
  for (int n = 0; n < NUMCOMPS;n++)
    {
      a_U(n) = 0.;
      for (int term = 0; term < RKORDER; term++)
        {
          a_U(n) = a_U(n)*a_chi + a_UTaylor(RKORDER-term-1 + n*RKORDER);
        }
    }
}
PROTO_KERNEL_END(f_RKTimeInterpF,f_RKTimeInterp)

PROTO_KERNEL_START
void f_stageUpdateF(
                 Var<double,NUMCOMPS>& a_UStage,
                 Var<double,NUMCOMPS>& a_rhsTot,
                 const Var<double,NUMCOMPS>& a_U0,
                 const Var<double,NUMCOMPS>& a_kStage,
                 const double& a_dtStageNext,
                 const double& a_stageWeight)
{
  for (int comp = 0; comp < NUMCOMPS; comp++)
    {
      a_rhsTot(comp) += a_stageWeight*a_kStage(comp);
      a_UStage(comp) = a_U0(comp) + a_dtStageNext*a_kStage(comp); 
    }
}
PROTO_KERNEL_END(f_stageUpdateF,f_stageUpdate)

template <class OpType, typename T,unsigned int C,MemType MEM>
AMRRK4<OpType,T,C,MEM>::AMRRK4(
        AMRData<T, C, MEM>& a_data,
        T                   a_dt)
{
    define(a_data, a_dt, 0);
}

template <class OpType, typename T,unsigned int C,MemType MEM>
void AMRRK4<OpType,T,C,MEM>::define(
        AMRData<T, C, MEM>& a_data,
        std::array<T, DIM>& a_dx,
        int                 a_level)
{
    int numLevels = a_data.grid().numLevels();
    
    m_data = &a_data;
    
    m_ops.resize(numLevels);
    m_UStage.resize(numLevels);
    m_RHSTotal.resize(numLevels);
    m_timeRefRatios.resize(numLevels - 1);
    m_spaceRefRatios.resize(numLevels - 1);
    m_cfInterp.resize(numLevels - 1);
    m_UTaylor.resize(numLevels - 1);
    m_UTaylorC.resize(numLevels - 1);
    m_fluxRegisters.resize(numLevels - 1);

    std::array<T, DIM> dx = a_dx;
    for (int lvl = a_level; lvl < numLevels; lvl++)
    {
        auto& layout = a_data.grid()[lvl].layout();
        
        m_ops[lvl].define(dx);
        m_UStage[lvl].define(layout, OpType::ghost());
        m_RHSTotal[lvl].define(layout, Point::Zeros());
        if (lvl < numLevels-1)
        {
            auto& fineLayout = a_data.grid()[lvl+1].layout();
            Point spaceRefRatio = a_data.grid().refRatio(lvl);
            auto cfLayout = fineLayout.coarsen(spaceRefRatio);

            m_spaceRefRatios[lvl] = spaceRefRatio;
            m_timeRefRatios[lvl] = spaceRefRatio.max();
            m_cfInterp[lvl].define(OpType::order()+1, spaceRefRatio);
            m_UTaylor[lvl].define(layout, Point::Zeros());
            m_fluxRegisters[lvl].define(layout, fineLayout, spaceRefRatio);
            
            // Point::Ones() are here to make the division a ceil instead of a floor
            Point interpGhost = OpType::ghost() / spaceRefRatio + Point::Ones();
            interpGhost += m_cfInterp[lvl].ghost();
            
            m_UTayolorC[lvl].define(cfLayout, interpGhost);
            for (int dir = 0; dir < DIM; dir++)
            {
                dx[dir] /= spaceRefRatio[dir];
            }
        }
    }
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














