
template <
    template<typename, template<MemType> typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, template<MemType> typename, MemType, Centering> class BCType,
    MemType MEM>
MBAMR<OPType, MAP, T, BCType, MEM>::MBAMR(
            const MBAMRGrid&            a_grid,
            Point                       a_refRatio)
{
    std::vector<Point> refRatios(a_grid.numBlocks(), a_refRatio);
    define(a_grid, refRatios);
}

template <
    template<typename, template<MemType> typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, template<MemType> typename, MemType, Centering> class BCType,
    MemType MEM>
MBAMR<OPType, MAP, T, BCType, MEM>::MBAMR(
            const MBAMRGrid&            a_grid,
            const std::vector<Point>&   a_refRatios)
{
    define(a_grid, a_refRatios);
}

template <
    template<typename, template<MemType> typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, template<MemType> typename, MemType, Centering> class BCType,
    MemType MEM>
void MBAMR<OPType, MAP, T, BCType, MEM>::define(
            const MBAMRGrid&            a_grid,
            const std::vector<Point>&   a_refRatios)
{
    PROTO_ASSERT(a_refRatios.size() == a_grid.numBlocks(),
            "MBAMR::define | Error: must provide a refinement ratio for each block.");

    m_numLevels = a_grid.numLevels(); 
    m_mg.clear();

    m_residual.clear();
    m_state.clear();
    m_state_0.clear();
    m_force.clear();
    m_local.clear();

    //m_boundRegisters.clear();
    m_constInterp.clear();
    m_BCInterp.clear();
    m_average.clear();
    m_localBlockInterp.clear();

    //m_boundRegisters.resize(m_numLevels);
    m_residual.resize(m_numLevels);
    m_state.resize(m_numLevels-1);
    m_state_0.resize(m_numLevels-1);
    m_force.resize(m_numLevels-1);
    m_local.resize(m_numLevels-1);
    m_localBlockInterp.resize(m_numLevels-1);

    m_map.define(a_grid, OP::ghost());

    for (int bi = 0; bi < a_grid.numBlocks(); bi++)
    {
        m_average.push_back(Stencil<T>::AvgDown(a_refRatios[bi]));
        m_constInterp.push_back(InterpStencil<T>::Constant(a_refRatios[bi]));
        m_BCInterp.push_back(InterpStencil<T>::FiniteVolume(a_refRatios[bi], OP::order()+1));
    }
   
    Point localGhost = OP::ghost() + m_BCInterp[0].ghost();

    for (int li = 0; li < m_numLevels; li++)
    {
        // build multigrid operators
        if (li == 0)
        {
            //FIXME: assumes uniform block size
            Box B0 = a_grid.getLevel(0).getBlock(0).domain().box();
            double L0 = B0.sizes().min();
            int N0 = log(L0)/log(2.0);
            m_mg.push_back(MBMultigrid<OPType, MAP, T, BCType, MEM>(
                        a_grid[li], Point::Ones(2), N0));
        } else {
            //FIXME: assumes no subcycling
            m_mg.push_back(MBMultigrid<OPType, MAP, T, BCType, MEM>(
                        a_grid[li], Point::Ones(), 1));
        }

        // build data holders
        m_residual[li] = std::make_shared<LevelStateData>(a_grid.getLevel(li), Point::Zeros());
        //m_boundRegisters[li].define(a_grid.getLevel(li), 1, Point::Zeros());
        if (li < m_numLevels-1)
        {
            m_state[li] = std::make_shared<LevelStateData>(a_grid.getLevel(li), OP::ghost());
            m_state_0[li] = std::make_shared<LevelStateData>(a_grid.getLevel(li),Point::Zeros());
            m_force[li] = std::make_shared<LevelStateData>(a_grid.getLevel(li), Point::Zeros());

            auto crseFineLayout = a_grid.getLevel(li+1).coarsen(a_refRatios);
            m_local[li] = std::make_shared<LevelStateData>(crseFineLayout, localGhost);
            
            MBLevelMap<MAP,MEM> localMap(crseFineLayout, localGhost);
            m_localBlockInterp[li] = std::make_shared<MBInterpOp>(localMap, OP::order()+1);
        }
    }
}

template <
    template<typename, template<MemType> typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, template<MemType> typename, MemType, Centering> class BCType,
    MemType MEM>
T MBAMR<OPType, MAP, T, BCType, MEM>::solve(
        AMRStateData&   a_state,
        AMRStateData&   a_force,
        int             a_maxIter,
        T               a_tolerance)
{
    PROTO_ASSERT(a_state.numLevels() == m_numLevels,
            "MBAMR::solve | Error: input does not have the right number of levels");
    PROTO_ASSERT(a_force.numLevels() == m_numLevels,
            "MBAMR::solve | Error: input does not have the right number of levels");

    // average down state
    // compute initial residual
    auto& crseState = a_state[0];
    const auto& crseForce = a_force[0];
    const auto& fineForce = a_force[m_numLevels-1];
    LevelStateData crseRes(crseState.layout(), Point::Zeros());
    T r0 = m_mg[0].resnorm(crseRes, crseState, crseForce);
    T r = r0;
    for (int ii = 0; ii < a_maxIter; ii++)
    {
        vCycle(a_state, a_force, fineForce, m_numLevels-1);

        //average down state
        r = m_mg[0].resnorm(crseRes, crseState, crseForce);
#if PR_VERBOSE > 0
        std::cout << "ii = " << ii << " | residual: " << r << std::endl;
#endif
        if (r/r0 < a_tolerance) {
            return r;
        }
    }
    return r;
}

template <
    template<typename, template<MemType> typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, template<MemType> typename, MemType, Centering> class BCType,
    MemType MEM>
void MBAMR<OPType, MAP, T, BCType, MEM>::vCycle(
        AMRStateData&         a_state,
        const AMRStateData&   a_force,
        const LevelStateData& a_levelForce,
        int                   a_level)
{
    auto& state = a_state[a_level];
    const auto& force = a_force[a_level]; 
    if (a_level == 0)
    {
        m_mg[a_level].solve(state, a_levelForce, 1, 1e-12);
    } else {
        auto& crseState = a_state[a_level-1];
        const auto& crseForce = a_force[a_level-1]; 
        
        interpBounds(crseState, state, a_level);
        m_mg[a_level].solve(state, force, 1, 1e-12);
        averageDown(crseState, state, a_level);
        crseState.copyTo(*(m_state_0[a_level-1]));
        coarseResidual(*m_residual[a_level-1],
                a_state, crseForce, force, a_level);
        op(a_level-1)(*m_force[a_level-1], crseState); //LC
        m_force[a_level-1]->increment(*m_residual[a_level-1]);
        vCycle(a_state, a_force, *m_force[a_level-1], a_level-1);
        fineCorrect(state, crseState, *m_state_0[a_level-1], a_level);
        interpBounds(crseState, state, a_level);
        m_mg[a_level].solve(state, force, 1, 1e-12);
    }
}

template <
    template<typename, template<MemType> typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, template<MemType> typename, MemType, Centering> class BCType,
    MemType MEM>
void MBAMR<OPType, MAP, T, BCType, MEM>::interpBounds(
        LevelStateData&     a_crseData,
        LevelStateData&     a_fineData,
        int                 a_fineLevel)
{
    HDF5Handler h5;

    auto& crseTmp = *m_local[a_fineLevel-1];
    a_crseData.copyTo(crseTmp);
    crseTmp.exchange();
    m_localBlockInterp[a_fineLevel-1]->apply(crseTmp, crseTmp);

    auto& layout = a_fineData.layout();
    for (auto iter : layout)
    {
        int block = layout.block(iter);
        Box box = layout[iter];
        auto& fine_i = a_fineData[iter];
        auto& crse_i = crseTmp[iter];
        StateData fine_0(box);
        fine_i.copyTo(fine_0);
        fine_i |= m_BCInterp[block](crse_i);
        fine_0.copyTo(fine_i);
    }
}

template <
    template<typename, template<MemType> typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, template<MemType> typename, MemType, Centering> class BCType,
    MemType MEM>
void MBAMR<OPType, MAP, T, BCType, MEM>::averageDown(
        LevelStateData&     a_crseData,
        LevelStateData&     a_fineData,
        int                 a_fineLevel)
{
    auto& crseTmp = *m_local[a_fineLevel-1];
    auto& layout = a_fineData.layout();
    for (auto iter : layout)
    {
        int block = layout.block(iter);
        auto& fine_i = a_fineData[iter];
        auto& crse_i = crseTmp[iter];
        crse_i |= m_average[block](fine_i);
    }
    crseTmp.copyTo(a_crseData);
}

template <
    template<typename, template<MemType> typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, template<MemType> typename, MemType, Centering> class BCType,
    MemType MEM>
void MBAMR<OPType, MAP, T, BCType, MEM>::coarseResidual(
        LevelStateData&     a_crseResidual,
        AMRStateData&       a_state,
        const LevelStateData&     a_crseForce,
        const LevelStateData&     a_fineForce,
        int                 a_fineLevel)
{
    PROTO_ASSERT(a_fineLevel > 0,
            "MBAMR::coarseResidual | Error: Fine Level must be greater than 0");
    PROTO_ASSERT(a_fineLevel < a_state.numLevels(),
            "MBAMR::coarseResidual | Error: Fine Level index is out of bounds");
    auto& crseState = a_state[a_fineLevel-1];
    auto& fineState = a_state[a_fineLevel];
    auto& crseOp = op(a_fineLevel-1);
    auto& fineOp = op(a_fineLevel);
    
    //interpolate onto file level
    //TODO: potentially redundant
    interpBounds(crseState, fineState, a_fineLevel);
    fineState.exchange();
    m_localBlockInterp[a_fineLevel]->apply(fineState, fineState);
    
    // compute coarse residual
    crseOp(a_crseResidual, crseState, -1.0); //-LC
    a_crseResidual.increment(a_crseForce);   //GC - LC
    
    for (auto iter : fineState.layout())
    {
        auto block = fineState.layout().block(iter);
        auto& fineState_i = fineState[iter];
        auto& fineForce_i = a_fineForce[iter];
        auto& crseTmp_i = (*m_local[a_fineLevel-1])[iter];
        
        // compute fine level residual
        BoxData<T, OP::numState(), MEM> fineTmp(fineState.layout()[iter]);
        fineOp[iter](fineTmp, fineState_i, -1.0);
        fineTmp += fineForce_i;
        // average onto crse temp
        crseTmp_i |= m_average[block](fineTmp);
    }
    // copy crse temp into crseResidual, overwriting the refined region
    m_local[a_fineLevel-1]->copyTo(a_crseResidual);
}

template <
    template<typename, template<MemType> typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, template<MemType> typename, MemType, Centering> class BCType, MemType MEM>
void MBAMR<OPType, MAP, T, BCType, MEM>::fineCorrect(
        LevelStateData&         a_fineState,
        const LevelStateData&   a_crseState,
        LevelStateData&         a_crseState_0,
        int                     a_fineLevel)
{
    for (auto iter : a_crseState.layout())
    {
        auto& crse_0 = a_crseState_0[iter];
        const auto& crse   = a_crseState[iter];
        crse_0 -= crse;
        crse_0 *= -1;
    }
    a_crseState_0.copyTo(*m_local[a_fineLevel]);
    for (auto iter : a_fineState.layout())
    {
        auto block = a_fineState.layout().block(iter);
        auto& delta = (*m_local[a_fineLevel])[iter];
        auto& state = a_fineState[iter];

        state += m_constInterp[block](delta);
    }
}

template <
    template<typename, template<MemType> typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, template<MemType> typename, MemType, Centering> class BCType,
    MemType MEM>
MBLevelMap<MAP>& MBAMR<OPType, MAP, T, BCType, MEM>::map(int a_level)
{
    PROTO_ASSERT(a_level < m_numLevels,
            "MBAMR::map | Error: Level index %i is out of bounds", a_level);
    return m_map[a_level];
}

template <
    template<typename, template<MemType> typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, template<MemType> typename, MemType, Centering> class BCType,
    MemType MEM>
MBInterpOp& MBAMR<OPType, MAP, T, BCType, MEM>::interpOp(int a_level)
{
    PROTO_ASSERT(a_level < m_numLevels,
            "MBAMR::interpOp | Error: Level index %i is out of bounds", a_level);
    return m_mg[a_level].interpOp();
}

template <
    template<typename, template<MemType> typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, template<MemType> typename, MemType, Centering> class BCType,
    MemType MEM>
MBLevelOp<OPType, MAP, double, BCType, MEM>& MBAMR<OPType, MAP, T, BCType, MEM>::op(int a_level)
{
    PROTO_ASSERT(a_level < m_numLevels,
            "MBAMR::op | Error: Level index %i is out of bounds", a_level);
    return m_mg[a_level].op();
}
