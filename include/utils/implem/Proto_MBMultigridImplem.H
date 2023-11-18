template <
    template<typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
MBMultigrid<OPType, MAP, T, BCType, MEM>::MBMultigrid(
        const MBDisjointBoxLayout&  a_fineLayout,
        Point                       a_refRatio,
        int                         a_numLevels)
{
    std::vector<Point> refRatios;
    for (int bi = 0; bi < a_fineLayout.numBlocks(); bi++)
    {
        refRatios.push_back(a_refRatio);
    }
    define(a_fineLayout, refRatios, a_numLevels);
}

template <
    template<typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
MBMultigrid<OPType, MAP, T, BCType, MEM>::MBMultigrid(
        const MBDisjointBoxLayout&  a_fineLayout,
        const std::vector<Point>&   a_refRatios,
        int                         a_numLevels)
{
    define(a_fineLayout, a_refRatios, a_numLevels);
}

template <
    template<typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBMultigrid<OPType, MAP, T, BCType, MEM>::define(
        const MBDisjointBoxLayout&  a_fineLayout,
        const std::vector<Point>&   a_refRatios,
        int                         a_numLevels)
{
    PR_TIME("MBMultigrid::define");

    // delete any existing data
    m_state_0.clear();
    m_state.clear();
    m_force.clear();
    m_residual.clear();
    m_local.clear();
    
    m_maps.clear();
    m_blockInterp.clear();
    
    m_ops.clear();
    
    m_increment.clear();
    m_average.clear();
    m_CFInterp.clear();

    // create coarse domains
    std::vector<MBProblemDomain>      problemDomains;
    std::vector<std::vector<Point>>   boxSizes;
    std::vector<MBDisjointBoxLayout> layouts;
    {
        PR_TIME("MBMultigrid::define::build_layouts");
        problemDomains.resize(a_numLevels);
        boxSizes.resize(a_numLevels);
        MBProblemDomain nextDomain = a_fineLayout.domain();
        std::vector<Point> nextBoxSizes = a_fineLayout.boxSizes();
        int level = a_numLevels-1;
        while (level >= 0)
        {
            /*
            std::cout << "Level: " << level << std::endl;
            for (int bi = 0; bi < a_fineLayout.numBlocks(); bi++)
            {
                std::cout << "\tBlock: " << bi << std::endl;
                std::cout << "\tDomain Box: " << nextDomain.getBlock(bi).box() << std::endl;
                std::cout << "\tBox Size: " << nextBoxSizes[bi] << std::endl;
            }
            */

            problemDomains[level] = nextDomain;
            boxSizes[level] = nextBoxSizes;
            if (level > 0)
            {
                // create the next coarser domain
                nextDomain = nextDomain.coarsen(a_refRatios);
                // make sure boxSize <= domainSize
                for (int bi = 0; bi < a_fineLayout.numBlocks(); bi++)
                {
                    Box domainBox = nextDomain.getBlock(bi).box();
                    if (domainBox.sizes() < nextBoxSizes[bi])
                    {
                        nextBoxSizes[bi] = domainBox.sizes();
                    }
                }
            }
            level--;
        }

        // Build the layouts
        for (int li = 0; li < a_numLevels; li++)
        {
            layouts.push_back(MBDisjointBoxLayout(problemDomains[li], boxSizes[li]));
        }
    }
    // initialize all of the state data
    Array<Point,DIM+1> ghost;
    Array<Point,DIM+1> noGhost;
    ghost.fill(OP::ghost());
    noGhost.fill(Point::Zeros());

    m_residual.resize(a_numLevels);
    m_state.resize(a_numLevels-1);
    m_state_0.resize(a_numLevels-1);
    m_force.resize(a_numLevels-1);
    m_local.resize(a_numLevels-1);
    
    for (int li = 0; li < a_numLevels; li++)
    {
        PR_TIME("MBMultigrid::define::allocate_data");
        m_residual[li] = std::make_shared<LevelStateData>(layouts[li], noGhost);
        
        // we don't need these quantities on the finest level
        if (li < a_numLevels-1)
        {
            m_state[li] = std::make_shared<LevelStateData>(layouts[li], ghost);
            m_state_0[li] = std::make_shared<LevelStateData>(layouts[li], ghost);
            m_force[li] = std::make_shared<LevelStateData>(layouts[li], noGhost);
           
            // fine index space, coarse grid. Used for CF operations.
            auto crseFineLayout = layouts[li+1].coarsen(a_refRatios);
            m_local[li] = std::make_shared<LevelStateData>(crseFineLayout, noGhost);
        }
    }

    // build maps, operators, lambda, and increment stencils
    m_maps.resize(a_numLevels);
    m_ops.resize(a_numLevels);
    m_lambda.resize(a_numLevels);
    m_increment.resize(a_numLevels);
    for (int li = 0; li < a_numLevels; li++)
    {
        PR_TIME("MBMultigrid::define::build_maps");
        m_maps[li].define(layouts[li], ghost);
        m_ops[li].define(m_maps[li]);
        m_lambda[li] = 0.5 / m_ops[li].spectralRadius();
        m_increment[li] = (-m_lambda[li])*Shift::Zeros();
    }

    // allocate space for block interpolation operators
    m_blockInterp.resize(a_numLevels);
    for (int li = 0; li < a_numLevels; li++)
    {
        PR_TIME("MBMultigrid::define::build_interp_ops");
        m_blockInterp[li] = std::make_shared<MBInterpOp>(m_maps[li]);
    }

    // define the coarse-fine stencil operators
    for (int bi = 0; bi < a_fineLayout.numBlocks(); bi++)
    {
        PR_TIME("MBMultigrid::define::build_stencils");
        m_average.push_back(Stencil<T>::AvgDown(a_refRatios[bi]));
        m_CFInterp.push_back(InterpStencil<T>::Constant(a_refRatios[bi]));
    }
    
    // define iteration limits
    m_numPreRelax = 4;
    m_numPostRelax = 4;
    m_numBottomRelax = 20;
    m_numLevels = a_numLevels;
}

template <
    template<typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
bool MBMultigrid<OPType, MAP, T, BCType, MEM>::validate(
        const MBDisjointBoxLayout&  a_fineLayout,
        const std::vector<Point>&   a_refRatios,
        int                         a_numLevels)
{
    if (m_state.size() != a_numLevels-1) {return false;}
    if (m_state_0.size() != a_numLevels-1) {return false;}
    if (m_force.size() != a_numLevels-1) {return false;}
    if (m_residual.size() != a_numLevels) {return false;}
    if (m_local.size() != a_numLevels-1) {return false;}

    for (int li = a_numLevels-1; li > 0; li--)
    {
        auto& layout = m_residual[li]->layout();
        if (!m_local[li-1]->layout().compatible(layout)) {return false;}
        if (li < a_numLevels-1)
        {
            if (!m_state[li]->layout().compatible(layout)) {return false;}
            if (!m_state_0[li]->layout().compatible(layout)) {return false;}
            if (!m_force[li]->layout().compatible(layout)) {return false;}
        }
    }
    return true;
}

template <
    template<typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
T MBMultigrid<OPType, MAP, T, BCType, MEM>::solve(
        LevelStateData& a_state,
        LevelStateData& a_force,
        int             a_maxIter,
        T               a_tolerance)
{
    m_solveIter = 0;
    int fineLevel = numLevels()-1;
    T res0 = resnorm(*m_residual[fineLevel], a_state, a_force, fineLevel);
    T res = res0;
#if PR_VERBOSE > 0
    pr_out() << "MBMultigrid::solve | initial residual: " << res0 << std::endl;
#endif
    for (int ii = 0; ii < a_maxIter; ii++)
    {
        vCycle(a_state, a_force, fineLevel);
        res = resnorm(*m_residual[fineLevel], a_state, a_force, fineLevel);
#if PR_VERBOSE > 0
        pr_out() << "MBMultigrid::solve | iter: " << ii << " | residual: " << res << std::endl;
#endif
        if (res < a_tolerance*res0) { break; }
        m_solveIter++;
    }
    return res;
}

template <
    template<typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBMultigrid<OPType, MAP, T, BCType, MEM>::fineCorrect(
        LevelStateData&         a_fine,
        const LevelStateData&   a_crse,
        LevelStateData&         a_crse0,
        int                     a_level)
{
    // Do coarse correction
    for (auto iter : a_crse.layout())
    {
        auto& crse_i  = a_crse[iter];
        auto& crse0_i = a_crse0[iter];

        crse0_i -= crse_i;
        crse0_i *= -1;
    }
    a_crse0.copyTo(*m_local[a_level-1]);

    // increment with the correction
    for (auto iter : a_fine.layout())
    {
        auto block = a_fine.layout().block(iter);
        auto& fine_i = a_fine[iter];
        auto& crse_i = (*m_local[a_level-1])[iter];
        fine_i += m_CFInterp[block](crse_i);
    }
}

template <
    template<typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBMultigrid<OPType, MAP, T, BCType, MEM>::averageDown(
        LevelStateData&         a_crse,
        const LevelStateData&   a_fine,
        int                     a_level)
{
    for (auto iter : a_fine.layout())
    {
        auto block = a_fine.layout().block(iter);
        auto& crse_i = (*m_local[a_level-1])[iter];
        auto& fine_i = a_fine[iter];
        crse_i |= m_average[block](fine_i);
    }
    m_local[a_level-1]->copyTo(a_crse);
}

template <
    template<typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBMultigrid<OPType, MAP, T, BCType, MEM>::crseForce(
        LevelStateData&         a_crseForce,
        LevelStateData&         a_state,
        const LevelStateData&   a_force,
        LevelStateData&         a_crseState,
        int                     a_level)
{
    // fill ghost cells on state data
    a_state.exchange();
    m_blockInterp[a_level]->apply(a_state, a_state);
    a_crseState.exchange();
    m_blockInterp[a_level-1]->apply(a_crseState, a_crseState);

    for (auto iter : a_state.layout())
    {
        auto block = a_state.layout().block(iter);
        auto& state_i = a_state[iter];
        auto& resC_i = (*m_local[a_level-1])[iter];
        auto& force_i = a_force[iter];
        auto res = m_ops[a_level][iter](state_i, a_state.layout()[iter]);
        res += force_i;
        resC_i |= m_average[block](res);
    }
    
    m_local[a_level-1]->copyTo(a_crseForce);
    for (auto iter : a_crseState.layout())
    {
        auto& state_i = a_crseState[iter];
        auto& force_i = a_crseForce[iter];
        auto lphi = m_ops[a_level-1][iter](state_i, a_crseState.layout()[iter]);
        force_i -= lphi;
    }
}

template <
    template<typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBMultigrid<OPType, MAP, T, BCType, MEM>::relax(
        LevelStateData&         a_state,
        const LevelStateData&   a_force,
        int                     a_numIter,
        int                     a_level)
{
    m_relaxIter = 0;
    HDF5Handler h5;
#if PR_VERBOSE > 1
    h5.writeMBLevel({"rhs"}, m_maps[a_level], a_force,
            "MBMultigrid_JF_S%i_L%i_R%i", m_solveIter, m_solveLevel, m_relaxIter);
#endif
    for (int ii = 0; ii < a_numIter; ii++)
    {
        a_state.exchange();
        m_blockInterp[a_level]->apply(a_state, a_state);
#if PR_VERBOSE > 1
        h5.writeMBLevel({"phi"}, m_maps[a_level], a_state,
                "MBMultigrid_Phi_S%i_L%i_R%i", m_solveIter, m_solveLevel, m_relaxIter);
#endif
        for (auto iter : a_state.layout())
        {
            auto& state_i = a_state[iter];
            auto& force_i = a_force[iter];
            auto& res_i = (*m_residual[a_level])[iter];
            m_ops[a_level][iter](res_i, state_i, -1);
            res_i += force_i;
            state_i += m_increment[a_level](res_i);
        }
#if PR_VERBOSE > 1
        h5.writeMBLevel({"residual"}, m_maps[a_level], *m_residual[a_level],
                "MBMultigrid_Res_S%i_L%i_R%i", m_solveIter, m_solveLevel, m_relaxIter);
    
        auto norm = resnorm(*m_residual[a_level], a_state, a_force, a_level);
        pr_out() << "\tMBMultigrid::relax | iter: " << ii << " | residual: " << norm << std::endl;
#endif
        m_relaxIter++;
    }
}

template <
    template<typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBMultigrid<OPType, MAP, T, BCType, MEM>::vCycle(
        LevelStateData&         a_state,
        const LevelStateData&   a_force,
        int                     a_level)
{
    HDF5Handler h5;
    m_solveLevel = a_level;
    //h5.writeMBLevel(m_maps[a_level], a_state, "MBMultigrid_vCycle_Phi_S%i_L%i_0", m_solveIter, m_solveLevel);
    if (a_level > 0)
    {
        //pre relax
        relax(a_state, a_force, m_numPreRelax, a_level);
        //average down state
        //h5.writeMBLevel(m_maps[a_level], a_state, "MBMultigrid_vCycle_Phi_S%i_L%i_1", m_solveIter, m_solveLevel);
        averageDown(*m_state[a_level-1], a_state, a_level);
        //h5.writeMBLevel(m_maps[a_level-1], *m_state[a_level-1], "MBMultigrid_vCycle_PhiC_S%i_L%i_0", m_solveIter, m_solveLevel);
        //save copy of coarse state
        m_state[a_level-1]->copyTo(*m_state_0[a_level-1]);
        //h5.writeMBLevel(m_maps[a_level-1], *m_state_0[a_level-1], "MBMultigrid_vCycle_PhiC0_S%i_L%i_0", m_solveIter, m_solveLevel);
        //compute coarse RHS
        crseForce(*m_force[a_level-1], a_state, a_force, *m_state[a_level-1], a_level);
        //h5.writeMBLevel(m_maps[a_level-1], *m_force[a_level-1], "MBMultigrid_vCycle_FC_S%i_L%i_0", m_solveIter, m_solveLevel);
        //recursive vcycle on coarse state/RHS
        vCycle(*m_state[a_level-1], *m_force[a_level-1], a_level-1);
        m_solveLevel = a_level;
        //h5.writeMBLevel(m_maps[a_level-1], *m_state[a_level-1], "MBMultigrid_vCycle_PhiC_S%i_L%i_1", m_solveIter, m_solveLevel);
        //correct this level's state
        fineCorrect(a_state, *m_state[a_level-1], *m_force[a_level-1], a_level);
        //h5.writeMBLevel(m_maps[a_level], a_state, "MBMultigrid_vCycle_Phi_S%i_L%i_2", m_solveIter, m_solveLevel);
        //post relax
        relax(a_state, a_force, m_numPostRelax, a_level);
        //h5.writeMBLevel(m_maps[a_level], a_state, "MBMultigrid_vCycle_Phi_S%i_L%i_3", m_solveIter, m_solveLevel);
    } else {
        // bottom relax
        relax(a_state, a_force, m_numBottomRelax, a_level);
        //h5.writeMBLevel(m_maps[a_level], a_state, "MBMultigrid_vCycle_Phi_S%i_L%i_1", m_solveIter, m_solveLevel);
    }
}

template <
    template<typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBMultigrid<OPType, MAP, T, BCType, MEM>::residual(
        LevelStateData&         a_residual,
        LevelStateData&         a_state,
        const LevelStateData&   a_force,
        int                     a_level)
{
    a_state.exchange();
    m_blockInterp[a_level]->apply(a_state, a_state);
    for (auto iter : a_residual.layout())
    {
        auto& res_i = a_residual[iter];
        auto& state_i = a_state[iter];
        auto& force_i = a_force[iter];

        m_ops[a_level][iter](res_i, state_i, -1);
        res_i += force_i;
    }
}

template <
    template<typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
T MBMultigrid<OPType, MAP, T, BCType, MEM>::resnorm(
        LevelStateData&         a_residual,
        LevelStateData&         a_state,
        const LevelStateData&   a_force,
        int                     a_level)
{
    a_state.exchange();
    m_blockInterp[a_level]->apply(a_state, a_state);
    Reduction<T,Abs> rxn;
    for (auto iter : a_residual.layout())
    {
        auto& res_i = a_residual[iter];
        auto& state_i = a_state[iter];
        auto& force_i = a_force[iter];
        m_ops[a_level][iter](res_i, state_i, -1);
        res_i += force_i;
        res_i.reduce(rxn);
    }
    return rxn.fetch();
}
