template <
    template<typename, typename, MemType> class OPType,
    class MAP,
    typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType,
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
    template<typename, typename, MemType> class OPType,
    class MAP,
    typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType,
    MemType MEM>
MBMultigrid<OPType, MAP, T, BCType, MEM>::MBMultigrid(
        const MBDisjointBoxLayout&  a_fineLayout,
        const std::vector<Point>&   a_refRatios,
        int                         a_numLevels)
{
    define(a_fineLayout, a_refRatios, a_numLevels);
}

template <
    template<typename, typename, MemType> class OPType,
    class MAP,
    typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType,
    MemType MEM>
void MBMultigrid<OPType, MAP, T, BCType, MEM>::define(
        const MBDisjointBoxLayout&  a_fineLayout,
        const std::vector<Point>&   a_refRatios,
        int                         a_numLevels)
{
    PR_TIME("MBMultigrid::define");

    m_numLevels = a_numLevels;

    // delete any existing data
    m_state_0.clear();
    m_state.clear();
    m_force.clear();
    m_residual.clear();
    m_local.clear();
    
    m_maps.clear();
    m_blockInterp.clear();
    
    m_ops.clear();
    
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
        m_residual[li] = std::make_shared<LevelStateData>(layouts[li], Point::Zeros());
        
        // we don't need these quantities on the finest level
        if (li < a_numLevels-1)
        {
            m_state[li] = std::make_shared<LevelStateData>(layouts[li], OP::ghost());
            m_state[li]->setVal(0);
            m_state_0[li] = std::make_shared<LevelStateData>(layouts[li], OP::ghost());
            m_state_0[li]->setVal(0);
            m_force[li] = std::make_shared<LevelStateData>(layouts[li], Point::Zeros());
           
            // fine index space, coarse grid. Used for CF operations.
            auto crseFineLayout = layouts[li+1].coarsen(a_refRatios);
            m_local[li] = std::make_shared<LevelStateData>(crseFineLayout, Point::Zeros());
        }
    }

    // build maps, operators, lambda, and increment stencils
    m_maps.resize(a_numLevels);
    m_ops.resize(a_numLevels);
    m_lambda.resize(a_numLevels);
    for (int li = 0; li < a_numLevels; li++)
    {
        PR_TIME("MBMultigrid::define::build_maps");
        m_maps[li].define(layouts[li], OP::ghost());
        
        m_ops[li].define(m_maps[li]);
        m_lambda[li] = 0.5 / m_ops[li].spectralRadius();
    }

    // allocate space for block interpolation operators
    m_blockInterp.resize(a_numLevels);
    for (int li = 0; li < a_numLevels; li++)
    {
        PR_TIME("MBMultigrid::define::build_interp_ops");
        m_blockInterp[li] = std::make_shared<MBInterpOp>(m_maps[li], 5);
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
    m_numBottomRelax = 2;

    m_solveIter = 0;
    m_solveLevel = m_numLevels-1;

    // print domains
    // for (int li = 0; li < a_numLevels; li++)
    // {
    //     pr_out() << std::endl << "LEVEL: " << li << std::endl;
    //     pr_out() << "Layout of state, state_0, force, residual: (size = " << layouts[li].numBoxes() << ")" << std::endl;
    //     layouts[li].print();
    //     if (li < a_numLevels-1)
    //     {
    //         pr_out() << "Layout of local: (size = " << m_local[li]->layout().numBoxes() << ")" << std::endl;
    //         m_local[li]->layout().print();
    //     }
    // }
    PROTO_ASSERT(validate(a_fineLayout, a_refRatios, a_numLevels), "MBMultigrid Object Failed Validation on Construction");
}

template <
    template<typename, typename, MemType> class OPType,
    class MAP,
    typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType,
    MemType MEM>
bool MBMultigrid<OPType, MAP, T, BCType, MEM>::validate(
        const MBDisjointBoxLayout&  a_fineLayout,
        const std::vector<Point>&   a_refRatios,
        int                         a_numLevels)
{
    if (a_numLevels != m_numLevels) {return false;}
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
    template<typename, typename, MemType> class OPType,
    class MAP,
    typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType,
    MemType MEM>
T MBMultigrid<OPType, MAP, T, BCType, MEM>::solve(
        LevelStateData& a_state,
        const LevelStateData& a_force,
        int             a_maxIter,
        T               a_tolerance)
{
    m_solveIter = 0;
    int fineLevel = numLevels()-1;
    T res0 = resnorm(*m_residual[fineLevel], a_state, a_force, fineLevel);
    T res = res0;
#if PR_VERBOSE > 0
    HDF5Handler h5;
    pr_out()  << "MBMultigrid::solve | initial residual: " << res0 << std::endl;
    std::cout << "MBMultigrid::solve | initial residual: " << res0 << std::endl;
    h5.writeMBLevel(m_maps[fineLevel], *m_residual[fineLevel],
            "MBMultigrid_Solve_Res_S%i", m_solveIter);
#endif
    for (int ii = 0; ii < a_maxIter; ii++)
    {
        vCycle(a_state, a_force, fineLevel);
        res = resnorm(*m_residual[fineLevel], a_state, a_force, fineLevel);
        m_solveIter++;
#if PR_VERBOSE > 0
        if (procID() == 0)
        {
            pr_out()  << "MBMultigrid::solve | after iter: " << ii << " | residual: " << res << std::endl;
            std::cout << "MBMultigrid::solve | after iter: " << ii << " | residual: " << res << std::endl;
        }
        
        h5.writeMBLevel(m_maps[fineLevel], *m_residual[fineLevel],
                "MBMultigrid_Solve_Res_S%i", m_solveIter);
        
#endif
        if (res < a_tolerance*res0) { break; }
    }
    return res;
}

template <
    template<typename, typename, MemType> class OPType,
    class MAP,
    typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType,
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
    template<typename, typename, MemType> class OPType,
    class MAP,
    typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType,
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
    template<typename, typename, MemType> class OPType,
    class MAP,
    typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType,
    MemType MEM>
void MBMultigrid<OPType, MAP, T, BCType, MEM>::crseForce(
        LevelStateData&         a_crseForce,
        LevelStateData&         a_state,
        const LevelStateData&   a_force,
        LevelStateData&         a_crseState,
        int                     a_level)
{
    interpBlockBoundaries(a_state, a_level);
    interpBlockBoundaries(a_crseState, a_level-1);

    for (auto iter : a_state.layout())
    {
        auto block = a_state.layout().block(iter);
        auto& state_i = a_state[iter];
        auto& resC_i = (*m_local[a_level-1])[iter];
        auto& force_i = a_force[iter];
        auto res = m_ops[a_level][iter](state_i, a_state.layout()[iter],-1);
        res += force_i;
        resC_i |= m_average[block](res);
    }
    m_local[a_level-1]->copyTo(a_crseForce);
    for (auto iter : a_crseState.layout())
    {
        auto& state_i = a_crseState[iter];
        auto& force_i = a_crseForce[iter];
        auto lphi = m_ops[a_level-1][iter](state_i, a_crseState.layout()[iter]);
        force_i += lphi;
    }
}

template <
    template<typename, typename, MemType> class OPType,
    class MAP,
    typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType,
    MemType MEM>
void MBMultigrid<OPType, MAP, T, BCType, MEM>::relax(
        LevelStateData&         a_state,
        const LevelStateData&   a_force,
        int                     a_numIter,
        int                     a_level)
{
    m_relaxIter = 0;
    // lambda ~ 0.5 / spectral radius
    Stencil<T> INCR = (-lambda(a_level))*Shift::Zeros();
    #if PR_VERBOSE > 1
    HDF5Handler h5;
    h5.writeMBLevel(m_maps[a_level], a_force, "MBMultigrid_Relax_Force_S%i_L%i", m_solveIter, m_solveLevel);
    #endif
    for (int ii = 0; ii < a_numIter; ii++)
    {
        #if PR_VERBOSE > 1
        h5.writeMBLevel(m_maps[a_level], a_state, "MBMultigrid_Relax_Phi_S%i_L%i_R%i_0", m_solveIter, m_solveLevel, m_relaxIter);
        h5.writeMBLevelBoundsUnified({"phi"}, a_state, "MBMultigrid_Relax_PhiBounds_S%i_L%i_R%i_0", m_solveIter, m_solveLevel, m_relaxIter);
        #endif
        interpBlockBoundaries(a_state, a_level);
        #if PR_VERBOSE > 1
        h5.writeMBLevel(m_maps[a_level], a_state, "MBMultigrid_Relax_Phi_S%i_L%i_R%i_1", m_solveIter, m_solveLevel, m_relaxIter);
        h5.writeMBLevelBoundsUnified({"phi"}, a_state, "MBMultigrid_Relax_PhiBounds_S%i_L%i_R%i_1", m_solveIter, m_solveLevel, m_relaxIter);
        #endif
        residual(*m_residual[a_level], a_state, a_force, a_level);
        #if PR_VERBOSE > 1
        h5.writeMBLevel(m_maps[a_level], a_state, "MBMultigrid_Relax_Phi_S%i_L%i_R%i_2", m_solveIter, m_solveLevel, m_relaxIter);
        h5.writeMBLevel(m_maps[a_level], *m_residual[a_level], "MBMultigrid_Relax_Residual_S%i_L%i_R%i", m_solveIter, m_solveLevel, m_relaxIter);
        std::cout << "Relax | Iter: " << m_relaxIter << " | Level: " << a_level << " | Residual: " << m_residual[a_level]->absMax() << std::endl;
        #endif
        for (auto iter : a_state.layout())
        {
            auto& state_i = a_state[iter];
            auto& res_i = (*m_residual[a_level])[iter];
            state_i += INCR(res_i);
        }
        #if PR_VERBOSE > 1
        h5.writeMBLevel(m_maps[a_level], a_state, "MBMultigrid_Relax_Phi_S%i_L%i_R%i_3", m_solveIter, m_solveLevel, m_relaxIter);
        #endif
        m_relaxIter++;
    }
}

template <
    template<typename, typename, MemType> class OPType,
    class MAP,
    typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType,
    MemType MEM>
void MBMultigrid<OPType, MAP, T, BCType, MEM>::vCycle(
        LevelStateData&         a_state,
        const LevelStateData&   a_force,
        int                     a_level)
{
    HDF5Handler h5;
    m_solveLevel = a_level;
    h5.writeMBLevel(m_maps[a_level], a_state, "MBMultigrid_vCycle_Phi_S%i_L%i_0", m_solveIter, m_solveLevel);
    if (a_level > 0)
    {
        //pre relax
        relax(a_state, a_force, m_numPreRelax, a_level);
        //average down state
        h5.writeMBLevel(m_maps[a_level], a_state, "MBMultigrid_vCycle_Phi_S%i_L%i_1", m_solveIter, m_solveLevel);
        averageDown(*m_state[a_level-1], a_state, a_level);
        h5.writeMBLevel(m_maps[a_level-1], *m_state[a_level-1], "MBMultigrid_vCycle_PhiC_S%i_L%i_0", m_solveIter, m_solveLevel);
        //save copy of coarse state
        m_state[a_level-1]->copyTo(*m_state_0[a_level-1]);
        h5.writeMBLevel(m_maps[a_level-1], *m_state_0[a_level-1], "MBMultigrid_vCycle_PhiC0_S%i_L%i_0", m_solveIter, m_solveLevel);
        //compute coarse RHS
        crseForce(*m_force[a_level-1], a_state, a_force, *m_state[a_level-1], a_level);
        h5.writeMBLevel(m_maps[a_level-1], *m_force[a_level-1], "MBMultigrid_vCycle_FC_S%i_L%i_0", m_solveIter, m_solveLevel);
        //recursive vcycle on coarse state/RHS
        vCycle(*m_state[a_level-1], *m_force[a_level-1], a_level-1);
        m_solveLevel = a_level;
        h5.writeMBLevel(m_maps[a_level-1], *m_state[a_level-1], "MBMultigrid_vCycle_PhiC_S%i_L%i_1", m_solveIter, m_solveLevel);
        //correct this level's state
        fineCorrect(a_state, *m_state[a_level-1], *m_state_0[a_level-1], a_level);
        h5.writeMBLevel(m_maps[a_level], a_state, "MBMultigrid_vCycle_Phi_S%i_L%i_2", m_solveIter, m_solveLevel);
        //post relax
        relax(a_state, a_force, m_numPostRelax, a_level);
        h5.writeMBLevel(m_maps[a_level], a_state, "MBMultigrid_vCycle_Phi_S%i_L%i_3", m_solveIter, m_solveLevel);
    } else {
        // bottom relax
        relax(a_state, a_force, m_numBottomRelax, a_level);
        h5.writeMBLevel(m_maps[a_level], a_state, "MBMultigrid_vCycle_Phi_S%i_L%i_1", m_solveIter, m_solveLevel);
    }
}

template <
    template<typename, typename, MemType> class OPType,
    class MAP,
    typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType,
    MemType MEM>
void MBMultigrid<OPType, MAP, T, BCType, MEM>::residual(
        LevelStateData&         a_residual,
        LevelStateData&         a_state,
        const LevelStateData&   a_force,
        int                     a_level)
{
    interpBlockBoundaries(a_state, a_level);
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
    template<typename, typename, MemType> class OPType,
    class MAP,
    typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType,
    MemType MEM>
T MBMultigrid<OPType, MAP, T, BCType, MEM>::resnorm(
        LevelStateData&         a_residual,
        LevelStateData&         a_state,
        const LevelStateData&   a_force,
        int                     a_level)
{
    LevelStateData LPhi(a_state.layout(), Point::Zeros());
    interpBlockBoundaries(a_state, a_level);
    Reduction<T,Abs> rxn;
    for (auto iter : a_residual.layout())
    {
        auto& res_i = a_residual[iter];
        auto& state_i = a_state[iter];
        auto& force_i = a_force[iter];
        auto& lphi_i = LPhi[iter];
        m_ops[a_level][iter](res_i, state_i, -1);
        res_i.copyTo(lphi_i);
        res_i += force_i;
        res_i.reduce(rxn);
    }
#if PR_VERBOSE > 0
    HDF5Handler h5;
    h5.writeMBLevel({"-lphi"}, m_maps[a_level], LPhi, "RESNORM_LPHI");
    h5.writeMBLevel({"phi"}, m_maps[a_level], a_state, "RESNORM_PHI");
    h5.writeMBLevel({"force"}, m_maps[a_level], a_force, "RESNORM_RHS");
#endif
    return rxn.fetch();
}
