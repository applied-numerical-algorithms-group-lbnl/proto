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
        m_increment[li] = m_lambda[li]*Shift::Zeros();
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
        T               a_tolerance){}

template <
    template<typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBMultigrid<OPType, MAP, T, BCType, MEM>::fineCorrect(
        LevelStateData&         a_fine,
        const LevelStateData&   a_crse,
        const LevelStateData&   a_crse0){}

template <
    template<typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBMultigrid<OPType, MAP, T, BCType, MEM>::crseForce(
        LevelStateData&         a_crseForce,
        const LevelStateData&   a_state,
        const LevelStateData&   a_force,
        const LevelStateData&   a_crseState){}

template <
    template<typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBMultigrid<OPType, MAP, T, BCType, MEM>::relax(
        LevelStateData&         a_state,
        const LevelStateData&   a_rhs,
        int                     a_numIter){}

template <
    template<typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBMultigrid<OPType, MAP, T, BCType, MEM>::vCycle(
        LevelStateData&         a_state,
        const LevelStateData&   a_force){}

template <
    template<typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBMultigrid<OPType, MAP, T, BCType, MEM>::residual(
        LevelStateData&         a_residual,
        const LevelStateData&   a_state,
        const LevelStateData&   a_force){}
