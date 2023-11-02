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
    std::vector<ProblemDomain>      problemDomains;
    std::vector<std::vector<Point>> boxSizes;
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
                    nextBoxSizes[bi] = domainBox;
                }
            }
            level--;
        }
    }

    // Build the layouts
    std::vector<MBDisjointBoxLayout> layouts;
    for (int li = 0; li < a_numLevels; li++)
    {
        layouts.push_back(MBDisjointBoxLayout(problemDomains[li], boxSizes[li]));
    }

    // initialize all of the state data
    Array<Point,DIM+1> ghost;
    Array<Point,DIM+1> noGhost;
    ghost.fill(OP::ghost());
    noGhost.fill(Point::Zeros());
    for (int li = 0; li < a_numLevels; li++)
    {
        m_residual.push_back(LevelStateData(layouts[li], noGhost));
        
        // we don't need these quantities on the finest level
        if (li < a_numLevels-1)
        {
            m_state.push_back(   LevelStateData(layouts[li], ghost));
            m_state_0.push_back( LevelStateData(layouts[li], ghost));
            m_force.push_back(   LevelStateData(layouts[li], noGhost));
           
            // fine index space, coarse grid. Used for CF operations.
            auto crseFineLayout = layouts[li+1].coarsen(a_refRatios);
            m_local.push_back(   LevelStateData(crseFineLayout, noGhost));
        }
    }

    // build maps, operators, lambda, and increment stencils
    m_maps.resize(a_numLevels);
    m_ops.resize(a_numLevels);
    m_lambda.resize(a_numLevels);
    m_increment.resize(a_numLevels);
    for (int li = 0; li < a_numLevels; li++)
    {
        m_maps[li].define(layouts[li], ghost);
        m_ops[li].define(m_maps[li]);
        m_lambda[li] = 0.5 / m_ops[li].spectralRadius();
        m_increment[li] = m_lambda[li]*Shift::Zeros();
    }


    // allocate space for block interpolation operators
    m_blockInterp.resize(a_numLevels);
   
    // define the coarse-fine stencil operators
    for (int bi = 0; bi < a_fineLayout.numBlocks(); bi++)
    {
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
