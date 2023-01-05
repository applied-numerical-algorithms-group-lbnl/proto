template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
AMRSolver_FASMultigrid<OpType, T, MEM>::AMRSolver_FASMultigrid(
    AMRGrid& a_grid,
    Array<T, DIM>&   a_dx)
{
    define(a_grid, a_dx);
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
AMRSolver_FASMultigrid<OpType, T, MEM>::define(
    AMRGrid& a_grid,
    Array<T, DIM>&   a_dx)
{
    m_residual.define(a_grid, Point::Zeros());
    m_amrOp = std::make_shared<AMROP>(a_grid, a_dx);
    m_amrLevel = std::make_shared<AMRLevel>(a_grid, *m_amrOp, a_grid.numLevels() - 1);
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
double
AMRSolver_FASMultigrid<OpType, T, MEM>::solve(
    AMRStateData& a_state,
    AMRStateData& a_force,
    AMRAuxData&   a_aux,
    int           a_maxIter,
    double        a_tolerance)
{
    return 0;
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
double
AMRSolver_FASMultigrid<OpType, T, MEM>::solve(
    AMRStateData& a_state,
    AMRStateData& a_force,
    int           a_maxIter,
    double        a_tolerance)
{
#ifdef PR_HDF5
    HDF5Handler h5;
#endif
    auto& crseLayout = a_state.grid()[0];
    double res0 = resnorm(m_residual, a_state, a_force);
#ifdef PR_HDF5
    h5.writeAMRData(m_amrOp->dx(0), m_residual, "AMR_RES_I%i", 0);
    h5.writeAMRData(m_amrOp->dx(0), a_state, "AMR_STATE_I%i", 0);
    h5.writeAMRData(m_amrOp->dx(0), a_force, "AMR_FORCE");
#endif
    //Proto::pout() << "AMR solve: initial residual = " << res0 << std::endl;
    if (procID() == 0) { 
        std::cout << "AMR solve: initial residual = " << res0 << std::endl;
    }
    double res;
    int finestLevel = a_state.grid().numLevels() - 1;
    for (int ii = 0; ii < a_maxIter; ii++)
    {
        m_amrLevel->vCycle(a_state, a_force, a_force[finestLevel]);
        res = resnorm(m_residual, a_state, a_force);
#ifdef PR_HDF5
        h5.writeAMRData(m_amrOp->dx(0), m_residual, "AMR_RES_I%i", ii+1);
        h5.writeAMRData(m_amrOp->dx(0), a_state, "AMR_STATE_I%i", ii+1);
#endif
        //Proto::pout() << "\tAMR solve: resnorm = " << res;
        //pout() << " | integral = " << m_residual.integrate(m_amrOp->dx(0)) << std::endl;
        if (procID() == 0)
        {
            std::cout << "\tAMR solve: resnorm = " << res;
            std::cout << " | integral = " << m_residual.integrate(m_amrOp->dx(0)) << std::endl;
        }
        if (res < a_tolerance*res0) { break; }
    }
    return res;
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
double
AMRSolver_FASMultigrid<OpType, T, MEM>::resnorm(
    AMRStateData& a_residual,
    AMRStateData& a_state,
    AMRStateData& a_force)
{
    a_state.averageDown();
    for (int lvl = 1; lvl < a_state.numLevels(); lvl++)
    {
        auto& interpBC = m_amrOp->interpStencil(lvl-1);
        interpBoundaries(a_state[lvl-1], a_state[lvl], interpBC);
    }
    a_state.exchange();
    m_amrOp->operator()(a_residual, a_state, -1); //apply with reflux
    a_residual.increment(a_force);
    a_residual.averageDown();
    return a_residual.absMax();
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
AMRSolver_FASMultigrid<OpType, T, MEM>::AMRLevel::AMRLevel(
    AMRGrid& a_grid,
    AMROP&   a_op,
    int      a_level)
{
    define(a_grid, a_op, a_level);
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
AMRSolver_FASMultigrid<OpType, T, MEM>::AMRLevel::define(
    AMRGrid& a_grid,
    AMROP&   a_op,
    int      a_level)
{
    m_level = a_level;

    m_numPreRelax = 2*DIM;
    m_numPostRelax = 2*DIM;
    m_numBottomRelax = 10;

    m_op = &a_op;
    
    m_residual.define(a_grid[a_level], Point::Zeros()); 
    if (a_level > 0)
    {
        Point refRatio = a_grid.refRatio(a_level-1);

        // define the level MG solver
        int numMGLevels = 1;
        m_mgSolver = std::make_shared<LevelSolver_FASMultigrid<OpType, T, MEM>>(
            a_grid[m_level], Point::Ones(2), numMGLevels, m_op->dx(a_level));
   
        // fine correction stencil 
        m_interpolate = InterpStencil<T>::PiecewiseConstant(refRatio);

        // define temporaries
        Point interpGhost = m_op->interpStencil(a_level-1).ghost();
        DisjointBoxLayout crseLayout = a_grid[a_level].coarsen(refRatio);
        // TODO: we probably don't need both of these
        m_crseLocalInterp.define(crseLayout, BOP::ghost() + interpGhost);
        m_crseLocal.define(crseLayout, Point::Zeros());

        m_crseRes.define(a_grid[a_level-1], Point::Zeros());
        m_crseState_0.define(a_grid[a_level-1], Point::Zeros());
        m_crseForce.define(a_grid[a_level-1], Point::Zeros());
        // define next level
        m_crseAMRLevel = std::make_shared<AMRLevel>(a_grid, a_op, a_level - 1);
         
    } else {
        // define the level MG solver
        Point domainSizeV = a_grid[0].domain().sizes();
        int minDomainSize = domainSizeV[0];
        for (int dir = 1; dir < DIM; dir++)
        {
            minDomainSize = min(domainSizeV[dir], minDomainSize);
        }
        int numMGLevels = log(minDomainSize*1.0) / log(2.0);
        m_mgSolver = std::make_shared<LevelSolver_FASMultigrid<OpType, T, MEM>>(
            a_grid[0], Point::Ones(2), numMGLevels, m_op->dx(a_level));
    }
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
AMRSolver_FASMultigrid<OpType, T, MEM>::AMRLevel::fineCorrect(
        LevelStateData& a_state,
        LevelStateData& a_crseState,
        LevelStateData& a_crseState_0)
{
    for (auto iter = a_crseState.begin(); iter.ok(); ++iter)
    {
        auto& crse_i  = a_crseState[*iter];
        auto& crse0_i = a_crseState_0[*iter];
        crse0_i -= crse_i;
        crse0_i *= -1; //should be crse - crse0;
    }
    a_crseState_0.copyTo(m_crseLocal);
    for (auto iter = a_state.begin(); iter.ok(); ++iter)
    {
        auto& state_i = a_state[*iter];
        auto& delta_i = m_crseLocal[*iter];

        state_i += m_interpolate(delta_i);
    }
}
template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
AMRSolver_FASMultigrid<OpType, T, MEM>::AMRLevel::crseResidual(
        LevelStateData& a_ResC,
        AMRStateData&   a_Phi,
        LevelStateData& a_GC,
        LevelStateData& a_G)
{
#ifdef PR_HDF5
    HDF5Handler h5;
#endif
    Point refRatio = a_Phi.grid().refRatio(m_level-1);

    auto& Phi  = a_Phi[m_level];
    auto& PhiC = a_Phi[m_level-1];

    auto& interpBC = m_op->interpStencil(m_level-1);
    interpBoundaries(PhiC, Phi, m_crseLocalInterp, interpBC); 
    Phi.exchange(); 
    m_op->levelApply(a_ResC, a_Phi, m_level-1, -1.0);   //-LC(PhiC) with refluxing
    a_ResC.increment(a_GC); // GC - LC(PhiC)

    auto AVG = Stencil<T>::AvgDown(refRatio);
    for (auto iter = Phi.begin(); iter.ok(); ++iter)
    {
        auto& phi_i  = Phi[*iter];
        auto& g_i    = a_G[*iter];
        auto& temp_i = m_crseLocal[*iter];
        
        auto res_i = m_op->boxOp(m_level)(phi_i, Phi.layout()[*iter], -1.0); // -L(Phi)
        res_i += g_i;
        temp_i |= AVG(res_i); // m_crseLocal now holds average of residual.
    }
    m_crseLocal.copyTo(a_ResC); // overwrite only the refined region. 
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
AMRSolver_FASMultigrid<OpType, T, MEM>::AMRLevel::vCycle(
        AMRStateData& a_state,
        AMRStateData& a_force,
        LevelStateData& a_levelForce)
{
    PR_TIMERS("AMRFAS::VCycle");
#ifdef PR_HDF5
    HDF5Handler h5; 
#endif
    auto dx = m_op->dx(m_level);
    auto& state = a_state[m_level];
    if (m_level == 0)
    {
        m_mgSolver->solve(state, a_levelForce, m_numBottomRelax, 0); 
    }
    else
    {
        auto& crseState = a_state[m_level - 1];
        auto& crseForce = a_force[m_level - 1];

        // interpolate boundary conditions
        auto& interpBC = m_op->interpStencil(m_level-1);
        interpBoundaries(crseState, state, m_crseLocalInterp, interpBC); 
        // smooth
        m_mgSolver->solve(state, a_levelForce, m_numPreRelax, 0); //0 tolerance = do max number of iterations 
        // average down
        averageDown(crseState, state, m_crseLocal, a_state.grid().refRatio(m_level-1));
        // save a copy of the coarse state
        crseState.copyTo(m_crseState_0);
        // compute the coarse force
        crseResidual(m_crseRes, a_state, crseForce, a_levelForce);
        //m_op->levelApply(m_crseForce, a_state, m_level-1);
        m_op->levelOp(m_level-1)(m_crseForce, crseState);
        m_crseForce.increment(m_crseRes); // LC(PhiC) + ResC
        // recursive call
        m_crseAMRLevel->vCycle(a_state, a_force, m_crseForce); 
        // fix up this level's state
        fineCorrect(state, crseState, m_crseState_0);
        // interpolate boundary conditions
        interpBoundaries(crseState, state, m_crseLocalInterp, interpBC); 
        // smooth
        m_mgSolver->solve(state, a_levelForce, m_numPostRelax, 0); // 0 tolerance = do max number of iterations
    }
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
AMRSolver_FASMultigrid<OpType, T, MEM>::AMRLevel::levelResidual(
        LevelStateData& a_res,
        AMRStateData& a_state,
        AMRStateData& a_force, 
        int a_level)
{
    a_state.exchange(); 
    m_op->levelApply(a_res, a_state, a_level, -1); //apply with refluxing
    for (auto iter = a_res.begin(); iter.ok(); ++iter)
    {
        auto& res_i = a_res[*iter];
        auto& force_i = a_force[a_level][*iter];
        res_i += force_i;
    }
}
