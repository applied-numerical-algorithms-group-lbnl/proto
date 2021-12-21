template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
AMRSolver_FASMultigrid<OpType, T, MEM>::AMRSolver_FASMultigrid(
    AMRGrid& a_grid,
    double   a_dx)
{
    define(a_grid, a_dx);
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
AMRSolver_FASMultigrid<OpType, T, MEM>::define(
    AMRGrid& a_grid,
    double   a_dx)
{
    m_residual.define(a_grid, Point::Zeros());
    m_amrOp = std::make_shared<AMROP>(a_grid, a_dx);
    m_amrOp-> setFluxScale(-1);
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
    HDF5Handler h5;
    auto& crseLayout = a_state.grid()[0];
    double res0 = resnorm(m_residual, a_state, a_force);
    h5.writeAMRData(m_amrOp->dx(0), m_residual, "AMR_RES_I%i", 0);
    h5.writeAMRData(m_amrOp->dx(0), a_state, "AMR_STATE_I%i", 0);
    h5.writeAMRData(m_amrOp->dx(0), a_force, "AMR_FORCE");
    std::cout << "AMR solve: initial residual = " << res0 << std::endl;
    double res;
    int finestLevel = a_state.grid().numLevels() - 1;
    for (int ii = 0; ii < a_maxIter; ii++)
    {
        AMR_SOLVE_ITER = ii;
        m_amrLevel->vCycle(a_state, a_force, a_force[finestLevel]);
        res = resnorm(m_residual, a_state, a_force);
        h5.writeAMRData(m_amrOp->dx(0), m_residual, "AMR_RES_I%i", ii+1);
        h5.writeAMRData(m_amrOp->dx(0), a_state, "AMR_STATE_I%i", ii+1);
        std::cout << "\tAMR solve: resnorm = " << res << std::endl;
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
    a_state.interp();
    a_state.exchange();
    m_amrOp->operator()(a_residual, a_state); //apply with reflux
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

    Point refRatio = Point::Ones(PR_AMR_REFRATIO);

    m_average = Stencil<T>::AvgDown(refRatio);
    m_interpolate = InterpStencil<T>::PiecewiseConstant(refRatio);
    
    m_op = &a_op;
    
    m_residual.define(a_grid[a_level], Point::Zeros()); 
    if (a_level > 0)
    {
        // define the level MG solver
        int numMGLevels = 1;
        //int numMGLevels = log(PR_AMR_REFRATIO*1.0) / log(PR_MG_REFRATIO*1.0);
        m_mgSolver = std::make_shared<LevelSolver_FASMultigrid<OpType, T, MEM>>(
            a_grid[m_level], numMGLevels, m_op->dx(a_level));
    
        // define temporaries
        DisjointBoxLayout crseLayout = a_grid[a_level].coarsen(refRatio);
        m_crseLocal.define(crseLayout, Point::Zeros());
        m_crseLocalState.define(crseLayout, m_op->ghost());
        m_crseState_0.define(a_grid[a_level-1], Point::Zeros());
        m_crseForce.define(a_grid[a_level-1], Point::Zeros());
        // define next level
        m_crseAMRLevel = std::make_shared<AMRLevel>(a_grid, a_op, a_level - 1);
         
    } else {
        // define the level MG solver
        Point domainSizeV = a_grid[0].domain().size();
        int minDomainSize = domainSizeV[0];
        for (int dir = 1; dir < DIM; dir++)
        {
            minDomainSize = min(domainSizeV[dir], minDomainSize);
        }
        int numMGLevels = log(minDomainSize*1.0) / log(2.0);
        m_mgSolver = std::make_shared<LevelSolver_FASMultigrid<OpType, T, MEM>>(
            a_grid[0], numMGLevels, m_op->dx(a_level));
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
AMRSolver_FASMultigrid<OpType, T, MEM>::AMRLevel::coarseForce(
        LevelStateData& a_crseForce,
        LevelStateData& a_fineForce,
        AMRStateData& a_state,
        AMRStateData& a_force)
{
    HDF5Handler h5;
    double dx = m_op->dx(m_level);
    double cdx = dx*PR_AMR_REFRATIO;

    auto& FC    = a_crseForce;
    auto& G     = a_fineForce;
    auto& GC    = a_force[m_level - 1];
    auto& Phi   = a_state[m_level];
    auto& PhiC  = a_state[m_level - 1];
   
    //a_state.interpLevel(m_level); //TODO: should do nothing 
    // TODO: potential non-periodic boundary condition issue here
    PhiC.copyTo(m_crseLocalState);
    m_crseLocalState.exchange(); //TODO: only really needed for corner case where refined region abuts a periodic boundary
    m_op->levelOp(m_level-1)(m_crseLocal, m_crseLocalState);
    double intLPhi_C = -m_crseLocal.integrate(cdx);
    for (auto iter = Phi.begin(); iter.ok(); ++iter)
    {
        auto& phi_i = Phi[*iter];
        auto& G_i = G[*iter];
        auto& F_i = m_crseLocal[*iter]; //currently holding -LPhi_C

        auto res_i = m_op->boxOp(m_level)(phi_i, iter.box());
        res_i += G_i;
        F_i *= -1;  //now it's +Lphi_C
        F_i += m_average(res_i); // now it's Lphi_C + <G_f - Lphi_f>
    }
    GC.copyTo(FC); // initialize non-refined area
    m_crseLocal.copyTo(FC); // overwrite refined area
    m_op->reflux(FC, PhiC, Phi, m_level-1, 1.0); //negative scaling because m_op is -L
}
template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
AMRSolver_FASMultigrid<OpType, T, MEM>::AMRLevel::vCycle(
        AMRStateData& a_state,
        AMRStateData& a_force,
        LevelStateData& a_levelForce)
{
    PR_TIMERS("LevelMultigrid::VCycle");
    HDF5Handler h5; 
    double dx = m_op->dx(m_level);
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
        a_state.interpLevel(m_level);
        // smooth
        m_mgSolver->solve(state, a_levelForce, m_numPreRelax, 0); 
        // average down
        state.coarsenTo(crseState, Point::Ones(PR_AMR_REFRATIO));
        // save a copy of the coarse state
        crseState.copyTo(m_crseState_0);
        // compute the coarse force
        coarseForce(m_crseForce, a_levelForce, a_state, a_force);
        // recursive call
        m_crseAMRLevel->vCycle(a_state, a_force, m_crseForce); 
        // fix up this level's state
        fineCorrect(state, crseState, m_crseState_0);
        // interpolate boundary conditions
        a_state.interpLevel(m_level);
        // smooth
        m_mgSolver->solve(state, a_levelForce, m_numPostRelax, 0); 
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
    m_op->levelApply(a_res, a_state, a_level); //apply with refluxing
    for (auto iter = a_res.begin(); iter.ok(); ++iter)
    {
        auto& res_i = a_res[*iter];
        auto& force_i = a_force[a_level][*iter];
        res_i += force_i;
    }
}
