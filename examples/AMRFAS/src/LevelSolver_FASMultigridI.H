template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
LevelSolver_FASMultigrid<OpType, T, MEM>::LevelSolver_FASMultigrid(
    DisjointBoxLayout& a_layout,
    int                a_numLevels,
    double             a_dx)
{
    define(a_layout, a_numLevels, a_dx);
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
LevelSolver_FASMultigrid<OpType, T, MEM>::define(
    DisjointBoxLayout& a_layout,
    int                a_numLevels,
    double             a_dx)
{
    m_residual.define(a_layout, Point::Zeros());
    m_mg = std::make_shared<MGLevel>(a_layout, a_dx, a_numLevels - 1);
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
LevelSolver_FASMultigrid<OpType, T, MEM>::defineAsSubcycler(
    DisjointBoxLayout& a_layout,
    int                a_numLevels,
    double             a_dx)
{
    m_mg = std::make_shared<MGLevel>(a_layout, a_dx, a_numLevels - 1);
    m_mg->defineAsSubcycler(a_layout, a_dx, a_numLevels - 1);
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
double
LevelSolver_FASMultigrid<OpType, T, MEM>::solve(
    LevelStateData& a_state,
    LevelStateData& a_force,
    LevelAuxData&   a_aux,
    int             a_maxIter,
    double          a_tolerance)
{
    // not implemented yet
    return 0;
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
double
LevelSolver_FASMultigrid<OpType, T, MEM>::solve(
    LevelStateData& a_state,
    LevelStateData& a_force,
    int             a_maxIter,
    double          a_tolerance)
{
    HDF5Handler h5;
    double res0 = resnorm(m_residual, a_state, a_force);
    //std::cout << "\tMG solve: initial residual = " << res0 << std::endl;
    double res;
    for (int ii = 0; ii < a_maxIter; ii++)
    {
        MG_SOLVE_ITER = ii;
        m_mg->vCycle(a_state, a_force);
        res = resnorm(m_residual, a_state, a_force);
        /*
        h5.writeLevel(m_mg->op().dx(), m_residual, "RES_AMRI%i_AMRL%i_AMRS%i_MGI%i", 
                AMR_SOLVE_ITER, AMR_LEVEL, AMR_RELAX_STAGE, MG_SOLVE_ITER); 
        */
        //std::cout << "\t\tMG solve: resnorm = " << res << std::endl;
        if (res < a_tolerance*res0) { break; }
    }
    return res;
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
double LevelSolver_FASMultigrid<OpType, T, MEM>::resnorm(
    LevelStateData& a_residual,
    LevelStateData& a_state,
    LevelStateData& a_force)
{
    a_state.exchange();
    Reduction<double, Abs> rxn;
    for (auto iter = a_state.begin(); iter.ok(); ++iter)
    {
        auto& state_i = a_state[*iter];
        auto& force_i = a_force[*iter];
        auto& res = a_residual[*iter];
        //typename BOP::StateData res(iter.box());
        //auto res = m_mg->op().op()(state_i, iter.box());
        m_mg->op().op()(res, state_i);
        res += force_i;
        res.reduce(rxn);
    }
    // TODO: MPI communication
    return rxn.fetch();
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void LevelSolver_FASMultigrid<OpType, T, MEM>::interpBC(
    LevelStateData& a_crseStateBC)
{
    m_mg->interpBC(a_crseStateBC);
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
LevelSolver_FASMultigrid<OpType, T, MEM>::MGLevel::MGLevel(
    DisjointBoxLayout& a_layout,
    double             a_dx,
    int                a_level)
{
    define(a_layout, a_dx, a_level);
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
LevelSolver_FASMultigrid<OpType, T, MEM>::MGLevel::define(
    DisjointBoxLayout& a_layout,
    double a_dx,
    int a_level)
{
    m_level = a_level;
    m_numPreRelax = 2*DIM;
    m_numPostRelax = 2*DIM;
    m_numBottomRelax = 10;

    m_levelOp.define(a_dx);
    m_levelOp.setFluxScale(-1);
    m_lambda = 0.5/m_levelOp.spectralRadius();
    Point refRatio = Point::Ones(PR_MG_REFRATIO);

    m_average = Stencil<T>::AvgDown(refRatio);
    m_interpolate = InterpStencil<T>::PiecewiseConstant(refRatio);
    m_increment = (-m_lambda)*Shift::Zeros();

    m_residual.define(a_layout, Point::Zeros());
    if (a_level > 0)
    {
        DisjointBoxLayout crseLocalLayout = a_layout.coarsen(refRatio);
        ProblemDomain crseDomain = a_layout.domain().coarsen(refRatio);
        DisjointBoxLayout crseLayout;
        Point boxSize = a_layout.boxSize();
        if (crseDomain.size() % boxSize == Point::Zeros())
        {
            crseLayout.define(crseDomain, boxSize);
        } else {
            crseLayout.define(crseDomain, crseDomain.size());
        }
        m_crseLocal.define(crseLocalLayout, Point::Zeros());
        m_crseForce.define(crseLayout,      Point::Zeros());
        m_crseState.define(crseLayout,      m_levelOp.ghost());
        m_crseState_0.define(crseLayout,      Point::Zeros());
        if (BOP::numAux() > 0)
        {
            m_crseAux.define(crseLayout,    m_levelOp.auxGhost());
        }

        m_crseMG = std::make_shared<MGLevel>(crseLayout,
            a_dx*PR_MG_REFRATIO, a_level - 1);
    }
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
LevelSolver_FASMultigrid<OpType, T, MEM>::MGLevel::defineAsSubcycler(
    DisjointBoxLayout& a_layout,
    double             a_dx,
    int                a_level)
{
    m_level = a_level;
    m_numPreRelax = 2*DIM;
    m_numPostRelax = 2*DIM;
    m_numBottomRelax = 10;

    m_levelOp.define(a_dx);
    m_levelOp.setFluxScale(-1);
    m_lambda = 0.5/m_levelOp.spectralRadius();
    Point refRatio = Point::Ones(PR_MG_REFRATIO);

    m_average = Stencil<T>::AvgDown(refRatio);
    m_interpolate = InterpStencil<T>::PiecewiseConstant(refRatio);
    m_increment = (-m_lambda)*Shift::Zeros();
    
    m_residual.define(a_layout, Point::Zeros());
    if (a_level > 0)
    {
        DisjointBoxLayout crseLayout = a_layout.coarsen(refRatio);
        // TODO: m_crseLocal isn't actually needed in this version
        m_crseLocal.define(crseLayout,   Point::Zeros());
        m_crseForce.define(crseLayout,   Point::Zeros());
        m_crseState.define(crseLayout,   m_levelOp.ghost());
        m_crseState_0.define(crseLayout, Point::Zeros());
        if (BOP::numAux() > 0)
        {
            m_crseAux.define(crseLayout,    m_levelOp.auxGhost());
        }

        // Boundary conditions
        int refRatioBC = pow(PR_MG_REFRATIO, a_level); // this looks wrong, but is correct
        m_interpBC = InterpStencil<T>::Build(4, Box::Kernel(2), 4, refRatioBC);
        m_crseMG = std::make_shared<MGLevel>(crseLayout,
            a_dx*PR_MG_REFRATIO, a_level - 1);
    }
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
LevelSolver_FASMultigrid<OpType, T, MEM>::MGLevel::interpBC(
    LevelStateData& a_crseStateBC)
{
    if (m_level > 0)
    {
        interpBoundaries(a_crseStateBC, m_crseState, m_interpBC);
    }
    if (m_level > 1)
    {
        m_crseMG->interpBC(a_crseStateBC);
    }
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
LevelSolver_FASMultigrid<OpType, T, MEM>::MGLevel::coarseForce(
        LevelStateData& a_crseForce,
        LevelStateData& a_state,
        LevelStateData& a_force,
        LevelStateData& a_crseState)
{
    PR_TIMERS("LevelMultigrid::coarseResidual");
    a_state.exchange();
    a_crseState.exchange();
    for (auto iter = a_state.begin(); iter.ok(); ++iter)
    {
        auto& state_i = a_state[*iter];
        auto& resC_i  = m_crseLocal[*iter];
        auto& force_i = a_force[*iter];
        //typename BOP::StateData res(iter.box());
        //m_levelOp.op()(res, state_i);
        auto res = m_levelOp.op()(state_i, iter.box());
        res += force_i;
        resC_i |= m_average(res);
    }
    m_crseLocal.copyTo(a_crseForce);
    for (auto iter = a_crseForce.begin(); iter.ok(); ++iter)
    {
        auto& state_i = a_crseState[*iter];
        auto& force_i = a_crseForce[*iter];
        //typename BOP::StateData lphi(iter.box());
        //m_crseMG->m_levelOp.op()(lphi, state_i);
        auto lphi = m_crseMG->m_levelOp.op()(state_i, iter.box());
        force_i -= lphi;
    }
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
LevelSolver_FASMultigrid<OpType, T, MEM>::MGLevel::residual(
        LevelStateData& a_residual,
        LevelStateData& a_state,
        LevelStateData& a_force)
{
    a_state.exchange();
    for (auto iter = a_residual.begin(); iter.ok(); ++iter)
    {
        auto& res_i   = a_residual[*iter];
        auto& state_i = a_state[*iter];
        auto& force_i = a_force[*iter];

        m_levelOp.op()(res_i, state_i);
        res_i += force_i;
    }
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
LevelSolver_FASMultigrid<OpType, T, MEM>::MGLevel::relax(
        LevelStateData& a_state,
        LevelStateData& a_force,
        int a_numIter)
{
    PR_TIMERS("LevelMultigrid::relax");
    for (int ii = 0; ii < a_numIter; ii++)
    {
        RELAX_ITER = ii;
        a_state.exchange();
        for (auto iter = a_state.begin(); iter.ok(); ++iter)
        {
            double dx = m_levelOp.dx();
            auto& state_i = a_state[*iter];
            auto& force_i = a_force[*iter];
            //typename BOP::StateData res(iter.box());
            //m_levelOp.op()(res, state_i);
            auto res = m_levelOp.op()(state_i, iter.box());
            res += force_i;
            state_i += m_increment(res);
        }
    }
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
LevelSolver_FASMultigrid<OpType, T, MEM>::MGLevel::fineCorrect(
        LevelStateData& a_fine,
        LevelStateData& a_crse,
        LevelStateData& a_crse0)
{
    for (auto iter = a_crse.begin(); iter.ok(); ++iter)
    {

        auto& crse_i = a_crse[*iter];
        auto& crse0_i = a_crse0[*iter];

        crse0_i -= crse_i;
        crse0_i *= -1;
    }
    a_crse0.copyTo(m_crseLocal);
    for (auto iter = a_fine.begin(); iter.ok(); ++iter)
    {
        auto& fine_i = a_fine[*iter];
        auto& crse_i = m_crseLocal[*iter];

        fine_i += m_interpolate(crse_i);
    }
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
LevelSolver_FASMultigrid<OpType, T, MEM>::MGLevel::vCycle(
        LevelStateData& a_state,
        LevelStateData& a_force)
{
    HDF5Handler h5;
    MG_LEVEL = m_level;
    PR_TIMERS("LevelMultigrid::VCycle");
    if (m_level > 0)
    {
        relax(a_state, a_force, m_numPreRelax);
        //FIXME: assumes periodic boundaries
        a_state.coarsenTo(m_crseState, Point::Ones(PR_MG_REFRATIO));
        m_crseState.copyTo(m_crseState_0);
        coarseForce(m_crseForce, a_state, a_force, m_crseState);
        m_crseMG->vCycle(m_crseState, m_crseForce);
        MG_LEVEL = m_level;
        fineCorrect(a_state, m_crseState, m_crseState_0);
        relax(a_state, a_force, m_numPostRelax);
        
    } else {
        relax(a_state, a_force, m_numBottomRelax);
    }
    /*
    residual(m_residual, a_state, a_force);
    h5.writeLevel(m_levelOp.dx(), m_residual, "RES_AMRI%i_AMRL%i_AMRS%i_MGI%i_MGL%i", 
        AMR_SOLVE_ITER, AMR_LEVEL, AMR_RELAX_STAGE, MG_SOLVE_ITER, MG_LEVEL); 
    */
}
