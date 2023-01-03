template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
LevelSolver_FASMultigrid<OpType, T, MEM>::LevelSolver_FASMultigrid(
    DisjointBoxLayout& a_layout,
    Point              a_refRatio,
    int                a_numLevels,
    Array<T, DIM> a_dx)
{
    define(a_layout, a_refRatio, a_numLevels, a_dx);
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
LevelSolver_FASMultigrid<OpType, T, MEM>::LevelSolver_FASMultigrid(
    DisjointBoxLayout& a_layout,
    Point              a_refRatio,
    int                a_numLevels,
    T                  a_dx)
{
    define(a_layout, a_refRatio, a_numLevels, a_dx);
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
LevelSolver_FASMultigrid<OpType, T, MEM>::define(
    DisjointBoxLayout& a_layout,
    Point              a_refRatio,
    int                a_numLevels,
    Array<T, DIM> a_dx)
{
    m_verbose = false;
    m_residual.define(a_layout, Point::Zeros());
    m_mg = std::make_shared<MGLevel>(a_layout, a_refRatio, a_dx, a_numLevels - 1);
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
LevelSolver_FASMultigrid<OpType, T, MEM>::define(
    DisjointBoxLayout& a_layout,
    Point              a_refRatio,
    int                a_numLevels,
    T                  a_dx)
{
    Array<T, DIM> dx;
    dx.fill(a_dx);
    define(a_layout, a_refRatio, a_numLevels, dx);
}

/*
template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
LevelSolver_FASMultigrid<OpType, T, MEM>::defineAsSubcycler(
    DisjointBoxLayout& a_layout,
    Point              a_refRatio,
    int                a_numLevels,
    T             a_dx)
{
    PROTO_ASSERT(a_refRatio.isotropic(),
        "LevelSolver_FASMultigrid::defineAsSubcycler | Error: Anisotropic refinement is not yet supported.");
    m_mg = std::make_shared<MGLevel>(a_layout, a_refRatio, a_dx, a_numLevels - 1);
    m_mg->defineAsSubcycler(a_layout, a_refRatio, a_dx, a_numLevels - 1);
}
*/

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
T
LevelSolver_FASMultigrid<OpType, T, MEM>::solve(
    LevelStateData& a_state,
    LevelStateData& a_force,
    LevelAuxData&   a_aux,
    int             a_maxIter,
    T          a_tolerance)
{
    // not implemented yet
    return 0;
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
T
LevelSolver_FASMultigrid<OpType, T, MEM>::solve(
    LevelStateData& a_state,
    LevelStateData& a_force,
    int             a_maxIter,
    T          a_tolerance)
{
#ifdef PR_HDF5
    HDF5Handler h5;
#endif
    using Proto::pout;

    T res0 = resnorm(m_residual, a_state, a_force);
    Array<T, DIM> dx = m_mg->op().dx();
    if (m_verbose)
    {
#ifdef PR_HDF5
        h5.writeLevel(dx, m_residual, "MG_RES_I%i", 0);
        h5.writeLevel(dx, a_state, "MG_STATE_I%i", 0);
        h5.writeLevel(dx, a_force, "MG_FORCE");
#endif
        pout() << "MG solve: initial residual = " << res0 << std::endl;
        if (procID() == 0)
        {
            std::cout << "MG solve: initial residual = " << res0 << std::endl;
        }
    }
    T res;
    for (int ii = 0; ii < a_maxIter; ii++)
    {
        m_mg->vCycle(a_state, a_force);
        res = resnorm(m_residual, a_state, a_force);
        if (m_verbose)
        {
#ifdef PR_HDF5
            h5.writeLevel(dx, m_residual, "MG_RES_I%i", ii+1);
            h5.writeLevel(dx, a_state, "MG_STATE_I%i", ii+1);
#endif
            pout() << "\tMG solve: resnorm = " << res;
            pout() << " | integral = " << m_residual.integrate(dx) << std::endl;
            if (procID() == 0)
            {
                std::cout << "\tMG solve: resnorm = " << res;
                std::cout << " | integral = " << m_residual.integrate(dx) << std::endl;
            }
        }
        if (res < a_tolerance*res0) { break; }
    }
    return res;
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
T LevelSolver_FASMultigrid<OpType, T, MEM>::resnorm(
    LevelStateData& a_residual,
    LevelStateData& a_state,
    LevelStateData& a_force)
{
    a_state.exchange();
    Reduction<T, Abs> rxn;
    for (auto iter = a_state.begin(); iter.ok(); ++iter)
    {
        auto& state_i = a_state[*iter];
        auto& force_i = a_force[*iter];
        auto& res = a_residual[*iter];
        m_mg->op().op()(res, state_i);
        res += force_i;
        res.reduce(rxn);
    }
    T localRes = rxn.fetch();
#ifdef PR_MPI
    T globalRes = 0;
    MPI_Allreduce(&localRes, &globalRes, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return globalRes;
#else
    return localRes;
#endif
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
void
LevelSolver_FASMultigrid<OpType, T, MEM>::setVerbose(
    bool a_flag)
{
    m_verbose = a_flag;
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
LevelSolver_FASMultigrid<OpType, T, MEM>::MGLevel::MGLevel(
    DisjointBoxLayout& a_layout,
    Point              a_refRatio,
    Array<T, DIM> a_dx,
    int                a_level)
{
    define(a_layout, a_refRatio, a_dx, a_level);
}

template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
LevelSolver_FASMultigrid<OpType, T, MEM>::MGLevel::define(
    DisjointBoxLayout& a_layout,
    Point              a_refRatio,
    Array<T, DIM> a_dx,
    int a_level)
{
    m_level = a_level;
    m_refRatio = a_refRatio;
    m_numPreRelax = 2*DIM;
    m_numPostRelax = 2*DIM;
    m_numBottomRelax = 10;

    m_levelOp.define(a_dx);
    m_levelOp.setFluxScale(-1);
    m_lambda = 0.5/m_levelOp.spectralRadius();

    m_average = Stencil<T>::AvgDown(a_refRatio);
    m_interpolate = InterpStencil<T>::PiecewiseConstant(a_refRatio);
    m_increment = (-m_lambda)*Shift::Zeros();

    m_residual.define(a_layout, Point::Zeros());
    if (a_level > 0)
    {
        DisjointBoxLayout crseLocalLayout = a_layout.coarsen(a_refRatio);
        ProblemDomain crseDomain = a_layout.domain().coarsen(a_refRatio);
        DisjointBoxLayout crseLayout;
        Point boxSize = a_layout.boxSize();
        if (crseDomain.sizes() % boxSize == Point::Zeros())
        {
            crseLayout.define(crseDomain, boxSize);
        } else {
            crseLayout.define(crseDomain, crseDomain.sizes());
        }
        m_crseLocal.define(crseLocalLayout, Point::Zeros());
        m_crseForce.define(crseLayout,      Point::Zeros());
        m_crseState.define(crseLayout,      m_levelOp.ghost());
        m_crseState_0.define(crseLayout,      Point::Zeros());
        if (BOP::numAux() > 0)
        {
            m_crseAux.define(crseLayout,    m_levelOp.auxGhost());
        }
        
        //FIXME: Assumes isotropic refinement
        Array<T, DIM> cdx;
        for (int dir = 0; dir < DIM; dir++) { cdx[dir] = a_dx[dir]*a_refRatio[dir]; }
        m_crseMG = std::make_shared<MGLevel>(crseLayout, a_refRatio, cdx, a_level - 1);
    }
}

/*
template<template<typename, MemType> class OpType,
    typename T, MemType MEM> 
void
LevelSolver_FASMultigrid<OpType, T, MEM>::MGLevel::defineAsSubcycler(
    DisjointBoxLayout& a_layout,
    Point              a_refRatio,
    T             a_dx,
    int                a_level)
{
    m_level = a_level;
    m_refRatio = a_refRatio;
    m_numPreRelax = 2*DIM;
    m_numPostRelax = 2*DIM;
    m_numBottomRelax = 10;

    m_levelOp.define(a_dx);
    m_levelOp.setFluxScale(-1);
    m_lambda = 0.5/m_levelOp.spectralRadius();

    m_average = Stencil<T>::AvgDown(a_refRatio);
    m_interpolate = InterpStencil<T>::PiecewiseConstant(a_refRatio);
    m_increment = (-m_lambda)*Shift::Zeros();
    
    m_residual.define(a_layout, Point::Zeros());
    if (a_level > 0)
    {
        DisjointBoxLayout crseLayout = a_layout.coarsen(a_refRatio);
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
        // FIXME: Assumes isotropic refinement
        int ref0 = a_refRatio[0];
        int refRatioBC = pow(ref0, a_level); // this looks wrong, but is correct
        m_interpBC = InterpStencil<T>::Build(BOP::order() + 1, a_refRatio);
        T cdx = ref0*a_dx;
        m_crseMG = std::make_shared<MGLevel>(crseLayout, a_refRatio, cdx, a_level - 1);
    }
}
*/

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
        auto res = m_levelOp.op()(state_i, a_state.layout()[*iter]);
        res += force_i;
        resC_i |= m_average(res);
    }
    m_crseLocal.copyTo(a_crseForce);
    for (auto iter = a_crseForce.begin(); iter.ok(); ++iter)
    {
        auto& state_i = a_crseState[*iter];
        auto& force_i = a_crseForce[*iter];
        auto lphi = m_crseMG->m_levelOp.op()(state_i, a_crseForce.layout()[*iter]);
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
        a_state.exchange();
        for (auto iter = a_state.begin(); iter.ok(); ++iter)
        {
            auto& state_i = a_state[*iter];
            auto& force_i = a_force[*iter];
            auto res = m_levelOp.op()(state_i, a_state.layout()[*iter]);
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
#ifdef PR_HDF5
    HDF5Handler h5;
#endif
    PR_TIMERS("LevelMultigrid::VCycle");
    if (m_level > 0)
    {
        relax(a_state, a_force, m_numPreRelax);
        //FIXME: assumes periodic boundaries
        averageDown(m_crseState, a_state, m_refRatio);
        m_crseState.copyTo(m_crseState_0);
        coarseForce(m_crseForce, a_state, a_force, m_crseState);
        m_crseMG->vCycle(m_crseState, m_crseForce);
        fineCorrect(a_state, m_crseState, m_crseState_0);
        relax(a_state, a_force, m_numPostRelax);
        
    } else {
        relax(a_state, a_force, m_numBottomRelax);
    }
}
