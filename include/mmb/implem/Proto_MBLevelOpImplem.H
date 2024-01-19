
template<template<typename, MemType> class OPType,
    template<MemType> typename MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
MBLevelOp<OPType,MAP,T,BCType,MEM>::MBLevelOp()
{
}

template<template<typename, MemType> class OPType,
    template<MemType> typename MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
MBLevelOp<OPType,MAP,T,BCType,MEM>::MBLevelOp(const MBLevelMap<MAP,MEM>& a_map)
{
    define(a_map);
}

template<template<typename, MemType> class OPType,
    template<MemType> typename MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBLevelOp<OPType,MAP,T,BCType,MEM>::define(const MBLevelMap<MAP,MEM>& a_map)
{
    PR_TIME("MBLevelOp::define");
    auto& layout = a_map.map().layout();
    m_map = &a_map;
    m_ops.clear();
    m_ops.resize(layout.localSize());
    m_boundRegister.define(layout,1,Point::Zeros());
    for (auto index : layout)
    {
        auto& x = a_map.map()[index];
        auto& J = a_map.jacobian()[index];
        auto block = layout.block(index);
        m_ops[index].define(layout[index], a_map.dx(block), x,J);
    }
}

template<template<typename, MemType> class OPType,
    template<MemType> typename MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBLevelOp<OPType,MAP,T,BCType,MEM>::operator()(
        LevelStateData&         a_output,
        const LevelStateData&   a_state,
        T                       a_scale)
{
    PR_TIME("MBLevelOp::operator()");
    auto& layout = m_map->map().layout();
    PROTO_ASSERT(a_output.compatible(layout),
            "MBLevelOp::operator() | Error: output layout is incompatible with map layout.");
    PROTO_ASSERT(a_state.compatible(layout),
            "MBLevelOp::operator() | Error: state layout is incompatible with map layout.");
    for (auto index : layout)
    {
        auto& op = m_ops[index];
        auto& state = a_state[index];
        auto& output = a_output[index];
        op(output, state, a_scale);
    }
    //matchFlux(a_output, a_state);
}

template<template<typename, MemType> class OPType,
    template<MemType> typename MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBLevelOp<OPType,MAP,T,BCType,MEM>::operator()(
        LevelStateData&         a_output,
        const LevelStateData&   a_state,
        const LevelAuxData&     a_aux,
        T                       a_scale)
{
    PR_TIME("MBLevelOp::operator()");
    auto& layout = m_map->map().layout();
    PROTO_ASSERT(a_output.compatible(layout),
            "MBLevelOp::operator() | Error: output layout is incompatible with map layout.");
    PROTO_ASSERT(a_state.compatible(layout),
            "MBLevelOp::operator() | Error: state layout is incompatible with map layout.");
    PROTO_ASSERT(a_aux.compatible(layout),
            "MBLevelOp::operator() | Error: aux data layout is incompatible with map layout.");
    for (auto index : layout)
    {
        auto& op = m_ops[index];
        auto& state = a_state[index];
        auto& output = a_output[index];
        auto& aux = a_aux[index];
        op(output, state, aux, a_scale);
    }
}

template<template<typename, MemType> class OPType,
    template<MemType> typename MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
OPType<T,MEM>& MBLevelOp<OPType,MAP,T,BCType,MEM>::operator[](MBIndex a_index)
{
    return m_ops[a_index];
}

template<template<typename, MemType> class OPType,
    template<MemType> typename MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
const OPType<T,MEM>& MBLevelOp<OPType,MAP,T,BCType,MEM>::operator[](MBIndex a_index) const
{
    return m_ops[a_index];
}

template<template<typename, MemType> class OPType,
    template<MemType> typename MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBLevelOp<OPType,MAP,T,BCType,MEM>::setDiagScale(T a_value)
{
    m_diagScale = a_value;
    for (auto iter : m_map->layout())
    {
        m_ops[iter].setDiagScale(a_value);
    }
}

template<template<typename, MemType> class OPType,
    template<MemType> typename MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBLevelOp<OPType,MAP,T,BCType,MEM>::setFluxScale(T a_value)
{
    m_fluxScale = a_value;
    for (auto iter : m_map->layout())
    {
        m_ops[iter].setFluxScale(a_value);
    }
}

template<template<typename, MemType> class OPType,
    template<MemType> typename MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBLevelOp<OPType,MAP,T,BCType,MEM>::setTime(T a_time)
{
    m_time = a_time;
    for (auto iter : m_map->layout())
    {
        m_ops[iter].setTime(a_time);
    }
}

template<template<typename, MemType> class OPType,
    template<MemType> typename MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBLevelOp<OPType,MAP,T,BCType,MEM>::setRKStage(unsigned int a_stage)
{
    m_rkStage = a_stage;
    for (auto iter : m_map->layout())
    {
        m_ops[iter].setRKStage(a_stage);
    }
}

template<template<typename, MemType> class OPType,
    template<MemType> typename MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType,
    MemType MEM>
void MBLevelOp<OPType,MAP,T,BCType,MEM>::matchFlux(
        LevelStateData& a_rhs,
        const LevelStateData& a_state)
{
    auto& layout = m_map->map().layout();
    PROTO_ASSERT(a_state.layout().compatible(layout),
            "MBLevelOp::matchFlux | Error: incompatible layouts");

    pr_out() << "MBLevelOp::matchFlux" << std::endl;
    // compute fluxes and copy them into the register
    for (auto iter : layout)
    {
        auto& op = (*this)[iter];
        const auto& state = a_state[iter];
        for (int dd = 0; dd < DIM; dd++)
        {
            SideIterator siter;
            for (siter.begin(); siter.ok(); ++siter)
            {
                Point dir = Point::Basis(dd, *siter);
                Box B0 = layout[iter];
                Box B1 = B0.grow((Centering)dd);
                StateData fluxData(B1); //conservative estimate
                // FIXME: compute local flux from adjacent coordinate system POV
                op.flux(fluxData, state, dd);
                pr_out() << "computing flux | base box: " << B0 << " | flux box: " << B1 << " | state box: " << state.box() << std::endl;
                if (layout.isBlockBoundary(iter, dir))
                {
                    auto& localFlux = m_boundRegister.local(iter, dd, *siter);
                    fluxData.copyTo(localFlux);
                    pr_out() << "copying flux from " << fluxData.box() << " into register " << localFlux.box() << std::endl;
                }
            }
        }
    }
    // fill adjacent registers
    m_boundRegister.exchange();

    for (auto iter : layout)
    {
        auto& op = (*this)[iter];
        auto& rhs = a_rhs[iter];
        const auto& state = a_state[iter];
        for (int dd = 0; dd < DIM; dd++)
        {
            SideIterator siter;
            for (siter.begin(); siter.ok(); ++siter)
            {
                Point dir = Point::Basis(dd, *siter);
                if (layout.isBlockBoundary(iter, dir))
                {
                    // FIXME: local register will have been computed from adjacent POV
                    // Need to recompute local flux and input to matchFlux function
                    auto& adjFlux = m_boundRegister.adjacent(iter, dd, *siter);
                    auto& locFlux = m_boundRegister.local(iter, dd, *siter);
                    pr_out() << "matching adj flux at " << adjFlux.box() << " with local flux at  " << locFlux.box() << std::endl;
                    op.matchFlux(rhs, locFlux, adjFlux, state, dir);
                }
            }
        }
    }
}
