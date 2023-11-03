
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

