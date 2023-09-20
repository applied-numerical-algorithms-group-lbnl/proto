template<template<typename, MemType> class OPType,
    template<MemType> typename MAP,
    typename T,  MemType MEM>
MBLevelOp<OPType, MAP, T,  MEM>::MBLevelOp()
{
}

template<template<typename, MemType> class OPType,
    template<MemType> typename MAP,
    typename T, MemType MEM>
MBLevelOp<OPType, MAP, T,  MEM>::MBLevelOp(MBLevelMap<MAP,MEM>& a_map)
{
    define(a_map);
}

template<template<typename, MemType> class OPType,
    template<MemType> typename MAP,
    typename T,  MemType MEM>
void MBLevelOp<OPType, MAP, T,  MEM>::define(MBLevelMap<MAP,MEM>& a_map)
{
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
    typename T,  MemType MEM>
void MBLevelOp<OPType, MAP, T,  MEM>::operator()(
        LevelStateData&         a_output,
        const LevelStateData&   a_state,
        T                       a_scale)
{
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
    typename T,  MemType MEM>
void MBLevelOp<OPType, MAP, T,  MEM>::operator()(
        LevelStateData&         a_output,
        const LevelStateData&   a_state,
        const LevelAuxData&     a_aux,
        T                       a_scale)
{
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

template <template<typename, MemType> class OpType,
    template<MemType> typename MAP,
         typename T,
         MemType MEM>
void MBLevelOp<OpType,MAP,T,MEM>::setDiagScale(T a_value)
{
    m_diagScale = a_value;
    for (auto iter : m_map->layout())
    {
        m_ops[iter].setDiagScale(a_value);
    }
}

template <template<typename, MemType> class OpType,
    template<MemType> typename MAP,
         typename T,
         MemType MEM>
void MBLevelOp<OpType,MAP,T,MEM>::setFluxScale(T a_value)
{
    m_fluxScale = a_value;
    for (auto iter : m_map->layout())
    {
        m_ops[iter].setFluxScale(a_value);
    }
}

template <template<typename, MemType> class OpType,
    template<MemType> typename MAP,
         typename T,
         MemType MEM>
void MBLevelOp<OpType,MAP,T,MEM>::setTime(T a_time)
{
    m_time = a_time;
    for (auto iter : m_map->layout())
    {
        m_ops[iter].setTime(a_time);
    }
}

template <template<typename, MemType> class OpType,
    template<MemType> typename MAP,
         typename T,
         MemType MEM>
void MBLevelOp<OpType,MAP,T,MEM>::setRKStage(unsigned int a_stage)
{
    m_rkStage = a_stage;
    for (auto iter : m_map->layout())
    {
        m_ops[iter].setRKStage(a_stage);
    }
}

