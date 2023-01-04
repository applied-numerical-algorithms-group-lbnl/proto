template<template<typename, MemType> class OPType,
    typename T, typename MAP, MemType MEM>
MBLevelOp<OPType, T, MAP, MEM>::MBLevelOp()
{
    m_defined = false;
}

template<template<typename, MemType> class OPType,
    typename T, typename MAP, MemType MEM>
MBLevelOp<OPType, T, MAP, MEM>::MBLevelOp(MBMap<MAP, MEM>& a_map)
{
    define(a_map);
}

template<template<typename, MemType> class OPType,
    typename T, typename MAP, MemType MEM>
void MBLevelOp<OPType, T, MAP, MEM>::define(MBMap<MAP, MEM>& a_map)
{
    auto& layout = a_map.map().layout();
    m_map = &a_map;
    m_ops.clear();
    for (auto index : layout)
    {
        auto& x = a_map.map()[index];
        auto& J = a_map.jacobian()[index];
        m_ops.push_back(BOP(x,J));
    }
    m_defined = true;
}

template<template<typename, MemType> class OPType,
    typename T, typename MAP, MemType MEM>
void MBLevelOp<OPType, T, MAP, MEM>::operator()(
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
    typename T, typename MAP, MemType MEM>
void MBLevelOp<OPType, T, MAP, MEM>::operator()(
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

