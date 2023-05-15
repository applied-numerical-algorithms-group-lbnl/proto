
template <template<typename, MemType> class OpType, typename T, MemType MEM>
LevelOp<OpType, T, MEM>::LevelOp() {}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
LevelOp<OpType, T, MEM>::LevelOp(DisjointBoxLayout& a_layout, T a_dx)
{
    define(a_layout, a_dx);
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
LevelOp<OpType, T, MEM>::LevelOp(DisjointBoxLayout& a_layout, Array<T, DIM> a_dx)
{
    define(a_layout, a_dx);
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::define(DisjointBoxLayout& a_layout, T a_dx)
{
    Array<T, DIM> dx;
    dx.fill(a_dx);
    define(a_layout, dx);
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::define(DisjointBoxLayout& a_layout, Array<T, DIM> a_dx)
{
    m_layout = a_layout;
    m_ops.clear();
    for (auto iter : a_layout)
    {
        m_ops.push_back(OpType(a_layout[iter], a_dx));
    }
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::operator()(
        LevelStateData& a_output,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux,
        T                     a_scale) const
{
    PROTO_ASSERT(a_output.layout().compatible(m_layout),
            "LevelOp::operator() | Error: Output has an incompatible layout.");
    PROTO_ASSERT(a_state.layout().compatible(m_layout),
            "LevelOp::operator() | Error: State has an incompatible layout.");
    PROTO_ASSERT(a_aux.layout().compatible(m_layout),
            "LevelOp::operator() | Error: Aux has an incompatible layout.");
     
    for (auto iter : m_layout)
    {
        auto& out_i = a_output[iter];
        const auto& state_i = a_state[iter];
        const auto& aux_i   = a_aux[iter];
        m_ops[iter](out_i, state_i, aux_i, a_scale);
    }
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::operator()(
        LevelStateData&       a_output,
        const LevelStateData& a_state,
        T                     a_scale) const
{
    PROTO_ASSERT(a_output.layout().compatible(m_layout),
            "LevelOp::operator() | Error: Output has an incompatible layout.");
    PROTO_ASSERT(a_state.layout().compatible(m_layout),
            "LevelOp::operator() | Error: State has an incompatible layout.");
    
    for (auto iter : m_layout)
    {
        auto& out_i = a_output[iter];
        const auto& state_i = a_state[iter];
        Array<StateData, DIM> fluxes;
        for (int dir = 0; dir < DIM; dir++)
        {
            fluxes[dir].define(out_i.box().grow(dir, Side::Hi, 1));
        }
        m_ops[iter](out_i, fluxes, state_i, a_scale);
    }
}
template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::setDiagScale(T a_value)
{
    m_diagScale = a_value;
    for (auto iter : m_layout)
    {
        m_ops[iter].setDiagScale(a_value);
    }
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::setFluxScale(T a_value)
{
    m_fluxScale = a_value;
    for (auto iter : m_layout)
    {
        m_ops[iter].setFluxScale(a_value);
    }
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::setTime(T a_time)
{
    m_time = a_time;
    for (auto iter : m_layout)
    {
        m_ops[iter].setTime(a_time);
    }
}
