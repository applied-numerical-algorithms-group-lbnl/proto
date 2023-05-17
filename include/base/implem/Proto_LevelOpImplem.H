
template <template<typename, MemType> class OpType,
         typename T,
         template<typename, unsigned int, MemType, Centering> class BCType,
         MemType MEM>
LevelOp<OpType, T, BCType, MEM>::LevelOp() {}

template <template<typename, MemType> class OpType,
         typename T,
         template<typename, unsigned int, MemType, Centering> class BCType,
         MemType MEM>
LevelOp<OpType, T, BCType, MEM>::LevelOp(DisjointBoxLayout& a_layout, T a_dx)
{
    define(a_layout, a_dx);
}

template <template<typename, MemType> class OpType,
         typename T,
         template<typename, unsigned int, MemType, Centering> class BCType,
         MemType MEM>
LevelOp<OpType, T, BCType, MEM>::LevelOp(DisjointBoxLayout& a_layout, Array<T, DIM> a_dx)
{
    define(a_layout, a_dx);
}

template <template<typename, MemType> class OpType,
         typename T,
         template<typename, unsigned int, MemType, Centering> class BCType,
         MemType MEM>
void LevelOp<OpType, T, BCType, MEM>::define(DisjointBoxLayout& a_layout, T a_dx)
{
    Array<T, DIM> dx;
    dx.fill(a_dx);
    define(a_layout, dx);
}

template <template<typename, MemType> class OpType,
         typename T,
         template<typename, unsigned int, MemType, Centering> class BCType,
         MemType MEM>
void LevelOp<OpType, T, BCType, MEM>::define(DisjointBoxLayout& a_layout, Array<T, DIM> a_dx)
{
    m_layout = a_layout;
    m_ops.clear();
    m_ops.resize(a_layout.localSize());
    for (auto iter : a_layout)
    {
        m_ops[iter].define(a_layout[iter], a_dx);
    }
}

template <template<typename, MemType> class OpType,
         typename T,
         template<typename, unsigned int, MemType, Centering> class BCType,
         MemType MEM>
void LevelOp<OpType, T, BCType, MEM>::operator()(
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

template <template<typename, MemType> class OpType,
         typename T,
         template<typename, unsigned int, MemType, Centering> class BCType,
         MemType MEM>
void LevelOp<OpType, T, BCType, MEM>::operator()(
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

template <template<typename, MemType> class OpType,
         typename T,
         template<typename, unsigned int, MemType, Centering> class BCType,
         MemType MEM>
void LevelOp<OpType, T, BCType, MEM>::applyBC(LevelStateData& a_state) const
{
    m_bc.apply(a_state);
}

template <template<typename, MemType> class OpType,
         typename T,
         template<typename, unsigned int, MemType, Centering> class BCType,
         MemType MEM>
void LevelOp<OpType, T, BCType, MEM>::setDiagScale(T a_value)
{
    m_diagScale = a_value;
    for (auto iter : m_layout)
    {
        m_ops[iter].setDiagScale(a_value);
    }
}

template <template<typename, MemType> class OpType,
         typename T,
         template<typename, unsigned int, MemType, Centering> class BCType,
         MemType MEM>
void LevelOp<OpType, T, BCType, MEM>::setFluxScale(T a_value)
{
    m_fluxScale = a_value;
    for (auto iter : m_layout)
    {
        m_ops[iter].setFluxScale(a_value);
    }
}

template <template<typename, MemType> class OpType,
         typename T,
         template<typename, unsigned int, MemType, Centering> class BCType,
         MemType MEM>
void LevelOp<OpType, T, BCType, MEM>::setTime(T a_time)
{
    m_time = a_time;
    for (auto iter : m_layout)
    {
        m_ops[iter].setTime(a_time);
    }
}

template <template<typename, MemType> class OpType,
         typename T,
         template<typename, unsigned int, MemType, Centering> class BCType,
         MemType MEM>
void LevelOp<OpType, T, BCType, MEM>::setRKStage(unsigned int a_stage)
{
    m_rkStage = a_stage;
    for (auto iter : m_layout)
    {
        m_ops[iter].setRKStage(a_stage);
    }
}
