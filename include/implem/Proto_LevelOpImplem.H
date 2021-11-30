
template <template<typename, MemType> class OPType, typename T, MemType MEM>
LevelOp<OPType, T, MEM>::LevelOp() {}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
LevelOp<OPType, T, MEM>::LevelOp(double a_dx)
{
    m_op.define(a_dx);
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void LevelOp<OPType, T, MEM>::operator()(
        LevelStateData& a_output,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux) const
{
    // TODO: check that all inputs have compatable layouts
    for (auto iter = a_state.begin(); iter.ok(); ++iter)
    {
        auto& out_i = a_output[*iter];
        const auto& state_i = a_state[*iter];
        const auto& aux_i   = a_aux[*iter];
        m_op(out_i, state_i, aux_i);
    }
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void LevelOp<OPType, T, MEM>::setDiagScale(double a_value)
{
    m_op.setDiagScale(a_value);
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void LevelOp<OPType, T, MEM>::setFluxScale(double a_value)
{
    m_op.setFluxScale(a_value);
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
double LevelOp<OPType, T, MEM>::dx() const
{
    return m_op.dx();
}
