
template <template<typename, MemType> class OPType, typename T, MemType MEM>
LevelOp<OPType, T, MEM>::LevelOp() {}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
LevelOp<OPType, T, MEM>::LevelOp(double a_dx)
{
    define(a_dx);
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void LevelOp<OPType, T, MEM>::define(double a_dx)
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
void LevelOp<OPType, T, MEM>::operator()(
        LevelStateData& a_output,
        const LevelStateData& a_state) const
{
    // TODO: check that all inputs have compatable layouts
    for (auto iter = a_state.begin(); iter.ok(); ++iter)
    {
        auto& out_i = a_output[*iter];
        const auto& state_i = a_state[*iter];
        m_op(out_i, state_i);
    }
}


template <template<typename, MemType> class OPType, typename T, MemType MEM>
void LevelOp<OPType, T, MEM>::coarseReflux(
        LevelFluxRegister<T, OP::numState(), MEM>&    a_fluxRegister,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux) const
{
    // TODO: check that all inputs have compatable layouts
    for (auto iter = a_state.begin(); iter.ok(); ++iter)
    {
        const auto& state_i = a_state[*iter];
        const auto& aux_i   = a_aux[*iter];
        for (int dir = 0; dir < DIM; dir++)
        {
            Box fluxBox = iter.box().grow(dir, Side::Hi, 1);
            BoxData<T, numState(), MEM> flux_id(fluxBox);
            m_op.flux(flux_id, state_i, aux_i, dir);
            a_fluxRegister.incrementCoarse(flux_id, *iter, 1.0, dir);
        }
    }
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void LevelOp<OPType, T, MEM>::coarseReflux(
        LevelFluxRegister<T, OP::numState(), MEM>&    a_fluxRegister,
        const LevelStateData& a_state) const
{
    // TODO: check that all inputs have compatable layouts
    for (auto iter = a_state.begin(); iter.ok(); ++iter)
    {
        const auto& state_i = a_state[*iter];
        for (int dir = 0; dir < DIM; dir++)
        {
            Box fluxBox = iter.box().grow(dir, Side::Hi, 1);
            BoxData<T, numState(), MEM> flux_id(fluxBox);
            m_op.flux(flux_id, state_i, dir);
            a_fluxRegister.incrementCoarse(flux_id, *iter, 1.0, dir);
        }
    }
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void LevelOp<OPType, T, MEM>::fineReflux(
        LevelFluxRegister<T, OP::numState(), MEM>&    a_fluxRegister,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux) const
{
    // TODO: check that all inputs have compatable layouts
    for (auto iter = a_state.begin(); iter.ok(); ++iter)
    {
        const auto& state_i = a_state[*iter];
        const auto& aux_i   = a_aux[*iter];
        for (int dir = 0; dir < DIM; dir++)
        {
            Box fluxBox = iter.box().grow(dir, Side::Hi, 1);
            BoxData<T, numState(), MEM> flux_id(fluxBox);
            m_op.flux(flux_id, state_i, aux_i, dir);
            a_fluxRegister.incrementFine(flux_id, *iter, 1.0, dir);
        }
    }
}

template <template<typename, MemType> class OPType, typename T, MemType MEM>
void LevelOp<OPType, T, MEM>::fineReflux(
        LevelFluxRegister<T, OP::numState(), MEM>&    a_fluxRegister,
        const LevelStateData& a_state) const
{
    // TODO: check that all inputs have compatable layouts
    for (auto iter = a_state.begin(); iter.ok(); ++iter)
    {
        const auto& state_i = a_state[*iter];
        for (int dir = 0; dir < DIM; dir++)
        {
            Box fluxBox = iter.box().grow(dir, Side::Hi, 1);
            BoxData<T, numState(), MEM> flux_id(fluxBox);
            m_op.flux(flux_id, state_i, dir);
            a_fluxRegister.incrementFine(flux_id, *iter, 1.0, dir);
        }
    }
}
/* TODO: Premature Optimization
template <template<typename, MemType> class OPType, typename T, MemType MEM>
void LevelOp<OPType, T, MEM>::refluxApply_CF(
        LevelStateData& a_output,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux,
        LevelFluxRegister&    a_frCoarse,
        LevelFluxRegister&    a_frFine) const
{
}
    
template <template<typename, MemType> class OPType, typename T, MemType MEM>
void LevelOp<OPType, T, MEM>::refluxApply_C(
        LevelStateData& a_output,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux,
        LevelFluxRegister&    a_frCoarse) const
{
}
    
template <template<typename, MemType> class OPType, typename T, MemType MEM>
void LevelOp<OPType, T, MEM>::refluxApply_F(
        LevelStateData& a_output,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux,
        LevelFluxRegister&    a_frFine) const
{
}
*/
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
