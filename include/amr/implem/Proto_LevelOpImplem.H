
template <template<typename, MemType> class OpType, typename T, MemType MEM>
LevelOp<OpType, T, MEM>::LevelOp() {}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
LevelOp<OpType, T, MEM>::LevelOp(T a_dx)
{
    define(a_dx);
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
LevelOp<OpType, T, MEM>::LevelOp(std::array<T, DIM> a_dx)
{
    define(a_dx);
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::define(T a_dx)
{
    std::array<T, DIM> dx;
    dx.fill(a_dx);
    define(dx);
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::define(std::array<T, DIM> a_dx)
{
    m_op.define(a_dx);
    m_crseFluxRegister = nullptr;
    m_fineFluxRegister = nullptr;
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::operator()(
        LevelStateData& a_output,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux,
        T                     a_scale) const
{
    PROTO_ASSERT(a_output.layout().compatible(a_state.layout()),
            "LevelOp::operator() | Error: Output and state have incompatible layouts.");
    PROTO_ASSERT(a_output.layout().compatible(a_aux.layout()),
            "LevelOp::operator() | Error: Output and aux have incompatible layouts.");
    
    for (auto iter = a_state.begin(); iter.ok(); ++iter)
    {
        auto& out_i = a_output[*iter];
        const auto& state_i = a_state[*iter];
        const auto& aux_i   = a_aux[*iter];
        std::array<StateData, DIM> fluxes;
        for (int dir = 0; dir < DIM; dir++)
        {
            Box fluxBox = iter.box().grow(dir, Side::Hi, 1);
            fluxes[dir].define(fluxBox);
        }
        m_op(out_i, fluxes, state_i, aux_i, a_scale);
        updateFluxRegisters(fluxes, *iter);
    }
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::operator()(
        LevelStateData&       a_output,
        const LevelStateData& a_state,
        T                     a_scale) const
{
    PROTO_ASSERT(a_output.layout().compatible(a_state.layout()),
            "LevelOp::operator() | Error: Output and state have incompatible layouts.");
    
    for (auto iter = a_state.begin(); iter.ok(); ++iter)
    {
        auto& out_i = a_output[*iter];
        const auto& state_i = a_state[*iter];
        std::array<StateData, DIM> fluxes;
        for (int dir = 0; dir < DIM; dir++)
        {
            Box fluxBox = iter.box().grow(dir, Side::Hi, 1);
            fluxes[dir].define(fluxBox);
        }
        m_op(out_i, fluxes, state_i, a_scale);
        updateFluxRegisters(fluxes, *iter);
    }
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::updateFluxRegisters(
        std::array<StateData, DIM>& a_fluxes,
        DataIndex& a_index)
{
    if (m_fineFluxRegister != nullptr)
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            m_fineFluxRegister->incrementFine(a_fluxes[dir], a_index, 1.0, dir);
        }
    }
    if (m_crseFluxRegister != nullptr)
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            m_crseFluxRegister->incrementCoarse(a_fluxes[dir], a_index, 1.0, dir);
        }
    }
}

/*
template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::coarseReflux(
        LevelFluxRegister<T, OP::numState(), MEM>&    a_fluxRegister,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux,
        T a_scale) const
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
            a_fluxRegister.incrementCoarse(flux_id, *iter, a_scale, dir);
        }
    }
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::coarseReflux(
        LevelFluxRegister<T, OP::numState(), MEM>&    a_fluxRegister,
        const LevelStateData& a_state,
        T a_scale) const
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
            a_fluxRegister.incrementCoarse(flux_id, *iter, a_scale, dir);
        }
    }
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::fineReflux(
        LevelFluxRegister<T, OP::numState(), MEM>&    a_fluxRegister,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux,
        T a_scale) const
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
            a_fluxRegister.incrementFine(flux_id, *iter, a_scale, dir);
        }
    }
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::fineReflux(
        LevelFluxRegister<T, OP::numState(), MEM>&    a_fluxRegister,
        const LevelStateData& a_state,
        T a_scale) const
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
            a_fluxRegister.incrementFine(flux_id, *iter, a_scale, dir);
        }
    }
}
*/
template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::setDiagScale(T a_value)
{
    m_op.setDiagScale(a_value);
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::setFluxScale(T a_value)
{
    m_op.setFluxScale(a_value);
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::setFineFluxRegister(
        LevelFluxRegister<T, OP::numState(), MEM>& a_register)
{
    m_fineFluxRegister = &a_register;
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::setCrseFluxRegister(
        LevelFluxRegister<T, OP::numState(), MEM>& a_register)
{
    m_crseFluxRegister = &a_register;
}

template <template<typename, MemType> class OpType, typename T, MemType MEM>
void LevelOp<OpType, T, MEM>::clearFluxRegisters()
{
    m_fineFluxRegister = nullptr;
    m_crseFluxRegister = nullptr;
}


template <template<typename, MemType> class OpType, typename T, MemType MEM>
std::array<T, DIM> LevelOp<OpType, T, MEM>::dx() const
{
    return m_op.dx();
}
