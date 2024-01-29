
template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::MBBoxOp(
        const MBLevelMap<MAP, MEM>& a_map,
        const MBIndex& a_index) : m_map(a_map)
{
    m_index = a_index;
    
    this->m_scaleDiag = 0;
    this->m_scaleFlux = 1;
    this->m_time = 0;
        
    this->m_definedFlux = true;
    this->m_definedSrce = true;
    init(); //user initialization
}

template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
void
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::operator()(
        StateData& a_output,
        const StateData& a_state,
        const AuxData&   a_aux,
        T                a_scale)  const
{
    Array<StateData, DIM> fluxes;
    for (int dir = 0; dir < DIM; dir++)
    {
        fluxes[dir].define(a_output.box().grow(dir, Side::Hi, 1));
    }
    (*this)(a_output, fluxes, a_state, a_aux, a_scale);
}

template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
void
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::operator()(
        StateData& a_output,
        const StateData& a_state,
        T                a_scale)  const
{
    Array<StateData, DIM> fluxes;
    for (int dir = 0; dir < DIM; dir++)
    {
        fluxes[dir].define(a_output.box().grow(dir, Side::Hi, 1));
    }
    (*this)(a_output, fluxes, a_state, a_scale);
}

template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
void
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::operator()(
        StateData& a_output,
        Array<StateData, DIM>& a_fluxes,
        const StateData& a_state,
        const AuxData&   a_aux,
        T                a_scale)  const
{
    //TODO: this can be modified later if the base class operation is not sufficient
    BoxOp<T,C_STATE,C_AUX,MEM>::operator()(a_output, a_fluxes, a_state, a_aux, a_scale);
}

template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
void
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::operator()(
        StateData& a_output,
        Array<StateData, DIM>& a_fluxes,
        const StateData& a_state,
        T                a_scale)  const
{
    //TODO: this can be modified later if the base class operation is not sufficient
    BoxOp<T,C_STATE,C_AUX,MEM>::operator()(a_output, a_fluxes, a_state, a_scale);
}

template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
typename MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::StateData
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::operator()(
        const StateData& a_state,
        const AuxData&   a_aux,
        Box              a_range,
        T                a_scale)  const
{
    StateData out(a_range);
    this->operator()(out, a_state, a_aux, a_scale);
    return out; 
}

template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
typename MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::StateData
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::operator()(
        const StateData& a_state,
        Box              a_range,
        T                a_scale)  const
{
    StateData out(a_range);
    this->operator()(out, a_state, a_scale);
    return out; 
}
#if 0
template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
void
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::flux(
        StateData& a_flux,
        const StateData& a_state,
        const AuxData&   a_aux, 
        int a_dir)  const
{
    // signals the operator that there is no user defined
    // function and that contributions should be ignored
    m_definedFlux = false;
}

template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
void
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::flux(
        StateData& a_flux,
        const StateData& a_state,
        int a_dir)  const
{
    // signals the operator that there is no user defined
    // function and that contributions should be ignored
    m_definedFlux = false;
}

template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
void
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::source(
        StateData&       a_source,
        const StateData& a_state,
        const AuxData&   a_aux)  const
{
    // signals the operator that there is no user defined
    // function and that contributions should be ignored
    m_definedSrce = false;
}

template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
void
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::source(
        StateData&       a_source,
        const StateData& a_state)  const
{
    // signals the operator that there is no user defined
    // function and that contributions should be ignored
    m_definedSrce = false;
}

template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
void
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::diag(
        StateData&       a_diag,
        const StateData& a_state,
        const AuxData&   a_aux)  const 
{
    a_state.copyTo(a_diag);
}

template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
void
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::diag(
        StateData&       a_diag,
        const StateData& a_state)  const
{
    a_state.copyTo(a_diag);
}
#endif

template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
const MBDisjointBoxLayout&
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::layout() const
{
    return m_map.layout();
}

template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
const MBLevelMap<MAP,MEM>&
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::map() const
{
    return m_map;
}

template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
MBIndex
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::index() const
{
    return m_index;
}

template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
unsigned int
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::block() const
{
    return layout().block(m_index);
}

template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
Box
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::box() const
{
    return layout()[m_index];
}

template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
Array<T, DIM>
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::dx() const
{
    return m_map.dx(block());
}

template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
const BoxData<double, DIM, MEM>&
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::x() const
{
    return m_map.map()[m_index];
}

template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM>
const BoxData<double, 1, MEM>&
MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>::jacobian() const
{
    return m_map.jacobian()[m_index];
}
