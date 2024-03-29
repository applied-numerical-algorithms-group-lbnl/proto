template <typename T,
        unsigned int C_STATE,
        unsigned int C_AUX,
        MemType MEM>
BoxOp<T, C_STATE, C_AUX, MEM>::BoxOp()
{
}

template <typename T,
        unsigned int C_STATE,
        unsigned int C_AUX,
        MemType MEM>
BoxOp<T, C_STATE, C_AUX, MEM>::BoxOp(const Box& a_box, T a_dx)
{
    define(a_box, a_dx);
}

template <typename T,
        unsigned int C_STATE,
        unsigned int C_AUX,
        MemType MEM>
BoxOp<T, C_STATE, C_AUX, MEM>::BoxOp(const Box& a_box, Array<T, DIM> a_dx)
{
    define(a_box, a_dx);
}


template <typename T,
        unsigned int C_STATE,
        unsigned int C_AUX,
        MemType MEM>
void
BoxOp<T, C_STATE, C_AUX, MEM>::define(const Box& a_box, T a_dx)
{
    Array<T, DIM> dx;
    dx.fill(a_dx);
    define(a_box, dx);
}

template <typename T,
        unsigned int C_STATE,
        unsigned int C_AUX,
        MemType MEM>
void
BoxOp<T, C_STATE, C_AUX, MEM>::define(const Box& a_box, Array<T, DIM> a_dx)
{
    m_scaleDiag = 0;
    m_scaleFlux = 1;
    m_dx = a_dx;
    m_box = a_box;
    m_time = 0;
        
    m_definedFlux = true;
    m_definedSrce = true;
    init(); //user initialization
}

template <typename T,
        unsigned int C_STATE,
        unsigned int C_AUX,
        MemType MEM>
void
BoxOp<T, C_STATE, C_AUX, MEM>::operator()(
        StateData& a_output,
        Array<StateData, DIM>& a_fluxes,
        const StateData& a_state,
        const AuxData&   a_aux,
        T                a_scale) const
{
    a_output.setVal(0);
    if (m_scaleFlux != 0)
    {
        StateData divFlux(a_output.box());
        for (int dir = 0; dir < DIM; dir++)
        {
            Stencil<T> DIV = 1.0*Shift::Basis(dir, 1)
                           - 1.0*Shift::Zeros();
            auto& flux_d = a_fluxes[dir];
            flux(flux_d, a_state, a_aux, dir);
            if (!m_definedFlux) {break; }
            if (dir == 0)
            {
                divFlux |= DIV(flux_d, 1.0/m_dx[dir]);
            } else {
                divFlux += DIV(flux_d, 1.0/m_dx[dir]);
            }
        }
        if (m_definedFlux)
        {
            a_output += divFlux;
        }
        StateData sourceData(a_output.box());
        source(sourceData, a_state, a_aux);
        if (m_definedSrce)
        {
            a_output += sourceData;
        }
        if (m_definedFlux || m_definedSrce)
        {
            a_output *= m_scaleFlux;
        }
    }
    if (m_scaleDiag != 0)
    {
        StateData diagData(a_output.box());
        diag(diagData, a_state, a_aux);
        diagData *= m_scaleDiag;
        a_output += diagData;
    }
    if (a_scale != 1.0)
    {
        a_output *= a_scale;
    }
}

template <typename T,
        unsigned int C_STATE,
        unsigned int C_AUX,
        MemType MEM>
void
BoxOp<T, C_STATE, C_AUX, MEM>::operator()(
        StateData& a_output,
        Array<StateData, DIM>& a_fluxes,
        const StateData& a_state,
        T                a_scale) const
{
    a_output.setVal(0);
    if (m_scaleFlux != 0)
    {
        StateData divFlux(a_output.box());
        divFlux.setVal(0);
        for (int dir = 0; dir < DIM; dir++)
        {
            Stencil<T> DIV = 1.0*Shift::Basis(dir, 1)
                           - 1.0*Shift::Zeros();
            auto& flux_d = a_fluxes[dir];
            flux(flux_d, a_state, dir);
            if (!m_definedFlux) {break; }
            divFlux += DIV(flux_d, 1.0/m_dx[dir]);
        }
        if (m_definedFlux)
        {
            a_output += divFlux;
        }
        StateData sourceData(a_output.box());
        source(sourceData, a_state);
        if (m_definedSrce)
        {
            a_output += sourceData;
        }
        if (m_definedFlux || m_definedSrce)
        {
            a_output *= m_scaleFlux;
        }
    }
    if (m_scaleDiag != 0)
    {
        StateData diagData(a_output.box());
        diag(diagData, a_state);
        diagData *= m_scaleDiag;
        a_output += diagData;
    }
    if (a_scale != 1.0)
    {
        a_output *= a_scale;
    }
}

template <typename T,
        unsigned int C_STATE,
        unsigned int C_AUX,
        MemType MEM>
void
BoxOp<T, C_STATE, C_AUX, MEM>::operator()(
        StateData& a_output,
        const StateData& a_state,
        const AuxData&   a_aux,
        T                a_scale) const
{
    Array<StateData, DIM> fluxes;
    for (int dir = 0; dir < DIM; dir++)
    {
        fluxes[dir].define(a_output.box().grow(dir, Side::Hi, 1));
    }
    (*this)(a_output, fluxes, a_state, a_aux, a_scale);
}

template <typename T,
        unsigned int C_STATE,
        unsigned int C_AUX,
        MemType MEM>
void
BoxOp<T, C_STATE, C_AUX, MEM>::operator()(
        StateData& a_output,
        const StateData& a_state,
        T                a_scale) const
{
    Array<StateData, DIM> fluxes;
    for (int dir = 0; dir < DIM; dir++)
    {
        fluxes[dir].define(a_output.box().grow(dir, Side::Hi, 1));
    }
    (*this)(a_output, fluxes, a_state, a_scale);
}

template <typename T,
        unsigned int C_STATE,
        unsigned int C_AUX,
        MemType MEM>
typename BoxOp<T, C_STATE, C_AUX, MEM>::StateData
BoxOp<T, C_STATE, C_AUX, MEM>::operator()(
        const StateData& a_state,
        Box              a_range,
        T                a_scale) const
{
    StateData out(a_range);
    this->operator()(out, a_state, a_scale);
    return out; 
}

template <typename T,
        unsigned int C_STATE,
        unsigned int C_AUX,
        MemType MEM>
typename BoxOp<T, C_STATE, C_AUX, MEM>::StateData
BoxOp<T, C_STATE, C_AUX, MEM>::operator()(
        const StateData& a_state,
        const AuxData&   a_aux,
        Box              a_range,
        T                a_scale) const
{
    StateData out(a_range);
    this->operator()(out, a_state, a_aux, a_scale);
    return out; 
}

/// User Defined Flux
template <typename T,
    unsigned int C_STATE,
    unsigned int C_AUX,
    MemType MEM>
void BoxOp<T, C_STATE, C_AUX, MEM>::flux(
        StateData& a_flux,
        const StateData& a_state,
        const AuxData&   a_aux, 
        int a_dir) const
{
    // signals the operator that there is no user defined
    // function and that contributions should be ignored
    m_definedFlux = false;
}

/// User Defined Flux
template <typename T,
    unsigned int C_STATE,
    unsigned int C_AUX,
    MemType MEM>
void BoxOp<T, C_STATE, C_AUX, MEM>::flux(
        StateData& a_flux,
        const StateData& a_state,
        int a_dir) const
{
    // signals the operator that there is no user defined
    // function and that contributions should be ignored
    m_definedFlux = false;
}

/// User Defined Source
template <typename T,
    unsigned int C_STATE,
    unsigned int C_AUX,
    MemType MEM>
void BoxOp<T, C_STATE, C_AUX, MEM>::source(
        StateData&       a_source,
        const StateData& a_state,
        const AuxData&   a_aux) const
{
    // signals the operator that there is no user defined
    // function and that contributions should be ignored
    m_definedSrce = false;
}

/// User Defined Source
template <typename T,
    unsigned int C_STATE,
    unsigned int C_AUX,
    MemType MEM>
void BoxOp<T, C_STATE, C_AUX, MEM>::source(
        StateData&       a_source,
        const StateData& a_state) const
{
    // signals the operator that there is no user defined
    // function and that contributions should be ignored
    m_definedSrce = false;
}

/// User Defined Diagonal
template <typename T,
    unsigned int C_STATE,
    unsigned int C_AUX,
    MemType MEM>
void BoxOp<T, C_STATE, C_AUX, MEM>::diag(
        StateData&       a_diag,
        const StateData& a_state,
        const AuxData&   a_aux) const 
{
    a_state.copyTo(a_diag);
}

/// User Defined Diagonal
template <typename T,
    unsigned int C_STATE,
    unsigned int C_AUX,
    MemType MEM>
void BoxOp<T, C_STATE, C_AUX, MEM>::diag(
        StateData&       a_diag,
        const StateData& a_state) const
{
    a_state.copyTo(a_diag);
}


/// Set Diagonal Term Scaling
template <typename T,
    unsigned int C_STATE,
    unsigned int C_AUX,
    MemType MEM>
void BoxOp<T, C_STATE, C_AUX, MEM>::setDiagScale(T a_value)
{
    m_scaleDiag = a_value;
}

/// Set Flux Term Scaling
template <typename T,
    unsigned int C_STATE,
    unsigned int C_AUX,
    MemType MEM>
void BoxOp<T, C_STATE, C_AUX, MEM>::setFluxScale(T a_value)
{
    m_scaleFlux = a_value;
}

/// Set Time
template <typename T,
    unsigned int C_STATE,
    unsigned int C_AUX,
    MemType MEM>
void BoxOp<T, C_STATE, C_AUX, MEM>::setTime(T a_time)
{
    m_time = a_time;
}

/// Set Runge Kutta Stage
template <typename T,
    unsigned int C_STATE,
    unsigned int C_AUX,
    MemType MEM>
void BoxOp<T, C_STATE, C_AUX, MEM>::setRKStage(unsigned int a_stage)
{
    m_RKStage = a_stage;
}

/// Get Grid Spacing
template <typename T,
    unsigned int C_STATE,
    unsigned int C_AUX,
    MemType MEM>
Array<T, DIM> BoxOp<T, C_STATE, C_AUX, MEM>::dx() const
{
    return m_dx;
}

/// Get Min Grid Spacing
template <typename T,
    unsigned int C_STATE,
    unsigned int C_AUX,
    MemType MEM>
T BoxOp<T, C_STATE, C_AUX, MEM>::dxMin() const
{
    T minValue = m_dx[0];
    for (int dir = 1; dir < DIM; dir++)
    {
        minValue = std::min(minValue, m_dx[dir]);
    }
    return minValue;
}

/// Get Grid Spacing
template <typename T,
    unsigned int C_STATE,
    unsigned int C_AUX,
    MemType MEM>
T BoxOp<T, C_STATE, C_AUX, MEM>::dxMax() const
{
    T maxValue = m_dx[0];
    for (int dir = 1; dir < DIM; dir++)
    {
        maxValue = std::max(maxValue, m_dx[dir]);
    }
    return maxValue;
}

#ifdef PR_MMB
template <typename T, unsigned int C_STATE,
    unsigned int C_AUX, MemType MEM>
BoxOp<T, C_STATE, C_AUX, MEM>::BoxOp(
    const Box& a_box,
    const Array<double, DIM>& a_dx,
    const BoxData<double, DIM, MEM>& a_x,
    const BoxData<double, 1,   MEM>& a_J)
{
    define(a_box, a_dx, a_x, a_J);
}

template <typename T, unsigned int C_STATE,
    unsigned int C_AUX, MemType MEM>
void BoxOp<T, C_STATE, C_AUX, MEM>::define(
    const Box& a_box,
    const Array<double, DIM>& a_dx,
    const BoxData<double, DIM, MEM>& a_x,
    const BoxData<double, 1,   MEM>& a_J)
{
    m_x = &a_x;
    m_J = &a_J;
    define(a_box, a_dx);
}
template <typename T, unsigned int C_STATE,
    unsigned int C_AUX, MemType MEM>
const BoxData<double, DIM, MEM>& BoxOp<T, C_STATE, C_AUX, MEM>::x() const
{
    return *m_x;
}
template <typename T, unsigned int C_STATE,
    unsigned int C_AUX, MemType MEM>
const BoxData<double, 1, MEM>& BoxOp<T, C_STATE, C_AUX, MEM>::jacobian() const
{
    return *m_J;
}
#endif
