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
BoxOp<T, C_STATE, C_AUX, MEM>::BoxOp(
        const DisjointBoxLayout& a_layout,
        const LevelIndex& a_index, T a_dx)
{
    define(a_layout, a_index, a_dx);
}

template <typename T,
        unsigned int C_STATE,
        unsigned int C_AUX,
        MemType MEM>
BoxOp<T, C_STATE, C_AUX, MEM>::BoxOp(
        const DisjointBoxLayout& a_layout,
        const LevelIndex& a_index, Array<T,DIM> a_dx)
{
    define(a_layout, a_index, a_dx);
}


template <typename T,
        unsigned int C_STATE,
        unsigned int C_AUX,
        MemType MEM>
void
BoxOp<T, C_STATE, C_AUX, MEM>::define(
        const DisjointBoxLayout& a_layout, 
        const LevelIndex& a_index, T a_dx)
{
    Array<T, DIM> dx;
    dx.fill(a_dx);
    define(a_layout, a_index, dx);
}

template <typename T,
        unsigned int C_STATE,
        unsigned int C_AUX,
        MemType MEM>
void
BoxOp<T, C_STATE, C_AUX, MEM>::define(
        const DisjointBoxLayout& a_layout,
        const LevelIndex& a_index, Array<T, DIM> a_dx)
{
    m_scaleDiag = 0;
    m_scaleFlux = 1;
    m_dx = a_dx;
    m_layout = a_layout;
    m_index = a_index;
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

template <typename T,
        unsigned int C_STATE,
        unsigned int C_AUX,
        MemType MEM>
void BoxOp<T, C_STATE, C_AUX, MEM>::applyBC(
    Array<StateData, DIM>& a_fluxes,
    const StateData& a_state) const
{
    for (int dd = 0; dd < DIM; dd++)
    {
        if (this->layout().domain().isPeriodic(dd)) { continue; }
        this->applyBC(a_fluxes, a_state, Face(dd,Side::Lo));
        this->applyBC(a_fluxes, a_state, Face(dd,Side::Hi));
    }
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

#if 0  // MOVED TO MBBoxOp
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
    m_isMapped = true;
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

template <typename T, unsigned int C_STATE,
    unsigned int C_AUX, MemType MEM>
void BoxOp<T, C_STATE, C_AUX, MEM>::matchFlux(
            StateData&       a_rhs,
            const StateData& a_locFlux,
            const StateData& a_adjFlux,
            const StateData& a_state,
            Point            a_localDir) const
{
    matchFluxAverage(a_rhs, a_locFlux, a_adjFlux, a_localDir);
}

template <typename T, unsigned int C_STATE,
    unsigned int C_AUX, MemType MEM>
void BoxOp<T, C_STATE, C_AUX, MEM>::matchFluxAverage(
            StateData&       a_rhs,
            const StateData& a_locFlux,
            const StateData& a_adjFlux,
            Point            a_localDir) const
{
    PROTO_ASSERT(a_localDir.codim() == 1,
            "BoxOp::matchFluxAverage | Error: direction should be codimension 1");
    for (int dd = 0; dd < DIM; dd++)
    {
        if (a_localDir[dd] == 0) {continue; }
        bool side = a_localDir[dd] > 0;
        //Point growDirs = Point::Ones() - Point::Basis(dd,1);
        Box B0 = this->box();//.grow(growDirs*this->ghost());
        Box flxBox = side ? B0.adjacent(a_localDir, 1) : B0.edge(a_localDir, 1);
        Box slnBox = B0.edge(a_localDir, 1);

        PROTO_ASSERT((slnBox & a_rhs.box()).size() > 0,
                "BoxOp::matchFluxAverage | Error: RHS not defined on block boundary");

        Shift shift = side ? Shift::Basis(dd,1) : Shift::Zeros();
        T coef = side*0.5/dx()[dd];
        Stencil<T> S = coef*shift;
        StateData tmp(slnBox);
        tmp |= S(a_adjFlux);
        tmp -= S(a_locFlux);
        
        a_rhs += tmp;
    }
}
#endif
#endif
