template <typename T,
        unsigned int C_STATE,
        unsigned int C_AUX,
        MemType MEM>
BoxOp<T, C_STATE, C_AUX, MEM>::BoxOp()
{
    m_scaleDiag = 0;
    m_scaleFlux = 1;
    m_dx = 1;

    m_definedFlux = true;
    m_definedSrce = true;
}

template <typename T,
        unsigned int C_STATE,
        unsigned int C_AUX,
        MemType MEM>
BoxOp<T, C_STATE, C_AUX, MEM>::BoxOp(double a_dx)
{
    m_scaleDiag = 0;
    m_scaleFlux = 1;
    m_dx = a_dx;
    
    m_definedFlux = true;
    m_definedSrce = true;
}

/// Apply (In Place)
template <typename T,
        unsigned int C_STATE,
        unsigned int C_AUX,
        MemType MEM>
void
BoxOp<T, C_STATE, C_AUX, MEM>::operator()(
        StateData& a_output,
        const StateData& a_state,
        const AuxData&   a_aux) const
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
            auto fluxBox = a_output.box().grow(dir, Side::Hi, 1);
            StateData flux_d(fluxBox);
            flux(flux_d, a_state, a_aux, dir);
            if (!m_definedFlux) {break; }
            divFlux += DIV(flux_d, 1.0/m_dx);
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

/// Set Diagonal Term Scaling
template <typename T,
    unsigned int C_STATE,
    unsigned int C_AUX,
    MemType MEM>
void BoxOp<T, C_STATE, C_AUX, MEM>::setDiagScale(double a_value)
{
    m_scaleDiag = a_value;
}

/// Set Flux Term Scaling
template <typename T,
    unsigned int C_STATE,
    unsigned int C_AUX,
    MemType MEM>
void BoxOp<T, C_STATE, C_AUX, MEM>::setFluxScale(double a_value)
{
    m_scaleFlux = a_value;
}

/// Get Grid Spacing
template <typename T,
    unsigned int C_STATE,
    unsigned int C_AUX,
    MemType MEM>
double BoxOp<T, C_STATE, C_AUX, MEM>::dx() const
{
    return m_dx;
}
