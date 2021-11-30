#pragma once
#ifndef _PROTO_BOX_OP_
#define _PROTO_BOX_OP_

#include "Proto_BoxData.H"
namespace Proto {

/// Abstract Box-Scope Operator
/**
    PatchOp is the class from which all operators on AMR
    hierarchies are derived. Representable operators are of the form:

    L(phi, rho) = d*D(phi, rho)*I + f*[Div(F(phi, rho)) + S(phi, rho)]
    L:      Full Operator
    phi:    State Variables
    rho:    Auxiliary Variables
    D:      Diagonal Component of L
    F:      Flux Component of L
    S:      Source Component of L
    I:      Identity
    d:      Diagonal Term Scaling
    f:      Flux and Source Term Scaling
*/
template <typename T, unsigned int C_STATE, unsigned int C_AUX, MemType MEM = MEMTYPE_DEFAULT>
class BoxOp
{
    public:
    typedef BoxData<T, C_STATE, MEM> StateData;
    typedef BoxData<T, C_AUX,   MEM> AuxData;

    inline static constexpr unsigned int numState() { return C_STATE; }
    inline static constexpr unsigned int numAux()   { return C_AUX; }

    inline BoxOp();
    inline BoxOp(double a_dx);
    inline void define(double a_dx);
    /// Apply (In Place)
    inline void
    operator()(
        StateData& a_output,
        const StateData& a_state,
        const AuxData&   a_aux) const;
    
    /// Apply (Out of Place, possibly unnecessary)
    inline StateData
    operator()(
        const StateData& a_state,
        const AuxData&   a_aux) const;

    /// User Defined Flux
    inline virtual void flux(
        StateData& a_flux,
        const StateData& a_state,
        const AuxData&   a_aux, 
        int a_dir) const;
    
    /// User Defined Source
    inline virtual void source(
        StateData&       a_source,
        const StateData& a_state,
        const AuxData&   a_aux) const; 
    
    /// User Defined Diagonal
    inline virtual void diag(
        StateData&       a_diag,
        const StateData& a_state,
        const AuxData&   a_aux) const; 
    
    /// Set Diagonal Term Scaling
    inline void setDiagScale(double a_value);
    
    /// Set Flux Term Scaling
    inline void setFluxScale(double a_value);
    
    /// Get Grid Spacing
    inline double dx() const;

    protected:

    private:

    double m_scaleDiag;
    double m_scaleFlux;

    mutable bool m_definedSrce;
    mutable bool m_definedFlux;

    double m_dx;
};

#include "implem/Proto_BoxOpImplem.H"
} // end namespace Proto

#endif // end include guard
