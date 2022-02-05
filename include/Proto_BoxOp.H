#pragma once
#ifndef _PROTO_BOX_OP_
#define _PROTO_BOX_OP_

#include "Proto_BoxData.H"
namespace Proto {

/// Abstract Box-Scope Operator
/**
    BoxOp is the class from which all operators on AMR
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

    TODO: potentially template the op on the data-holder itself (e.g. BoxData) for
    additional flexibility.  

    \tparam T       Datatype of the data holder (e.g. int, double, etc.)
    \tparam C_STATE Number of state variables
    \tparam C_AUX   Number of auxiliary (non-state) variables
    \tparam MEM     Proto::MemType of the data holder
*/
template <typename T, unsigned int C_STATE, unsigned int C_AUX, MemType MEM = MEMTYPE_DEFAULT>
class BoxOp
{
    public:
    typedef BoxData<T, C_STATE, MEM> StateData;
    typedef BoxData<T, C_AUX,   MEM> AuxData;

    /// Get Number of State Variables
    /**
        Allows C_STATE to be known after template parameter binding
    */
    inline static constexpr unsigned int numState() { return C_STATE; }

    /// Get Number of Auxiliary Variables
    /**
        Allows C_AUX to be known after template parameter binding
    */
    inline static constexpr unsigned int numAux()   { return C_AUX; }

    //TODO: In addition to the following, the interface requires these
    // static functions. There is no good way to specify necessary static
    // functions in an interface which is extremely awkward...
    //
    // static Point ghost() const;      (Needed ghost data in the state variables)
    // static Point auxGhost() const;   (Needed ghost data in the aux variables)

    /// Default Constructor
    inline BoxOp();

    /// Non-Trivial Constructor
    /**
        \param dx   Grid spacing
    */
    inline BoxOp(double a_dx);

    /// Lazy Constructor
    /**
        \param dx   Grid spacing
    */
    inline void define(double a_dx);

    /// Spectral Radius (User Defined)
    /**
        Approximate spectral radius of the operator. Not required, but 
        Needed for iterative solvers. 
    */
    inline virtual double spectralRadius() const { return 1; }

    /// Apply (In Place)
    /**
        Compute L(phi, rho)

        \param output   Evaluated operator (output)
        \param state    State variables
        \param aux      Auxiliary variables
    */
    inline virtual void
    operator()(
        StateData& a_output,
        const StateData& a_state,
        const AuxData&   a_aux,
        T                a_scale = 1.0) const;
    
    /// Apply (In Place)
    /**
        Compute L(phi)

        \param output   Evaluated operator (output)
        \param state    State variables
    */
    inline virtual void
    operator()(
        StateData& a_output,
        const StateData& a_state,
        T                a_scale = 1.0) const;
    
    /// Apply (Out of Place)
    /**
        Compute L(phi, rho).

        \param state    State variables
        \param aux      Auxiliary variables
        \param range    Range of the output data
    */
    inline virtual StateData
    operator()(
        const StateData& a_state,
        const AuxData&   a_aux,
        Box              a_range,
        T                a_scale = 1.0) const;
    
    /// Apply (Out of Place)
    /**
        Compute L(phi)

        \param state    State variables
        \param range    Range of the output data
    */
    inline virtual StateData
    operator()(
        const StateData& a_state,
        Box              a_range,
        T                a_scale = 1.0) const;

    /// User Defined Flux
    /**
        Compute F(phi, rho) in the direction dir
       
        \param flux     Computed flux (output) 
        \param state    State variables
        \param aux      Auxiliary variables
        \param dir      Direction in [0, DIM)
    */
    inline virtual void flux(
        StateData& a_flux,
        const StateData& a_state,
        const AuxData&   a_aux, 
        int a_dir) const;
    
    /// User Defined Flux
    /**
        Compute F(phi) in the direction dir
        
        \param flux     Computed flux (output) 
        \param state    State variables
        \param dir      Direction in [0, DIM)
    */
    inline virtual void flux(
        StateData& a_flux,
        const StateData& a_state,
        int a_dir) const;
    
    /// User Defined Source
    /**
        Compute S(phi, rho)
       
        \param source     Computed flux (output) 
        \param state    State variables
        \param aux      Auxiliary variables
    */
    inline virtual void source(
        StateData&       a_source,
        const StateData& a_state,
        const AuxData&   a_aux) const; 
    
    /// User Defined Source
    /**
        Compute S(phi)
       
        \param source   Computed flux (output) 
        \param state    State variables
    */
    inline virtual void source(
        StateData&       a_source,
        const StateData& a_state) const;
    
    /// User Defined Diagonal
    /**
        Compute D(phi, rho)
       
        \param diag     Computed flux (output) 
        \param state    State variables
        \param aux      Auxiliary variables
    */
    inline virtual void diag(
        StateData&       a_diag,
        const StateData& a_state,
        const AuxData&   a_aux) const; 
    
    /// User Defined Diagonal
    /**
        Compute D(phi)
       
        \param diag     Computed flux (output) 
        \param state    State variables
    */
    inline virtual void diag(
        StateData&       a_diag,
        const StateData& a_state) const;
    
    /// Set Diagonal Term Scaling
    /**
        Sets the value of d
    */
    inline void setDiagScale(double a_value);
    
    /// Set Flux Term Scaling
    /**
        Sets the value of f
    */
    inline void setFluxScale(double a_value);
    
    /// Read Diagonal Term Scaling
    inline double diagScale() const { return m_scaleDiag; }
    
    /// Read Flux Term Scaling
    inline double fluxScale() const { return m_scaleFlux; }
    
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
