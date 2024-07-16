#pragma once
#ifndef _PROTO_MB_BOX_OP_
#define _PROTO_MB_BOX_OP_

#include "Proto_BoxOp.H"
#include "Proto_MBLevelMap.H"
namespace Proto {

/// Abstract Box-Scope Operator for Mapped Multiblock
/**
    MBBoxOp is the class from which all mapped multiblock operators are
    derived. It is a special case of BoxOp with additional built in
    capabilities for providing mapped coordinate context.

    //FIXME: update this for mapped
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
    \tparam MAP     Datatype of the MBMapOp (this is the same as the MBLevelMap's MAP parameter)
    \tparam C_STATE Number of state variables
    \tparam C_AUX   Number of auxiliary (non-state) variables
    \tparam MEM     Proto::MemType of the data holder
*/
template <typename T, template<MemType> typename MAP,
         unsigned int C_STATE, unsigned int C_AUX, MemType MEM = MEMTYPE_DEFAULT>
class MBBoxOp : public BoxOp<T, C_STATE, C_AUX, MEM>
{
    public:
    typedef BoxData<T, C_STATE, MEM> StateData;
    typedef BoxData<T, C_AUX,   MEM> AuxData;
    
    //TODO: In addition to the following, the interface requires these
    // static functions. There is no good way to specify necessary static
    // functions in an interface which is extremely awkward...
    //
    // static Point ghost() const;      (Needed ghost data in the state variables)
    // static Point auxGhost() const;   (Needed ghost data in the aux variables)

    /// Non-Trivial Constructor
    inline MBBoxOp(const MBLevelMap<MAP, MEM>& a_map, const MBIndex& a_index);
    inline MBBoxOp(MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>&& a_bop) = default;
    inline MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>&
        operator=(MBBoxOp<T,MAP,C_STATE,C_AUX,MEM>&& a_bop)=default;
    
    /// Initialize
    /**
        User defined function that is called inside of define.
        Useful for caching things like Stencil instances.
    */
    inline virtual void init() {}

    /// Spectral Radius (User Defined)
    /**
        Approximate spectral radius of the operator. Not required, but 
        Needed for iterative solvers. 
    */
    inline virtual T spectralRadius() const  { return 1; }

#if 0 //CURRENTLY UNUSED

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
        T                a_scale = 1.0) const ;
#endif
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
        T                a_scale = 1.0) const ;
#if 0 // CURRENTLY UNUSED
      
    /// Apply (In Place, Flux Output)
    /**
        Compute L(phi, rho)

        \param output   Evaluated operator (output)
        \param state    State variables
        \param aux      Auxiliary variables
    */
    inline virtual void
    operator()(
        StateData& a_output,
        Array<StateData, DIM>& a_fluxes,
        const StateData& a_state,
        const AuxData&   a_aux,
        T                a_scale = 1.0) const ;
#endif
    /// Apply (In Place, Flux Output)
    /**
        Compute L(phi)

        \param output   Evaluated operator (output)
        \param state    State variables
    */
    inline virtual void
    operator()(
        StateData& a_output,
        Array<StateData, DIM>& a_fluxes,
        const StateData& a_state,
        T                a_scale = 1.0) const ;
    
#if 0 //CURRENTLY UNUSED

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
        T                a_scale = 1.0) const ;
#endif

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
        T                a_scale = 1.0) const ;

    /// Apply All Boundary Conditions
    /** 
     *  Executes logic for determining which Face of a patch is on a domain boundary.
     *  Calls the user defined boundary condition function as needed
     */
    inline void
    applyBC(
            Array<StateData, DIM>& a_fluxes,
            const StateData& a_state) const;


    /// Apply Boundary Condition (User Specified)
    /** 
     *  Use this function to specify problem domain boundary conditions.
     *  The default behavior is a null-op, leaving the flux argument unchanged.
     */
    inline virtual void
    applyBC(
            Array<StateData, DIM>& a_fluxes,
            const StateData& a_state,
            Face             a_face) const {}

#if 0
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
        int a_dir) const ;
    
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
        int a_dir) const ;
    
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
        const AuxData&   a_aux) const ; 
    
    /// User Defined Source
    /**
        Compute S(phi)
       
        \param source   Computed flux (output) 
        \param state    State variables
    */
    inline virtual void source(
        StateData&       a_source,
        const StateData& a_state) const ;
    
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
        const AuxData&   a_aux) const ; 
    
    /// User Defined Diagonal
    /**
        Compute D(phi)
       
        \param diag     Computed flux (output) 
        \param state    State variables
    */
    inline virtual void diag(
        StateData&       a_diag,
        const StateData& a_state) const ;
#endif 
    /// Get Map
    inline const MBLevelMap<MAP,MEM>& map() const;
    
    /// Get Index
    inline MBIndex index() const;

    /// Get Layout
    inline const MBDisjointBoxLayout& layout() const;
    
    /// Get Block
    inline unsigned int block() const;
    
    /// Get Box
    inline Box box() const;
    
    /// Get Grid Spacing (Mapped Coordinates)
    inline Array<T, DIM> dx() const;

    /// Get Coordinates
    /** Output is immutable and node centered */
    inline const BoxData<double, DIM, MEM>& x() const;
    
    /// Get Coordinates
    /** Output is immutable and cell averaged */
    inline const BoxData<double, 1, MEM>& jacobian() const;

#if 0
    /// Enforce Block Boundary Flux Matching condition
    inline virtual void matchFlux(
            StateData&       a_rhs,
            const StateData& a_locFlux,
            const StateData& a_adjFlux,
            const StateData& a_state,
            Point            a_localDir) const;

    /// Match Block Boundary Fluxes using the average of fluxes
    inline void matchFluxAverage(
            StateData&       a_rhs,
            const StateData& a_locFlux,
            const StateData& a_adjFlux,
            Point            a_localDir) const;
#endif
    protected:

    private:
    
    MBIndex m_index;
    const MBLevelMap<MAP, MEM>& m_map;

};

#include "implem/Proto_MBBoxOpImplem.H"
} // end namespace Proto

#endif // end include guard
