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

    /// Non-Trivial Isotropic Constructor
    /**
        \param dx   Grid spacing
    */
    inline BoxOp(const DisjointBoxLayout& a_layout, const LevelIndex& a_index, T a_dx);
    
    /// Non-Trivial Anisotropic Constructor
    /**
        \param dx   Grid spacing
    */
    inline BoxOp(const DisjointBoxLayout& a_layout, const LevelIndex& a_index, Array<T, DIM> a_dx);

    inline BoxOp(BoxOp<T,C_STATE,C_AUX,MEM>&& a_bop) = default;
    inline BoxOp<T,C_STATE,C_AUX,MEM>& operator=(BoxOp<T,C_STATE,C_AUX,MEM>&& a_bop) = default;
    /// Lazy Isotropic Constructor
    /**
        \param dx   Grid spacing
    */
    inline void define(const DisjointBoxLayout& a_layout, const LevelIndex& a_index, T a_dx);
    
    /// Lazy Anisotropic Constructor
    /**
        \param dx   Grid spacing
    */
    inline void define(
            const DisjointBoxLayout& a_layout, const LevelIndex& a_index, Array<T, DIM> a_dx);

    /// Initialize
    /**
        User defined function that is called inside of define.
        Useful for caching things like Stencil instances.
    */
    inline virtual void init(){}

    /// Spectral Radius (User Defined)
    /**
        Approximate spectral radius of the operator. Not required, but 
        Needed for iterative solvers. 
    */
    inline virtual T spectralRadius() const { return 1; }

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
        T                a_scale = 1.0) const;
    
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

    /// Apply All Boundary Conditions
    inline void
    applyBC(
        Array<StateData, DIM>& a_fluxes,
        const StateData& a_state) const;

    /// Apply User Specified Boundary Condition
    /**
    */
    inline virtual void
    applyBC(
        Array<StateData, DIM>& a_fluxes,
        const StateData& a_state,
        Face             a_face) const {}

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
   
    ///TODO: not implemented (see MMB version which has this implemented)
#if 0
    /// Apply All BCs relevant to a Box
    /** Automates the logic for determining if the patch
     *  is on a domain boundary and calls the user defined
     *  applyBC function with the correct Face argument
     */
    inline void applyBC(
        Array<StateData,DIM>& a_fluxes,
        const StateData&      a_state) const;

    /// Apply Specific BC (User Defined)
    inline virtual void applyBC(
            Array<StateData,DIM>& a_fluxes,
            const StateData&      a_state,
            Face                  a_face);

#endif
    /// Set Diagonal Term Scaling
    /**
        Sets the value of d
    */
    inline void setDiagScale(T a_value);
    
    /// Set Flux Term Scaling
    /**
        Sets the value of f
    */
    inline void setFluxScale(T a_value);
   
    /// Set Time
    inline void setTime(T a_time);

    /// Set Runge Kutta Stage
    inline void setRKStage(unsigned int a_stage);
    
    /// Read Diagonal Term Scaling
    inline T diagScale() const { return m_scaleDiag; }
    
    /// Read Flux Term Scaling
    inline T fluxScale() const { return m_scaleFlux; }
    
    /// Read Flux Term Scaling
    inline unsigned int RKStage() const { return m_RKStage; }
    
    /// Get Time
    inline T time() const { return m_time; }

    /// Get Grid Spacing
    inline Array<T, DIM> dx() const;

    /// Get Min Grid Spacing
    inline T dxMin() const;
    
    /// Get Max Grid Spacing
    inline T dxMax() const;
    
    /// Get Index
    inline const LevelIndex& index() const { return m_index; }

    /// Get Layout
    inline const DisjointBoxLayout& layout() const { return m_layout; }

    /// Get Box
    inline const Box& box() const {return m_layout[m_index]; }
    
    /// DEPRECATED
    /// Mapped Multiblock Utilities
#if 0

    /// Map Constructor
    /** a_x is node centered, a_J is cell averaged */
    inline BoxOp(
            const Box& a_box,
            const Array<double, DIM>& a_dx,
            const BoxData<double, DIM, MEM>& a_x, 
            const BoxData<double, 1, MEM>& a_J);
    
    /// Map Define
    /** a_x is node centered, a_J is cell averaged */
    inline void define(
            const Box& a_box,
            const Array<double, DIM>& a_dx,
            const BoxData<double, DIM, MEM>& a_x, 
            const BoxData<double, 1, MEM>& a_J);
    
    /// Get Coordinates
    /** Output is immutable and node centered */
    inline const BoxData<double, DIM, MEM>& x() const;
    
    /// Get Coordinates
    /** Output is immutable and cell averaged */
    inline const BoxData<double, 1, MEM>& jacobian() const;

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
    T m_scaleDiag;
    T m_scaleFlux;
    T m_time;
    unsigned int m_RKStage;
    LevelIndex m_index;
    DisjointBoxLayout m_layout;

    mutable bool m_definedSrce;
    mutable bool m_definedFlux;

    Array<T, DIM> m_dx;

    private:

};

#include "implem/Proto_BoxOpImplem.H"
} // end namespace Proto

#endif // end include guard
