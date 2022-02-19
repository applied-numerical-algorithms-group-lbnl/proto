#pragma once
#ifndef _PROTO_LEVEL_OP_
#define _PROTO_LEVEL_OP_

#include "Proto_LevelBoxData.H"
#include "Proto_LevelFluxRegister.H"

namespace Proto {

/// Level-Scope Operator
/**
    LevelOp contains the necessary tools to apply a user-defined descendent of 
    BoxOp on a level.

    TODO: This level of abstraction might be unnecessary.

    \tparam OpType  A BoxOp class
    \tparam T       Datatype of the data holder (e.g. int, double, etc.)
    \tparam MEM     Proto::MemType of the data holder
*/

template <template<typename, MemType> class OPType,
    typename T, MemType MEM = MEMTYPE_DEFAULT>
class LevelOp
{
    public:

    typedef OPType<T, MEM> OP;
    typedef LevelBoxData<T, OP::numState(), MEM, PR_CELL> LevelStateData;
    typedef LevelBoxData<T, OP::numAux(),   MEM, PR_CELL> LevelAuxData;
    
    static constexpr unsigned int numState() { return OP::numState(); }
    static constexpr unsigned int numAux()   { return OP::numAux(); }

    /// Ghost Region Size (State)
    /**
        Amount of data needed in ghost regions of state variables to apply the operator
    */ 
    static Point ghost() {return OP::ghost(); }
    
    /// Ghost Region Size (Auxiliary)
    /**
        Amount of data needed in ghost regions of auxiliary variables to apply the operator
    */ 
    static Point auxGhost()
    {
        if (numAux() > 0) { return OP::auxGhost(); }
        return Point::Zeros();
    }
  
    static int order()
    {
        return OP::order();
    }
   
    /// Default Constructor 
    inline LevelOp();

    /// Non-Trivial Constructor
    /**
        \param dx   Grid spacing
    */
    inline LevelOp(double a_dx);
    
    /// Lazy Constructor
    /**
        \param dx   Grid spacing
    */
    inline void define(double a_dx);

    /// Spectral Radius
    inline double spectralRadius() const {return m_op.spectralRadius(); }
            
    /// Apply
    /**
        Compute L(phi, rho)

        \param output   Evaluated operator (output)
        \param state    State variables
        \param aux      Auxiliary variables
    */
    inline void
    operator()(
        LevelStateData& a_output,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux,
        T                     a_scale = 1.0) const;

    /// Apply
    /**
        Compute L(phi)

        \param output   Evaluated operator (output)
        \param state    State variables
    */
    inline void
    operator()(
        LevelStateData& a_output,
        const LevelStateData& a_stat,
        T                     a_scale = 1.0) const;
    
    inline void setDiagScale(double a_value);
    inline void setFluxScale(double a_value);
    inline double diagScale() const { return m_op.diagScale(); }
    inline double fluxScale() const { return m_op.fluxScale(); }
    inline double dx() const;
    inline const OP& op() const {return m_op; }
   
    /// =======================================================================
    // Private Interface
     
    /// Helper functions for AMR refluxing
    inline void
    coarseReflux(
        LevelFluxRegister<T, OP::numState(), MEM>&    a_fluxRegister,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux) const;
    
    /// Helper functions for AMR refluxing
    inline void
    coarseReflux(
        LevelFluxRegister<T, OP::numState(), MEM>&    a_fluxRegister,
        const LevelStateData& a_state) const;
    
    /// Helper functions for AMR refluxing
    inline void
    fineReflux(
        LevelFluxRegister<T, OP::numState(), MEM>&    a_fluxRegister,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux) const;
    
    /// Helper functions for AMR refluxing
    inline void
    fineReflux(
        LevelFluxRegister<T, OP::numState(), MEM>&    a_fluxRegister,
        const LevelStateData& a_state) const;
    
    /* TODO: Premature optimization 
    inline void
    refluxApply_CF(
        LevelStateData& a_output,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux,
        LevelFluxRegister&    a_frCoarse,
        LevelFluxRegister&    a_frFine) const;
    
    inline void
    refluxApply_C(
        LevelStateData& a_output,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux,
        LevelFluxRegister&    a_frCoarse) const;
    
    inline void
    refluxApply_F(
        LevelStateData& a_output,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux,
        LevelFluxRegister&    a_frFine) const;
    */

    private:

    OP m_op;
};

#include "implem/Proto_LevelOpImplem.H"
} // end proto namespace
#endif // end include guard
