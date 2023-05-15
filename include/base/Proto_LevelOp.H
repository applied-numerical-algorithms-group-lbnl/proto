#pragma once
#ifndef _PROTO_LEVEL_OP_
#define _PROTO_LEVEL_OP_

#include "Proto_LevelBoxData.H"

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
    typedef BoxData<T, OP::numState(), MEM> StateData;
    typedef BoxData<T, OP::numAux(), MEM> AuxData;
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
    inline LevelOp(DisjointBoxLayout& a_layout, T a_dx);
    inline LevelOp(DisjointBoxLayout& a_layout, Array<T, DIM> a_dx);
    
    /// Lazy Constructor
    /**
        \param dx   Grid spacing
    */
    inline void define(DisjointBoxLayout& a_layout, T a_dx);
    inline void define(DisjointBoxLayout& a_layout, Array<T, DIM> a_dx);

    /// Spectral Radius
    inline T spectralRadius() const {return m_ops[0].spectralRadius(); }
            
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

    inline void setDiagScale(T a_value);
    inline void setFluxScale(T a_value);
    inline void setTime(T a_time);
    inline T diagScale() const { return m_diagScale; }
    inline T fluxScale() const { return m_fluxScale; }
    inline T time() const {return m_time; }
    inline Array<T, DIM> dx() const { return m_dx; }
    inline const OP& op(LevelIndex& a_index) const { return m_ops[a_index]; }
   
    private:
   
    T m_diagScale;
    T m_fluxScale;
    T m_time;
    Array<T, DIM> m_dx;

    DisjointBoxLayout m_layout;
    std::vector<OP> m_ops;
};

#include "implem/Proto_LevelOpImplem.H"
} // end proto namespace
#endif // end include guard
