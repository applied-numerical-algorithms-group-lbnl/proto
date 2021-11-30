#pragma once
#ifndef _PROTO_LEVEL_OP_
#define _PROTO_LEVEL_OP_

#include "Proto_LevelBoxData.H"

namespace Proto {

template <template<typename, MemType> class OPType,
    typename T, MemType MEM = MEMTYPE_DEFAULT>
class LevelOp
{
    public:
    typedef OPType<T, MEM> OP;
    
    static constexpr unsigned int numState() { return OP::numState(); }
    static constexpr unsigned int numAux()   { return OP::numAux(); }
    
    typedef LevelBoxData<T, OP::numState(), MEM, CELL> LevelStateData;
    typedef LevelBoxData<T, OP::numAux(),   MEM, CELL> LevelAuxData;
    
    inline LevelOp();
    inline LevelOp(double a_dx);
    
    inline void
    operator()(
        LevelStateData& a_output,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux) const;

    inline void setDiagScale(double a_value);
    inline void setFluxScale(double a_value);
    inline double dx() const;
    inline OP& op() {return m_op; }
    private:
    OP m_op;

};

#include "implem/Proto_LevelOpImplem.H"
} // end proto namespace
#endif // end include guard
