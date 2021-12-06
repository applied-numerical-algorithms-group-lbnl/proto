#pragma once
#ifndef _PROTO_LEVEL_SOLVER_
#define _PROTO_LEVEL_SOLVER_

#include "Proto_LevelOp.H"

namespace Proto {

template <template<typename, MemType> class OpType, typename T, MemType MEM = MEMTYPE_DEFAULT>
class LevelSolver
{
    public:
    typedef OpType<T, MEM> BOP;
    typedef LevelOp<OpType, T, MEM> LOP;
    typedef LevelBoxData<T, BOP::numState(), MEM, CELL> LevelStateData;
    typedef LevelBoxData<T, BOP::numAux(),   MEM, CELL> LevelAuxData;
    static constexpr unsigned int numState() { return BOP::numState(); }
    static constexpr unsigned int numAux()   { return BOP::numAux(); }
    
    inline virtual void
    solve(
        LevelStateData& a_state,
        LevelStateData& a_forcing,
        LevelAuxData&   a_aux,
        int             a_maxIter,
        double          a_tolerance) = 0;

};
} // end namespace Proto

#endif // end include guard
