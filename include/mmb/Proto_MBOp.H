#pragma once
#ifndef _PROTO_MB_OP_
#define _PROTO_MB_OP_

#include "Proto_MBLevelBoxData.H"

namespace Proto {

template<template<typename, MemType> class OPType,
    typename T, MemType MEM, typename MAP>
class MBOp
{
    public:
    typedef OPType<T, MEM> BOP;
    typedef BoxData<T, BOP::numState(), MEM> StateData;
    typedef BoxData<T, BOP::numAux(),   MEM>   AuxData;
    typedef MBLevelBoxData<T, BOP::numState(), MEM, PR_CELL> LevelStateData;
    typedef MBLevelBoxData<T, BOP::numAux(),   MEM, PR_CELL> LevelAuxData;

    static constexpr unsigned int numState() { return BOP::numState(); }
    static constexpr unsigned int numAux()   { return BOP::numAux(); }
    
    static Point ghost() { return BOP::ghost(); }
    static int order() { return BOP::order(); }
    
    inline MBOp();
    inline MBOp(MBMap<MAP, MEM>& a_map);
    inline void define(MBMap<MAP, MEM>& a_map);
    
    inline void operator()(
        LevelStateData&         a_output,
        const LevelStateData&   a_state,
        T                       a_scale = 1.0);
    inline void operator()(
        LevelStateData&         a_output,
        const LevelStateData&   a_state,
        const LevelAuxData&     a_aux,
        T                       a_scale = 1.0);

    private:
    
    bool m_defined;
    MBMap<MAP, MEM>* m_map;
    
}
#include "implem/Proto_MBOpImplem.H"
} // end namespace Proto
#endif //end include guard
