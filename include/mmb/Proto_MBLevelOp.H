#pragma once
#ifndef _PROTO_MB_OP_
#define _PROTO_MB_OP_

#include "Proto_MBLevelBoxData.H"
#include "Proto_MBLevelMap.H"

namespace Proto {

/// Level-Scope Operator
/**
    MBLevelOp contains the necessary tools to apply a user-defined descendent of 
    BoxOp on a mapped multiblock level

    \tparam OpType  A BoxOp class
    \tparam T       Datatype of the data holder (e.g. int, double, etc.)
    \tparam MEM     Proto::MemType of the data holder
*/
template<template<typename, MemType> class OPType,
    typename T,
    MemType MEM = MEMTYPE_DEFAULT>
class MBLevelOp
{
    public:
    typedef OPType<T, MEM> OP;
    typedef BoxData<T, OP::numState(), MEM> StateData;
    typedef BoxData<T, OP::numAux(),   MEM>   AuxData;
    typedef MBLevelBoxData<T, OP::numState(), MEM, PR_CELL> LevelStateData;
    typedef MBLevelBoxData<T, OP::numAux(),   MEM, PR_CELL> LevelAuxData;

    static constexpr unsigned int numState() { return OP::numState(); }
    static constexpr unsigned int numAux()   { return OP::numAux(); }
    
    static Point ghost() { return OP::ghost(); }
    static int order() { return OP::order(); }
    /// Spectral Radius
    inline T spectralRadius() const {return m_ops[0].spectralRadius(); }
    
    inline MBLevelOp();
    inline MBLevelOp(MBLevelMap<MEM>& a_map);
    inline MBLevelOp(MBLevelOp<OPType, T, MEM>&& a_op) = default;
    inline MBLevelOp<OPType, T, MEM>&
        operator=(MBLevelOp<OPType, T, MEM>&& a_op) = default;
    inline void define(MBLevelMap<MEM>& a_map);
    
    inline void operator()(
        LevelStateData&         a_output,
        const LevelStateData&   a_state,
        T                       a_scale = 1.0);
    inline void operator()(
        LevelStateData&         a_output,
        const LevelStateData&   a_state,
        const LevelAuxData&     a_aux,
        T                       a_scale = 1.0);

    inline void setDiagScale(T a_value);
    inline void setFluxScale(T a_value);
    inline void setTime(T a_time);
    inline void setRKStage(unsigned int a_stage);
    inline T diagScale() const { return m_diagScale; }
    inline T fluxScale() const { return m_fluxScale; }
    inline T time() const {return m_time; }
    inline unsigned int RKStage() const {return m_rkStage; }
    inline const MBDisjointBoxLayout& layout() const {return m_map->map().layout(); }
    inline const OP& operator[](const LevelIndex& a_index) const { return m_ops[a_index]; }
    private:
    
    T m_diagScale;
    T m_fluxScale;
    T m_time;
    unsigned int m_rkStage;
    
    MBLevelMap<MEM>* m_map;
    std::vector<OP> m_ops;
};

#include "implem/Proto_MBLevelOpImplem.H"
} // end namespace Proto
#endif //end include guard
