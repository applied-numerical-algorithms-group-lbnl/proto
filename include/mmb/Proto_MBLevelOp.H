#pragma once
#ifndef _PROTO_MB_OP_
#define _PROTO_MB_OP_

#include "Proto_MBLevelBoxData.H"
#include "Proto_MBLevelMap.H"
#include "Proto_MBBoxOp.H"

namespace Proto {

/// Level-Scope Operator
/**
    MBLevelOp contains the necessary tools to apply a user-defined descendent of 
    MBBoxOp on a mapped multiblock level

    \tparam OpType  A MBBoxOp class
    \tparam MAP     A MBMapOp class
    \tparam T       Datatype of the data holder (e.g. int, double, etc.)
    \tparam BCType  A LevelBC class
    \tparam MEM     Proto::MemType of the data holder
*/
template<template<typename, typename, MemType> class OPType,
    typename MAP,
    typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType = NullBC,
    MemType MEM = MEMTYPE_DEFAULT>
class MBLevelOp
{
    public:
    typedef OPType<T, MAP, MEM> OP;
    typedef BCType<T,OP::numState(), MAP, MEM, PR_CELL> BC;
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
    inline MBLevelOp(const MBLevelMap<MAP,MEM>& a_map);
    inline MBLevelOp(MBLevelOp<OPType, MAP, T, BCType, MEM>&& a_op) = default;
    inline MBLevelOp<OPType, MAP, T, BCType, MEM>&
        operator=(MBLevelOp<OPType, MAP, T, BCType, MEM>&& a_op) = default;
    inline void define(const MBLevelMap<MAP,MEM>& a_map);
    
    inline void operator()(
        LevelStateData&         a_output,
        const LevelStateData&   a_state,
        T                       a_scale = 1.0);
    inline void operator()(
        LevelStateData&         a_output,
        const LevelStateData&   a_state,
        const LevelAuxData&     a_aux,
        T                       a_scale = 1.0);
    inline OP& operator[](MBIndex a_index);
    inline const OP& operator[](MBIndex a_index) const;
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
    inline void matchFlux(LevelStateData& a_rhs, const LevelStateData& a_state);
    private:
    
    T m_diagScale;
    T m_fluxScale;
    T m_time;
    unsigned int m_rkStage;
       
    const MBLevelMap<MAP,MEM>* m_map;
   // mutable MBBoundaryRegister<T, OP::numState(), MEM>  m_boundRegister;
    std::vector<OP> m_ops;
};

#include "implem/Proto_MBLevelOpImplem.H"
} // end namespace Proto
#endif //end include guard
