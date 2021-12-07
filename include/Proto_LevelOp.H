#pragma once
#ifndef _PROTO_LEVEL_OP_
#define _PROTO_LEVEL_OP_

#include "Proto_LevelBoxData.H"
#include "Proto_LevelFluxRegister.H"

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
 
    static constexpr Point ghost() {return OP::ghost(); }
    static constexpr Point auxGhost()
    {
        if (numAux() > 0)
        {
            return OP::auxGhost();
        }
        return Point::Zeros();
    }
    
    inline LevelOp();
    inline LevelOp(double a_dx);
    inline void define(double a_dx);

    inline double spectralRadius() const {return m_op.spectralRadius(); }
            
    inline void
    operator()(
        LevelStateData& a_output,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux) const;

    inline void
    operator()(
        LevelStateData& a_output,
        const LevelStateData& a_stat) const;

    inline void
    coarseReflux(
        LevelFluxRegister<T, OP::numState(), MEM>&    a_fluxRegister,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux) const;
    
    inline void
    coarseReflux(
        LevelFluxRegister<T, OP::numState(), MEM>&    a_fluxRegister,
        const LevelStateData& a_state) const;
    
    inline void
    fineReflux(
        LevelFluxRegister<T, OP::numState(), MEM>&    a_fluxRegister,
        const LevelStateData& a_state,
        const LevelAuxData&   a_aux) const;
    
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

    inline void setDiagScale(double a_value);
    inline void setFluxScale(double a_value);
    inline double diagScale() const { return m_op.diagScale(); }
    inline double fluxScale() const { return m_op.fluxScale(); }
    inline double dx() const;
    inline const OP& op() {return m_op; }
    private:
    OP m_op;

};

#include "implem/Proto_LevelOpImplem.H"
} // end proto namespace
#endif // end include guard
