#pragma once
#ifndef _PROTO_AMR_OP_
#define _PROTO_AMR_OP_

#include "Proto_AMRData.H"

namespace Proto {

template <template<typename, MemType> class OPType,
    typename T, MemType MEM = MEMTYPE_DEFAULT>
class AMROp
{
    public:

    typedef OPType<T, MEM> BOP;
    typedef LevelOp<OPType, T, MEM> LOP;

    static constexpr unsigned int numState() { return BOP::numState(); }
    static constexpr unsigned int numAux()   { return BOP::numAux(); }
    static constexpr Point ghost()      { return BOP::ghost(); }
    static constexpr Point auxGhost()
    {
        if (numAux() > 0)
        {
            return BOP::auxGhost();
        }
        return Point::Zeros();
    }
    
    typedef AMRData<T, BOP::numState(), MEM, CELL> AMRStateData;
    typedef AMRData<T, BOP::numAux(),   MEM, CELL> AMRAuxData;
    
    inline AMROp();
    inline AMROp(AMRGrid& a_grid, double a_cdx);
    inline void define(AMRGrid& a_grid, double a_cdx);
    
    inline void
    operator()(
        AMRStateData& a_output,
        const AMRStateData& a_state,
        const AMRAuxData&   a_aux) const;
    
    inline void
    operator()(
        AMRStateData& a_output,
        const AMRStateData& a_state) const;

    inline void setDiagScale(double a_value, int a_gridFactor = 0);
    inline void setFluxScale(double a_value, int a_gridFactor = 0);
    inline double dx(int a_level) const { return m_dx[a_level]; }
    inline const LOP& levelOp(int a_level) const { return m_levelOps[a_level]; }
    inline const BOP& boxOp(int a_level) const { return levelOp(a_level).op(); }
    private:
    bool m_defined;    
    AMRGrid m_grid; //this is intentionally a deep copy.
    std::vector<double>      m_dx;
    std::vector<LOP>         m_levelOps;
    std::vector<std::shared_ptr<LevelFluxRegister<T>>> m_fluxRegisters;

};

#include "implem/Proto_AMROpImplem.H"

} // end namespace Proto
#endif // end include guard

