#pragma once
#ifndef _PROTO_MB_AMR_H_
#define _PROTO_MB__H_

#include "Proto.H"

namespace Proto {

template <
    template<typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType = PeriodicBC,
    MemType MEM = MEMTYPE_DEFAULT>
class MBAMR
{
public:
    typedef OPType<T, MEM> OP; 
    typedef BCType<T,OP::numState(), MEM, PR_CELL> BC;
    typedef MBLevelOp<OPType, MAP, double, BCType, MEM> LOP;
    typedef MBLevelBoxData<T, OP::numState(), MEM, PR_CELL> LevelStateData;
    typedef MBAMRData<T, OP::numState(), MEM, PR_CELL> AMRStateData;

    inline MBAMR(
            const MBAMRGrid&            a_grid,
            Point                       a_refRatio);

    inline MBAMR(
            const MBAMRGrid&            a_grid,
            const std::vector<Point>&   a_refRatios);

    inline void define(
            const MBAMRGrid&            a_grid,
            const std::vector<Point>&   a_refRatios);
    
    inline T solve(
            AMRStateData&   a_state,
            AMRStateData&   a_force,
            int             a_maxIter,
            T               a_tolerance);

    inline bool validate(
            const MBAMRGrid&            a_grid,
            const std::vector<Point>&   a_refRatios);

    inline MBLevelMap<MAP>& map(int a_level);
    inline MBInterpOp& interpOp(int a_level);
    inline MBLevelOp<OPType, MAP, double, BCType, MEM>& op(int a_level);
private:

    inline void vCycle(
            AMRStateData&         a_state,
            const AMRStateData&   a_force,
            const LevelStateData& a_levelForce,
            int                   a_level);
    
    std::vector<MBMultigrid<OPType, MAP, T, BCType, MEM>> m_mg;
    
    std::vector<std::shared_ptr<LevelStateData>> m_residual;
    std::vector<std::shared_ptr<LevelStateData>> m_state;
    std::vector<std::shared_ptr<LevelStateData>> m_state_0;
    std::vector<std::shared_ptr<LevelStateData>> m_force;
    std::vector<std::shared_ptr<LevelStateData>> m_local;

    std::vector<std::shared_ptr<MBInterpOp>>     m_localBlockInterp;
    std::vector<Stencil<T>> m_average; 
    std::vector<InterpStencil<T>> m_constInterp; 
    std::vector<InterpStencil<T>> m_BCInterp; 
    int m_numLevels;
};

#include "implem/Proto_MBAMRImplem.H"
} // end Proto namespace
#endif //end include guard
