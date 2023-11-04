#pragma once
#ifndef _PROTO_MB_MULTIGRID_H_
#define _PROTO_MB_MULTIGRID_H_

#include "Proto.H"

namespace Proto {

template <
    template<typename, MemType> class OPType,
    template<MemType> class MAP,
    typename T,
    template<typename, unsigned int, MemType, Centering> class BCType = PeriodicBC,
    MemType MEM = MEMTYPE_DEFAULT>
class MBMultigrid
{
public:
    typedef OPType<T, MEM> OP; 
    typedef BCType<T,OP::numState(), MEM, PR_CELL> BC;
    typedef MBLevelOp<OPType, MAP, double, BCType, MEM> LOP;
    typedef MBLevelBoxData<T, OP::numState(), MEM, PR_CELL> LevelStateData;
    
    inline MBMultigrid(
            const MBDisjointBoxLayout&  a_fineLayout,
            Point                       a_refRatio,
            int                         a_numLevels);
    
    inline MBMultigrid(
            const MBDisjointBoxLayout&  a_fineLayout,
            const std::vector<Point>&   a_refRatios,
            int                         a_numLevels);
    
    inline void define(
            const MBDisjointBoxLayout&  a_fineLayout,
            const std::vector<Point>&   a_refRatio,
            int                         a_numLevels);
    
    inline T solve(
            LevelStateData& a_state,
            LevelStateData& a_force,
            int             a_maxIter,
            T               a_tolerance);
    
    inline bool validate(
            const MBDisjointBoxLayout&  a_fineLayout,
            Point                       a_refRatio,
            int                         a_numLevels);

private:

    // V Cycle subroutines
    inline void fineCorrect(
            LevelStateData&         a_fine,
            const LevelStateData&   a_crse,
            const LevelStateData&   a_crse0);

    inline void crseForce(
            LevelStateData&         a_crseForce,
            const LevelStateData&   a_state,
            const LevelStateData&   a_force,
            const LevelStateData&   a_crseState);

    inline void relax(
            LevelStateData&         a_state,
            const LevelStateData&   a_rhs,
            int                     a_numIter);

    inline void vCycle(
            LevelStateData&         a_state,
            const LevelStateData&   a_force);

    inline void residual(
            LevelStateData&         a_residual,
            const LevelStateData&   a_state,
            const LevelStateData&   a_force);

    // solver temporaries
    std::vector<std::shared_ptr<LevelStateData>>  m_state;    
    std::vector<std::shared_ptr<LevelStateData>>   m_state_0;    
    std::vector<std::shared_ptr<LevelStateData>>   m_force;
    std::vector<std::shared_ptr<LevelStateData>>   m_residual; 
    std::vector<std::shared_ptr<LevelStateData>>   m_local; 
   
    std::vector<MBLevelMap<MAP>>  m_maps;
    std::vector<std::shared_ptr<MBInterpOp>>       m_blockInterp;

    std::vector<LOP>              m_ops;

    std::vector<Stencil<T>>       m_increment;    
    std::vector<Stencil<T>>       m_average;      
    std::vector<InterpStencil<T>> m_CFInterp;     

    std::vector<T>                m_lambda;
    int m_numPreRelax;
    int m_numPostRelax;
    int m_numBottomRelax;

};

#include "implem/Proto_MBMultigridImplem.H"
} // end namespace Proto
#endif // end include guard
