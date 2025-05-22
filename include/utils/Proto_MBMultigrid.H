#pragma once
#ifndef _PROTO_MB_MULTIGRID_H_
#define _PROTO_MB_MULTIGRID_H_

#include "Proto.H"

namespace Proto {

template <
    template<typename, typename, MemType> class OPType,
    class MAP,
    typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType = NullBC,
    MemType MEM = MEMTYPE_DEFAULT>
class MBMultigrid
{
public:
    typedef OPType<T, MAP, MEM> OP; 
    typedef BCType<T,OP::numState(), MAP, MEM, PR_CELL> BC;
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
            const LevelStateData& a_force,
            int             a_maxIter,
            T               a_tolerance);

    inline bool validate(
            const MBDisjointBoxLayout&  a_fineLayout,
            const std::vector<Point>&   a_refRatio,
            int                         a_numLevels);

    inline int numLevels() const { return m_numLevels; }
    
    inline void residual(
            LevelStateData&         a_residual,
            LevelStateData&         a_state,
            const LevelStateData&   a_force,
            int                     a_level);
    
    inline void residual(
            LevelStateData&         a_residual,
            LevelStateData&         a_state,
            const LevelStateData&   a_force)
    {
        residual(a_residual, a_state, a_force, m_numLevels-1);
    }
    
    inline T resnorm(
            LevelStateData&         a_residual,
            LevelStateData&         a_state,
            const LevelStateData&   a_force,
            int                     a_level);
    
    inline T resnorm(
            LevelStateData&         a_residual,
            LevelStateData&         a_state,
            const LevelStateData&   a_force)
    {
        return resnorm(a_residual, a_state, a_force, m_numLevels-1);
    }

    inline LOP& op(int a_level) { return m_ops[a_level]; }
    inline LOP& op() { return op(m_numLevels-1); }
    inline MBInterpOp& interpOp(int a_level) { return *m_blockInterp[a_level]; }
    inline MBInterpOp& interpOp() { return interpOp(m_numLevels-1); }
    inline double lambda(int a_level) const { return m_lambda[a_level]; }
    inline double lambda() const { return lambda(m_numLevels-1); }
    inline MBLevelMap<MAP>& map(int a_level) { return m_maps[a_level]; }
    inline MBLevelMap<MAP>& map() { return map(m_numLevels-1); }
    inline void relax(
            LevelStateData&         a_state,
            const LevelStateData&   a_force,
            int                     a_numIter,
            int                     a_level);
    inline void relax(
            LevelStateData&         a_state,
            const LevelStateData&   a_force,
            int                     a_numIter)
    {
        relax(a_state, a_force, a_numIter, m_numLevels-1);
    }
    inline void interpBlockBoundaries(LevelStateData& a_data, int a_level)
    {
        a_data.exchange();
        m_blockInterp[a_level]->apply(a_data, a_data);
    }

    inline void applyHomogeneousBC(LevelStateData& a_data);

    // V Cycle subroutines
    inline void fineCorrect(
            LevelStateData&         a_fine,
            const LevelStateData&   a_crse,
            LevelStateData&         a_crse0,
            int                     a_level);
    
    inline void averageDown(
            LevelStateData&         a_crse,
            const LevelStateData&   a_fine,
            int                     a_level);

    inline void crseForce(
            LevelStateData&         a_crseForce,
            LevelStateData&         a_state,
            const LevelStateData&   a_force,
            LevelStateData&         a_crseState,
            int                     a_level);

    inline void vCycle(
            LevelStateData&         a_state,
            const LevelStateData&   a_force,
            int                     a_level);
    
    inline void vCycle(
            LevelStateData&         a_state,
            const LevelStateData&   a_force)
    {
        vCycle(a_state, a_force, m_numLevels-1);
    }

    private:

    std::vector<MBDisjointBoxLayout> computeLayouts(const MBDisjointBoxLayout& fineLayout) const;
    bool canCombineCoarsePatches(const MBDisjointBoxLayout& fineLayout, Point refRatio) const;
    void initializeDataHolders(const std::vector<MBDisjointBoxLayout>& layouts);

    // solver temporaries
    std::vector<std::shared_ptr<LevelStateData>>   m_state;    
    std::vector<std::shared_ptr<LevelStateData>>   m_state_0;    
    std::vector<std::shared_ptr<LevelStateData>>   m_force;
    std::vector<std::shared_ptr<LevelStateData>>   m_residual; 
    std::vector<std::shared_ptr<LevelStateData>>   m_local; 
   
    std::vector<MBLevelMap<MAP>>                        m_maps;
    std::vector<std::shared_ptr<MBInterpOp>>   m_blockInterp;

    std::vector<LOP>              m_ops;

    std::vector<Stencil<T>>       m_increment;    
    std::vector<Stencil<T>>       m_average;      
    std::vector<InterpStencil<T>> m_CFInterp;     

    std::vector<T>                m_lambda;
    int m_numPreRelax;
    int m_numPostRelax;
    int m_numBottomRelax;
    int m_numLevels;
    std::vector<Point> m_refRatios;

    int m_solveIter;
    int m_relaxIter;
    int m_solveLevel;
};

#include "implem/Proto_MBMultigridImplem.H"
} // end namespace Proto
#endif // end include guard
