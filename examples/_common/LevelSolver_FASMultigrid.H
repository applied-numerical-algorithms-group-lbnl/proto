#pragma once
#ifndef _PROTO_LEVEL_SOLVER_FAS_MULTIGRID_
#define _PROTO_LEVEL_SOLVER_FAS_MULTIGRID_

#include "Proto.H"
#include "Proto_LevelSolver.H"

using namespace Proto;

template <template<typename, MemType> class OpType, typename T, MemType MEM = MEMTYPE_DEFAULT>
class LevelSolver_FASMultigrid : public LevelSolver<OpType, T, MEM>
{
    public:
   
    typedef OpType<T, MEM> BOP;
    typedef LevelOp<OpType, T, MEM> LOP;
    typedef LevelBoxData<T, BOP::numState(), MEM> LevelStateData;
    typedef LevelBoxData<T, BOP::numAux(), MEM> LevelAuxData;

    inline LevelSolver_FASMultigrid(
        DisjointBoxLayout& a_layout,
        Point              a_refRatio,
        int                a_numLevels,
        T                  a_dx);
    
    inline LevelSolver_FASMultigrid(
        DisjointBoxLayout& a_layout,
        Point              a_refRatio,
        int                a_numLevels,
        Array<T, DIM> a_dx);
    
    inline void define(
        DisjointBoxLayout& a_layout,
        Point              a_refRatio,
        int                a_numLevels,
        T                  a_dx);
    
    inline void define(
        DisjointBoxLayout& a_layout,
        Point              a_refRatio,
        int                a_numLevels,
        Array<T, DIM> a_dx);

    /*
    inline void defineAsSubcycler(
        DisjointBoxLayout& a_layout,
        Point              a_refRatio,
        int                a_numLevels,
        T             a_dx);
    */

    inline T solve(
        LevelStateData& a_state,
        LevelStateData& a_forcing,
        LevelAuxData&   a_aux,
        int             a_maxIter,
        T               a_tolerance);
    
    inline T solve(
        LevelStateData& a_state,
        LevelStateData& a_forcing,
        int             a_maxIter,
        T               a_tolerance);
    
    inline T resnorm(
        LevelStateData& a_residual,
        LevelStateData& a_state,
        LevelStateData& a_forcing);

    inline LevelStateData& res() {return m_residual; }

    inline void interpBC(LevelStateData& a_crseStateBC);

    inline void setVerbose(bool a_flag);

    private:
    
    // local class for implementing recursive Multigrid structure
    class MGLevel
    {
        public:
                
        MGLevel(DisjointBoxLayout&      a_layout,
                Point                   a_refRatio,
                Array<T, DIM>      a_dx,
                int                     a_level);
        
        void define(DisjointBoxLayout&  a_layout,
                Point                   a_refRatio,
                Array<T, DIM>      a_dx,
                int                     a_level);
        /*
        void defineAsSubcycler(DisjointBoxLayout&  a_layout,
                Point                   a_refRatio,
                T                       a_dx,
                int                     a_level);
        */

        void interpBC(LevelStateData& a_crseStateBC);

        void fineCorrect(
                LevelStateData& a_fine,
                LevelStateData& a_crse,
                LevelStateData& a_crse0);
        
        void coarseForce(
                LevelStateData& a_crseForce,
                LevelStateData& a_state,
                LevelStateData& a_force,
                LevelStateData& a_crseState);
        
        void relax(
                LevelStateData& a_state,
                LevelStateData& a_force,
                int a_numIter);
        
        void vCycle(
                LevelStateData& a_state,
                LevelStateData& a_force);

        void residual(
                LevelStateData& a_residual,
                LevelStateData& a_state,
                LevelStateData& a_force);
                
        LOP& op() { return m_levelOp; }
        
        private:
        
        std::shared_ptr<MGLevel> m_crseMG;
        LOP m_levelOp;
        int m_level;
        Point m_refRatio;
        T m_lambda;
         
        int m_numPreRelax;
        int m_numPostRelax;
        int m_numBottomRelax;

        LevelStateData m_crseState_0;
        LevelStateData m_crseState;
        LevelStateData m_crseForce;
        LevelStateData m_residual;
        LevelAuxData m_crseAux;
        LevelStateData m_crseLocal; // for averaging down
    
        Stencil<T> m_average;
        Stencil<T> m_increment;
        InterpStencil<T> m_interpolate;
        InterpStencil<T> m_interpBC;
    };
  
    bool m_verbose; 
    LevelStateData m_residual; 
    std::shared_ptr<MGLevel> m_mg;
};
#include "LevelSolver_FASMultigridI.H"
#endif // end include guard
