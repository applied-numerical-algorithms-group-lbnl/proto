#pragma once
#ifndef _PROTO_AMR_SOLVER_MULTIGRID_
#define _PROTO_AMR_SOLVER_MULTIGRID_

#include "Proto.H"
#include "LevelSolver_FASMultigrid.H"
#ifndef PR_MG_REFRATIO
#define PR_MG_REFRATIO 2
#endif

using namespace Proto;

extern int AMR_LEVEL;
extern int AMR_SOLVE_ITER;
extern int AMR_RELAX_STAGE;
template <template<typename, MemType> class OpType, typename T, MemType MEM = MEMTYPE_DEFAULT>
class AMRSolver_FASMultigrid : public AMRSolver<OpType, T, MEM>
{
    public:
    
    typedef OpType<T, MEM> BOP;
    typedef LevelOp<OpType, T, MEM> LOP;
    typedef AMROp<OpType, T, MEM> AMROP;
    typedef LevelBoxData<T, BOP::numState(), MEM> LevelStateData;
    typedef LevelBoxData<T, BOP::numAux(), MEM>   LevelAuxData;
    typedef AMRData<T, BOP::numState(), MEM>      AMRStateData;
    typedef AMRData<T, BOP::numAux(), MEM>        AMRAuxData;

    inline AMRSolver_FASMultigrid(
        AMRGrid& a_grid,
        double   a_dx);
    
    inline void define(
        AMRGrid& a_layout,
        double   a_dx);
    
    inline double solve(
        AMRStateData& a_state,
        AMRStateData& a_force,
        AMRAuxData&   a_aux,
        int           a_maxIter,
        double        a_tolerance);
    
    inline double solve(
        AMRStateData& a_state,
        AMRStateData& a_force,
        int           a_maxIter,
        double        a_tolerance);

    inline double resnorm(
        AMRStateData& a_residual,
        AMRStateData& a_state,
        AMRStateData& a_force);

    private:
    
    // local class for implementing recursive Multigrid structure
    class AMRLevel
    {
        public:
                
        AMRLevel(AMRGrid& a_grid,
                 AMROP&   a_op,
                 int      a_level);
        
        void define(AMRGrid& a_grid,
                AMROP&       a_op,
                int          a_level);
        
        void fineCorrect(
                LevelStateData& a_state,
                LevelStateData& a_crseState,
                LevelStateData& a_crseState_0);
       
        void crseResidual(
                LevelStateData& a_crseRes,
                LevelStateData& a_crseState,
                LevelStateData& a_crseForce,
                LevelStateData& a_fineState,
                LevelStateData& a_fineForce);
        
        void vCycle(
                AMRStateData& a_state,
                AMRStateData& a_force,
                LevelStateData& a_crseForce);

        //TODO: purely for testing 
        void levelResidual(
                LevelStateData& a_res,
                AMRStateData& a_state,
                AMRStateData& a_force,
                int a_level);

        private:
        
        int m_level;
         
        int m_numPreRelax;
        int m_numPostRelax;
        int m_numBottomRelax;

        Stencil<T> m_average;
        InterpStencil<T> m_interpolate;
        InterpStencil<T> m_interpBC;
         
        std::shared_ptr<LevelSolver_FASMultigrid<OpType, T, MEM>> m_mgSolver;
        AMROP* m_op;

        LevelStateData m_crseLocal;
        LevelStateData m_crseRes;
        LevelStateData m_crseState_0;
        LevelStateData m_crseForce;
        LevelStateData m_residual;  //TODO: for testing, remove later.
         
        std::shared_ptr<AMRLevel> m_crseAMRLevel;
    };
   
    AMRStateData            m_residual; 
    std::shared_ptr<AMRLevel> m_amrLevel;
    std::shared_ptr<AMROP>    m_amrOp;
    InterpStencil<T> m_interpBC;
};
#include "AMRSolver_FASMultigridI.H"
#endif // end include guard
