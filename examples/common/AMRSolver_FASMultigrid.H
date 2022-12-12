#pragma once
#ifndef _PROTO_AMR_SOLVER_MULTIGRID_
#define _PROTO_AMR_SOLVER_MULTIGRID_

#include "Proto.H"
#include "LevelSolver_FASMultigrid.H"

using namespace Proto;

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
        std::array<T, DIM>&   a_dx);
    
    inline void define(
        AMRGrid& a_layout,
        std::array<T, DIM>&   a_dx);
   
    /// Solve with Aux Data 
    /** Not implemented */ 
    inline double solve(
        AMRStateData& a_state,
        AMRStateData& a_force,
        AMRAuxData&   a_aux,
        int           a_maxIter,
        double        a_tolerance);
    
    /// Solve
    /**
        Execute AMR solve

        \param a_state      The state (Phi) on all levels
        \param a_force      The force (G) on all levels
        \param a_maxIter    The max number of AMR V-Cycles
        \param a_tolerance  Prescribed residual tolerance 
    */
    inline double solve(
        AMRStateData& a_state,
        AMRStateData& a_force,
        int           a_maxIter,
        double        a_tolerance);
    
    /// Compute Norm of Residual
    /**
        Compute the residual and it's norm simultaneously.

        R = G - L(Phi)

        \param a_residual   The computed residual (R) [output]
        \param a_state      The state (Phi) on all levels
        \param a_force      The force (G) on all levels 
    */
    inline double resnorm(
        AMRStateData& a_residual,
        AMRStateData& a_state,
        AMRStateData& a_force);

    private:
    
    // local class for implementing recursive Multigrid structure
    // TODO: refactor this into an iterative structure.
    //       the recursive design is awkward to work with. 
    class AMRLevel
    {
        public:
                
        AMRLevel(AMRGrid& a_grid,
                 AMROP&   a_op,
                 int      a_level);
        
        void define(AMRGrid& a_grid,
                AMROP&       a_op,
                int          a_level);
       
        // Compute Fine Correction
        /**
            from two states at level L-1, compute the correction to level L.
            fineState += interp(crseState - crseState_0)
            \param a_crseState      A corrected state at L-1            [input]
            \param a_crseState_0    An initial state at L-1             [input]
            \param a_state          A state at level L to be corrected  [input/output]
        */ 
        void fineCorrect(
                LevelStateData& a_state,
                LevelStateData& a_crseState,
                LevelStateData& a_crseState_0);
       
        // Compute the Coarse Residual
        /**
            Compute the residual on level L-1.
            1.) Interpolate boundary conditions for Phi on L-1
            2.) Compute L(Phi) on L-1 (reflux accounted for)
            3.) R = G - L(Phi) on L-1
            4.) r = G - L(Phi) on L
            5.) Overwrite R with Avg(r) on C(L) in L-1

            \param a_crseRes        The residual (R) on level L-1
            \param a_state          The state (Phi) on all levels
            \param a_crseForce      The forcing (G) (right-hand side) on level L-1
            \param a_fineForce      The forcing (G) (right-hand side) on level L
        */
        void crseResidual(
                LevelStateData& a_crseRes,
                AMRStateData&   a_state,
                LevelStateData& a_crseForce,
                LevelStateData& a_fineForce);
        
        // Execute V-Cycle
        /**
            Compute a V-Cycle on this level (L) and 
            recursively calls this function on L-1 if L > 0. 
            
            IF (L == 0)
            1.) Relax Phi with FASMultigrid on level L
            ELSE
            1.) Interpolate BCs for PhiC on L-1
            2.) Relax Phi with FASMultigrid on level L 
            3.) PhiC = Avg(Phi) on C(L) in L-1 (e.g. average down)
            4.) PhiC_0 = Phi on L-1
            5.) Compute residual RC on L-1
            6.) GC = LC(PhiC) + RC on L-1
            7.) vCycle (PhiC, GC)
            8.) Phi += Interp(PhiC - PhiC_0) on L
            9.) Interpolate BCs for PhiC on L-1 (***)
            10.) Relax Phi with FASMultigrid on L

            \param a_state      The state (Phi) on all levels
            \param a_force      The force (G) on all levels
            \param a_levelForce  The modified force (G) on L
        */
        void vCycle(
                AMRStateData& a_state,
                AMRStateData& a_force,
                LevelStateData& a_levelForce);
        
        /// Compute the residual on a level
        /**
            1.) exchange(Phi)
            2.) R = G - L(Phi)

            \param a_res    The residual (R) to be computed
            \param a_state  The state (Phi) on all levels
            \param a_force  The forcing (G) on all levels
            \param a_level  The level (L) corresponding to a_res
        */
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
         
        std::shared_ptr<LevelSolver_FASMultigrid<OpType, T, MEM>> m_mgSolver;
        AMROP* m_op;

        LevelStateData m_crseLocal;
        LevelStateData m_crseLocalInterp;
        LevelStateData m_crseRes;
        LevelStateData m_crseState_0;
        LevelStateData m_crseForce;
        LevelStateData m_residual;  //TODO: for testing, remove later.
         
        std::shared_ptr<AMRLevel> m_crseAMRLevel;
    };
   
    AMRStateData            m_residual; 
    std::shared_ptr<AMRLevel> m_amrLevel;
    std::shared_ptr<AMROP>    m_amrOp;
};
#include "AMRSolver_FASMultigridI.H"
#endif // end include guard
