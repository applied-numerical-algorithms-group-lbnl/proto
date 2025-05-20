#pragma once
#ifndef _PROTO_MB_AMR_H_
#define _PROTO_MB_AMR_H_

#include "Proto.H"
namespace Proto {

template <
    template<typename, typename, MemType> class OPType,
    class MAP,
    typename T,
    template<typename, unsigned int, typename, MemType, Centering> class BCType = NullBC,
    MemType MEM = MEMTYPE_DEFAULT>
class MBAMR
{
public:
    // typedefs for cleaner reading
    typedef OPType<T, MAP, MEM> OP; 
    typedef BCType<T,OP::numState(), MAP, MEM, PR_CELL> BC;
    typedef MBLevelOp<OPType, MAP, double, BCType, MEM> LOP;
    typedef BoxData<T, OP::numState(), MEM> StateData;
    typedef MBLevelBoxData<T, OP::numState(), MEM, PR_CELL> LevelStateData;
    typedef MBAMRData<T, OP::numState(), MEM, PR_CELL> AMRStateData;
    typedef MBLevelFluxRegister<T, OP::numState(), MEM> FluxRegister;

    /// Constructor
    /**
     * uniform refinement ratio
     */
    inline MBAMR(
            const MBAMRLayout&            a_grid,
            Point                       a_refRatio);

    /// Constructor
    /**
     * block-dependent refinement ratio
     */
    inline MBAMR(
            const MBAMRLayout&            a_grid,
            const std::vector<Point>&   a_refRatios); // one refinement ratio per block

    /// Define
    inline void define(
            const MBAMRLayout&            a_grid,
            const std::vector<Point>&   a_refRatios);
    
    /// Solve
    /**
     *  Solve F = OP(S)
     */
    inline T solve(
            AMRStateData&   a_state,
            AMRStateData&   a_force,
            int             a_maxIter,
            T               a_tolerance);

    /// Validate
    /**
     *  Debugging tool.
     *  Run internal checks to make sure that all data holders were initialized correctly
     */
    inline bool validate(
            const MBAMRLayout&            a_grid,
            const std::vector<Point>&   a_refRatios);
    
    /// Interp Boundaries
    /**
     *  Interpolate onto the ghost region of a fine data set F on level Lf from a coarse data set C
     *   - Assumes C is exactly one refinement ratio coarser than F
     *   - Assumes C and F obey proper nesting
     *   - Exchange and block boundary interpolation are first computed on a COPY of C
     *   - Yields an error if Lf <= 0
     *
     *   \param a_crseData      The coarse array C
     *   \param a_fineData      The fine array F
     *   \param a_fineLevel     The fine level Lf
     */
    inline void interpBounds(
        const LevelStateData&   a_crseData,   
        LevelStateData&         a_fineData,   
        int                     a_fineLevel); 

    /// @brief Exchange data and interpolate block boundary info
    /// @param a_data 
    /// @param a_level 
    inline void interpBlockBounds(
        LevelStateData& a_data,
        int             a_level);

    /// Average Down
    /**
     *  Averages a fine data set Sf on level Lf onto the region of a coarse data set Sc.
     *   - Only the region of Sc corresponding to the coarsened domain of Sf is affected
     *   - Assumes C and F obey proper nesting
     *   - No exchange or block boundary interpolation 
     *   - Yields an error if Lf <= 0
     *
     *   \param a_crseData      The coarse array C
     *   \param a_fineData      The fine array F
     *   \param a_fineLevel     The fine level Lf
     */
    inline void averageDown(
        LevelStateData&         a_crseData,  
        const LevelStateData&   a_fineData, 
        int                     a_fineLevel); 
    
    /// @brief Average Down Recursively starting from startLevel
    /// @param data         The data to be averaged down
    /// @param startLevel   Finest level contributing to the operation. If -1 (default), all levels are used
    inline void averageDown(AMRStateData& data, int startLevel = -1);

    /// @brief Compute L(phi) on a single level WITHOUT refluxing
    /// @param a_rhs    The value of L(phi) on the specified level
    /// @param a_state  The input state (phi) on the specified level
    /// @param a_level  The specified level
    inline void levelApply(
        LevelStateData&         a_rhs,
        LevelStateData&         a_state,
        int                     a_level,
        T                       a_scale = 1.0); 

    /// @brief Compute L(phi) on a single level with refluxing at coarse-fine boundaries if necessary
    /// @param a_rhs        The value of L(phi) on the specified level
    /// @param a_state      The input state (phi) on the specified level
    /// @param a_fineState  The input state (phi) on the next finer level (used for refluxing)
    /// @param a_level      The specfied level
    inline void levelApply(
        LevelStateData&         a_rhs,
        LevelStateData&         a_state,
        LevelStateData&         a_fineState,
        int                     a_level,
        T                       a_scale = 1.0); 

    /// Compute Coarse Residual
    /**
     *  Computes the coarse residual Rc on Bc. e.g.:
     *  Rc = Fc - OP(Sc) on Bc - C(Bf)
     *  Rc = Ff - OP(Sf) on C(Bf)
     *  Where:
     *   - Bc and Bf are coarse and fine domains respectively
     *   - Rc = coarse residual on Bc
     *   - Fc = coarse force on Bc
     *   - Ff = fine force on Bf
     *   - Sc = coarse state on Bc
     *   - Sf = fine state on Bf
     *   - Fine data are on level Lf, coarse data are on Level Lc
     *   - C  = coarsening operator between level Lf and Lf-1 == Lc
     *  Assumptions / Details:
     *   - Lf and Lc obey proper nesting constraints
     *   - Lf > 0
     *   - First calls interpBounds(Sc, Sf, Lf)
     *   - Executes exchange and block boundary interpolation on Sf
     *
     *  \param a_crseResidual   the output coarse residual Rc
     *  \param a_state          AMR hierarchy of state data. Contains Sf and Sc
     *  \param a_crseForce      the force Fc on level Lc
     *  \param a_fineForce      the force Ff on level Lf
     *  \param a_fineLevel      the level index Lf > 0
     */
    inline void coarseResidual(
            LevelStateData&         a_crseResidual,
            AMRStateData&           a_state,
            const LevelStateData&   a_crseForce,
            const LevelStateData&   a_fineForce,
            int                     a_fineLevel);

    /// Apply Correction to Fine Level
    /**
     *  Update State Sf on level Lf using data on the next coarser level Lc
     *  Before (Sc0) and after (Sc1) a V-Cycle. e.g.:
     *  Sf += I(Sc1 - Sc0)
     *  Where:
     *   - Sf is on domain Bf
     *   - Sc0 and Sc1 are on domain Bc
     *   - C is the coarsening operator between level Lf and Lc
     *   - I is a piecewise constant interpolation operator
     *   - Coarse data in Bc outside of C(Bf) do not participate
     *   - Sc0 is overwritten with Sc1 - Sc0
     *
     *  \param a_fineState      the state being corrected, Sf
     *  \param a_crseState      the updated coarse state, Sc1
     *  \param a_crseState_0    the previous coarse state, Sc0
     *  \param a_fineLevel      the fine level Lf > 0
     */
    inline void fineCorrect(
            LevelStateData&         a_fineState,
            const LevelStateData&   a_crseState,
            LevelStateData&         a_crseState_0,
            int                     a_fineLevel);
    
    /// Access Map
    inline MBLevelMap<MAP>& map(int a_level);
    
    /// Access Block Boundary Interpolation Operator
    inline MBInterpOp& interpOp(int a_level);

    /// Access Level Operator
    inline MBLevelOp<OPType, MAP, double, BCType, MEM>& op(int a_level);
    
    /// Compute V Cycle
    /**
     *  Implements the AMR algorithm.
     *  The "level force" is computed recursively
     *  from the finer level as: Fc = OP(Sc) + Rc
     *  Where Rc is the coarse residual computed as shown in coarseResidual
     *
     *  \param a_state      AMR hierarchy of all state data
     *  \param a_force      AMR hierarchy of all force data
     *  \param a_levelForce Force from finer level
     *  \param a_level      Level associated with the level force
     */
    inline void vCycle(
            AMRStateData&         a_state,
            const AMRStateData&   a_force,
            const LevelStateData& a_levelForce,
            int                   a_level);
private:

    int m_solveIter;
    MBAMRMap<MAP, MEM> m_map; ///< AMR Map handle

    std::vector<MBMultigrid<OPType, MAP, T, BCType, MEM>> m_mg; ///< MG solvers for each level
    
    std::vector<std::shared_ptr<LevelStateData>> m_residual; ///< temporary which holds Rc
    std::vector<std::shared_ptr<LevelStateData>> m_state_0; ///< temporary which holds Sc0
    std::vector<std::shared_ptr<LevelStateData>> m_force;   ///< temporary holding level forces
    std::vector<std::shared_ptr<LevelStateData>> m_local;   ///< temporary data holders for coarsened fine data (e.g. on domains C(Bf) for each fine domain Bf).
    std::vector<std::unique_ptr<FluxRegister>> m_fluxRegisters;

    std::vector<std::shared_ptr<MBInterpOp>>     m_localBlockInterp; ///< block boundary interp ops
    std::vector<Stencil<T>> m_average;          ///< average down ops
    std::vector<InterpStencil<T>> m_constInterp;///< piecewise constant interp ops
    std::vector<InterpStencil<T>> m_BCInterp;   ///< coarse-fine boundary interp ops (high order) 
    int m_numLevels;
};

#include "implem/Proto_MBAMRImplem.H"
} // end Proto namespace
#endif //end include guard
