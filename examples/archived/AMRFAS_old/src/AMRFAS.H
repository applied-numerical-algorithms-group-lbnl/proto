#pragma once

#include "Multigrid.H"

/// Adaptive Mesh Refinement Full Approximation Scheme Solver
/**
    Implements a framework for solving operators using Adaptive Mesh Refinement (AMR). 
    This particular implementation uses the Full Approximation Scheme (FAS) version of AMR wherein full right-hand side fields are
    computed and passed into each recusive level of the vcycle instead of working with residuals. Multigrid with sub-level cycling
    is used for each level solve.

    \tparam AMR_OP      A class derived from BaseOp which defines the operator that this framework will solve
*/
template<class AMR_OP>
class AMRFAS
{
    public:

        typedef AMRData<AMR_OP::numcomps()> AMR_DATA;

        /// Default Constructor
        /**
            Creates an undefined AMRFAS Solver
        */
        AMRFAS(){m_defined = false;}

        /// Simple Constructor
        /**
            Initialize an AMRFAS solver using an AMRLayout, grid spacing, and bottom solve depth.
            This is the recommended public user constructor.
            
            \param a_layout         Hierarchy of grid layouts assumed to satisfy proper nesting requirements 
            \param a_dx             Grid spacing at this level (in all spatial directions)
            \param a_bottomDepth    Number of levels of the Multigrid operator at the bottom of the AMR vcycle.
        */
        AMRFAS(const AMRLayout& a_layout,
                double a_dx,
                int a_bottomDepth)
        {
            define(a_layout, a_dx, a_layout.size()-1, a_bottomDepth);
        }
        
        /// General Constructor
        /**
            Creates a level of an AMRFAS hierarchy and all coarser levels recursively. 
            Used internally; not recommended for public use.

            \param a_layout         Hierarchy of grid layouts assumed to satisfy proper nesting requirements
            \param a_dx             Grid spacing at this level
            \param a_level          This levels in the hierarchy
            \param a_bottomDepth    Number of levels of the Multigrid operator at the bottom of the AMR vcycle.
        */
        AMRFAS(const AMRLayout& a_layout,
                double a_dx,
                int    a_level,
                int    a_bottomDepth)
        {
            define(a_layout, a_dx, a_level, a_bottomDepth);
        }
        
        /// Define
        /**
            Recursively define an AMRFAS. Not recommended for public use. 

            \param a_layout         Hierarchy of grid layouts assumed to satisfy proper nesting requirements
            \param a_dx             Grid spacing at this level
            \param a_level          This level in the hierarchy
            \param a_bottomDepth    Number of levels of the Multigrid operator at the bottom of the AMR vcycle.
        */ 
        inline void define(const AMRLayout& a_layout,
                double a_dx,
                int    a_level,
                int    a_bottomDepth);

        /// AMRVCycle
        /**
            AMRFAS VCycle algorithm with explicit residuals.  

            \param  a_phi       AMR Hierarchy of the independent variable data
            \param  a_rho       AMR Hierarchy of the right-hand side data
            \param  a_res       AMR Hierarchy of the residual field
            \param  a_n         Current iteration. Currently used for debugging. 

        */
        inline void vcycle(AMR_DATA& a_phi, const AMR_DATA& a_rho, AMR_DATA& a_res, int a_n = 0);
        
        /// Compute AMR Residual
        /**
            Initializes an AMR Hierarchy residual field using corresponding data for the independent variable and right-hand side.
            Not used inside of VCycle. This function is mostly just for utility and may be removed in a future release.
           
            \param  a_res       AMR Hierarchy of the residual field
            \param  a_phi       AMR Hierarchy of the independent variable data
            \param  a_rho       AMR Hierarchy of the right-hand side data
        */
        inline void residual(AMR_DATA& a_res, AMR_DATA& a_phi, const AMR_DATA& a_rho);
        
        /// AMR Apply
        /**
            Evaluate the operator on an AMR Hierarchy. This is the full AMR version of AMR_OP::apply. 
            
            \param  a_rhs       AMR Hierarchy of the independent variable data
            \param  a_phi       AMR Hierarchy of the right-hand side data
        */
        inline void operator()(AMR_DATA& a_rhs, AMR_DATA& a_phi, Real a_scale = 1.0);
        
    private:
       
        inline void vcycle(AMR_DATA& a_phi, const AMR_DATA& a_rho, AMR_DATA& a_res,
                const LevelData<FArrayBox>& a_R, int a_n = 0);
        
        int m_level; ///< This level

        Multigrid<AMR_OP>    m_mg;    ///< Multigrid instance used for this level
        AMR_OP               m_op;    ///< AMR operator instance for this level
        LevelData<FArrayBox> m_phiC0; ///< Temperary data storage for phi on the next coarser level (before recursive call)
        LevelData<FArrayBox> m_RC;    ///< Storage for the right-hand side computed for the next coarser level

        std::shared_ptr<AMRFAS<AMR_OP>> m_coarser;    ///< The next coarser AMRFAS level
        LevelFluxRegister    m_reflux;                ///< The flux register used for correcting the next coarser level's flux
        bool                 m_defined;               ///< True if define was called on this (false for default constructed objects)
};

#include "AMRFASImplem.H"
