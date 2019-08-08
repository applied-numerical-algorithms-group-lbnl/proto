#pragma once

#include "Multigrid.H"

/// Adaptive Mesh Refinement: Full Approximation Scheme (AMRFAS)
/**
    Implements a framework for solving operators using Adaptive Mesh Refinement (AMR). 
    This particular implementation uses the Full Approximation Scheme (FAS) version of AMR wherein full right-hand side fields are
    computed and passed into each recusive level of the vcycle instead of working with residuals.

    \tparam AMR_OP      A class derived from BaseOp which defines the operator that this framework will solve
    \tparam FArrayBox        The data holder for each patch. For now this must be a valid template parameter for Chombo's LevelData class
*/
template<class AMR_OP>
class AMRFAS
{
    public:

        typedef AMRData<AMR_OP::numcomps()> AMR_DATA;

        /// Default Constructor
        /**
            Creates an undefined AMRFAS
        */
        AMRFAS(){m_defined = false;}

        /// Simple Constructor
        /**
            Create an AMRFAS hierarchy from a vector of grids.
            
            \param a_grids          std::vector of Chombo DisjointBoxLayout objects that defines the grid at each level
            \param a_dx             Grid spacing at this level
            \param a_bottomDepth    Number of levels of the Multigrid operator at the bottom of the AMR vcycle.
        */
        //AMRFAS(const std::vector<DisjointBoxLayout>& a_grids,
        AMRFAS(const AMRLayout& a_layout,
                double a_dx,
                int a_bottomDepth)
        {
            define(a_layout, a_dx, a_layout.size()-1, a_bottomDepth);
        }
        
        /// General Constructor
        /**
            Creates a level of an AMRFAS hierarchy and all coarser levels. 

            \param a_grids          std::vector of Chombo DisjointBoxLayout objects that defines the grid at each level
            \param a_dx             Grid spacing at this level
            \param a_level          This levels in the hierarchy
            \param a_bottomDepth    Number of levels of the Multigrid operator at the bottom of the AMR vcycle.
        */
        //AMRFAS(const std::vector<DisjointBoxLayout>& a_grids,
        AMRFAS(const AMRLayout& a_layout,
                double a_dx,
                int    a_level,
                int    a_bottomDepth)
        {
            define(a_layout, a_dx, a_level, a_bottomDepth);
        }
        
        /// Define
        /**
            Recursively define an AMRFAS

            \param a_grids          std::vector of Chombo DisjointBoxLayout objects that defines the grid at each level
            \param a_dx             Grid spacing at this level
            \param a_level          This level in the hierarchy
            \param a_bottomDepth    Number of levels of the Multigrid operator at the bottom of the AMR vcycle.
        */ 
        //inline void define(const std::vector<DisjointBoxLayout>& a_grids,
        inline void define(const AMRLayout& a_layout,
                double a_dx,
                int    a_level,
                int    a_bottomDepth);

        /// AMRVCycle
        /**
            
            AMRFAS VCycle algorithm. All inputs are AMR Hierarchies implemented as std::vector of std::shared_ptr to Chombo LevelData<FArrayBox>.

            \param  a_phi       AMR Hierarchy of the independent variable data
            \param  a_rho       AMR Hierarchy of the right-hand side data
            \param  a_res       AMR Hierarchy of the residual field
            \param  a_n         Current iteration. Used for debugging. 

        */
        inline void vcycle(AMR_DATA& a_phi, const AMR_DATA& a_rho, AMR_DATA& a_res, int a_n = 0);
        /*
        inline void vcycle(std::vector<std::shared_ptr<LevelData<FArrayBox>>> a_phi,
                const std::vector<std::shared_ptr<LevelData<FArrayBox>>> a_rho,
                const std::vector<std::shared_ptr<LevelData<FArrayBox>>> a_res,
                int a_n = 0);
        */
        /// Compute AMR Residual
        /**
            Initializes an AMR Hierarchy residual field using corresponding data for the independent variable and right-hand side.
            Not used inside of VCycle. This function is mostly just for utility.
            
            \param  a_res       AMR Hierarchy of the residual field
            \param  a_phi       AMR Hierarchy of the independent variable data
            \param  a_rho       AMR Hierarchy of the right-hand side data
        */
        inline void residual(AMR_DATA& a_res, AMR_DATA& a_phi, const AMR_DATA& a_rho);
        /*
        inline void residual(
            std::vector<std::shared_ptr<LevelData<FArrayBox>>> a_res,
            const std::vector<std::shared_ptr<LevelData<FArrayBox>>> a_phi,
            const std::vector<std::shared_ptr<LevelData<FArrayBox>>> a_rho);
        */
        /// AMR Apply
        /**
            Evaluate the operator on an AMR Hierarchy. This is the full AMR version of AMR_OP::apply. 
            
            \param  a_rhs       AMR Hierarchy of the independent variable data
            \param  a_phi       AMR Hierarchy of the right-hand side data
        */
        inline void operator()(AMR_DATA& a_rhs, AMR_DATA& a_phi, Real a_scale = 1.0);
        /// Get Operator
        /**
            Return the AMR_OP for this level
        */
        inline AMR_OP& op(){return m_op;}
        
        /// Get All Operators
        /**
            Return all AMR_OP instances associated with this. 

            \param a_ops    std::vector of AMR_OP* which will be populated
        */
        inline void ops(std::vector<AMR_OP*>& a_ops)
        {
            if (m_level > 0){m_coarser->ops(a_ops);}
            a_ops.push_back(&m_op);
        }

        /// Write AMR Hierarchy
        /**
            Utility function for writing an AMR Hierarchy to HDF5. Useful for debugging. 
            Will likely be moved to an as-of-yet unwritten AMRHierarchy class.
            
            \param a_data   AMR Hierarchy to be written
            \param a_fname  Desired filename. May include printf-style formatting.
            \param a_params Optional printf-style arguments for a_fname
        */
        /*
        template<typename... Args>
        inline void write(std::vector<std::shared_ptr<LevelData<FArrayBox>>> a_data,
                const char* a_fname,
                Args... a_params);
        */
    private:
       
        inline void vcycle(AMR_DATA& a_phi, const AMR_DATA& a_rho, AMR_DATA& a_res,
                const LevelData<FArrayBox>& a_R, int a_n = 0);
        /* 
        inline void vcycle(std::vector<std::shared_ptr<LevelData<FArrayBox>>> a_phi,
                const std::vector<std::shared_ptr<LevelData<FArrayBox>>> a_rho,
                const std::vector<std::shared_ptr<LevelData<FArrayBox>>> a_res,
                const LevelData<FArrayBox>& a_R,
                int a_n = 0);
        */
        
        int                        m_level; ///< This level

        Multigrid<AMR_OP>    m_mg;    ///< Multigrid instance used for this level
        AMR_OP                     m_op;    ///< AMR operator instance for this level
        LevelData<FArrayBox>            m_phiC0; ///< Temperary data storage for phi on the next coarser level (before recursive call)
        LevelData<FArrayBox>            m_RC;    ///< Storage for the right-hand side computed for the next coarser level

        std::shared_ptr<AMRFAS<AMR_OP>> m_coarser;    ///< The next coarser AMRFAS level
        LevelFluxRegister          m_reflux;                ///< The flux register used for correcting the next coarser level's flux
        bool m_defined;                                     ///< True if this is a well define was called on this
};

#include "AMRFASImplem.H"
