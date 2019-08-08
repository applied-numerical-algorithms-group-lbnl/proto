#pragma once

#include "Proto.H"
#include "BaseOp.H"

#include <vector>
#include <iostream>

/// Multigrid Algorithm
/**
   Implements the Multigrid algorithm for standalone use or as part of an AMR hierarchy. 
*/

template <class OP>
class Multigrid
{
    public:
        typedef LevelData<FArrayBox> LD;

        /// Default Constructor
        /**
            Builds an undefined Multigrid. Call define to finish initialization.
        */
        inline Multigrid(){m_defined = false;}

        /// General Constructor
        /**
            Builds a general Multigrid.

            \param a_layout     Chombo DisjointBoxLayout of the input data
            \param a_dx         Finest grid spacing
            \param a_numLevels  Number of legels for the Multigrid vcycle
            \param a_amr        AMR flag. If true, extra space is allocated for computing boundary conditions
            \param a_amrLevel   Level of the AMR hierarchy. Unused if a_amr is false.
        */
        inline Multigrid(const DisjointBoxLayout& a_layout,
                double     a_dx,
                int        a_numLevels,
                bool       a_amr = false,
                int        a_amrLevel = 0) 
        {
            define(a_layout, a_dx, a_numLevels-1, a_amr, a_amrLevel);
        };

        /// Define
        /**
            Initializes an undefined Multigrid

            \param a_layout     Chombo DisjointBoxLayout of this level
            \param a_dx         Grid spacing at this level
            \param a_level      Level of the Multigrid vcycle this Multigrid represents
            \param a_amr        AMR flag. If true, extra space is allocated for computing boundary conditions
            \param a_amrLevel   Level of the AMR hierarchy. Unused if a_amr is false.
        */
        inline void define(const DisjointBoxLayout& a_layout,
                double     a_dx,
                int        a_level,
                bool       a_amr = false,
                int        a_amrLevel = 0);

        /// VCycle (Standalone)
        /**
            Canonical VCycle algorithm for multigrid. Used as a standalone solver and as the bottom solver in AMR algorithms.

            \param a_phi        Chombo LevelData<FArrayBox> reference containing independent variable data
            \param a_phi        Chombo LevelData<FArrayBox> reference containing right-hand side data
        */ 
        inline void vcycle(LD& a_phi, const LD& a_R);

        /// VCycle (AMR Version)
        /**
            VCycle algorithm for use in an AMR hierarchy.
            Includes additional code for computing boundary conditions at intermediate Multigrid levels.
            This functions behavior is undefined if called from a Multigrid initialized with a_amr = false.

            \param a_phi        Chombo LevelData<FArrayBox> reference containing independent variable data
            \param a_phiC       Chombo LevelData<FArrayBox> reference containing independent variable data on the next coarser AMR level
            \param a_phi        Chombo LevelData<FArrayBox> reference containing right-hand side data
        */ 
        inline void vcycle(LD& a_phi, LD& a_phiC, const LD& a_R);

        /// Get Operator
        inline OP& op(){return m_op;}
    private:

        int m_level;    ///< This level in the Multigrid hierarchy
        OP  m_op;       ///< Instance of OP associated with this level

        LevelData<FArrayBox> m_phiC0;    ///< Temporary storage for phi at the next coarser level (before recursive call)
        LevelData<FArrayBox> m_phiC;     ///< Storage for phi at the next coarser level
        LevelData<FArrayBox> m_RC;       ///< Storage for the right-hand side at the next coarser level
        std::shared_ptr<Multigrid<OP>> m_coarser;  ///< Multigrid instance at the next coarser level
        
        LevelData<FArrayBox> m_phiCAMR;  ///< Temporary storage used for regridding in AMR. Unused with standalone Multigrid
        bool m_defined;             ///< True if define(...) has been called.
};

#include "MultigridImplem.H"
