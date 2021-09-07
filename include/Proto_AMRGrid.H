#pragma once
#ifndef _PROTO_AMR_GRID_
#define _PROTO_AMR_GRID_

#include "Proto_FinitePointSet.H"
#include "Proto_DisjointBoxLayout.H"
#include "Proto_LevelBoxData.H"
#include "Proto_HDF5.H"

#define PR_AMR_REFRATIO 4

namespace Proto
{

    typedef LevelBoxData<char, 1, MEMTYPE_DEFAULT> LevelTagData;
    typedef BoxData<char, 1, MEMTYPE_DEFAULT> TagData;
    
    /// AMR Grid
    /**
        AMRGrid is a collection of DisjointBoxLayout objects that define
        properly nested unions of logically rectangular space. A non-trivially
        constructed AMRGrid is guaranteed to satisfy the proper nesting
        requirement of AMR algorithms.
    */
    class AMRGrid
    {
        public:

        /// Buffer Tags
        /**
            Grows a level of tag data by bufferSize.
        */
        inline static void buffer(LevelTagData& a_tags, int a_bufferSize);
        
        /// Trivial Constructor
        AMRGrid() {m_defined = false;}

        /// Multi-Level Constructor
        /**
            maxLevels must be >= the size of layouts
        */
        AMRGrid(std::vector<DisjointBoxLayout>& a_layouts, int a_maxLevels)
        { define(a_layouts, a_maxLevels); }
        
        /// Single-Level Constructor
        /**
            maxLevels must be positive
        */
        AMRGrid(DisjointBoxLayout& a_layout, int a_maxLevels)
        { define(a_layout, a_maxLevels); }
        
        /// Multi-Level Define
        /**
            maxLevels must be >= the size of layouts
        */
        inline void
        define(std::vector<DisjointBoxLayout>& a_layouts, int a_maxLevels);        
        
        /// Single-Level Define
        /**
            maxLevels must be positive
        */
        inline void
        define(DisjointBoxLayout& a_layout, int a_maxLevels);        
        
        /// Layout Access (Const)
        inline const DisjointBoxLayout&
        operator[](unsigned int a_level) const;
        
        /// Layout Access (Non-Const)
        inline DisjointBoxLayout&
        operator[](unsigned int a_level);

        /// Regrid
        /**
            Given tag data on a level, modify or generate a grid on the next finer level.
            Input level must correspond to the grid on which input tags are defined. 
            Input level must be coarser than the finest permissible level in the grid.
            For grids with fewer than their predefined max number of levels, this
            function can add at most one level og grid.
            Proper nesting is NOT enforced by this function (see enforceNesting).
        */
        inline void regrid(LevelTagData& a_tags, unsigned int a_level);
        
        /// Enforce Nesting
        /**
            Enforce the proper nesting criteria on a given level.
            Input level must be an "internal" level. That is to say, greater than 0 
            and finer than the finest level in the array.
        */
        inline void enforceNesting(unsigned int a_level);

        /// Query Number of Levels
        inline int numLevels() { return m_layouts.size(); }

        /// Query Max Number of Levels
        inline int maxLevels() { return m_maxLevels; }

        private:

        bool m_defined;
        std::vector<DisjointBoxLayout> m_layouts;
        int m_maxLevels;
    };


#include "implem/Proto_AMRGridImplem.H"
}
#endif //end include guard
