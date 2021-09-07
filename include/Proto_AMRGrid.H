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

        inline static void buffer(LevelTagData& a_tags, int a_bufferSize);
        
        /// Trivial Constructor
        AMRGrid() {m_defined = false;}

        /// Multi-Level Constructor
        AMRGrid(std::vector<DisjointBoxLayout>& a_layouts) { define(a_layouts); }
        
        /// Single-Level Constructor
        AMRGrid(DisjointBoxLayout& a_layout) { define(a_layout); }
        
        /// Multi-Level Constructor
        inline void
        define(std::vector<DisjointBoxLayout>& a_layouts);        
        
        /// Single-Level Constructor
        inline void
        define(DisjointBoxLayout& a_layout);        
        
        /// Layout Access (Const)
        inline const DisjointBoxLayout&
        operator[](unsigned int a_level) const;
        
        /// Layout Access (Non-Const)
        inline DisjointBoxLayout&
        operator[](unsigned int a_level);

        ///
        inline void regrid(LevelTagData& a_tags, unsigned int a_level);
        
        inline void enforceNesting(unsigned int a_level = 0);

        private:

        bool m_defined;
        std::vector<DisjointBoxLayout> m_layouts;
    };


#include "implem/Proto_AMRGridImplem.H"
}
#endif //end include guard
