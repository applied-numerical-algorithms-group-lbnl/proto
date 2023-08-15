#pragma once
#ifndef _PROTO_MB_AMR_GRID_
#define _PROTO_MB_AMR_GRID_

#include "Proto_AMRGrid.H"
#include "Proto_MBDisjointBoxLayout.H"

namespace Proto
{
    /// Multiblock AMR Grid
    class MBAMRGrid
    {
        public:

        inline void define(

        private:
        
        std::vector<AMRGrid> m_grids;  ///< One grid per block
        std::vector<MBDisjointBoxLayout> m_layouts; ///< One layout per level
        std::vector<std::shared_ptr<MBBoxPartition>> m_partitions; ///< One per level
    };

#include "implem/Proto_MBAMRDataImplem.H"
} // end namespace Proto
#endif //end include guard
