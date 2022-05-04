#pragma once
#ifndef _PROTO_BLOCK_BOUNDARY_
#define _PROTO_BLOCK_BOUNDARY_

namespace Proto
{
    /// Block Boundary
    /**
        Contains the abstraction for a boundary of codimension CoDim
        Between two blocks. 
    */
    template<unsigned int CoDim>
    class BlockBoundary
    {
        public:

        /// Dimensionality
        inline constexpr int dim() { return DIM - CoDim; }

        /// Get Neighbor
        /**
            Get the index associated with the adjacent block so long
            as index is in m_blocks.
        */
        inline BlockIndex neighbor(BlockIndex index);
        inline void print() const;
        private:

        std::array<BlockIndex, 2> m_blocks;
        // TODO: something here to represent the coordinate rotation between blocks
    }; // end class BlockBoundary
   
    /// Stream Operator 
    inline std::ostream& operator<<(
        std::ostream&           a_os,
        const BlockBoundary&    a_boundary);

}
#include "implem/Proto_BlockBoundaryImplem.H"
#endif // end include guard
