#pragma once

#ifndef _PROTO_BLOCK_GRAPH_
#define _PROTO_BLOCK_GRAPH_

#include "Proto.H"
#include "Proto_CoordPermutation.H"

namespace Proto
{
    unsigned int codim(Point a_dir);
    
    class BlockGraphNode
    {
        public:
        BlockGraphNode(unsigned int a_index);
    
        inline bool addBoundary(
                Point                   a_dir,
                std::shared_ptr<BlockGraphNode>  a_block,
                CoordPermutation&       a_rotation);
        
        inline bool addBoundary(
                Point                   a_dir,
                std::shared_ptr<BlockGraphNode>  a_block);
        
        inline unsigned int id() const { return m_index; }
        
        inline const std::vector<std::shared_ptr<BlockGraphNode>>& boundaries(Point& a_dir) const;
        
        inline unsigned int numBoundaries(Point& a_dir) const;
        
        private:
      
        unsigned int m_index;
        std::map<Point, std::vector<std::shared_ptr<BlockGraphNode>>> m_boundaries;
        std::map<Point, CoordPermutation> m_rotations;
    };
   
    class BlockGraph
    {
        public:

        inline BlockGraph(unsigned int a_numBlocks);
       
        inline void addBoundary(
                unsigned int        a_srcBlock,
                unsigned int        a_dstBlock,
                unsigned int        a_dir,
                Side::LoHiSide      a_side,
                CoordPermutation&   a_rotation);
        
        inline unsigned int size() const {return m_blocks.size(); }

        inline const BlockGraphNode& operator[](unsigned int a_index) const;

        //std::vector<std::shared_ptr<Block>>::iterator begin() { return m_blocks.begin(); }
        //std::vector<std::shared_ptr<Block>>::iterator end() { return m_blocks.end(); }
    
        private:

        std::vector<std::shared_ptr<BlockGraphNode>> m_blocks;
    };
#include "implem/Proto_BlockGraphImplem.H"
} // end namespace Proto
#endif // end include guard
