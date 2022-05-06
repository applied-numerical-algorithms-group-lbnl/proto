#pragma once

#ifndef _PROTO_BLOCK_GRAPH_
#define _PROTO_BLOCK_GRAPH_

#include "Proto.H"
#include "Proto_CoordPermutation.H"

namespace Proto
{
    unsigned int codim(Point a_dir);
    
    class BlockGraph;

    class Block
    {
        friend class BlockGraph;
        public:
        
        inline unsigned int id() { return m_index; }

        inline const std::vector<std::shared_ptr<Block>>& operator[](Point& a_dir) const;

        inline unsigned int numBoundaries(Point& a_dir) const;
        
        private:
       
        Block(unsigned int a_index);
        
        bool addBoundary(
                Point                   a_dir,
                std::shared_ptr<Block>  a_block,
                CoordPermutation&       a_rotation);
        
        bool addBoundary(
                Point                   a_dir,
                std::shared_ptr<Block>  a_block);
        
        unsigned int m_index;
        std::unordered_map<Point, std::vector<std::shared_ptr<Block>>> m_boundaries;
        std::unordered_map<Point, CoordPermutation> m_rotations;
    };

    class BlockGraph
    {
        public:

        BlockGraph(unsigned int a_numBlocks);
       
        unsigned int size() const {return m_blocks.size(); }

        void addBoundary(
                unsigned int        a_srcBlock,
                unsigned int        a_dstBlock,
                unsigned int        a_dir,
                Side::LoHiSide      a_side,
                CoordPermutation&   a_rotation);
        
        const Block& operator[](unsigned int a_index) const { return *m_blocks[a_index]; }

        //std::vector<std::shared_ptr<Block>>::iterator begin() { return m_blocks.begin(); }
        //std::vector<std::shared_ptr<Block>>::iterator end() { return m_blocks.end(); }
    
        private:

        std::vector<std::shared_ptr<Block>> m_blocks;
    };
    
#include "implem/Proto_BlockGraphImplem.H"
} // end namespace Proto
