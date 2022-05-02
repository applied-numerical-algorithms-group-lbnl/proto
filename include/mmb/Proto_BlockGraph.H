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
        

        private:
        
        Block(unsigned int a_index, Point a_blockSize);
        
        
        bool addBoundary(
                Point                   a_dir,
                std::shared_ptr<Block>  a_block,
                CoordPermutation&       a_rotation);
        
        bool addBoundary(
                Point                   a_dir,
                std::shared_ptr<Block>  a_block);
        
        ProblemDomain m_domain;
        unsigned int m_index;
        std::unordered_map<Point, std::vector<std::shared_ptr<Block>>> m_boundaries;
        std::unordered_map<Point, CoordPermutation> m_rotations;
    };

    class BlockGraph
    {
        public:

        BlockGraph(unsigned int a_numBlocks, Point a_blockSize);
       
        unsigned int size() const {return m_blocks.size(); }
        Point blockSize() const { return m_blockSize; }

        void addBoundary(
                unsigned int        a_srcBlock,
                unsigned int        a_dstBlock,
                unsigned int        a_dir,
                Side::LoHiSide      a_side,
                CoordPermutation&   a_rotation);

        private:

        
        std::vector<std::shared_ptr<Block>> m_blocks;
        Point m_blockSize;
    };
    
#include "implem/Proto_BlockGraphImplem.H"
} // end namespace Proto
