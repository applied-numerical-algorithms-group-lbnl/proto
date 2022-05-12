#pragma once

#ifndef _PROTO_BLOCK_GRAPH_
#define _PROTO_BLOCK_GRAPH_

#include "Proto.H"
#include "Proto_CoordPermutation.H"

namespace Proto
{
    unsigned int codim(Point a_dir);
   
    std::vector<Point> codimDirs(unsigned int a_codim);

    struct BlockGraphNode;
    typedef std::pair<std::shared_ptr<BlockGraphNode>, Point> ArcPair;
    
    struct BlockGraphNode
    {
        BlockGraphNode(unsigned int a_index);
     
        void addArc(Point a_dir, ArcPair& a_arc)
        {
            // ignore self arcs
            if (a_arc.first->index == this->index) { return; }
            for (auto arc : boundaries[a_dir])
            {
                // ignore duplicate arcs
                if (arc == a_arc) { return; }
            }
            boundaries[a_dir].push_back(a_arc);
        }

        unsigned int index;
        std::map<Point, std::vector<ArcPair>> boundaries;
        std::map<Point, CoordPermutation> rotations;
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

        inline void closeCircuit(
                unsigned int a_srcBlock,
                unsigned int a_dstBlock);

        inline unsigned int size() const {return m_blocks.size(); }

        //inline const BlockGraphNode& operator[](unsigned int a_index) const;

        inline unsigned int numBoundaries(unsigned int a_blockIndex, Point& a_dir) const;

        inline std::set<Point> connectivity(unsigned int a_srcBlock, unsigned int a_dstBlock) const;

        inline unsigned int adjacent(unsigned int a_srcBlock, unsigned int a_dir, Side::LoHiSide a_side) const;
        
        inline void print() const;
    
        private:

        std::vector<std::shared_ptr<BlockGraphNode>> m_blocks;
    };
#include "implem/Proto_BlockGraphImplem.H"
} // end namespace Proto
#endif // end include guard
