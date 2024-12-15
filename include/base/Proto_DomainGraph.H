#pragma once
#ifndef _PROTO_DOMAIN_GRAPH_
#define _PROTO_DOMAIN_GRAPH_

#include "Proto.H"

namespace Proto {
    
    typedef int BlockIndex;
    
    struct GraphArc {
        
        inline GraphArc(
            BlockIndex a_srcBlock, BlockIndex a_dstBlock, Point a_srcToDst, Point a_dstToSrc) : 
            sourceBlock(a_srcBlock), destBlock(a_dstBlock), sourceToDest(a_srcToDst), destToSource(a_dstToSrc)
            {
                rotation = CoordPermutation::identity();
            }
        
        inline bool operator==(const GraphArc& a_arc) const;

        BlockIndex sourceBlock;
        BlockIndex destBlock;
        Point sourceToDest;
        Point destToSource;
        CoordPermutation rotation;
    };

    struct GraphNode {

        inline GraphNode(BlockIndex a_index);

        inline void addArc(
            BlockIndex a_dstBlock,
            Point a_toDst,
            Point a_fromSrc);

        BlockIndex block;

        private: 
        std::map<Point, std::vector<GraphArc>> arcs;
    };

    class DomainGraph {

        public:
        inline DomainGraph(int a_numBlocks);

        inline void defineBoundary(
            BlockIndex a_srcBlock,
            BlockIndex a_dstBlock,
            Point a_srcToDst,
            CoordPermutation a_rotation);

        inline Point connectivity(
            BlockIndex a_srcBlock,
            BlockIndex a_dstBlock);

        inline CoordPermutation rotation(
            BlockIndex a_srcBlock,
            BlockIndex a_dstBlock);

        inline Point mapDirection(
            Point a_dir,
            BlockIndex a_srcBlock,
            BlockIndex a_dstBlock);

        inline Point mirrorDirection(
            Point a_dir,
            BlockIndex a_srcBlock,
            BlockIndex a_dstBlock);
        
        inline std::vector<BlockIndex> adjacent(
            BlockIndex a_srcBlock,
            Point a_srcToDst);
            
        inline void close();
        private:
        std::vector<std::shared_ptr<GraphNode>> m_nodes;
    };
}
#endif