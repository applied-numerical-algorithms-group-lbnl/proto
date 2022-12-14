#pragma once

#ifndef _PROTO_MB_GRAPH_
#define _PROTO_MB_GRAPH_

//#include "Proto.H"

namespace Proto
{
    inline unsigned int codim(Point a_dir);
   
    inline std::vector<Point> codimDirs(unsigned int a_codim);

    inline std::vector<Point> getBasisVectors(const Point& a_pt);

    struct MBGraphArc
    {
        inline MBGraphArc(
                unsigned int a_srcBlock,
                unsigned int a_dstBlock,
                Point a_srcToDst,
                Point a_dstToSrc)
        {
            srcBlock = a_srcBlock;
            dstBlock = a_dstBlock;
            srcToDst = a_srcToDst;
            dstToSrc = a_dstToSrc;
        }

        inline void print() {
            std::cout << "srcBlock: " << srcBlock;
            std::cout << " | dstBlock: " << dstBlock;
            std::cout << " | srcToDst: " << srcToDst;
            std::cout << " | dstToSrc: " << dstToSrc << std::endl;
        }

        unsigned int srcBlock;
        unsigned int dstBlock;
        Point srcToDst;
        Point dstToSrc;
        CoordPermutation R;
    };

    struct MBGraphNode
    {
        inline MBGraphNode(unsigned int a_index);
     
        inline void addArc(unsigned int a_block, Point a_toDst, Point a_fromDst);
        inline CoordPermutation& rotation(Point a_dir, unsigned int a_block);

        unsigned int index;
        mutable std::map<Point, std::vector<MBGraphArc>> arcs;
    };
  
    class MBProblemDomain;

    class MBGraph
    {
        friend class MBProblemDomain;
        public:

        inline MBGraph(unsigned int a_numBlocks);

        inline void defineBoundary(
                unsigned int        a_srcBlock,
                unsigned int        a_dstBlock,
                unsigned int        a_dir,
                Side::LoHiSide      a_side,
                CoordPermutation&   a_rotation);

        inline void defineBoundary(
                unsigned int        a_srcBlock,
                unsigned int        a_dstBlock,
                Point               a_dir,
                CoordPermutation&   a_rotation);

        inline void fixRotations();
        inline void fixRotationsOld();

        inline unsigned int size() const {return m_blocks.size(); }

        inline unsigned int numBoundaries(unsigned int a_blockIndex, Point& a_dir) const;

        inline std::set<Point> fullConnectivity(
                unsigned int a_srcBlock,
                unsigned int a_dstBlock) const;
        
        inline Point connectivity(unsigned int a_srcBlock, unsigned int a_dstBlock) const;

        inline Point reverseArc(
                unsigned int a_srcBlock,
                unsigned int a_dstBlock,
                Point a_dir) const;

        inline Point reverseDir(
                unsigned int a_srcBlock,
                unsigned int a_dstBlock,
                Point a_dir) const;

        inline CoordPermutation& rotation(
                unsigned int a_block,
                Point a_dir) const;
        
        inline CoordPermutation& rotation(
                unsigned int a_srcBlock,
                Point a_dir,
                unsigned int a_dstBlock) const;

        inline unsigned int adjacent(
                unsigned int    a_srcBlock,
                unsigned int    a_dir,
                Side::LoHiSide  a_side) const;

        inline unsigned int adjacent(
                unsigned int    a_srcBlock,
                Point           a_dir) const;

        inline std::vector<MBGraphArc> boundaries(
                unsigned int a_srcBlock,
                Point        a_dir) const;

        inline bool operator==(const MBGraph& a_rhs) const;

        inline void print() const;
        inline void close();
        private:
        bool m_closed;
        std::vector<std::shared_ptr<MBGraphNode>> m_blocks;
    };
#include "implem/Proto_MBGraphImplem.H"
} // end namespace Proto
#endif // end include guard
