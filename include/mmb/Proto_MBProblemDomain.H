#pragma once
#ifndef _PROTO_MB_PROBLEM_DOMAIN_
#define _PROTO_MB_PROBLEM_DOMAIN_

#include "Proto.H"
#include "Proto_MBGraph.H"

namespace Proto
{
    struct MBPoint
    {
        PatchID point;
        BlockIndex block;

        inline MBPoint(PatchID a_patch, BlockIndex a_block) : point(a_patch), block(a_block) {}
        inline MBPoint(const MBPoint& a_rhs) : point(a_rhs.point), block(a_rhs.block) {}
        MBPoint& operator=(const MBPoint& a_rhs) = default;
        inline bool operator<(const MBPoint& a_rhs) const
        {
            if (block == a_rhs.block)
            {
                return (point < a_rhs.point);
            } else {
                return (block < a_rhs.block);
            }
        }
        
    };
    /// Mapped Multi-Block Problem Domain
    /**
        The ProblemDomain equivalent for a MMB domain. Contains all of the
        normal ProblemDomain objects for each individual block as well
        as the connectivity information between blocks.
    */
    class MBProblemDomain
    {
        public:
       
        inline MBProblemDomain();

        inline MBProblemDomain(unsigned int a_numBlocks);

        inline MBProblemDomain(std::shared_ptr<MBGraph> a_graph);

        // Set Size of Specified Block
        inline void defineDomain(BlockIndex a_blockID, Point a_domainSize);
        
        // Add a Codimension-1 Boundary to the Graph
        inline void defineBoundary(
                BlockIndex        a_srcBlock,
                BlockIndex        a_dstBlock,
                unsigned int        a_dir,
                Side::LoHiSide      a_side,
                CoordPermutation&   a_rotation);
        
        // Add a Codimension-1 Boundary to the Graph
        inline void defineBoundary(
                BlockIndex        a_srcBlock,
                BlockIndex        a_dstBlock,
                Point               a_dir,
                CoordPermutation&   a_rotation);

        // Number of Blocks
        inline unsigned int size() const { return m_graph->size(); }
        
        // Number of Blocks
        inline unsigned int numBlocks() const { return size(); }
        
        // Get Graph Reference
        inline MBGraph& graph() const { return *m_graph; }
       
        // Get Graph Reference
        inline std::shared_ptr<MBGraph> graphPtr() const { return m_graph; }
        
        // Get Local Problem Domain
        inline const ProblemDomain& getBlock(BlockIndex a_blockID) const;
        
        // Get Size of Block
        inline Point blockSize(BlockIndex a_blockID) const;
        
        // Get Differential Volume
        inline double dv(BlockIndex a_blockID) const;


        inline std::vector<MBGraphArc> boundaries(
                BlockIndex a_srcBlock,
                Point        a_dir) const;
        
        // Convert a Point Between Block Coordinate Systems
        inline Point convertNode(
                Point a_point,
                BlockIndex a_srcBlock,
                BlockIndex a_dstBlock) const;

        // Convert a Point Between Block Coordinate Systems
        inline Point convertPoint(
                Point a_point,
                BlockIndex a_srcBlock,
                BlockIndex a_dstBlock,
                Centering a_ctr) const;

        // Convert a Box Between Block Coordinate Systems
        inline Box convertBox(
                Box a_box,
                BlockIndex a_srcBlock,
                BlockIndex a_dstBlock,
                Centering a_ctr = PR_CELL) const;

        inline MBProblemDomain coarsen(Point a_refRatio) const;
        inline MBProblemDomain coarsen(const std::vector<Point>& a_refRatios) const;
        inline MBProblemDomain refine(Point a_refRatio) const;
        inline MBProblemDomain refine(const std::vector<Point>& a_refRatios) const;

        inline bool coarsenable(Point a_refRatio) const;
        inline bool coarsenable(BlockIndex a_block, Point a_refRatio) const;
        inline bool defined() const;
        inline bool operator==(const MBProblemDomain& a_rhs) const;
        inline void close() const;
        
        inline std::vector<MBPoint> patches(int a_boxSize) const;
        inline std::vector<MBPoint> patches(Point a_boxSize) const;
        inline std::vector<MBPoint> patches(const std::vector<Point>& a_boxSizes) const;

        inline Point getDirectionFromPoint(Point a_point, BlockIndex a_block) const;
        inline bool inTriplePoint(Point a_point, BlockIndex a_block) const;
        inline bool inBlockBoundary(Point a_point, BlockIndex a_block) const;
        inline bool inDomainBoundary(Point a_point, BlockIndex a_block) const;
        inline bool inInterior(Point a_point, BlockIndex a_block) const;
        inline bool onBlockBoundary(Point a_patch, BlockIndex a_block) const;
        inline bool onBlockBoundary(Point a_patch, BlockIndex a_block, Point a_dir) const;
        inline bool onDomainBoundary(Point a_patch, BlockIndex a_block) const;
        inline bool onDomainBoundary(Point a_patch, BlockIndex a_block, Point a_dir) const;

        private:

        mutable bool m_closed;

        std::vector<ProblemDomain>  m_domains;
        std::shared_ptr<MBGraph> m_graph;

        std::vector<bool> m_domainsDefined;

    }; // end class MBProblemDomain
#include "implem/Proto_MBProblemDomainImplem.H"
}

#endif // end include guard
