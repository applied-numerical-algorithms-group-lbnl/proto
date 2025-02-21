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
        inline bool operator==(const MBPoint a_rhs) const{
            return point == a_rhs.point && block == a_rhs.block;
        }

        struct Hash
        {
            std::size_t operator()(const MBPoint &key) const
            {
                std::size_t seed = 0;

                // Hash combine function
                auto hash_combine = [](std::size_t& seed, std::size_t value) {
                    seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                };
        
                // Hash the vector elements
                for (int dd = 0; dd < DIM; dd++) {
                    hash_combine(seed, std::hash<int>{}(key.point[dd]));
                }
        
                // Hash the unsigned int
                hash_combine(seed, std::hash<int>{}(key.block));
        
                return seed;
            }
        };
        
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

        inline bool isPointInTriplePoint(Point a_point, BlockIndex a_block) const;
        inline bool isPointInTriplePoint(Point a_point, BlockIndex a_block, Point a_dir) const;
        
        inline bool isPointInBlockBoundary(Point a_point, BlockIndex a_block) const;
        inline bool isPointInBlockBoundary(Point a_point, BlockIndex a_block, Point a_dir) const;
        
        inline bool isPointInDomainBoundary(Point a_point, BlockIndex a_block) const;
        inline bool isPointInDomainBoundary(Point a_point, BlockIndex a_block, Point a_dir) const;
        
        inline bool isPointInInterior(Point a_point, BlockIndex a_block) const;
        
        inline bool isPointOnBlockBoundary(Point a_point, BlockIndex a_block) const;
        inline bool isPointOnBlockBoundary(Point a_point, BlockIndex a_block, Point a_dir) const;
        
        inline bool isPointOnDomainBoundary(Point a_point, BlockIndex a_block) const;
        inline bool isPointOnDomainBoundary(Point a_point, BlockIndex a_block, Point a_dir) const;

        private:

        mutable bool m_closed;

        std::vector<ProblemDomain>  m_domains;
        std::shared_ptr<MBGraph> m_graph;

        std::vector<bool> m_domainsDefined;

    }; // end class MBProblemDomain
#include "implem/Proto_MBProblemDomainImplem.H"
}

#endif // end include guard
