#pragma once
#ifndef _PROTO_MB_PROBLEM_DOMAIN_
#define _PROTO_MB_PROBLEM_DOMAIN_

#include "Proto.H"
#include "Proto_MBGraph.H"

namespace Proto
{
    struct MBPatchID
    {
        PatchID patch;
        BlockIndex block;

        inline MBPatchID(PatchID a_patch, BlockIndex a_block) : patch(a_patch), block(a_block) {}
        inline MBPatchID(const MBPatchID& a_rhs) : patch(a_rhs.patch), block(a_rhs.block) {}
        MBPatchID& operator=(const MBPatchID& a_rhs) = default;
        inline bool operator<(const MBPatchID& a_rhs) const
        {
            if (block == a_rhs.block)
            {
                return (patch < a_rhs.patch);
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
        
        inline std::vector<MBPatchID> patches(int a_boxSize) const;
        inline std::vector<MBPatchID> patches(Point a_boxSize) const;
        inline std::vector<MBPatchID> patches(const std::vector<Point>& a_boxSizes) const;

        inline Point getDirectionFromPatch(MBPatchID a_patch) const;
        inline bool inTriplePoint(MBPatchID a_patch) const;
        inline bool inBlockBoundary(MBPatchID a_patch) const;
        inline bool inDomainBoundary(MBPatchID a_patch) const;
        inline bool inInterior(MBPatchID a_patch) const;
        inline bool onBlockBoundary(PatchID a_patch, Point a_dir) const;
        inline bool onDomainBoundary(PatchID a_patch, Point a_dir) const;

        private:

        mutable bool m_closed;

        std::vector<ProblemDomain>  m_domains;
        std::shared_ptr<MBGraph> m_graph;

        std::vector<bool> m_domainsDefined;

    }; // end class MBProblemDomain
#include "implem/Proto_MBProblemDomainImplem.H"
}

#endif // end include guard
