#pragma once
#ifndef _PROTO_MB_PROBLEM_DOMAIN_
#define _PROTO_MB_PROBLEM_DOMAIN_

#include "Proto.H"
#include "Proto_MBGraph.H"

namespace Proto
{
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
       
        // Set Size of Specified Block
        inline void defineDomain(unsigned int a_blockID, Point a_domainSize);
        
        // Add a Codimension-1 Boundary to the Graph
        inline void defineBoundary(
                unsigned int        a_srcBlock,
                unsigned int        a_dstBlock,
                unsigned int        a_dir,
                Side::LoHiSide      a_side,
                CoordPermutation&   a_rotation);
        
        // Add a Codimension-1 Boundary to the Graph
        inline void defineBoundary(
                unsigned int        a_srcBlock,
                unsigned int        a_dstBlock,
                Point               a_dir,
                CoordPermutation&   a_rotation);

        // Number of Blocks
        inline unsigned int size() const { return m_graph->size(); }
        
        // Number of Blocks
        inline unsigned int numBlocks() const { return size(); }
        
        // Get Graph Reference
        inline MBGraph& graph() const { return *m_graph; }
        
        // Get Local Problem Domain
        inline const ProblemDomain& blockDomain(unsigned int a_blockID) const;
        
        // Get Size of Block
        inline Point blockSize(unsigned int a_blockID) const;
        inline std::vector<MBGraphArc> boundaries(
                unsigned int a_srcBlock,
                Point        a_dir) const;
        
        // Convert a Point Between Block Coordinate Systems
        inline Point convert(
                Point a_point,
                unsigned int a_srcBlock,
                unsigned int a_dstBlock) const;

        // Convert a Box Between Block Coordinate Systems
        inline Box convert(
                Box a_box,
                unsigned int a_srcBlock,
                unsigned int a_dstBlock) const;

        inline MBProblemDomain coarsen(Point a_refRatio) const;
        inline MBProblemDomain coarsen(std::vector<Point> a_refRatios) const;
        inline MBProblemDomain refine(Point a_refRatio) const;

        inline bool coarsenable(Point a_refRatio) const;
        inline bool coarsenable(unsigned int a_block, Point a_refRatio) const;
        inline bool defined();
        inline bool operator==(const MBProblemDomain& a_rhs) const;
        inline void close();
        private:

        mutable bool m_closed;

        std::vector<ProblemDomain>  m_domains;
        std::shared_ptr<MBGraph> m_graph;

        std::vector<bool> m_domainsDefined;

    }; // end class MBProblemDomain
#include "implem/Proto_MBProblemDomainImplem.H"
}

#endif // end include guard
