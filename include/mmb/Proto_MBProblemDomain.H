#pragma once
#ifndef _PROTO_MB_PROBLEM_DOMAIN_
#define _PROTO_MB_PROBLEM_DOMAIN_

#include "Proto_ProbemDomain.H"
#include "Proto_BlockGraph.H"

namespace Proto
{
    /// Mapped Multi-Block Problem Domain
    /**
        The ProblemDomain equivalent for a MMB domain. Contains all of the
        normal ProblemDomain objects for each individual block as well
        as the connectivity information between blocks.
    */
    template<typename T>
    class MBProblemDomain
    {
        public:
        inline MBProblemDomain();
        inline MBProblemDomain(unsigned int a_numBlocks);
        inline void define(unsigned int a_numBlocks);
        inline void defineBlock(unsigned int a_blockID, Point a_domainSize, BlockMap<T> a_map);
        
        inline BlockGraph& graph() { return *m_graph; }
        inline const ProblemDomain& domain(unsigned int a_blockID) const { return m_domains[a_blockID]; }
        inline const BlockMap<T>& map(unsigned int a_blockID) const { return m_maps[a_blockID]; }
        
        inline MBProblemDomain coarsen(Point a_refRatio) const;
        inline MBProblemDomain refine(Point a_refRatio) const;

        private:
        inline bool coarsenable(Point a_refRatio);
        inline bool allBlocksDefined();

        std::vector<ProblemDomain>  m_domains;
        std::vector<BlockMap<T>>       m_maps;
        std::shared_ptr<BlockGraph> m_graph;
        

    }; // end class MBProblemDomain
}

#include "implem/Proto_MBProblemDomainImplem.H"
#endif // end include guard
