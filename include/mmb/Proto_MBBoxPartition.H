#pragma once
#ifndef _PROTO_MB_BOX_PARTITION_
#define _PROTO_MB_BOX_PARTITION_

#include "Proto_MBProblemDomain.H"
#include <unordered_map>
namespace Proto
{
    struct MBBoxPartition
    {
        public:

        inline MBBoxPartition(const MBProblemDomain& a_domain);
        inline bool defineBlock(unsigned int a_block, std::vector<Point>& a_patches);
        inline void loadBalance();
        inline bool compatible(const MBBoxPartition& a_rhs);
        inline unsigned int numProcs() const;
        inline unsigned int numBoxes() const;
        inline unsigned int numBoxes(unsigned int a_proc) const;
        inline unsigned int procStartIndex(unsigned int a_proc) const;    
        inline unsigned int procEndIndex(unsigned int a_proc) const;    
        inline unsigned int find(unsigned int a_block, Point a_patch);
        inline void print() const;
    
        private:

        inline void assign(
                std::vector<Point>& a_patches,
                unsigned int        a_globalIndex,
                int                 a_proc,
                unsigned int        a_num);

        MBProblemDomain m_patchDomains;
        std::unordered_map<uint64_t, uint64_t> m_indexMap; ///< Maps local to global
        std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int>> m_procMap; ///< Maps proc to global indices
        std::vector<std::shared_ptr<BoxPartition>> m_blockPartitions;
        std::vector<std::tuple<unsigned int, Point, unsigned int>> m_partition;
    };
} // end namespace Proto

#endif //end include guard
