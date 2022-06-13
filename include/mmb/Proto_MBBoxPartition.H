#pragma once
#ifndef _PROTO_MB_BOX_PARTITION_
#define _PROTO_MB_BOX_PARTITION_

#include "Proto_MBProblemDomain.H"
#include <unordered_map>
namespace Proto
{
    typedef std::pair<patchID_t, blockID_t> MBPatchID_t;
    typedef unsigned int blockID_t;
    
    class MBBoxPartition
    {
        typedef Point patchID_t;
        typedef unsigned int procID_t;
        typedef uint64_t globalIndex_t;
        
        public:

        inline MBBoxPartition(
                const MBProblemDomain& a_domain,
                const std::vector<MBPatchID_t>& a_patches);
        inline void loadBalance(const std::vector<MBPatchID_t>& a_patches);
        inline bool compatible(const MBBoxPartition& a_rhs);
        inline unsigned int numProcs() const;
        inline unsigned int numBoxes() const;
        inline unsigned int numBoxes(procID_t a_proc) const;
        inline unsigned int procStartIndex(procID_t a_proc) const;    
        inline unsigned int procEndIndex(procID_t a_proc) const;    
        inline unsigned int find(blockID_t a_block, patchID_t a_patch);
        inline const std::vector<std::pair<MBPatchID_t, procID_t>>& partition() const { return m_partition; }
        inline std::shared_ptr<BoxPartition> blockPartition(blockID_t a_block) const;
        inline void print() const;
    
        private:

        inline void assign(
                std::vector<Point>& a_patches,
                unsigned int        a_globalIndex,
                int                 a_proc,
                unsigned int        a_num);

        inline void buildLocalMaps();
        inline void buildGlobalMaps();

        MBProblemDomain m_patchDomains;
        std::unordered_map<MBPatchID, globalIndex_t> m_indexMap; ///< Maps local to global
        std::unordered_map<procID_t, std::pair<globalIndex_t, globalIndex_t>> m_procMap; ///< Maps proc to global indices
        std::vector<std::shared_ptr<BoxPartition>> m_blockPartitions;
        std::vector<std::pair<MBPatchID_t, procID_t>> m_partition;
    };
} // end namespace Proto

#endif //end include guard
