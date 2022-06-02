#pragma once
#ifndef _PROTO_BOX_PARTITION_
#define _PROTO_BOX_PARTITION_

#include "Proto_Morton.H"

namespace Proto
{
    class BoxPartition
    {
        
        inline BoxPartition(
                const ProblemDomain& a_patchDomain,
                const std::vector<Point>& a_patches);
        
        inline BoxPartition(
                const ProblemDomain& a_patchDomain,
                const std::vector<Point>& a_patches,
                unsigned int a_startProc,
                unsigned int a_endProc);
        
        inline void define(
                const ProblemDomain& a_patchDomain,
                const std::vector<Point>& a_patches,
                unsigned int a_startProc,
                unsigned int a_endProc);
        
        inline void loadBalance(
                std::vector<Point>& a_patches,
                unsigned int a_startProc,
                unsigned int a_endProc);
        
        inline bool compatible(const BoxPartition& a_rhs);
   
        inline unsigned int startProc() { return m_startProc; }
        inline const std::vector<std::pair<Point, unsigned int>>& partition() { return m_partition; }

        private:
        ProblemDomain m_patchDomain; ///< Domain in patch space
        std::unordered_map<uint64_t, int> m_indexMap; ///< Maps Morton index to global index
        unsigned int m_startProc; ///< The first proc on which this partition has data
        std::vector<int> m_procIndices; ///< Maps processor number to global index space
        std::vector<std::pair<Point, unsigned int>> m_partition; ///< Maps each patch to a proc
    };

#include "implem/Proto_BoxPartitionImplem.H"
} // end namespace Proto
#endif //end include guard
