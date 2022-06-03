#pragma once
#ifndef _PROTO_BOX_PARTITION_
#define _PROTO_BOX_PARTITION_

#include "Proto_Morton.H"
#include <unordered_map>
namespace Proto
{
    class BoxPartition
    {
        public: 
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
   
        inline unsigned int numProcs() const;
        inline unsigned int numBoxes() const;
        inline unsigned int numBoxes(unsigned int a_proc) const;
        inline unsigned int startProc() const { return m_startProc; }
        inline unsigned int endProc() const { return m_startProc + numProcs(); }
        inline const std::vector<std::pair<Point, unsigned int>>& partition() const
        { return m_partition; }
        inline unsigned int procStartIndex(unsigned int a_proc) const;
        inline unsigned int procEndIndex(unsigned int a_proc) const;
        inline unsigned int patchIndex(Point a_pt);
        inline unsigned int find(Point a_pt);
        inline void print() const;
        private:
        ProblemDomain m_patchDomain; ///< Domain in patch space
        std::unordered_map<uint64_t, int> m_indexMap; ///< Maps Morton index to global index
        unsigned int m_startProc; ///< The first proc on which this partition has data
        std::vector<unsigned int> m_procIndices; ///< Maps processor number to global index space
        std::vector<std::pair<Point, unsigned int>> m_partition; ///< Maps each patch to a proc
    };

#include "implem/Proto_BoxPartitionImplem.H"
} // end namespace Proto
#endif //end include guard
