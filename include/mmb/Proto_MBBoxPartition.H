#pragma once
#ifndef _PROTO_MB_BOX_PARTITION_
#define _PROTO_MB_BOX_PARTITION_

#include "Proto_MBProblemDomain.H"
#include <unordered_map>
#include <unordered_set>
namespace Proto
{
   
    struct PatchBound
    {
        uint64_t index;
        Point dir;

        bool operator==(const PatchBound& a_rhs) const
        {
            if (index != a_rhs.index) { return false; }
            if (dir != a_rhs.dir) { return false; }
            return true;
        }
    
        struct Hash
        {
            std::size_t operator()(const PatchBound& key) const {

                std::size_t hash = std::hash<uint64_t>{}(key.index);
                for (int ii = 0; ii < DIM; ++ii) {
                    hash ^= std::hash<int>{}(key.dir[ii]) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
                }
                return hash;
            }
        };
    };

    class MBBoxPartition
    {
        public:

        //inline MBBoxPartition(){}
        inline MBBoxPartition(
                const MBProblemDomain& a_domain);
        inline MBBoxPartition(
                const MBProblemDomain& a_domain,
                const std::vector<MBPatchID>& a_patches);
        inline void define(
                const MBProblemDomain& a_domain,
                const std::vector<MBPatchID>& a_patches);
        inline void define(
                const std::vector<MBPatchID>& a_patches);
        inline void loadBalance();
        inline void findBoundaries();
        inline bool compatible(const MBBoxPartition& a_rhs);
        const MBProblemDomain& domain() const {return m_patchDomain; }
        inline unsigned int numProcs() const;
        inline unsigned int numBoxes() const;
        inline unsigned int numBoxes(unsigned int a_proc);
        inline BlockIndex block(uint64_t a_globalIndex) const;
        inline BlockIndex numBlocks() const;
        inline uint64_t procStartIndex(unsigned int a_proc) const;    
        inline uint64_t procEndIndex(unsigned int a_proc) const;    
        inline uint64_t localIndex(uint64_t a_globalIndex, unsigned int a_proc) const;
        inline uint64_t globalIndex(uint64_t a_localIndex, unsigned int a_proc) const;
        inline uint64_t find(Point a_patch, BlockIndex a_block) const;
        inline bool contains(Point a_patch, BlockIndex a_block) const;
        
        inline bool isBlockBoundary(
                uint64_t  a_index,
                Point           a_dir) const;
        inline bool isBlockBoundary(
                uint64_t  a_index,
                Point           a_dir,
                BlockIndex    a_block) const;
        inline bool isDomainBoundary(
                uint64_t  a_index,
                Point           a_dir) const;
        inline bool isInteriorBoundary(
                uint64_t  a_index,
                Point           a_dir) const;

        inline const std::vector<std::pair<MBPatchID, unsigned int>>& partition() const
        { return m_partition; }
        inline std::shared_ptr<BoxPartition> blockPartition(BlockIndex a_block) const;
        inline void print();
    
        private:

        inline void assign(
                std::vector<Point>& a_patches,
                unsigned int        a_globalIndex,
                int                 a_proc,
                unsigned int        a_num);

        inline void buildLocalMaps();
        inline void buildGlobalMaps();

        MBProblemDomain m_patchDomain;
        mutable std::map<MBPatchID, uint64_t> m_indexMap; ///< Maps local to global
        mutable std::unordered_map<unsigned int, std::pair<uint64_t, uint64_t>> m_procMap; ///< Maps proc to global indices
        std::vector<std::shared_ptr<BoxPartition>> m_blockPartitions;
        std::vector<std::pair<MBPatchID, unsigned int>> m_partition;
        mutable std::unordered_map<PatchBound, std::vector<uint64_t>, PatchBound::Hash> m_blockBoundaries;
        mutable std::unordered_set<PatchBound, PatchBound::Hash> m_domainBoundaries;
    };
#include "implem/Proto_MBBoxPartitionImplem.H"
} // end namespace Proto

#endif //end include guard
