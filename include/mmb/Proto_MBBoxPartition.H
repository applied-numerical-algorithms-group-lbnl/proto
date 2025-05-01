#pragma once
#ifndef _PROTO_MB_BOX_PARTITION_
#define _PROTO_MB_BOX_PARTITION_

#include "Proto_MBProblemDomain.H"
#include <unordered_map>
#include <unordered_set>

namespace Proto
{
    namespace
    {
        struct PatchBound
        {
            uint64_t index;
            Point dir;

            bool operator==(const PatchBound &a_rhs) const
            {
                if (index != a_rhs.index)
                {
                    return false;
                }
                if (dir != a_rhs.dir)
                {
                    return false;
                }
                return true;
            }

            struct Hash
            {
                std::size_t operator()(const PatchBound &key) const
                {

                    std::size_t hash = std::hash<uint64_t>{}(key.index);
                    for (int ii = 0; ii < DIM; ++ii)
                    {
                        hash ^= std::hash<int>{}(key.dir[ii]) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
                    }
                    return hash;
                }
            };
        };
        
        enum BoundaryType 
        {
            PR_DOMAIN_BOUNDARY,
            PR_BLOCK_BOUNDARY,
            PR_INTERIOR_BOUNDARY
        };
        
        struct MBPatchBoundInfo {
            MBPatchBoundInfo()
            {
                adjacentPatches.clear();
                boundaryType = PR_INTERIOR_BOUNDARY;
                isRefinementBoundary = false;
            }
            std::set<MBPoint> adjacentPatches;
            BoundaryType boundaryType;
            bool isRefinementBoundary;
        };
        struct MBPatchInfo {
            MBPatchInfo() : patch(Point::Zeros(), 0) {}
            MBPatchInfo(const MBPoint& _patch, bool _exists)
            : patch(_patch), exists(_exists)
            {
                boundaries.clear();
                for (auto dir : Point::Directions())
                {
                    boundaries[dir] = MBPatchBoundInfo();
                }
            }

            MBPoint patch;
            bool exists;
            std::map<Point, MBPatchBoundInfo> boundaries;
        };

    }
    class MBBoxPartition
    {
        public:

        inline MBBoxPartition(
                const MBProblemDomain& a_domain);
        inline MBBoxPartition(
                const MBProblemDomain& a_domain,
                const std::vector<MBPoint>& a_patches);
        inline void define(
                const MBProblemDomain& a_domain,
                const std::vector<MBPoint>& a_patches);
        inline void define(
                const std::vector<MBPoint>& a_patches);
        inline void loadBalance();
        inline bool compatible(const MBBoxPartition& a_rhs);
        const MBProblemDomain& domain() const {return m_patchDomain; }
        inline unsigned int numProcs() const;
        inline unsigned int numBoxes() const;
        inline unsigned int numBoxes(unsigned int a_proc);
        inline BlockIndex numBlocks() const;
        inline uint64_t procStartIndex(unsigned int a_proc) const;    
        inline uint64_t procEndIndex(unsigned int a_proc) const;    
        inline uint64_t localIndex(uint64_t a_globalIndex, unsigned int a_proc) const;
        inline uint64_t globalIndex(uint64_t a_localIndex, unsigned int a_proc) const;
        inline uint64_t find(Point a_patch, BlockIndex a_block) const;
        inline bool contains(Point a_patch, BlockIndex a_block) const;
        
        inline void validatePatch(const MBPoint& patch, std::string methodName) const;
        inline void validateIndex(uint64_t index, std::string methodName) const;
        inline Point connectivity(const MBPoint& src, const MBPoint& dst) const;
        unsigned int getProc(uint64_t globalIndex) const;
        inline MBPoint getPatch(uint64_t globalIndex) const;
        inline BlockIndex getBlock(uint64_t globalIndex) const;
        inline Point getPoint(uint64_t globalIndex) const; 
        inline std::set<MBPoint> adjacentPatches(const MBPoint& patch) const;
        inline bool isBlockBoundary(
                MBPoint patch,
                Point dir) const;
        inline bool isBlockBoundary(
                MBPoint patch,
                Point dir,
                BlockIndex adjBlock) const;
        inline bool isDomainBoundary(
                MBPoint patch,
                Point dir) const;
        inline bool isInteriorBoundary(
                MBPoint patch,
                Point dir) const;
        inline bool isRefinementBoundary(
                MBPoint patch,
                Point dir) const;

        inline const std::vector<std::pair<MBPoint, unsigned int>>& partition() const
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
        inline void findAdjacencies();
        inline void getBoundaryInfo(MBPatchInfo& patchInfo);
        inline void updateRefinementBoundaries();

        MBProblemDomain m_patchDomain;
        mutable std::map<MBPoint, uint64_t> m_indexMap; ///< Maps local to global
        mutable std::unordered_map<unsigned int, std::pair<uint64_t, uint64_t>> m_procMap; ///< Maps proc to global indices
        std::vector<std::shared_ptr<BoxPartition>> m_blockPartitions;
        std::vector<std::pair<MBPoint, unsigned int>> m_partition;
        mutable std::unordered_map<MBPoint, MBPatchInfo, MBPoint::Hash> m_patchAdjacencies;

    };
#include "implem/Proto_MBBoxPartitionImplem.H"
} // end namespace Proto

#endif //end include guard
