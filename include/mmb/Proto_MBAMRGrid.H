#pragma once
#ifndef _PROTO_MB_AMR_GRID_
#define _PROTO_MB_AMR_GRID_

#include "Proto.H"
#include "Proto_MBAMRBoxPartition.H"
#include "Proto_MBDisjointBoxLayout.H"

namespace Proto
{
    /// Multiblock AMR Grid
    class MBAMRGrid
    {
        public:

        inline MBAMRGrid(
            const MBProblemDomain&          a_coarseDomain,
            const std::vector<MBPatchID_t>& a_coarsePatches,
            const std::vector<Point>&       a_coarseBoxSizes,
            const std::vector<Point>&       a_refRatios);
        
        inline MBAMRGrid(
            const MBProblemDomain&          a_coarseDomain,
            const std::vector<Point>&       a_coarseBoxSizes,
            const std::vector<Point>&       a_refRatios);

        inline MBAMRGrid(
                const MBProblemDomain&      a_coarseDomain,
                Point                       a_coarseBoxSize,
                Point                       a_refRatio,
                int                         a_numLevels);

        inline void define(
            const MBProblemDomain&          a_coarseDomain,
            const std::vector<MBPatchID_t>& a_coarsePatches,
            const std::vector<Point>&       a_coarseBoxSizes,
            const std::vector<Point>&       a_refRatios);
        
        inline void define(
            const MBProblemDomain&          a_coarseDomain,
            const std::vector<Point>&       a_coarseBoxSizes,
            const std::vector<Point>&       a_refRatios);

        inline void setPatches(
            const std::vector<MBPatchID_t>& a_coarsePatches);
        
        inline void setPatches(
            int a_level,
            const std::vector<MBPatchID_t>& a_coarsePatches);
        
        inline const AMRGrid& getBlock(BlockIndex a_block) const;
        inline AMRGrid& getBlock(BlockIndex a_block);
        inline const MBDisjointBoxLayout& getLevel(int a_level) const;
        inline MBDisjointBoxLayout& getLevel(int a_level);
        inline const MBDisjointBoxLayout& operator[](int a_level) const {return getLevel(a_level); }
        inline MBDisjointBoxLayout& operator[](int a_level) {return getLevel(a_level); }
        inline unsigned int numLevels() const { return m_layouts.size(); }
        inline unsigned int maxLevels() const { return numLevels(); }
        inline unsigned int numBlocks() const { return m_grids.size(); }
        inline std::vector<Point> refRatios() const { return m_refRatios; } 
        inline MBIndex find(Point a_point, BlockIndex a_block, int a_level) const;
        inline MBIterator begin(int a_level) const;
        inline MBIterator end(int a_level) const;
        private:
        
        std::vector<Point>                  m_refRatios;
        std::shared_ptr<MBAMRBoxPartition>  m_partition; 
        std::vector<AMRGrid>                m_grids;  ///< One grid per block
        std::vector<MBDisjointBoxLayout>    m_layouts; ///< One layout per level
    };

#include "implem/Proto_MBAMRGridImplem.H"
} // end namespace Proto
#endif //end include guard
