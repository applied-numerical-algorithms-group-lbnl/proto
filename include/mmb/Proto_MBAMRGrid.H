#pragma once
#ifndef _PROTO_MB_AMR_GRID_
#define _PROTO_MB_AMR_GRID_

#include "Proto_AMRGrid.H"
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

        inline void define(
            const MBProblemDomain&          a_coarseDomain,
            const std::vector<MBPatchID_t>& a_coarsePatches,
            const std::vector<Point>&       a_coarseBoxSizes,
            const std::vector<Point>&       a_refRatios);
        
        inline void define(
            const MBProblemDomain&          a_coarseDomain,
            const std::vector<Point>&       a_coarseBoxSizes,
            const std::vector<Point>&       a_refRatios);

        inline const AMRGrid& getBlock(unsigned int a_block) const;
        inline const MBDisjointBoxLayout& getLevel(unsigned int a_level) const;
        inline const unsigned int numLevels() const { return m_layouts.size(); }
        inline const unsigned int maxLevels() const { return numLevels(); }
        inline const unsigned int numBlocks() const { return m_grids.size(); }
        private:
        
        std::vector<Point>   m_refRatios;
        MBAMRBoxPartition    m_partition; ///
        std::vector<AMRGrid> m_grids;  ///< One grid per block
        std::vector<MBDisjointBoxLayout> m_layouts; ///< One layout per level
    };

#include "implem/Proto_MBAMRGridImplem.H"
} // end namespace Proto
#endif //end include guard
