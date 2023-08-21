#pragma once
#ifndef _PROTO_MB_AMR_BOX_PARTITION_
#define _PROTO_MB_AMR_BOX_PARTITION_

#include "Proto_MBBoxPartition.H"

namespace Proto
{
    /// Multiblock AMR Box Partition
    class MBAMRBoxPartition
    {
        public:

        inline MBAMRBoxPartition(
                const MBProblemDomain& a_coarseDomain,
                std::vector<Point>&    a_refRatios);
        inline void define(
                const MBProblemDomain& a_coarseDomain,
                std::vector<Point>&    a_refRatios);
        inline void loadBalance();
        inline bool compatible(const MBAMRBoxPartition& a_rhs) const;
        inline std::shared_ptr<BoxPartition> blockPartition(
                unsigned int a_level,
                unsigned int a_block) const;
        inline std::shared_ptr<MBBoxPartition> levelPartition(
                unsigned int a_level) const;
        inline unsigned int numLevels() { return m_levelPartitions.size(); }
        inline void print();

        private:
        
        std::vector<std::shared_ptr<MBBoxPartition>> m_levelPartitions;
    };

#include "implem/Proto_MBAMRBoxPartitionImplem.H"
} // end namespace Proto
#endif //end include guard
