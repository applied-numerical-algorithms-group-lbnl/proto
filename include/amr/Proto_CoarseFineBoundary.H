#pragma once
#ifndef _PROTO_COARSE_FINE_BOUNDARY_
#define _PROTO_COARSE_FINE_BOUNDARY_

#include "Proto.H"

namespace Proto
{

    template < typename P>
    class CoarseFineBoundary
    {
        public:
            CoarseFineBoundary(std::shared_ptr<P> coarsePartition, std::shared_ptr<P> finePartition);
            
            bool onCoarseBoundary(DataIndex<P> coarseIndex) const;
            bool onFineBoundary(DataIndex<P> fineIndex) const;
            std::vector<DataIndex<P>> getCoarseFromFine(DataIndex<P> fineIndex) const;
            std::vector<DataIndex<P>> getFineFromCoarse(DataIndex<P> coarseIndex) const;
        private:
            std::shared_ptr<P> m_coarsePartition;
            std::shared_ptr<P> m_finePartition;
            std::map<DataIndex<P>, std::vector<DataIndex<P>>> m_coarseMap;
            std::map<DataIndex<P>, std::vector<DataIndex<P>>> m_fineMap;

    };

    template <typename P>
    CoarseFineBoundary<P>::CoarseFineBoundary(
        std::shared_ptr<P> coarsePartition,
        std::shared_ptr<P> finePartition)
    {
        m_coarsePartition = coarsePartition;
        m_finePartition = finePartition;

        DataIterator<P> coarseIter(coarsePartition);
        DataIterator<P> fineIter(finePartition);

        for (auto iter = fineIter.begin(); iter != fineIter.end(); ++fineIter)
        {
            DataIndex<P> fineIndex = *iter;
            DataIndex<P> coarseIndex = finePartition->getCoarseIndexFromFine(fineIndex, coarsePartition);
            PROTO_ASSERT(coarseIndex != *coarseIter.end(),
                "CoarseFineBoundary::constructor | Error: Proper Nesting Violation");
            
            for (Point dir : Point::DirectionsOfCodim(1))
            {
                if (finePartition->isInteriorBoundary(fineIndex, dir))
                {
                    DataIndex<P> adjFineIndex = finePartition->getAdjacentIndex(fineIndex);
                    if (adjFineIndex == *fineIter.end())
                    {
                        DataIndex<P> adjCoarseIndex = finePartition->getCoarseIndexFromFine(fineIndex, coarsePartition, dir);
                        PROTO_ASSERT(adjCoarseIndex != *coarseIter.end(),
                            "CoarseFineBoundary::constructor | Error: Proper Nesting Violation");
                        m_fineMap[fineIndex].push_back(adjCoarseIndex);
                        m_coarseMap[adjCoarseIndex].push_back(fineIndex);
                    }
                }
            }

        }
    }

    template<typename P>
    bool CoarseFineBoundary<P>::onCoarseBoundary(DataIndex<P> coarseIndex) const
    {
        PROTO_ASSERT(coarseIndex.compatible(*m_coarsePartition),
            "CoarseFineBoundary::onCoarseBoundary | Error: Incompatible Index");
        bool exists = (m_coarseMap.find(coarseIndex) != m_coarseMap.end());
        return exists;
    }

    template<typename P>
    bool CoarseFineBoundary<P>::onFineBoundary(DataIndex<P> fineIndex) const
    {
        PROTO_ASSERT(fineIndex.compatible(*m_finePartition),
            "CoarseFineBoundary::onFineBoundary | Error: Incompatible Index");
        bool exists = (m_fineMap.find(fineIndex) != m_fineMap.end());
        return exists;
    }

    template<typename P>
    std::vector<DataIndex<P>> CoarseFineBoundary<P>::getCoarseFromFine(DataIndex<P> fineIndex) const
    {
        if (onFineBoundary(fineIndex))
        {
            return m_fineMap[fineIndex];
        }
        return std::vector<DataIndex<P>>();
    }
    
    template<typename P>
    std::vector<DataIndex<P>> CoarseFineBoundary<P>::getFineFromCoarse(DataIndex<P> coarseIndex) const
    {
        if (onCoarseBoundary(coarseIndex))
        {
            return m_coarseMap[coarseIndex];
        }
        return std::vector<DataIndex<P>>();
    }
}
#endif