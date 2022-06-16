#pragma once
#ifndef _PROTO_MB_BOX_LAYOUT_
#define _PROTO_MB_BOX_LAYOUT_

#include "Proto.H"
#include "Proto_MBProblemDomain.H"
#include "Proto_MBBoxPartition.H"

namespace Proto
{
    class MBDisjointBoxLayout
    {
        typedef DataIterator<MBBoxPartition> iterator_t;
        typedef DataIndex<MBBoxPartition> index_t;

        public:
        inline MBDisjointBoxLayout(
                const MBProblemDomain& a_domain,
                Point a_boxSize);
        inline MBDisjointBoxLayout(
                const MBProblemDomain& a_domain,
                const std::vector<Point>& a_boxSizes);
        inline MBDisjointBoxLayout(
                const MBProblemDomain& a_domain,
                const std::vector<MBPatchID_t>& a_patches,
                const std::vector<Point>& a_boxSizes);
        
        inline void define(
                const MBProblemDomain& a_domain,
                const std::vector<Point>& a_boxSizes);
        inline void define(
                const MBProblemDomain& a_domain,
                const std::vector<MBPatchID_t>& a_patches,
                const std::vector<Point>& a_boxSizes);
        
        inline Box operator[](const index_t& a_index);
        inline iterator_t begin();
        inline iterator_t end();

        private:
        std::vector<DisjointBoxLayout>  m_layouts;
        std::shared_ptr<MBBoxPartition> m_partition;
    };
#include "implem/Proto_MBDisjointBoxLayoutImplem.H"
} // end namespace Proto
#endif //end include guard
