#pragma once
#ifndef _PROTO_MB_BOX_LAYOUT_
#define _PROTO_MB_BOX_LAYOUT_

#include "Proto.H"
#include "Proto_MBProblemDomain.H"
#include "Proto_MBDataIterator.H"
#include "Proto_MBBoxPartition.H"

namespace Proto
{
    class MBDisjointBoxLayout
    {
        public:
        MBDisjointBoxLayout(MBProblemDomain& a_domain, Point a_boxSize);
        MBDisjointBoxLayout(MBProblemDomain& a_domain, std::vector<Point> a_boxSize);

        inline Box operator[](const MBDataIndex& a_index);
        inline MBDataIterator& begin();
        inline MBDataIterator& end();

        private:
        std::vector<DisjointBoxLayout>  m_layouts;
        std::shared_ptr<MBBoxPartition> m_partition;
    };
#include "implem/Proto_MBDisjointBoxLayoutImplem.H"
} // end namespace Proto
#endif //end include guard
