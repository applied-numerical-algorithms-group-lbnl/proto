#pragma once
#ifndef _PROTO_MB_BOX_LAYOUT_
#define _PROTO_MB_BOX_LAYOUT_

#include "Proto.H"
#include "Proto_MBProblemDomain.H"
#include "Proto_MBBoxPartition.H"

namespace Proto
{
    typedef DataIndex<MBBoxPartition> MBIndex;
    typedef DataIterator<MBBoxPartition> MBIterator;

    class MBDisjointBoxLayout
    {
        public:
        typedef DataIndex<MBBoxPartition> MBIndex;

        inline MBDisjointBoxLayout(){}

        inline MBDisjointBoxLayout(
                MBProblemDomain& a_domain,
                Point a_boxSize);

        inline MBDisjointBoxLayout(
                MBProblemDomain& a_domain,
                const std::vector<Point>& a_boxSizes);
        
        inline MBDisjointBoxLayout(
                MBProblemDomain& a_domain,
                const std::vector<MBPatchID_t>& a_patches,
                const std::vector<Point>& a_boxSizes);
        
        inline MBDisjointBoxLayout(
                std::shared_ptr<MBBoxPartition> a_partition,
                const std::vector<Point>& a_boxSizes);
        
        inline void define(
                MBProblemDomain& a_domain,
                const std::vector<Point>& a_boxSizes);
        
        inline void define(
                MBProblemDomain& a_domain,
                const std::vector<MBPatchID_t>& a_patches,
                const std::vector<Point>& a_boxSizes);

        inline void define(
                std::shared_ptr<MBBoxPartition> a_partition,
                const std::vector<Point>& a_boxSizes);

        inline const MBProblemDomain& domain() const; 
        inline unsigned int procID(const MBIndex& a_index) const;
        inline unsigned int numBlocks() const;
        inline unsigned int numBoxes() const;
        inline unsigned int numBoxes(unsigned int a_proc) const;
        inline unsigned int numProcs() const;
        inline std::vector<MBPatchID_t> patches() const;
        inline DataIndex<BoxPartition> blockIndex(const MBIndex& a_index) const;
        inline Box operator[](const MBIndex& a_index) const;
        inline MBIterator begin() const;
        inline MBIterator end() const;
        inline Point point(const MBIndex& a_index) const; 
        inline Box box(const MBIndex& a_index) const; 
        inline unsigned int block(const MBIndex& a_index) const; 
        inline DisjointBoxLayout& getBlock(unsigned int a_block);
        inline const DisjointBoxLayout& getBlock(unsigned int a_block) const;
        inline DisjointBoxLayout levelLayout() const;
        inline MBIndex find(Point a_point, unsigned int a_block) const;
        inline const MBBoxPartition& partition() const { return *m_partition; }
        inline size_t offset(unsigned int a_proc) const;
        inline bool compatible(const MBDisjointBoxLayout& a_rhs) const;
        inline Point connectivity(
                const MBIndex&  a_i1,
                const MBIndex&  a_i2) const;
        inline bool isBlockBoundary(
                const MBIndex&  a_i1,
                const MBIndex&  a_i2) const;
        inline bool isBlockBoundary(
                const MBIndex&  a_index,
                Point           a_dir) const;
        inline bool isBlockBoundary(
                const MBIndex&  a_index,
                Point           a_dir,
                unsigned int    a_block) const;
        inline bool isDomainBoundary(
                const MBIndex&  a_index,
                Point           a_dir) const;
        inline bool isInteriorBoundary(
                const MBIndex&  a_index,
                Point           a_dir) const;
        inline std::vector<MBIndex> bounds(
                const MBIndex&  a_index) const;
        inline unsigned int size() const;
        inline unsigned int localSize() const;
        inline std::vector<Point> boxSizes() const;
        inline MBDisjointBoxLayout coarsen(const Point& a_refRatio) const;
        inline MBDisjointBoxLayout coarsen(const std::vector<Point>& a_refRatios) const;
        inline void print() const;
        private:
        MBProblemDomain                 m_domain;
        std::vector<DisjointBoxLayout>  m_layouts;
        std::shared_ptr<MBBoxPartition> m_partition;
    };

    inline std::ostream& operator<< (std::ostream& os, const DataIndex<MBBoxPartition>& a_index);
#include "implem/Proto_MBDisjointBoxLayoutImplem.H"
} // end namespace Proto
#endif //end include guard
