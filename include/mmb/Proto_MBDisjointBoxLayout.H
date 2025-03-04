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

        inline MBDisjointBoxLayout(){m_defined = false;}

        inline MBDisjointBoxLayout(
                MBProblemDomain& a_domain,
                Point a_boxSize);

        inline MBDisjointBoxLayout(
                MBProblemDomain& a_domain,
                std::vector<Point>               a_boxSizes);
        
        inline MBDisjointBoxLayout(
                MBProblemDomain& a_domain,
                std::vector<MBPoint>         a_patches,
                std::vector<Point>               a_boxSizes);
        
        inline MBDisjointBoxLayout(
                std::shared_ptr<MBBoxPartition>  a_partition,
                std::vector<Point>               a_boxSizes);
        
        inline void define(
                MBProblemDomain& a_domain,
                std::vector<Point>               a_boxSizes);
        
        inline void define(
                MBProblemDomain& a_domain,
                std::vector<MBPoint>         a_patches,
                std::vector<Point>               a_boxSizes);

        inline void define(
                std::shared_ptr<MBBoxPartition>  a_partition,
                std::vector<Point>               a_boxSizes);

        inline const MBProblemDomain& domainPtr() const; 
        inline const MBProblemDomain& domain() const; 
        inline unsigned int procID(const MBIndex& a_index) const;
        inline unsigned int numBlocks() const;
        inline unsigned int numBoxes() const;
        inline unsigned int numBoxes(unsigned int a_proc) const;
        inline unsigned int numProcs() const;
        inline std::vector<MBPoint> patches() const;
        inline DataIndex<BoxPartition> blockIndex(const MBIndex& a_index) const;
        inline Box operator[](const MBIndex& a_index) const;
        inline MBIterator begin() const;
        inline MBIterator end() const;
        inline Point point(const MBIndex& a_index) const; 
        inline Box box(const MBIndex& a_index) const; 
        inline Array<double, DIM> dx(BlockIndex a_block) const;
        inline Array<double, DIM> dx(MBIndex a_index) const;
        inline Box blockDomainBox(BlockIndex a_block) const { return getBlock(a_block).domain().box(); }
        inline unsigned int block(const MBIndex& a_index) const; 
        inline MBPoint patch(const MBIndex& a_index) const;
        inline DisjointBoxLayout& getBlock(BlockIndex a_block);
        inline const DisjointBoxLayout& getBlock(BlockIndex a_block) const;
        inline MBIndex find(Point a_point, BlockIndex a_block) const;
        inline const MBBoxPartition& partition() const { return *m_partition; }
        inline const MBProblemDomain& patchDomain() const { return partition().domain(); }
        inline std::shared_ptr<MBBoxPartition> partitionPtr() const { return m_partition; }
        inline size_t offset(unsigned int a_proc) const;
        inline bool compatible(const MBDisjointBoxLayout& a_rhs) const;
        inline bool compatible(const MBIndex& a_index) const;
        inline Point connectivity(
            const MBPoint&  patch_1,
            const MBPoint&  patch_2) const;
        inline Point connectivity(
            const MBIndex& i1,
            const MBIndex& i2) const;
        inline bool isBlockBoundary(
            const MBPoint&  patch_1,
            const MBPoint&  patch_2) const;
        inline bool isBlockBoundary(
            const MBIndex&  patch_1,
            const MBIndex&  patch_2) const;
        inline bool isBlockBoundary(
            const MBPoint&  patch,
            Point           dir) const;
        inline bool isBlockBoundary(
            const MBIndex&  patch,
            Point           dir) const;
        inline bool isBlockBoundary(
            const MBPoint&  patch,
            Point           dir,
            BlockIndex      adjBlock) const;
        inline bool isPatchOnBlockBoundary( const MBPoint& patch) const;
        inline bool isPatchOnBlockBoundary(
            const MBPoint& patch,
            Point dir) const;
        inline bool isPatchInBlockBoundary( const MBPoint& patch) const;
        inline bool isPatchInBlockBoundary( const MBPoint& patch, Point dir) const;
        inline bool isDomainBoundary(
            const MBPoint&  patch,
            Point           dir) const;
        inline bool isDomainBoundary(
            const MBIndex&  patch,
            Point           dir) const;
        inline bool isPatchOnDomainBoundary( const MBPoint& patch) const;
        inline bool isPatchOnDomainBoundary(
            const MBPoint& patch,
            Point dir) const;
        inline bool isPatchInDomainBoundary( const MBPoint& patch) const;
        inline bool isPatchInDomainBoundary(
            const MBPoint& patch,
            Point dir) const;
        inline bool isRefinementBoundary(
            const MBPoint&  patch,
            Point           dir) const;
        inline bool isPatchOnRefinementBoundary( const MBPoint& patch) const;
        inline bool isInteriorBoundary(
            const MBPoint&  patch,
            Point           dir) const;
        inline bool isPatchInTriplePointRegion(const MBPoint& patch) const;
        inline bool isPatchInTriplePointRegion(
            const MBPoint& patch,
            Point dir) const;
        inline std::vector<MBIndex> adjacentIndices( const MBPoint& patch) const;
        inline std::vector<MBIndex> adjacentIndices(
            const MBPoint& patch,
            Point dir) const;
        inline unsigned int size() const;
        inline unsigned int localSize() const;
        inline std::vector<Point> boxSizes() const;
        inline Point boxSize(BlockIndex a_block) const;
        MBDisjointBoxLayout coarsen(const Point& a_refRatio) const;
        MBDisjointBoxLayout coarsen(const std::vector<Point>& a_refRatios) const;
        void validateSelf(std::string methodName) const;
        void validateIndex(const MBIndex& index, std::string methodName) const;
        void validateBlock(const BlockIndex& block, std::string methodName) const; 
        inline void print() const;
        inline bool defined() const {return m_defined; }
        private:
        bool                            m_defined;
        MBProblemDomain                 m_domain;
        std::vector<DisjointBoxLayout>  m_layouts;
        std::shared_ptr<MBBoxPartition> m_partition;
    };

    inline std::ostream& operator<< (std::ostream& os, const DataIndex<MBBoxPartition>& a_index);
#include "implem/Proto_MBDisjointBoxLayoutImplem.H"
} // end namespace Proto
#endif //end include guard
