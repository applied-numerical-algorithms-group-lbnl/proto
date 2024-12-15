#pragma once
#ifndef _PROTO_DOMAIN_
#define _PROTO_DOMAIN_

#include "Proto.H"

namespace Proto
{
    struct PointInfo {
        inline PointInfo(Point a_point, BlockIndex a_block = 0, int a_proc = 0) 
        : m_point(a_point), m_block(a_block), m_proc(a_proc) {}
        
        inline BlockIndex block() const { return m_block; }
        inline Point point() const { return m_point; }
        inline int proc() const { return m_proc; }

        inline bool operator<(const PointInfo& info)
        {
            if (m_proc != info.m_proc)
            {
                return (m_proc < info.m_proc);
            } else if(m_block != info.m_block) {
                return (m_block < info.m_block);
            } else {
                return m_point < info.m_point;
            }
        }

        private: 
        Point m_point;
        BlockIndex m_block;
        int m_proc;
    };
    
    typedef PointInfo PatchID;

    class Domain {

        public:

        inline Domain(std::vector<Box> domainBoxes);
        inline Domain(Box domainBox);
        inline Domain(Box domainBox, int numBlocks);

        inline void defineBoundary(
            BlockIndex srcBlock,
            BlockIndex dstBlock,
            Point srcToDst,
            CoordPermutation rotation);

        inline void close();
        
        /** Periodicity is defined in a direction if a
         * subdomain has a boundary with itself with 
         * a coordinate rotation of identity */
        inline bool isPeriodic(
            Point dir,
            BlockIndex block = 0);

        inline bool isValidTileSize(
            Point size,
            BlockIndex block = 0);

        inline Box box(BlockIndex block = 0);

        inline int numBlocks();

        inline Domain coarsen(int refRatio);
        inline Domain coarsen(Point refRatio);
        inline Domain coarsen(std::vector<Point> refRatios);

        inline Domain refine(int refRatio);
        inline Domain refine(Point refRatio);
        inline Domain refine(std::vector<Point> refRatios);

        inline std::vector<PatchID> patchesOfSize(int boxSize);
        inline std::vector<PatchID> patchesOfSize(Point boxSize);
        inline std::vector<PatchID> patchesOfSize(std::vector<Point> boxSizes);

        private:

        bool m_defined;
        std::vector<Box> m_subdomains;
        std::shared_ptr<DomainGraph> m_graph;
    };
}

#endif