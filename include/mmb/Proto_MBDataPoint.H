#pragma once
#ifndef _PROTO_MB_DATA_POINT_
#define _PROTO_MB_DATA_POINT_

#include "Proto.H"
#include "Proto_MBLevelBoxData.H"

namespace Proto
{
    class MBDataPoint
    {
        public:
        inline MBDataPoint() = default;
        inline MBDataPoint(const MBDataPoint& a_rhs) = default;
        inline MBDataPoint(MBIndex a_index, Point a_point,
                const MBDisjointBoxLayout& a_layout);
        inline MBDataPoint(MBIndex a_index, Point a_point,
                const MBDisjointBoxLayout& a_layout, 
                Point a_boundaryDir, BlockIndex a_block);


        inline bool operator==(const MBDataPoint& a_rhs) const;

        inline bool operator<(const MBDataPoint& a_rhs) const;
        inline bool inBoundary() const {return m_inBoundary; }
        inline const Point& point() const { return m_point; }
        inline MBIndex index() const { return MBIndex(m_partition, m_index); }
        inline const BlockIndex& block() const { return m_block; }
        inline const Point& boundaryDir() const { return m_boundaryDir; }
        inline BlockIndex srcBlock() const;
        inline BlockIndex dstBlock() const;
        inline MBDataPoint shift(Point a_shift) const; 
        template<typename T, unsigned int C, MemType MEM, Centering CTR>
            inline BoxData<T, C, MEM>&
            patch(MBLevelBoxData<T, C, MEM, CTR>& a_data) const;
        template<typename T, unsigned int C, MemType MEM, Centering CTR>
            inline BoxData<T, C, MEM>&
            srcPatch(MBLevelBoxData<T, C, MEM, CTR>& a_data) const;

        // MEMBER DATA
        private:
        unsigned int    m_index;
        Point       m_boundaryDir;
        Point       m_point;
        BlockIndex  m_block;
        bool        m_inBoundary;
        std::shared_ptr<MBBoxPartition> m_partition;
    };

   inline std::ostream& operator<<(std::ostream& a_os, const MBDataPoint& a_point)
   {
       a_os << "point: " << a_point.point() << "| srcBlock: " << a_point.srcBlock() << " dstBlock: " << a_point.dstBlock() << " | bdry: " << a_point.boundaryDir();
       return a_os;
   }

#include "implem/Proto_MBDataPointImplem.H"
} //end namespace Proto
#endif // end include guard
