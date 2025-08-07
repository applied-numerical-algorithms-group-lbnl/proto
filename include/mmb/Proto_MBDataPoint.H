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
        inline MBDataPoint(const MBDataPoint& a_rhs);
        inline MBDataPoint(MBIndex a_dstIndex, Point a_dstPoint,
            const MBDisjointBoxLayout& a_layout);
        inline MBDataPoint(MBIndex a_dstIndex, Point a_dstPoint,
            MBIndex a_srcIndex, Point a_srcPoint,
            const MBDisjointBoxLayout& a_layout);
        
        // Comparitors
        inline MBDataPoint& operator=(const MBDataPoint& a_rhs) = default;
        inline bool operator==(const MBDataPoint& a_rhs) const;
        inline bool operator<(const MBDataPoint& a_rhs) const;

        // Accessors
        inline bool inBoundary() const {return m_boundaryDir != Point::Zeros(); }
        inline const Point& dstPoint() const { return m_dstPoint; }
        inline const MBIndex& dstIndex() const { return m_dstIndex; }
        inline const BlockIndex& dstBlock() const { return m_dstBlock; }
        inline const Point& srcPoint() const { return m_srcPoint; }
        inline const MBIndex& srcIndex() const { return m_srcIndex; }
        inline const BlockIndex& srcBlock() const { return m_srcBlock; }
        inline const Point& boundaryDir() const { return m_boundaryDir; }

        template<typename T, unsigned int C, MemType MEM, Centering CTR>
        inline BoxData<T, C, MEM>&
        dstPatch(MBLevelBoxData<T, C, MEM, CTR>& a_data) const;

        template<typename T, unsigned int C, MemType MEM, Centering CTR>
        inline BoxData<T, C, MEM>&
        srcPatch(MBLevelBoxData<T, C, MEM, CTR>& a_data) const;

        inline MBDataPoint shift(Point a_dstShift, const MBDisjointBoxLayout& a_layout) const; 

        // MEMBER DATA
        private:
        MBIndex     m_dstIndex;
        MBIndex     m_srcIndex;
        Point       m_srcPoint;
        Point       m_dstPoint;

        BlockIndex  m_srcBlock;
        BlockIndex  m_dstBlock;
        Point       m_boundaryDir;
    };

   inline std::ostream& operator<<(std::ostream& a_os, const MBDataPoint& a_point)
   {
        a_os << "dstPoint: " << a_point.dstPoint() << "| dstBlock: " << a_point.dstBlock() << " | srcPoint: " << a_point.srcPoint() << " | srcBlock: " << a_point.srcBlock() << " | boundaryDir: " << a_point.boundaryDir();
        return a_os;
   }

#include "implem/Proto_MBDataPointImplem.H"
} //end namespace Proto
#endif // end include guard
