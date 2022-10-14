#pragma once
#ifndef _PROTO_MB_INTERP_OP_
#define _PROTO_MB_INTERP_OP_

#include "Proto.H"
#include "Proto_MBLevelBoxData.H" //for MBDataPoint definition
#include "Proto_MBMap.H"
#include "Proto_Matrix.H"

namespace Proto
{
    template<typename T>
    class MBInterpOp
    {
        public:

        inline MBInterpOp(
            MBDataPoint a_point,
            MBMap<T>& a_map,
            const std::vector<Point>& a_footprint);

        template<unsigned int C, MemType MEM>
        inline void execute(
                MBLevelBoxData<T, C, MEM>& a_dst,
                MBLevelBoxData<T, C, MEM>& a_src);
        
        inline unsigned int size() const;

        private:

        MBDataPoint m_dst;
        std::set<MBDataPoint> m_srcs;
        std::vector<T> m_coefs;
    };
#include "implem/Proto_MBInterpOpImplem.H"
} //end namespace Proto
#endif // end include guard
