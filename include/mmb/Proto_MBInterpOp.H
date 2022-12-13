#pragma once
#ifndef _PROTO_MB_INTERP_OP_
#define _PROTO_MB_INTERP_OP_

#include "Proto.H"
#include "Proto_MBLevelBoxData.H" //for MBDataPoint definition
#include "Proto_MBMap.H"
#include "Proto_Matrix.H"
#include "Proto_Operator.H"

namespace Proto
{
    class MBPointInterpOp
    {
        public:

        //TODO: Should be able to deduce the order from the footprint
        template<typename Func>
        inline MBPointInterpOp(
            MBDataPoint                 a_dst,
            Point                       a_ghost,
            MBMap<Func>&                a_map,
            const std::vector<Point>&   a_footprint,
            unsigned int                a_order);

        template<typename T, unsigned int C, MemType MEM>
        inline void apply(
                MBLevelBoxData<T, C, MEM>& a_dst,
                MBLevelBoxData<T, C, MEM>& a_src);
        
        inline unsigned int size() const;

        inline const MBDataPoint& target() const {return m_dst;}

        private:

        MBDataPoint m_dst;
        std::vector<MBDataPoint> m_srcs;
        std::vector<double> m_coefs;
    };

    class MBInterpOp
    {
        public:
        
        template<typename Func>
        inline MBInterpOp(
            MBMap<Func>&                a_map,
            const std::vector<Point>&   a_footprint,
            Point                       a_ghost,
            unsigned int                a_order);

        template<typename T, unsigned int C, MemType MEM>
        inline void apply(
                MBLevelBoxData<T, C, MEM>& a_dst,
                MBLevelBoxData<T, C, MEM>& a_src);
        
        private:
        
        Point m_ghost;
        int m_order;
        std::vector<Point> m_footprint;
        std::vector<MBPointInterpOp> m_ops;
    };
#include "implem/Proto_MBInterpOpImplem.H"
} //end namespace Proto
#endif // end include guard
