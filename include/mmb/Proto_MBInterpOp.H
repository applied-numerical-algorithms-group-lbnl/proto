#pragma once
#ifndef _PROTO_MB_INTERP_OP_
#define _PROTO_MB_INTERP_OP_

#include "Proto.H"
#include "Proto_MBLevelBoxData.H" //for MBDataPoint definition
#include "Proto_MBLevelMap.H"
#include "Proto_Matrix.H"
#include "Proto_Operator.H"

namespace Proto
{
    class MBPointInterpOp
    {
        public:

        //TODO: Should be able to deduce the order from the footprint
        template<MemType MEM>
        inline MBPointInterpOp(
            MBDataPoint                 a_dst,
            Point                       a_ghost,
            MBLevelMap<MEM>&            a_map,
            const std::vector<Point>&   a_footprint,
            unsigned int                a_order);

        template<typename T, unsigned int C, MemType MEM>
        inline void apply(
                MBLevelBoxData<T, C, MEM>& a_dst,
                MBLevelBoxData<T, C, MEM>& a_src);
        
        inline unsigned int size() const;
        
        template<typename T, unsigned int C, MemType MEM>
        inline Matrix<double> coefs(
                MBLevelBoxData<T, C, MEM>& a_src,
                int a_srcComponent = 0);
        
        inline const MBDataPoint& target() const {return m_dst;}
        inline const std::vector<MBDataPoint>& sources() const {return m_srcs;}
        inline std::vector<Point> exponents() const {return m_exponents; }
        inline Matrix<double>& CMatrix() {return m_C; }
        inline Matrix<double>& DMatrix() {return m_D; }
        inline Matrix<double>& SMatrix() {return m_S; }
     
        private:
       
        std::vector<Point> m_exponents;
        MBDataPoint m_dst;
        Matrix<double> m_C;
        Matrix<double> m_D;
        Matrix<double> m_S;
        Matrix<double> m_Cinv;
        std::vector<MBDataPoint> m_srcs;
        std::vector<double> m_coefs;
    };

    class MBInterpOp
    {
        public:
        
        inline MBInterpOp(Point a_ghost, unsigned int a_order);
        
        template<MemType MEM>
        inline void define(
                MBLevelMap<MEM>&            a_map,
                const std::vector<Point>&   a_footprint,
                int block = -1);

        template<typename T, unsigned int C, MemType MEM>
        inline void apply(
                MBLevelBoxData<T, C, MEM>& a_dst,
                MBLevelBoxData<T, C, MEM>& a_src);
        
        template<typename T, unsigned int C, unsigned int P, MemType MEM>
        inline void coefs(
                MBLevelBoxData<T, P, MEM>& a_coefs,
                MBLevelBoxData<T, C, MEM>& a_src,
                int a_srcComponent = 0);

        inline MBPointInterpOp& operator()(MBDataPoint& a_target)
        {
            for (auto& op : m_ops)
            {
                if (a_target == op.target()) {return op;}
            }
            MayDay<void>::Abort("Failed to find op at target");
        }
        private:
       
        Point m_ghost;
        int m_order;
        std::vector<MBPointInterpOp> m_ops;
    };
#include "implem/Proto_MBInterpOpImplem.H"
} //end namespace Proto
#endif // end include guard
