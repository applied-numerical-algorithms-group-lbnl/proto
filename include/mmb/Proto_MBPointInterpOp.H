#pragma once
#ifndef _PROTO_MB_POINT_INTERP_OP_
#define _PROTO_MB_POINT_INTERP_OP_
#include "Proto.H"

namespace Proto
{
    /// Mapped Multiblock Block Boundary Interpolation Operator
    /** MBPointInterpOp interpolates to exactly one point in a
     *  MBLevelBoxData. This object is used by MBInterpOp which executes the
     *  interpolation across the entire MBLevelBoxData */
    class MBPointInterpOp
    {
        public:

        /// Constructor
        /** Creates a PointInterpOp to the destination point dst which
         *  Which defines a specific block boundary location in an
         *  MBLevelBoxData defined on the specified map with the specified
         *  number of ghost cells. The operator itself is generated using
         *  a list of relative Points which correspond to a traditional
         *  Stencil operation on a Cartesion coordinate system. The order
         *  specifies the intended order of accuracy of the interpolation.
         *
         *  TODO: order can probably be inferred from the footprint or vice versa
         */
        template<typename MAP, MemType MEM>
        inline MBPointInterpOp(
            MBDataPoint                 a_dst,
            Point                       a_ghost,
            const MBLevelMap<MAP, MEM>& a_map,
            std::set<Point>          a_footprint,
            unsigned int                a_order);
        
        inline MBPointInterpOp(
            MBDataPoint a_dst,
            std::vector<MBDataPoint> a_srcs,
            std::vector<double> a_coefs);
        
        /// Apply
        /** computes the interpolation */
        template<typename T, unsigned int C, MemType MEM>
        inline void apply(
                MBLevelBoxData<T, C, MEM>& a_dst,
                MBLevelBoxData<T, C, MEM>& a_src);
        
        /// Size
        /** Returns the number of coeffcients used in the interpolation */
        inline unsigned int size() const;
        
        /// Coefficients
        /** Returns the coefficients of the interpolating polynomial for a
         *  specified component in a specified data set.
         *
         *  Mostly used for debugging. */
        //template<typename T, unsigned int C, MemType MEM>
        //inline Matrix<double> coefs(
        //        MBLevelBoxData<T, C, MEM>& a_src,
        //        int a_srcComponent = 0);
        
        /// Target
        /** Returns the MBDataPoint corresponding to the point of interpolation */
        inline const MBDataPoint& target() const {return m_dst;}

        /// Sources
        /** Returns the list of locations used for input data to the interpolation
         *  as a list of MBDataPoint objects. The order of the returnd points corresponds
         *  to the order of the coefficients in the operator matrix (outpout of SMatrix) */
        inline std::vector<MBDataPoint> sources() const {return m_srcs;}
        inline std::vector<double> coefs() const {return m_coefs; }
        inline void print() const; 
        inline void writeFootprint(std::string filename) const;
        
        private:

        inline void computeSourcePoints(
            const std::set<Point>& footprint,
            const MBDisjointBoxLayout& layout,
            Point ghostSizes);
        inline void computeExponents();
        inline void validateMatrixDimensions(const std::set<Point> &footprint) const;
        template <typename MAP, MemType MEM>
        double computeSourceDataRadius(
            const MBLevelMap<MAP, MEM> &map,
            Array<double, DIM>& xCenter) const;
        template<typename MAP, MemType MEM>
        std::unordered_map<BlockIndex, std::vector<BoxData<double>>> computeMoments(
            const MBLevelMap<MAP, MEM>& map) const;
        void computeCoefficients(std::unordered_map<BlockIndex, std::vector<BoxData<double>>>& momentData);

        MBDataPoint m_dst;
        std::vector<MBDataPoint> m_srcs;
        std::vector<Point> m_exponents;
        std::vector<double> m_coefs;
        int m_order;
    };
    #include "implem/Proto_MBPointInterpOpImplem.H"
} // end proto namespace
#endif // end include guard