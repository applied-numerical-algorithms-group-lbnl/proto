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
        template<MemType MEM>
        inline MBPointInterpOp(
            MBDataPoint                 a_dst,
            Point                       a_ghost,
            MBLevelMap<MEM>&            a_map,
            const std::vector<Point>&   a_footprint,
            unsigned int                a_order);
        
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
        template<typename T, unsigned int C, MemType MEM>
        inline Matrix<double> coefs(
                MBLevelBoxData<T, C, MEM>& a_src,
                int a_srcComponent = 0);
        
        /// Target
        /** Returns the MBDataPoint corresponding to the point of interpolation */
        inline const MBDataPoint& target() const {return m_dst;}

        /// Sources
        /** Returns the list of locations used for input data to the interpolation
         *  as a list of MBDataPoint objects. The order of the returnd points corresponds
         *  to the order of the coefficients in the operator matrix (outpout of SMatrix) */
        inline const std::vector<MBDataPoint>& sources() const {return m_srcs;}

        /// Exponents
        /** Returns the list of multidimensional exponents used in the interpolating
         *  polynomial. The order of these exponents corresponds with the columns of the
         *  C and D matrices (outputs of CMatrix and DMatrix) */
        inline std::vector<Point> exponents() const {return m_exponents; }

        /// C Matrix
        /** Returns the matrix of source moments used to generate the interpolation operator.
         *  The number of rows corresponds to the number of source points and the number of
         *  columns correspond to the polynomial terms (see exponents)*/
        inline Matrix<double>& CMatrix() {return m_C; }
        
        /// D Matrix
        /** Returns the matrix of destination moments used to generate the operator.
         *  The output has a single row associated with the destination point and the
         *  columns correspond to the polynomial terms (see exponents)*/
        inline Matrix<double>& DMatrix() {return m_D; }

        /// S Matrix
        /** Matrix of coefficients used to weight the source data and compute the interpolation.
         *  S has a single row with entries equal to size(). Each entry is the coeffcient that
         *  scales the data point from the associated source */
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

    /// Mapped Multiblock Block Boundary Interpolation Operator
    /** MBInterpOp interpolates data to all block boundary ghost cells in a 
     *  MBLevelBoxData. Ghost cells associated with domain boundaries are not
     *  affected. */
    class MBInterpOp
    {
        public:
        
        /// Initialize 
        /** This function does not create any operators. The user must call define separately */
        inline MBInterpOp(Point a_ghost, unsigned int a_order);
        
        /// Define Block
        /** Builds the necessary MBPointInterpOp objects for each block boundary ghost cell
         *  implied by the input MBLevelMap and the number of ghost cells this was initialized
         *  with (see constructor). This function only defines the operators used to interpolate
         *  into the specified block. A user may choose to specify a different version of
         *  physical space when building the operators for each block */
        template<MemType MEM>
        inline void define(
                MBLevelMap<MEM>&            a_map,
                const std::vector<Point>&   a_footprint,
                int block = -1);
        
        /// Apply
        /** Interpolates from source data to all block boundary ghost cells in the destination.
         *  It is assumed that MBLevelBoxData::fillBoundaries() has been called on the source
         *  data beforehand. */
        template<typename T, unsigned int C, MemType MEM>
        inline void apply(
                MBLevelBoxData<T, C, MEM>& a_dst,
                MBLevelBoxData<T, C, MEM>& a_src);
        
        /// Coefficients
        /** Writes the first P coefficients of the interpolating polynomial to a specified
         *  component of a specified source data array. Useful for debugging */
        template<typename T, unsigned int C, unsigned int P, MemType MEM>
        inline void coefs(
                MBLevelBoxData<T, P, MEM>& a_coefs,
                MBLevelBoxData<T, C, MEM>& a_src,
                int a_srcComponent = 0);

        /// Get Point Operator
        /** Returns the MBPointInterpOp associated with a specified target cell. */
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
