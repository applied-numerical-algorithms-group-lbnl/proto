#pragma once
#ifndef _PROTO_MB_INTERP_OP_
#define _PROTO_MB_INTERP_OP_

#include "Proto.H"
#include "Proto_MBLevelBoxData.H" //for MBDataPoint definition
#include "Proto_MBLevelMap.H"
#include "Proto_Matrix.H"
#include "Proto_Operator.H"
#include <unordered_set>
#include <algorithm>

namespace Proto
{
    // used for the hash table
    class MBPointID
    {
        public:

        inline MBPointID(Point a_pt, int a_index) : m_point(a_pt), m_index(a_index) {}

        inline MBPointID(const MBPointID& a_rhs) : m_point(a_rhs.m_point), m_index(a_rhs.m_index) {}
        inline bool operator<(const MBPointID& a_rhs) const
        {
            if (m_index != a_rhs.m_index) { return m_index < a_rhs.m_index; }
            if (m_point != a_rhs.m_point) { return m_point < a_rhs.m_point; }
            return false;
        }

        inline bool operator==(const MBPointID& a_rhs) const
        {
            return ((m_index == a_rhs.m_index) && (m_point == a_rhs.m_point));
        }
        
        inline Point point() const { return m_point; }
        inline BlockIndex index() const {return m_index; }

        struct Hash {
            std::size_t operator()(const MBPointID& a_pt) const {
                std::size_t hashValue = 0;
                hash_combine(hashValue, a_pt.index());
                for (int ii = 0; ii < DIM; ii++)
                {
                    hash_combine(hashValue, a_pt.point()[ii]);
                }
                return hashValue;
            }
            
            void hash_combine(std::size_t& seed, int value) const {
                seed ^= std::hash<int>{}(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
        };
        
        private: 
        
        const Point m_point;
        const int   m_index;
    };



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
        template<template<MemType> typename MAP, MemType MEM>
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
        
        private:

        inline void computeSourcePoints(
            const std::set<Point>& footprint,
            const MBDisjointBoxLayout& layout,
            Point ghostSizes);
        inline void computeExponents();
        inline void validateMatrixDimensions(const std::set<Point> &footprint) const;
        template <template <MemType> typename MAP, MemType MEM>
        double computeSourceDataRadius(
            const MBLevelMap<MAP, MEM> &map,
            Array<double, DIM>& xCenter) const;
        template<template<MemType> typename MAP, MemType MEM>
        std::unordered_map<BlockIndex, std::vector<BoxData<double>>> computeMoments(
            const MBLevelMap<MAP, MEM>& map) const;
        void computeCoefficients(std::unordered_map<BlockIndex, std::vector<BoxData<double>>>& momentData);
        
        MBDataPoint m_dst;
        std::vector<MBDataPoint> m_srcs;
        std::vector<Point> m_exponents;
        std::vector<double> m_coefs;
        int m_order;
    };

    /// Mapped Multiblock Block Boundary Interpolation Operator
    /** MBInterpOp interpolates data to all block boundary ghost cells in a 
     *  MBLevelBoxData. Ghost cells associated with domain boundaries are not
     *  affected. */
    class MBInterpOp
    {
        public:
       
        inline MBInterpOp();

        /// Initialize 
        /** This function does not create any operators. The user must call define separately */
        inline MBInterpOp(Point a_ghost, unsigned int a_order = 4);

        template<template<MemType> typename MAP, MemType MEM>
        inline MBInterpOp(
                const MBLevelMap<MAP, MEM>& a_map,
                int                         a_order = 4);

        /// Define Block
        /** Builds the necessary MBPointInterpOp objects for each block boundary ghost cell
         *  implied by the input MBLevelMap and the number of ghost cells this was initialized
         *  with (see constructor). This function only defines the operators used to interpolate
         *  into the specified block. A user may choose to specify a different version of
         *  physical space when building the operators for each block */
        template<template<MemType> typename MAP, MemType MEM>
        inline void define(
                const MBLevelMap<MAP, MEM>& a_map,
                int                         a_order = 4);
        
        /// Define Block
        /** Builds the necessary MBPointInterpOp objects for each block boundary ghost cell
         *  implied by the input MBLevelMap and the number of ghost cells this was initialized
         *  with (see constructor). This function only defines the operators used to interpolate
         *  into the specified block. A user may choose to specify a different version of
         *  physical space when building the operators for each block */
        template<template<MemType> typename MAP, MemType MEM>
        inline void define(
                const MBLevelMap<MAP, MEM>& a_map,
                MBInterpLayout              a_interpLayout,
                int                         a_order = 4);

        template<template<MemType> typename MAP, MemType MEM>
        inline void define(
                const std::vector<std::shared_ptr<MBLevelMap<MAP, MEM>>>& a_maps,
                MBInterpLayout              a_interpLayout,
                int                         a_order = 4);
        
        /// Manual Operator Insertion
        /** Allows a user to build an MBInterpOp by manually constructing the point ops
         */
        inline void insert(std::shared_ptr<MBPointInterpOp> a_op);

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
        /*
        template<typename T, unsigned int C, unsigned int P, MemType MEM>
        inline void coefs(
                MBLevelBoxData<T, P, MEM>& a_coefs,
                MBLevelBoxData<T, C, MEM>& a_src,
                int a_srcComponent = 0);
        */

        /// Get Point Operator
        /** Returns the MBPointInterpOp associated with a specified target cell. */
        inline MBPointInterpOp& operator()(MBDataPoint& a_target)
        {
            for (auto& op : m_ops)
            {
                if (a_target == op->target()) {return *op;}
            }
            MayDay<void>::Abort("Failed to find op at target");
        }

        inline void copyAxis(unsigned int a_axis, bool a_flag);

        void print() const;
       
        template<typename T, unsigned int C, MemType MEM>
        void printErrorPoints(
                MBLevelBoxData<T,C,MEM>& a_err,
                T                        a_eps = 0);

        private:
     
        template<template<MemType> typename MAP, MemType MEM>
        void addPointOps(const MBIndex& index, Point dir, const MBLevelMap<MAP,MEM>& map);
        // template<template<MemType> typename MAP, MemType MEM>
        // void addTriplePointOps(const MBIndex& index, Point dir, const MBLevelMap<MAP,MEM>& map);
        // template<template<MemType> typename MAP, MemType MEM>
        // void addBlockBoundaryOps(const MBIndex& index, Point dir, const MBLevelMap<MAP,MEM>& map);
        // template<template<MemType> typename MAP, MemType MEM>
        // void addTriplePointAdjacentOps(const MBIndex& index, Point dir, const MBLevelMap<MAP,MEM>& map);
        // std::vector<Point> getTriplePointAdjacentFootprint(
        //     Point centerPoint,
        //     const MBIndex& index,
        //     Point boundaryDir,
        //     Point triplePointDir) const;
        bool containsOp(Point dataPoint, const MBIndex& index) const;
        // std::vector<Point> getTriplePointFootprint(
        //     Point centerPoint,
        //     const MBIndex& index,
        //     Point dir) const;
        template<template<MemType> typename MAP, MemType MEM>
        void addPointOp(
            Point dataPoint,
            const MBIndex& index,
            const std::set<Point>& footprint,
            const MBLevelMap<MAP,MEM>& map);
        void addShiftedOp(
            const MBPointInterpOp& baseOp,
            Point shift);
        void copyAlongSymmetryAxes(const MBPointInterpOp& pointOp);

        MBInterpLayout m_interpLayout;
        //mutable Array<bool, DIM> m_copyAxes;
        //MBDisjointBoxLayout m_layout;
        //std::vector<Point> m_baseFootprint;
        //Point m_footprintDims;
        Point m_ghost;
        int m_order;
        std::vector<std::shared_ptr<MBPointInterpOp>> m_ops;
        std::unordered_set<MBPointID, MBPointID::Hash> m_points;
    };

    // Interpolate Block Boundaries
    template<template<MemType> typename MAP, typename T, unsigned int C, MemType MEM, Centering CTR>
    inline void interpBoundaries(
            MBLevelBoxData<T, C, MEM, CTR>& a_data, unsigned int a_order = 4);
    
    // Interpolate Block Boundaries
    template<template<MemType> typename MAP, typename T, unsigned int C, MemType MEM, Centering CTR>
    inline void interpBoundaries(
            MBLevelBoxData<T, C, MEM, CTR>& a_data,
            MBLevelMap<MAP, MEM>& a_map,
            unsigned int a_order = 4);
    
#include "implem/Proto_MBInterpOpImplem.H"
} //end namespace Proto
#endif // end include guard
