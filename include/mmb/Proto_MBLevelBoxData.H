#pragma once
#ifndef _PROTO_MB_LEVEL_DATA_
#define _PROTO_MB_LEVEL_DATA_

#include "Proto.H"
#include "Proto_MBDisjointBoxLayout.H"
#include "Proto_MBBoundaryData.H"

namespace Proto
{
    // Forward Declaration
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class MBLevelBoxData;
    template<template<MemType> typename MAP, MemType MEM>
    class MBLevelMap;

    typedef MotionItem<MBBoxPartition, MBBoxPartition> MBMotionItem;
  
    /// Multiblock Data Point
    /** This structure is used to specify a data point in the context of mapped multiblock
     *  which may exist in either the interior of a block or in the boundary region
     *  which is allocated between each pair of blocks according to the supplied ghost sizing
     *  when a MBLevelBoxData is created.
     *
     *  This structure is used heavily in the implementation of block boundary interpolation,
     *  but it is otherwise not very useful and is not recommended for public use.
     */
    struct MBDataPoint
    {
        public:
            MBIndex index;
            Point   boundaryDir;
            Point   point;
            
            inline MBDataPoint(){}
            inline MBDataPoint(MBIndex a_index, Point a_point,
                    const MBDisjointBoxLayout& a_layout);
            inline MBDataPoint(MBIndex a_index, Point a_point,
                    const MBDisjointBoxLayout& a_layout, 
                    Point a_boundaryDir, int a_block);
            inline bool operator==(const MBDataPoint& a_rhs) const;

            inline bool operator<(const MBDataPoint& a_rhs) const;
            inline bool inBoundary() const {return m_inBoundary; }
            inline unsigned int srcBlock() const;
            inline unsigned int dstBlock() const;

            template<typename T, unsigned int C, MemType MEM, Centering CTR>
            inline BoxData<T, C, MEM>& patch(MBLevelBoxData<T, C, MEM, CTR>& a_data) const;
            template<typename T, unsigned int C, MemType MEM, Centering CTR>
            inline BoxData<T, C, MEM>& srcPatch(MBLevelBoxData<T, C, MEM, CTR>& a_data) const;
        private:
            int  block;
            bool m_inBoundary;
            const MBDisjointBoxLayout* layout;
    };

   inline std::ostream& operator<<(std::ostream& a_os, const MBDataPoint& a_point)
   {
       a_os << "point: " << a_point.point << "| srcBlock: " << a_point.srcBlock() << " dstBlock: " << a_point.dstBlock() << " | bdry: " << a_point.boundaryDir;
   }
 
    /// Level Copier Operator
    template<typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
    class MBLevelCopierOp
    {
        public:
        MBLevelBoxData<T, C, SRC_MEM, CTR>* m_src;
        MBLevelBoxData<T, C, DST_MEM, CTR>* m_dst;
        
        inline MBLevelCopierOp(){};
        inline MBLevelCopierOp(
                MBLevelBoxData<T, C, SRC_MEM, CTR>& a_srcData,
                MBLevelBoxData<T, C, DST_MEM, CTR>& a_dstData);
        inline uint64_t linearSize(const Box& a_bx, const MBIndex a_index) const;
        inline void linearOut(void* a_buffer, const MBMotionItem& a_info) const;
        inline void linearIn(void* a_buffer, const MBMotionItem& a_info) const;
        inline void localCopy(const MBMotionItem& a_info) const;
    };

    /// Level Copier
    template<typename T, unsigned int C, MemType SRC_MEM, MemType DST_MEM, Centering CTR>
    class MBLevelCopier
    : public Copier
      <MBLevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>, MBBoxPartition, MBBoxPartition, SRC_MEM, DST_MEM>
    {
        public:
        inline MBLevelCopier(){};
        inline void buildMotionPlans(
                MBLevelCopierOp<T, C, SRC_MEM, DST_MEM, CTR>& a_op);
    };
    
    /// Exchange Copier Operator
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class MBLevelExchangeCopierOp
    {
        public:
        MBLevelBoxData<T, C, MEM, CTR>* m_data;
        
        inline MBLevelExchangeCopierOp(){};
        inline MBLevelExchangeCopierOp(MBLevelBoxData<T, C, MEM, CTR>& a_src);
        inline uint64_t linearSize(
                const Box& a_bx,
                const MBIndex& a_index) const;
        inline void linearOut(void* a_buffer, const MBMotionItem& a_info) const;
        
        inline void linearIn(void* a_buffer, const MBMotionItem& a_info) const;
        
        inline void localCopy(const MBMotionItem& a_info) const;
    };
    
    /// Exchange Copier
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class MBLevelExchangeCopier
    : public Copier<MBLevelExchangeCopierOp<T, C, MEM, CTR>, MBBoxPartition, MBBoxPartition, MEM, MEM>
    {
        public:
        inline MBLevelExchangeCopier(){};
        inline void buildMotionPlans(MBLevelExchangeCopierOp<T, C, MEM, CTR>& a_op);
    };
    
    /// Multiblock Level Box Data
    /** This class is the analog of LevelBoxData in the context of mapped multiblock applications.
     *  The most critical differences between this object and the single block version are the
     *  presence of additional data buffers allocated at block boundaries which are utilized for
     *  block boundary interpolation.
     */
    template<typename T, unsigned int C, MemType MEM, Centering CTR=PR_CELL>
    class MBLevelBoxData
    {
        public:
        // PatchID, BlockID, Dir
        typedef std::tuple<Point, unsigned int, Point> boundsKey_t;

        inline MBLevelBoxData();

        /// Main Constructor
        /** Creates and allocates data for an MBLevelBoxData defined on a specified layout and
         *  allocating a specified number of ghost cells.
         *
         *  The ghost input has different semantics from the non-
         *  multiblock context. Here, ghost[0] refers to the number of ghost cells that each patch
         *  will allocate - this is the same as the ghost input used for LevelBoxData.
         *  The additional inputs ghost[n] for n > 0 define the size of the block boundary data
         *  buffers that will be allocated for a block boundary of codimension-n. For 4th order
         *  interpolation at block boundaries, it is recommended that ghost[n] for n>0 be greater
         *  than ghost[0] by at least 2 units in each direction. This recommendation holds even for
         *  cases where ghost[0] is equal to zero in one or more coordinate directions.
         *
         *  \param a_layout     A mapped multiblock layout
         *  \param a_ghost      An array of ghost sizes. See the description above.
         *  \param a_boundGhost Additional allocation for block boundary buffers.
         *                      Use this if there is a need for block boundary buffers to overlap
         *                      with the domain interior. This should remain defaulted unless you
         *                      really know what you are doing.
         */
        inline MBLevelBoxData(
                const MBDisjointBoxLayout& a_layout,
                const Array<Point, DIM+1>& a_ghost,
                Point                      a_boundGhost = Point::Zeros());
        
        /// Simplified Constructor
        /** The ghost input for this constructor corresponds to ghost[0] in the general constructor. The boundary
         *  data buffers are constructed conservatively assuming 4th order block-boundary interpolation is desired.
         */
        inline MBLevelBoxData(
                const MBDisjointBoxLayout& a_layout,
                Point                      a_ghost,
                Point                      a_boundGhost = Point::Zeros());

        
        /// Lazy Constructor
        inline void define(
                const MBDisjointBoxLayout&    a_layout,
                const Array<Point, DIM+1>&    a_ghost,
                Point                         a_boundGhost = Point::Zeros());
       
        /// Initialize FIXME (NOT CURRENTLY WORKING)
        /** Initialize the data in this using forall-style syntax. Func should
         *  be a valid forall kernel with the following signature:
         *
         *  void foo(Point& p, Var<T,C,MEM>& d, unsigned int b, ...);
         *
         *  p:    The point of evaluation
         *  d:    A data point in this
         *  b:    The block in which foo is being evaluated
         *  ...:  Other arguments. These may be MBLevelBoxData, MBLevelMap, or any kind of POD
         *        If MBLevelMap is used as an input, the argument will be converted to an
         *        MBLevelBoxData containing the coordinate variables with the same centering
         *        as *this.
         */
        template<typename Func, typename... Srcs>
        inline void initialize(Func& a_func, Srcs&... a_srcs);
        
        /// Constant Initialization
        inline void setVal(T a_value);

        /// Initialize with Random Noise
        inline void setRandom(T a_low, T a_high); 

        /// FIXME (NOT CURRENTLY WORKING)
        template<typename Func, typename... Srcs>
        inline void initConvolve(Func& a_func, Srcs&... a_srcs);

        /// Exchange Ghost Data
        /** Executes the parallel copy necessary to fill ghost regions.This includes exchanging data
         *  between patches inside of the same block as well as filling block boundary buffers. In
         *  order to fill ghost cells which lie on a block boundary, an MBInterpOp operation must be
         *  applied after this function is called. 
         */
        inline void exchange();

        /// Parallel Data Copy
        /** Executes a parallel copy between two compatible MBLevelBoxData. Only copies valid data;
         *  ghost data will need to be initialized afterwards by calling exchange().
         *  TODO: exchange might actually be copied as well, but I need to test this -CLG
         */
        template<MemType DST_MEM>
        inline void copyTo(MBLevelBoxData<T,C,DST_MEM,CTR>& a_rhs) const;

        /// Get Ghost Size Array
        Array<Point, DIM+1> ghost() const;

        /// Get Patch Reference
        inline BoxData<T, C, MEM>& operator[](MBIndex a_index);
        
        /// Get Patch Reference (Const Overload)
        inline const BoxData<T, C, MEM>& operator[](MBIndex a_index) const;

        /// Get Block Data
        /** Returns the normal LevelBoxData associated with a specified block */
        inline LevelBoxData<T, C, MEM, CTR>& blockData(unsigned int a_block);
        
        /// Get Block Data (Const Overload)
        inline const LevelBoxData<T, C, MEM, CTR>& blockData(unsigned int a_block) const;
        
        /// Get Layout
        inline const MBDisjointBoxLayout& layout() const {return *m_layout; }    
        
        /// Get Patch Box
        /** multiblock analog of LevelBoxData::patchBox. Returns the Box which defines the domain
         *  of a specified patch index, correcting for centering */
        inline Box patchBox(const MBIndex& a_index) const;
        
        /// Get Patch Size
        /** Multiblock analog of LevelBoxData::patchSize. Returns the number of data points in any
         *  patch within the specified block accounting for centering and ghost size.
         */
        inline unsigned int patchSize(unsigned int a_block) const;

        /// Linear Size of Patch
        inline size_t linearSize(const MBIndex& a_index) const;
        
        /// Linear Size
        inline size_t linearSize() const;

        /// Linear Out
        inline void linearOut(void* a_buffer, const MBIndex& a_index, const Box& a_box) const;
        
        /// Linear Out
        inline void linearOut(void* a_buffer, const MBIndex& a_index) const;
        
        /// Linear Out
        inline void linearOut(void* a_buffer) const;
        
        /// Linear In
        inline void linearIn(void* a_buffer);
        
        /// Linear Offset
        inline size_t offset(unsigned int a_proc) const;
        
        /// Get Boundary Buffers
        /** Returns the block boundary buffers for the block boundary specified by an MBIndex
         *  and a direction in the coordinate system of the block specified by the MBIndex.
         *
         *  This function is mostly used for debugging and is not recommended for public use.
         */
        inline std::vector<MBBoundaryData<T, C, MEM>>& bounds(
                MBIndex a_index, Point a_dir);

        /// Get Boundary Buffers (Const Overload)
        inline const std::vector<MBBoundaryData<T, C, MEM>>& bounds(
                MBIndex a_index, Point a_dir) const;
        
        /// Detect Block Boundary
        /** Determine if a specified index and direction correspond to a block boundary with a
         *  specified block index */
        inline bool isBlockBoundary(
                MBIndex a_index, Point a_dir, unsigned int a_adjBlock);
        
        /// Get Block Boundary
        /** Retrieve the block boundary buffer associated with the specified block boundary
         *  information. This function yields an error if isBlockBoundary is false for these
         *  inputs. 
         */
        inline MBBoundaryData<T, C, MEM>& bound(
                MBIndex a_index, Point a_dir, unsigned int a_adjBlock);
        
        /// Get Block Boundary (Const Overload)
        inline const MBBoundaryData<T, C, MEM>& bound(
                MBIndex a_index, Point a_dir, unsigned int a_adjBlock) const;

        /// Compatibility Check
        /** This is the same as this->layout().compatible(a_layout) */
        inline bool compatible(const MBDisjointBoxLayout& a_layout) const;

        /// Convert MBDataPoint to Var
        /** Used by interpolation API. Not recommended for public use. */
        inline Var<T, C, MEM> operator[](const MBDataPoint& a_point);
        
        /// Convert MBDataPoint to Var (Const Overload)
        /** Used by interpolation API. Not recommended for public use. */
        inline Var<T, C, MEM> operator[](const MBDataPoint& a_point) const;

        /// Domain Boundary Check
        /** Determines which domain boundaries - if any - a specified Point in a specified patch
         *  resides on. The output is a vector of directions which correspond to the domain
         *  boundaries the specified point abuts. If the specified point is not adjacent to a
         *  domain boundary, an empty vector is returned.
         */
        inline std::vector<Point> onDomainBoundary(
                const Point& a_point, const MBIndex& a_index) const;

        /// Compute Interpolation Footprint
        /** Computes the list of data points which will be used as the source data for an interpolation
         *  operator to the data point specified by a_center and a_index.
         *
         *  TODO: a_ghost appears here because of a quirk of organization. In the future, this function
         *  will live inside of MBInterpOp instead of here. 
         *
         *  \param a_center     The point to which interpolation will be computed
         *  \param a_ghost      The size of the ghost region on the block boundary
         *  \param a_footprint  A set of Points which represent the interpolation stencil footprint if
         *                      the stencil were to be used in a pure cartesian context
         *  \param a_index      The patch that contains a_center (in its ghost region)
         */
        inline std::vector<MBDataPoint> interpFootprint(
                const Point& a_center,
                const Point& a_ghost,
                const std::vector<Point>& a_footprint,
                const MBIndex& a_index);
        
        /// Print Boundary Buffer Data
        /** Prints the contents of the boundary buffers associated with a specified index.
         */
        inline void printBounds(const MBIndex& a_index, int a_prec = 4);
        
        /// Get Number of Blocks
        inline unsigned int numBlocks() const { return m_layout->numBlocks(); }
        inline MBIterator begin();
        inline MBIterator end();
        inline const MBIterator begin() const;
        inline const MBIterator end() const;
        private:
       
        MBLevelExchangeCopier<T, C, MEM, CTR> m_exchangeCopier;
        Array<Point, DIM+1> m_ghost;
        const MBDisjointBoxLayout*     m_layout;
        std::vector<std::shared_ptr<LevelBoxData<T, C, MEM, CTR>>> m_data;
        mutable std::map<boundsKey_t, std::vector<MBBoundaryData<T, C, MEM>>> m_bounds;
    };
#include "implem/Proto_MBLevelBoxDataImplem.H"
} // end namespace Proto
#endif //end include guard
