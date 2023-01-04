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
    template<typename Func, MemType MEM>
    class MBMap;

    typedef MotionItem<MBBoxPartition, MBBoxPartition> MBMotionItem;
    
    struct MBDataPoint
    {
        public:
            MBIndex index;
            int     block;
            Point   boundaryDir;
            Point   point;

            inline MBDataPoint(){};
            inline MBDataPoint(MBIndex a_index, Point a_point);
            inline MBDataPoint(MBIndex a_index, Point a_point, Point a_boundaryDir, int a_block);
            inline bool operator==(const MBDataPoint& a_rhs) const;

            inline bool operator<(const MBDataPoint& a_rhs) const;
            inline bool inBoundary() const {return m_inBoundary; }

            template<typename T, unsigned int C, MemType MEM, Centering CTR>
            inline BoxData<T, C, MEM>& patch(MBLevelBoxData<T, C, MEM, CTR>& a_data) const;
        private:
            bool m_inBoundary;
    };

   inline std::ostream& operator<<(std::ostream& a_os, const MBDataPoint& a_point)
   {
       a_os << "point: " << a_point.point << "| block: " << a_point.block << " | bdry: " << a_point.boundaryDir;
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
    template<typename T, unsigned int C, MemType MEM, Centering CTR=PR_CELL>
    class MBLevelBoxData
    {
        public:
        // PatchID, BlockID, Dir
        typedef std::tuple<Point, unsigned int, Point> boundsKey_t;

        inline MBLevelBoxData();

        inline MBLevelBoxData(
                const MBDisjointBoxLayout& a_layout,
                Point                      a_ghost,
                Point                      a_boundGhost = Point::Zeros());

        inline MBLevelBoxData(
                const MBDisjointBoxLayout& a_layout,
                const Array<Point, DIM+1>& a_ghost,
                Point                      a_boundGhost = Point::Zeros());
        
        inline void define(
                const MBDisjointBoxLayout&    a_layout,
                const Array<Point, DIM+1>&    a_ghost,
                Point                         a_boundGhost = Point::Zeros());
       
        /// Initialize
        /** Initialize the data in this using forall-style syntax. Func should
         *  be a valid forall kernel with the following signature:
         *
         *  void foo(Point& p, Var<T,C,MEM>& d, unsigned int b, ...);
         *
         *  p:    The point of evaluation
         *  d:    A data point in this
         *  b:    The block in which foo is being evaluated
         *  ...:  Other arguments. These may be MBLevelBoxData, MBMap, or any kind of POD
         *        If MBMap is used as an input, the argument will be converted to an
         *        MBLevelBoxData containing the coordinate variables with the same centering
         *        as *this.
         */
        template<typename Func, typename... Srcs>
        inline void initialize(Func& a_func, Srcs&... a_srcs);
        
        inline void setVal(T a_value);
        template<typename Func, typename... Srcs>
        inline void initConvolve(Func& a_func, Srcs... a_srcs);
        inline void fillBoundaries();
        inline void interpolate();
        inline void exchange();
        template<MemType DST_MEM>
        inline void copyTo(MBLevelBoxData<T,C,DST_MEM,CTR>& a_rhs) const;

        Array<Point, DIM+1> ghost() const;
        inline BoxData<T, C, MEM>& operator[](MBIndex a_index);
        inline const BoxData<T, C, MEM>& operator[](MBIndex a_index) const;
        inline LevelBoxData<T, C, MEM, CTR>& blockData(unsigned int a_block);
        inline const LevelBoxData<T, C, MEM, CTR>& blockData(unsigned int a_block) const;
        inline const MBDisjointBoxLayout& layout() const {return *m_layout; }    
        inline Box patchBox(const MBIndex& a_index) const;
        inline unsigned int patchSize(unsigned int a_block) const;
        inline size_t linearSize(const MBIndex& a_index) const;
        inline size_t linearSize() const;
        inline void linearOut(void* a_buffer, const MBIndex& a_index, const Box& a_box) const;
        inline void linearOut(void* a_buffer, const MBIndex& a_index) const;
        inline void linearOut(void* a_buffer) const;
        inline void linearIn(void* a_buffer);
        inline size_t offset(unsigned int a_proc) const;
        inline std::vector<MBBoundaryData<T, C, MEM>>& bounds(
                MBIndex a_index, Point a_dir);
        inline const std::vector<MBBoundaryData<T, C, MEM>>& bounds(
                MBIndex a_index, Point a_dir) const;
        inline bool isBlockBoundary(
                MBIndex a_index, Point a_dir, unsigned int a_adjBlock);
        inline MBBoundaryData<T, C, MEM>& bound(
                MBIndex a_index, Point a_dir, unsigned int a_adjBlock);
        inline const MBBoundaryData<T, C, MEM>& bound(
                MBIndex a_index, Point a_dir, unsigned int a_adjBlock) const;
        inline bool compatible(const MBDisjointBoxLayout& a_layout) const;

        // Interpolation Interface
        inline Var<T, C, MEM> operator[](const MBDataPoint& a_point);
        inline Var<T, C, MEM> operator[](const MBDataPoint& a_point) const;
        inline std::vector<MBDataPoint> interpFootprint(
                const Point& a_center,
                const Point& a_ghost,
                const std::vector<Point>& a_footprint,
                const MBIndex& a_index);
        inline void printBounds(const MBIndex& a_index, int a_prec = 4);
#ifdef PR_OPS
#endif
        inline unsigned int numBlocks() const { return m_layout->numBlocks(); }
        inline MBIterator begin();
        inline MBIterator end();
        private:
       
        MBLevelExchangeCopier<T, C, MEM, CTR> m_exchangeCopier;
        Array<Point, DIM+1> m_ghost;
        const MBDisjointBoxLayout*     m_layout;
        std::vector<std::shared_ptr<LevelBoxData<T, C, MEM, CTR>>> m_data;
        std::map<boundsKey_t, std::vector<MBBoundaryData<T, C, MEM>>> m_bounds;
    };
#include "implem/Proto_MBLevelBoxDataImplem.H"
} // end namespace Proto
#endif //end include guard
