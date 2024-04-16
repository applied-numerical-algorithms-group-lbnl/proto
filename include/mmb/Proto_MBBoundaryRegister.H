#pragma once
#ifndef _PROTO_MB_BOUNDARY_REGISTER_
#define _PROTO_MB_BOUNDARY_REGISTER_

#include "Proto.H"
#include "Proto_MBDisjointBoxLayout.H"
#include "Proto_MBBoundaryData.H"

namespace Proto
{
    typedef MotionItem<MBBoxPartition, MBBoxPartition> MBMotionItem;

    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class MBBoundaryRegister;
    
    /// Copier Operator
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class MBBoundaryCopierOp
    {
        public:
        MBBoundaryRegister<T,C,MEM,CTR>* m_data;
        
        inline MBBoundaryCopierOp(){};
        
        inline MBBoundaryCopierOp(
                MBBoundaryRegister<T, C, MEM, CTR>& a_src);
        
        inline uint64_t linearSize(
                const Box&          a_bx,
                const MBIndex&      a_index) const;
        
        inline void linearOut(
                void*               a_buffer,
                const MBMotionItem& a_info) const;
        
        inline void linearIn(
                void*               a_buffer,
                const MBMotionItem& a_info) const;
        
        inline void localCopy(
                const MBMotionItem& a_info) const;
    };
    
    /// Copier
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class MBBoundaryCopier
    : public Copier<MBBoundaryCopierOp<T, C, MEM, CTR>,
        MBBoxPartition, MBBoxPartition, MEM, MEM>
    {
        public:
        inline MBBoundaryCopier();

        inline void buildMotionPlans(
                MBBoundaryCopierOp<T, C, MEM, CTR>& a_op);
    };
    
    /// Boundary Register
    /** A specialized data holder that contains data in the vicinity of block boundaries
     */
    template<typename T, unsigned int C, MemType MEM=MEMTYPE_DEFAULT, Centering CTR = PR_CELL>
    class MBBoundaryRegister
    {
        public:
        // PatchID, BlockID, Dir
        typedef std::tuple<Point, unsigned int, Point> boundsKey_t;

        /// Trivial Constructor
        inline MBBoundaryRegister(){}

        /// Constructor
        /** Creates a register from a layout. Depth is the distance in cells from a boundary that
         *  will be allocated in each MBBoundaryData's base buffer. If depth is negative, the 
         *  register's base region is inside the problem domain and vice versa. Each of these 
         *  buffers can optionally be enlarged anisotropically using the ghost parameter.
         *  If "bothSides" is true, all registers will allocate space on both sides of the
         *  boundary. In this case, the sign of depth has no effect.
         */
        inline MBBoundaryRegister(
                MBDisjointBoxLayout         a_layout,           
                int                         a_depth,            
                Point                       a_ghost,            
                bool                        a_bothSides = true);
        
        /// Lazy Constructor
        inline void define(
                MBDisjointBoxLayout         a_layout,   
                int                         a_depth,    
                Point                       a_ghost,    
                bool                        a_bothSides = true);
        
        // Get Boundary Data
        inline MBBoundaryData<T,C,MEM>& bound(
                MBIndex         a_index,
                Point           a_dir,
                unsigned int    a_adjBlock);
        
        // Get Boundary Data
        inline MBBoundaryData<T,C,MEM>& bound(
                MBIndex         a_index,
                unsigned int    a_dir,
                Side::LoHiSide  a_side);
        
        // Get All Boundary Data
        inline std::vector<MBBoundaryData<T,C,MEM>> bounds(
                MBIndex a_index,
                Point   a_dir);
        
        // Get All Boundary Data
        inline std::vector<MBBoundaryData<T,C,MEM>> bounds(
                MBIndex a_index);

        // Get Local Data
        inline BoxData<T,C,MEM>& local(
                MBIndex         a_index,
                Point           a_dir,
                unsigned int    a_adjBlock);
        
        // Get Local Data
        inline BoxData<T,C,MEM>& local(
                MBIndex         a_index,
                unsigned int    a_dir,
                Side::LoHiSide  a_side);
        
        // Get Local Data
        inline BoxData<T,C,MEM>& local(
                MBIndex         a_localIndex,
                MBIndex         a_adjIndex);

        // Get Adjacent Data
        inline BoxData<T,C,MEM>& adjacent(
                MBIndex         a_index,
                unsigned int    a_dir,
                Side::LoHiSide  a_side);

        // Get Adjacent Data
        inline BoxData<T,C,MEM>& adjacent(
                MBIndex         a_index,
                Point           a_dir,
                unsigned int    a_adjBlock);
        
        // Get Adjacent Data
        inline BoxData<T,C,MEM>& adjacent(
                MBIndex         a_localIndex,
                MBIndex         a_adjIndex);
        
        // Get Rotation
        inline CoordPermutation toLocal(
                MBIndex         a_index,
                Point           a_dir,
                unsigned int    a_adjBlock);
        
        // Get Rotation
        inline CoordPermutation toLocal(
                MBIndex         a_index,
                unsigned int    a_dir,
                Side::LoHiSide  a_side);

        // Zero All Registers
        inline void clear();

        // Zero Local Registers
        inline void clearLocal();
        
        // Zero Adjacent Registers
        inline void clearAdjacent();

        // Fill Adjacent Registers
        inline void exchange();
        
        // Get Layout
        inline MBDisjointBoxLayout layout() {return m_layout; }
        
        // Get Ghost Size
        inline Point ghost() {return m_ghost; }

        private:

        Point m_ghost;
        MBDisjointBoxLayout m_layout;
        mutable std::map<boundsKey_t, std::vector<MBBoundaryData<T, C, MEM>>> m_bounds;
        MBBoundaryCopier<T,C,MEM,CTR> m_copier; 
    };
#include "implem/Proto_MBBoundaryRegisterImplem.H"
} // end namespace Proto
#endif // end include guard
