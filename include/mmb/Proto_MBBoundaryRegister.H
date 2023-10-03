#pragma once
#ifndef _PROTO_MB_BOUNDARY_REGISTER_
#define _PROTO_MB_BOUNDARY_REGISTER_

#include "Proto.H"
#include "Proto_MBDisjointBoxLayout.H"
#include "Proto_MBBoundaryData.H"

namespace Proto
{
    template<typename T, unsigned int C, MemType MEM=MEMTYPE_DEFAULT>
    class MBBoundaryRegister
    {
        public:
        // PatchID, BlockID, Dir
        typedef std::tuple<Point, unsigned int, Point> boundsKey_t;

        inline MBBoundaryRegister(
                const MBDisjointBoxLayout&  a_layout,   // layout
                int                         a_depth,    // distance from boundary
                Point                       a_ghost);   // additional storage in all directions

        inline void define(
                const MBDisjointBoxLayout&  a_layout,   // layout
                int                         a_depth,    // distance from boundary
                Point                       a_ghost);   // additional storage in all directions
     
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
        inline std::vector<MBBoundaryData<T,C,MEM>>>& bounds(
                MBIndex a_index,
                Point   a_dir);

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
                MBIndex         a_index,
                unsigned int    a_dir,
                Side::LoHiSide  a_side);

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

        // Set all data to zero
        inline void clear();

        // Exchange and fill adjacent data
        /** Local data is maintained
         */
        inline void exchange();
        
        // Exchange and fill local data
        /** Local data is replaced by adjacent (properly rotated)
         */
        inline void exchangeToLocal();

        private:

        MBDisjointBoxLayout m_layout;
        mutable std::map<boundsKey_t, std::vector<MBBoundaryData<T, C, MEM>>> m_bounds;
    }
} // end namespace Proto
#endif // end include guard
