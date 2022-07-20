#pragma once
#ifndef _PROTO_MB_LEVEL_DATA_
#define _PROTO_MB_LEVEL_DATA_

#include "Proto.H"
#include "Proto_MBDisjointBoxLayout.H"
#include "Proto_MBBoundaryData.H"

namespace Proto
{
    

    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class MBLevelBoxData
    {
        public:
        typedef DataIterator<MBBoxPartition> iterator_t;
        typedef DataIndex<MBBoxPartition> index_t;
        typedef std::pair<Point, index_t> boundsKey_t;

        inline MBLevelBoxData(
                MBDisjointBoxLayout& a_layout,
                unsigned int         a_ghost);

        inline MBLevelBoxData(
                MBDisjointBoxLayout&     a_layout,
                std::array<Point, DIM+1> a_ghost);
        
        inline define(
                MBDisjointBoxLayout&     a_layout,
                std::array<Point, DIM+1> a_ghost);
        
        template<typename Func, typename... Srcs>
        inline void initialize(Func& a_func, Srcs... a_srcs);
        
        template<typename Func, typename... Srcs>
        inline void initConvolve(Func& a_func, Srcs... a_srcs);

        inline BoxData<T, C, MEM>& operator[](index_t a_index);
        inline const BoxData<T, C, MEM>& operator[](index_t a_index) const;
        inline LevelBoxData<T, C, MEM, CTR>& block(unsigned int a_block);
        inline const LevelBoxData<T, C, MEM, CTR>& block(unsigned int a_block) const;
   
        inline std::vector<MBBoundaryData>& bounds(index_t a_patch, Point a_dir);
        inline const std::vector<MBBoundaryData>& bounds(index_t a_patch, Point a_dir) const;

        inline iterator_t begin();
        inline iterator_t end();
        private:
    
        std::array<Point, DIM+1> m_ghost;
        MBDisjointBoxLayout*     m_layout
        std::vector<std::shared_ptr<LevelBoxData<T, C, MEM, CTR>>> m_data;
        std::unordered_map<boundsKey_t, std::vector<MBBoundaryData>> m_bounds;
    };
#include "implem/Proto_MBLevelBoxDataImplem.H"
} // end namespace Proto
#endif //end include guard
