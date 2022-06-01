#pragma once
#ifndef _PROTO_MB_LEVEL_BOX_DATA_
#define _PROTO_MB_LEVEL_BOX_DATA_

#include "Proto.H"
#include "Proto_MBDisjointBoxLayout.H"
#include "Proto_MBDataIterator.H"

namespace Proto
{
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class MBLevelBoxData
    {
        public:
        MBLevelBoxData(const MBDisjointBoxLayout& a_layout, unsigned int a_ghost);
        MBLevelBoxData(const MBDisjointBoxLayout& a_layout, std::array<unsigned int, DIM> a_ghost);

        template<typename Func, typename... Srcs>
        inline void initialize(Func& a_func, Srcs... a_srcs);
        
        template<typename Func, typename... Srcs>
        inline void initConvolve(Func& a_func, Srcs... a_srcs);
        
        inline BoxData<T, C, MEM>& operator[](const MBDataIndex& a_index);
        inline const BoxData<T, C, MEM>& operator[](const MBDataIndex& a_index) const;
        inline LevelBoxData<T, C, MEM, CTR>& block(unsigned int a_blockID);
        inline const LevelBoxData<T, C, MEM, CTR>& block(unsigned int a_blockID) const;

        inline void copyTo(MBLevelBoxData<T, C, MEM, CTR>& a_dest);
        inline void exchange();     // normal exchange + fill m_boundData
        inline void interpolate();  // interpolate from m_boundData into m_data

        inline MBDataIterator& begin();
        inline MBDataIterator& end();

        private:
        MBLevelBoxData(const MBLevelBoxData& a_rhs); // copying may lead to unexpected behavior
        MBDisjointBoxLayout m_layout;
        std::array<unsigned int, DIM> m_ghost;
        std::vector<std::shared_ptr<LevelBoxData<T, C, MEM, CTR>>> m_data;
        std::vector<std::shared_ptr<MBLevelBoundaryData<T, C, MEM>>> m_boundData;
    };
#include "implem/Proto_MBLevelBoxDataImplem.H"
} // end namespace Proto
#endif //end include guard
