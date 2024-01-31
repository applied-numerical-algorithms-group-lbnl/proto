#pragma once
#ifndef __PROTO_LEVEL_BC_LIB__
#define __PROTO_LEVEL_BC_LIB__

#include "Proto_LevelBC.H"
namespace Proto {
    
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class PeriodicBC : public LevelBC<T,C,MEM,CTR>
    {
        public:
        using LevelBC<T,C,MEM,CTR>::LevelBC;
        inline void apply(LevelBoxData<T,C,MEM,CTR>& a_data) const;
    };

    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class ConstDirichletBC : public LevelBC<T,C,MEM,CTR>
    {
        public:
        using LevelBC<T,C,MEM,CTR>::LevelBC;
        inline void setVal(T a_value, int a_comp = -1);
        inline void apply(LevelBoxData<T,C,MEM,CTR>& a_data) const;
#if PR_MMB
        inline void apply(MBLevelBoxData<T,C,MEM,CTR>& a_data) const;
#endif
        private:
        std::vector<std::pair<int, T>> m_compValues;
        T m_globalValue = 0;
    };
   
#ifdef PR_MMB
#endif

    #include "implem/Proto_LevelBCLibImplem.H"
} // end namespace Proto
#endif //end include guard
