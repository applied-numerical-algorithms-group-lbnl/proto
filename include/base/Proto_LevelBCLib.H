#pragma once
#ifndef __PROTO_LEVEL_BC_LIB__
#define __PROTO_LEVEL_BC_LIB__

#include "Proto_LevelBC.H"
namespace Proto {

    /// Periodic Boundary Condition
    /** This is the default boundary condition for most of Proto. 
     *  All this does is call a_data.exchange()
     */
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class PeriodicBC : public LevelBC<T,C,MEM,CTR>
    {
        public:
        using LevelBC<T,C,MEM,CTR>::LevelBC;
        inline void apply(LevelBoxData<T,C,MEM,CTR>& a_data) const;
    };

    /// Not Really Dirichlet BC
    /** Test Boundary Condition that fills ghost cells with a constant
     */
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class ConstDirichletBC : public LevelBC<T,C,MEM,CTR>
    {
        public:
        using LevelBC<T,C,MEM,CTR>::LevelBC;
        inline void setVal(T a_value, int a_comp = -1);
        inline void apply(LevelBoxData<T,C,MEM,CTR>& a_data) const;
        private:
        std::vector<std::pair<int, T>> m_compValues;
        T m_globalValue = 0;
    };
    #include "implem/Proto_LevelBCLibImplem.H"
} // end namespace Proto
#endif //end include guard
