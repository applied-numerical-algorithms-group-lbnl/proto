#pragma once
#ifndef __MB_PROTO_LEVEL_BC_LIB__
#define __MB_PROTO_LEVEL_BC_LIB__

#include "Proto_MBLevelBC.H"
namespace Proto {

    /// No Boundary Condition
    /** This is the default boundary condition for MMB
     *  and it does nothing
     */
    template<typename T, unsigned int C, 
        template<MemType> typename MAP, 
        MemType MEM, Centering CTR>
    class NullBC : public MBLevelBC<T,C,MAP,MEM,CTR>
    {
        public:
        using MBLevelBC<T,C,MAP,MEM,CTR>::MBLevelBC;
        inline void apply(MBLevelBoxData<T,C,MEM,CTR>& a_data) const {}
    };

    /// Constant Ghost Value
    /** Fills ghost cells with a constant value which can be different
     *  for each component of the state. This is mostly just for testing
     *  purposes.
     */
    template<typename T, unsigned int C,
        template<MemType> typename MAP,
        MemType MEM, Centering CTR>
    class ConstBC : public MBLevelBC<T,C,MAP,MEM,CTR>
    {
        public:
        using MBLevelBC<T,C,MAP,MEM,CTR>::MBLevelBC;
        inline void setVal(T a_value, int a_comp = -1);
        inline void apply(MBLevelBoxData<T,C,MEM,CTR>& a_data) const;
        private:
        std::vector<std::pair<int, T>> m_compValues;
        T m_globalValue = 0;
    };
    #include "implem/Proto_MBLevelBCLibImplem.H"
} // end namespace Proto
#endif //end include guard
