#pragma once
#ifndef __PROTO_LEVEL_BC_LIB__
#define __PROTO_LEVEL_BC_LIB__

#include "Proto_LevelBoxData.H"

namespace Proto {
    
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class PeriodicBC : LevelBC<T,C,MEM,CTR>
    {
        public:
        inline void apply(LevelBoxData<T,C,MEM,CTR>& a_data) const;
    };

    #include "implem/Proto_LevelBCLibImplem.H"
} // end namespace Proto
#endif //end include guard
