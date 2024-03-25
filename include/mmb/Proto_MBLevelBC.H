#pragma once
#ifndef __MB_PROTO_LEVEL_BC__
#define __MB_PROTO_LEVEL_BC__

#include "Proto_MBLevelBoxData.H"
namespace Proto {
    
    template<typename T, unsigned int C, 
        template<MemType> typename MAP,
        MemType MEM=MEMTYPE_DEFAULT, Centering CTR=PR_CELL>
    class MBLevelBC : public LevelBC<T,C,MEM,CTR> {
        public:

        // Constructor
        inline MBLevelBC(const MBLevelMap<MAP, MEM>& a_map);
        
        // Basic interface
        inline virtual void init() {}
        inline virtual void apply(MBLevelBoxData<T,C,MEM,CTR>& a_data) const {}

        inline const MBLevelMap<MAP, MEM> map() { return *m_map; }

        private:

        const MBLevelMap<MAP, MEM>* m_map;
    };

    #include "implem/Proto_MBLevelBCImplem.H"
} // end namespace Proto
#include "Proto_MBLevelBCLib.H"
#endif //end include guard
