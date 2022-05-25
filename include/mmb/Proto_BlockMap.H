#pragma once
#ifndef _PROTO_BLOCK_MAP_
#define _PROTO_BLOCK_MAP_

namespace Proto
{
    template<typename T>
    class BlockMap
    {
        typedef T(*mapFunc_t)(T);
        public:
        
        private:
        
        mapFunc_t m_compToPhys;
        mapFunc_t m_physToComp;
        std::array<T, DIM> m_origin;
        std::array<T, DIM> m_dx;
    }
} // end namespace Proto

#include "implem/Proto_BlockMapImplem.H"
#endif // end include guard
