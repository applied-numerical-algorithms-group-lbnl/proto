#pragma once
#ifndef _PROTO_MB_LEVEL_BOUNDARY_DATA_
#define _PROTO_MB_LEVEL_BOUNDARY_DATA_

#include "Proto.H"
#include "Proto_MBDisjointBoxLayout.H"

namespace Proto
{
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class MBLevelBoundaryData
    {
        MBLevelBoundaryData(
                const MBDisjointBoxLayout& a_layout,
                unsigned int a_ghost);
        MBLevelBoundaryData(
                const MBDisjointBoxLayout& a_layout,
                std::array<unsigned int, DIM> a_ghost);
    };
#include "implem/Proto_MBLevelBoundaryDataImplem.H"
} // end namespace Proto
#endif //end include guard
