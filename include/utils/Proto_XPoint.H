#pragma once
#ifndef _PROTO_X_POINT_H_
#define _PROTO_X_POINT_H_
#include "Proto.H"
#include "MBMap_XPointRigid.H"

namespace Proto {

    template<int NBLOCK>
    class XPoint
    {
        public:

        static MBProblemDomain Domain(
            int a_domainSize,
            int a_thickness);

        template<MemType MEM = MEMTYPE_DEFAULT>
        static MBLevelMap<MBMap_XPointRigid<NBLOCK, MEM>, MEM> Map(
            const MBDisjointBoxLayout& a_layout,
            Point a_ghost);
        
        template<typename T, unsigned int C, MemType MEM, Centering CTR>
        static MBLevelMap<MBMap_XPointRigid<NBLOCK, MEM>, MEM> Map(
            const MBLevelBoxData<T,C,MEM,CTR>& a_data);
    };
#include "implem/Proto_XPointImplem.H"
} // end namespace Proto
#endif // end include guard
