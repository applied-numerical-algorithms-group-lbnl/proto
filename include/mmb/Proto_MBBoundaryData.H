#pragma once
#ifndef _PROTO_MB_BOUNDARY_DATA_
#define _PROTO_MB_BOUNDARY_DATA__

#include "Proto.H"
#include "Proto_MBDisjointBoxLayout.H"

namespace Proto
{
    template<typename T, unsigned int C, MemType MEM>
    struct MBBoundaryData
    {
        MBBoundaryData();

        MBBoundaryData(
            MBIndex a_localIndex,
            MBIndex a_adjIndex,
            Box a_localBox,
            Box a_adjBox,
            CoordPermutation a_adjToLocal);
        
        inline void define(
            MBIndex a_localIndex,
            MBIndex a_adjIndex,
            Box a_localBox,
            Box a_adjBox,
            CoordPermutation a_adjToLocal);

        inline void fill(const BoxData<T, C, MEM>& a_data);

        inline void print();

        MBIndex localIndex;
        MBIndex adjIndex;
        CoordPermutation adjToLocal;
        std::shared_ptr<BoxData<T, C, MEM>> localData;
        std::shared_ptr<BoxData<T, C, MEM>> adjData;
    };
    
    template<typename T, unsigned int C, MemType MEM>
    MBBoundaryData<T, C, MEM>::MBBoundaryData() {}

    template<typename T, unsigned int C, MemType MEM>
    MBBoundaryData<T, C, MEM>::MBBoundaryData(
            MBIndex a_localIndex,
            MBIndex a_adjIndex,
            Box a_localBox,
            Box a_adjBox,
            CoordPermutation a_adjToLocal)
    {
        define(a_localIndex, a_adjIndex, a_localBox, a_adjBox, a_adjToLocal);
    }
    
    template<typename T, unsigned int C, MemType MEM>
    void MBBoundaryData<T, C, MEM>::define(
            MBIndex a_localIndex,
            MBIndex a_adjIndex,
            Box a_localBox,
            Box a_adjBox,
            CoordPermutation a_adjToLocal)
    {
        localIndex = a_localIndex;
        adjIndex = a_adjIndex;
        localData = std::make_shared<BoxData<T, C, MEM>>(a_localBox, 7);
        adjData = std::make_shared<BoxData<T, C, MEM>>(a_adjBox, 7);
        adjToLocal = a_adjToLocal;
    }
    
    template<typename T, unsigned int C, MemType MEM>
    void MBBoundaryData<T, C, MEM>::fill(
            const BoxData<T, C, MEM>& a_data)
    {
        a_data.copyTo(*adjData);
        adjData->copyTo(*localData, adjToLocal);
    }
    
    template<typename T, unsigned int C, MemType MEM>
    void MBBoundaryData<T, C, MEM>::print()
    {
        std::cout << " | localBoundary: " << localData->box();
        std::cout << " | adjBoundary: " << adjData->box() << std::endl;
    }
} // end namespace Proto
#endif //end include guard
