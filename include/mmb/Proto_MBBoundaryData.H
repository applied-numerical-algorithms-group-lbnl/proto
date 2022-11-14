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
            CoordPermutation a_adjToLocal,
            Point a_ghost = Point::Zeros());
        
        inline void define(
            MBIndex a_localIndex,
            MBIndex a_adjIndex,
            Box a_localBox,
            Box a_adjBox,
            CoordPermutation a_adjToLocal,
            Point a_ghost = Point::Zeros());

        inline void fill(const BoxData<T, C, MEM>& a_data);

        inline Box localBox() {return localData->box().grow(-localGhost);}
        inline Box adjBox() {return adjData->box().grow(-adjGhost);}

        inline void print();

        MBIndex localIndex;
        MBIndex adjIndex;
        CoordPermutation adjToLocal;
        Point localGhost;
        Point adjGhost;
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
            CoordPermutation a_adjToLocal,
            Point a_ghost)
    {
        define(a_localIndex, a_adjIndex, a_localBox, a_adjBox, a_adjToLocal, a_ghost);
    }
    
    template<typename T, unsigned int C, MemType MEM>
    void MBBoundaryData<T, C, MEM>::define(
            MBIndex a_localIndex,
            MBIndex a_adjIndex,
            Box a_localBox,
            Box a_adjBox,
            CoordPermutation a_adjToLocal,
            Point a_ghost)
    {
        localGhost = a_ghost;
        localIndex = a_localIndex;
        adjIndex = a_adjIndex;
        Box localBox = a_localBox.grow(a_ghost);
        localData = std::make_shared<BoxData<T, C, MEM>>(localBox, 7);
        adjGhost = a_adjToLocal.inverse()(a_ghost).abs();
        Box adjBox = a_adjBox.grow(adjGhost);
        adjData = std::make_shared<BoxData<T, C, MEM>>(adjBox, 7);
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
