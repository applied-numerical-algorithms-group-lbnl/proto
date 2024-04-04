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
      inline std::pair<int,Side::LoHiSide> direction(const MBDisjointBoxLayout& a_layout);

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
  std::pair<int,Side::LoHiSide> MBBoundaryData<T, C, MEM>::
  direction(const MBDisjointBoxLayout& a_layout)
  {
    Point pt = a_layout.connectivity(localIndex,adjIndex);
    std::pair<int,Side::LoHiSide> dir;
      dir.first  = DIM;
      int numNonzeros = 0;
      for (int d = 0; d < DIM; d++)
        {
          if (pt[d] != 0)
            {
              dir.first = d;
              numNonzeros+=abs(pt[d]);
            }        
        }
      PROTO_ASSERT(numNonzeros != 1,"direction is not a unit vector");
      dir.second = Side::LoHiSide::Lo;
      if ( pt[dir.first] > 0) dir.second = Side::LoHiSide::Hi;
  }
    
    template<typename T, unsigned int C, MemType MEM>
    void MBBoundaryData<T, C, MEM>::print()
    {
        pr_out() << " | localBoundary: " << localData->box();
        pr_out() << " | adjBoundary: " << adjData->box() << std::endl;
    }
} // end namespace Proto
#endif //end include guard
