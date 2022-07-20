#pragma once
#ifndef _PROTO_MB_BOUNDARY_DATA_
#define _PROTO_MB_BOUNDARY_DATA__

#include "Proto.H"

namespace Proto
{
    template<typename T, unsigned int C, MemType MEM>
    struct MBBoundaryData
    {
        MBBoundaryData(
            unsigned int a_localBlock,
            unsigned int a_adjBlock,
            Box a_localBox,
            Box a_adjBox);
        
        inline void define(
            unsigned int a_localBlock,
            unsigned int a_adjBlock,
            Box a_localBox,
            Box a_adjBox);

        unsigned int localBlock;
        unsigned int adjBlock;
        BoxData<T, C, MEM> localData;
        BoxData<T, C, MEM> adjData;
    };

    template<typename T, unsigned int C, MemType MEM>
    MBBoundaryData<T, C, MEM>::MBBoundaryData(
            unsigned int a_localBlock,
            unsigned int a_adjBlock,
            Box a_localBox,
            Box a_adjBox)
    {
        define(a_localBlock, a_adjBlock, a_localBox, a_adjBox);
    }
    
    template<typename T, unsigned int C, MemType MEM>
    void MBBoundaryData<T, C, MEM>::define(
            unsigned int a_localBlock,
            unsigned int a_adjBlock,
            Box a_localBox,
            Box a_adjBox)
    {
        localBlock = a_localBlock;
        adjBlock = a_adjBlock;
        localData.define(a_localBox);
        adjData.define(a_adjBox);
    }
} // end namespace Proto

