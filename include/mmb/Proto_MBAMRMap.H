#pragma once
#ifndef _PROTO_MB_AMR_MAP_
#define _PROTO_MB_AMR_MAP_

#include "Proto_MBAMRData.H"

namespace Proto
{

template<template<MemType> typename MAP, MemType MEM = MEMTYPE_DEFAULT>
class MBAMRMap
{
    public:
        /// Trivial Constructor
        inline MBAMRMap(){}
        
        /// Non-Trivial Constructor
        inline MBAMRMap(
                const MBAMRGrid&            a_grid,
                const Array<Point, DIM+1>&  a_ghost,
                unsigned int                a_focalBlock = 0);
        
        /// Lazy Constructor
        inline void define(
                const MBAMRGrid&            a_grid,
                const Array<Point, DIM+1>&  a_ghost,
                unsigned int                a_focalBlock = 0);
        
        inline void define(
                const MBAMRGrid&            a_grid,
                const Point                 a_ghost,
                unsigned int                a_focalBlock = 0);
        
        /// Get Level Map
        inline const MBLevelMap<MAP,MEM>& getLevel(unsigned int a_level) const;

        /// Get Grid
        inline const MBAMRGrid& grid() const {return *m_grid;}

        /// Get Cached Coordinates
        inline const MBAMRData<double, DIM, MEM, PR_NODE>& map() const {return m_X; }

        /// Get Cached Jacobian
        inline const MBAMRData<double, 1, MEM, PR_NODE>& jacobian() const {return m_J; }

        inline MBLevelMap<MAP,MEM>& operator[](unsigned int a_level){ return *m_maps[a_level];}

    private:
        const MBAMRGrid* m_grid; 
        MBAMRData<double, DIM, MEM, PR_NODE> m_X;
        MBAMRData<double, 1, MEM, PR_CELL> m_J;
        std::vector<std::shared_ptr<MBLevelMap<MAP,MEM>>> m_maps;
};
#include "implem/Proto_MBAMRMapImplem.H"
} //end Proto namespace
#endif
