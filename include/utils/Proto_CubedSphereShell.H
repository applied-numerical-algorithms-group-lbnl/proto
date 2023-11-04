#pragma once
#include "Proto.H"

#define CUBED_SPHERE_SHELL_R0 1.0
#define CUBED_SPHERE_SHELL_R1 1.1
#define CUBED_SPHERE_SHELL_RADIAL_COORD 0

namespace Proto
{

template<MemType MEM>
class MBMap_CubedSphereShell : public MBMapOp<MEM>
{
public:

    using MBMapOp<MEM>::MBMapOp;

    inline void apply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J,
            FluxBoxData<double, DIM, MEM>& a_NT);
};

template<MemType MEM>
    class MBMap_CubedSphereShellPolar : public MBMapOp<MEM>
{
public:

    using MBMapOp<MEM>::MBMapOp;

    inline void apply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J,
            FluxBoxData<double, DIM, MEM>& a_NT);
};

class CubedSphereShell 
{
    public:

    /// Get Problem Domain
    /**
     *  \param a_domainSize - Size of each block in the angular directions
     *  \param a_thickness  - Size of each block in the radial direction
     *  \param a_radialDir  - Which coordinate corresponds to the radial direction (0 or 2)
     */
    static MBProblemDomain Domain(
        int a_domainSize,
        int a_thickness,
        int a_radialDir = 0);
    
    /// Get Map
    /**
     *  \param a_layout
     *  \param a_ghost
     */
    template<MemType MEM = MEMTYPE_DEFAULT>
    static MBLevelMap<MBMap_CubedSphereShell, MEM> Map(
            const MBDisjointBoxLayout& a_layout,
            const Array<Point, DIM+1> a_ghost);
    
    /// Get Map (Convenient Overload)
    /**
     *  \param a_data
     */
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    static MBLevelMap<MBMap_CubedSphereShell, MEM> Map(
            const MBLevelBoxData<T,C,MEM,CTR>& a_data);
    
    /// Get MB Level Operator
    /** 
     *  \param a_map - The map object produced by one of the Map functions above
     */
    template<template<typename, MemType> class OpType,
        typename T,
        MemType MEM = MEMTYPE_DEFAULT>
    static MBLevelOp<OpType, MBMap_CubedSphereShell, T>
    Operator(MBLevelMap<MBMap_CubedSphereShell, MEM>& a_map);
    
    /// Get MB Interpolation Operator
    /** 
     */
    template<MemType MEM>
    static MBInterpOp InterpOp(
            const MBDisjointBoxLayout& a_layout,
            const Array<Point, DIM+1>& a_ghost,
            unsigned int a_order = 4);
    
    /// Get MB Interpolation Operator
    /** 
     */
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    static MBInterpOp InterpOp(
            MBLevelBoxData<T,C,MEM,CTR>& a_data,
            unsigned int a_order = 4);
    
    /// Interpolate Block Boundary Ghost Cells
    /** Creates and throws away an MBInterpOp. This operator is expensive to build,
     *  so it's usually better to build the operator using InterpOp and keeping it cached.
     */
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    static void InterpBoundaries(
            MBLevelBoxData<T,C,MEM,CTR>& a_data,
            unsigned int a_order = 4);
};
#include "implem/Proto_CubedSphereShellImplem.H"
} // end namespace Proto
