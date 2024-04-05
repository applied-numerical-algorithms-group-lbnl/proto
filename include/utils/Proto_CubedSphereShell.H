#pragma once
#ifndef _PROTO_CUBED_SPHERE_SHELL_H_
#define _PROTO_CUBED_SPHERE_SHELL_H_
#include "Proto.H"
#include "MBMap_CubedSphereShell.H"
#include "Proto_CubedSphereKernels.H"

namespace Proto {

  /// Cubed Sphere Shell Utility
  /** This class serves as a catchall namespace for generating Proto data structures
   *  that use MBMap_CubedSphereShell
   */

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
                                                       Point a_ghost);
    
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
    template<template<typename, template<MemType> typename, MemType> class OpType,
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
                               Point a_ghost,
                               int a_order = 4);
    
    /// Get MB Interpolation Operator
    /** 
     */
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    static MBInterpOp InterpOp(
                               MBLevelBoxData<T,C,MEM,CTR>& a_data,
                               int a_order = 4);
    
    /// Interpolate Block Boundary Ghost Cells
    /** Creates and throws away an MBInterpOp. This operator is expensive to build,
     *  so it's usually better to build the operator using InterpOp and keeping it cached.
     */
  
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    static void InterpBoundaries(
                                 MBLevelBoxData<T,C,MEM,CTR>& a_data,
                                 int a_order = 4);
    
    template<typename T, unsigned int C, MemType MEM>
    static void consToSphNGEuler( 
                                 BoxData<T,C,MEM>& a_JU,
                                 BoxData<T,1,MEM>& a_dVolr,
                                 Box a_bx,
                                 Box a_blockBox,
                                 T a_dxiPerp,
                                 int a_block,                 
                                 unsigned int a_order = 4);
    
     template<typename T, unsigned int C, MemType MEM, Centering CTR>
     static void consToSphInterpEuler( 
                                      MBLevelBoxData<T,C,MEM,CTR>& a_JU,
                                      MBInterpOp& a_iop,
                                      MBLevelBoxData<T,1,MEM,CTR>& a_dVolr, 
                                      unsigned int a_order = 4);
     template<typename T, unsigned int C, MemType MEM, Centering CTR>
     static void consToSphInterpEuler( 
                                      MBLevelBoxData<T,C,MEM,CTR>& a_JU,
                                      MBLevelBoxData<T,1,MEM,CTR>& a_dVolr,
                                      Array<T,DIM> a_dx,
                                      unsigned int a_order = 4);
    
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    static void interpBoundariesMHD( 
                                    MBLevelBoxData<T,C,MEM,CTR>& a_data,
                                    const MBInterpOp& a_op,
                                    unsigned int a_order = 4);
    
  
    /// Average of n x \partial_{xi_2} n over the {xi_1 = const.} face.

    template<typename T,MemType MEM>
    static void nCrossD2n(
                          BoxData<T,DIM,MEM>& a_n_x_d2n,
                          Box a_bx,
                          T a_h,
                          uint a_block);

    /// Average of \partial_{xi_1} n x n over the {xi_2 = const.} face.

    template<typename T,MemType MEM>
    static void d1nCrossN(
                          BoxData<T,DIM,MEM>& a_d1n_x_n,
                          Box a_bx,
                          T a_h,
                          uint a_block);

    /// Average of A over faces. 
    template<typename T,MemType MEM>
    static BoxData<T,DIM,MEM,DIM> AMatrix(
                                          Box a_bx,
                                          T a_h,
                                          unsigned int a_block,
                                          unsigned int a_centering);

    /// Average of (det(A) A) over faces. 
    template<typename T,MemType MEM>
    static BoxData<T,DIM,MEM,DIM> detA_A(
                                         Box a_bx,
                                         T a_h,
                                         unsigned int a_block,
                                         unsigned int a_centering);

    /// Average of rows of adj(A) over faces.

    template<typename T,MemType MEM>
    static BoxData<T,DIM,MEM> adjA(
                                   Box a_bx,
                                   T a_h,
                                   int a_block,
                                   int a_centering);

    /// Average of determinant(A) over faces.

    template<typename T,MemType MEM>  
    static BoxData<T,1,MEM> detA(
                                 Box a_bx,
                                 T a_h,
                                 unsigned int a_centering);

    /// Get Array that encodes the permutation matrix relative to block 3.  
    inline static Array<uint,DIM> perm(uint a_block)
    {
      Array<Array<uint,DIM>,6> permvals
        =  {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
      return permvals[a_block];
    };
    /// Get Array that encodes the signs of the permuted coordinates, relative to block 3.
    inline static Array<int,DIM> permsign(uint a_block)
    {
      Array<Array<int,DIM>,6> signvals
        = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};
      return signvals[a_block];
    };

    template <typename T, MemType MEM>
    inline
    static void cellMatrixQuotientMatAndInv(
                                     BoxData<T,DIM,MEM,1>& a_quot,
                                     BoxData<T,DIM,MEM,1>& a_quotBar,
                                     const BoxData<T,DIM,MEM,DIM>& a_denominator,
                                     const BoxData<T,DIM,MEM,DIM>& a_denomInv,
                                     const BoxData<T,DIM,MEM,1>& a_numerator);
    template <typename T, MemType MEM>
    inline
    static void cellMatrixQuotientMatAndInv(
                                     BoxData<T,DIM,MEM,1>& a_quot,
                                     const BoxData<T,DIM,MEM,DIM>& a_denom,
                                     const BoxData<T,DIM,MEM,DIM>& a_denomInv,
                                     const BoxData<T,DIM,MEM,1>& a_numerator,
                                     Box a_quotBox,
                                     Box a_blockBox);
    template <typename T, MemType MEM>
    inline
    static void cellQuotientValAndInv(
                               BoxData<T,1,MEM,1>& a_quot,
                               BoxData<T,1,MEM,1>& a_quotBar,
                               const BoxData<T,1,MEM>& a_denominator,
                               const BoxData<T,1,MEM>& a_denomInv,
                               const BoxData<T,1,MEM>& a_numerator);
    template <typename T, MemType MEM>
    inline
    static void cellQuotientValAndInv(
                               BoxData<T,1,MEM>& a_quot,
                               const BoxData<T,1,MEM>& a_denominator,
                               const BoxData<T,1,MEM>& a_denomInv,
                               const BoxData<T,1,MEM>& a_numerator,
                               Box a_quotBox,
                               Box a_blockBox);
    
    template <typename T, MemType MEM>
    inline
    static void cellQuotientValAndInvRadial(
                                            BoxData<T,1,MEM>& a_quot,
                                            const BoxData<T,1,MEM>& a_denominator,
                                            const BoxData<T,1,MEM>& a_denomInv,
                                            const BoxData<T,1,MEM>& a_numerator,
                                            Box a_quotBox,
                                            Box a_blockBox);
    
    template <typename T, unsigned int C, MemType MEM>
    static BoxData<T,C,MEM> deriv3case(
                                       const BoxData<T,C,MEM>& a_input,
                                       Box a_outputBox,
                                       Box a_validBox,
                                       uint a_dir);
  };
#if 0
  template<typename T,MemType MEM>
  inline BoxData<T,1,MEM>
  detA(
       Box a_bx,
       T a_h,
       unsigned int a_centering);
#endif
#include "implem/Proto_CubedSphereShellImplem.H"
} // end namespace Proto
#endif // end include guard
