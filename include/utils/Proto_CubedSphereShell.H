#pragma once
#include "Proto.H"

#define CUBED_SPHERE_SHELL_R0 1.0
#define CUBED_SPHERE_SHELL_R1 1.1
#define CUBED_SPHERE_SHELL_RADIAL_COORD 0
typedef Proto::Array<Proto::Array<double,DIM>,DIM> matType;
typedef Proto::Array<double,DIM> vecType;
using namespace Proto;
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
  template <typename T, unsigned int C, MemType MEM>
  void cellMatrixQuotientMatAndInv(
                                   BoxData<T,C,MEM,1>& a_quot,
                                   BoxData<T,C,MEM,1>& a_quotBar,
                                   const BoxData<T,C,MEM,C>& a_denominator,
                                   const BoxData<T,C,MEM,C>& a_denomInv,
                                   const BoxData<T,C,MEM,1>& a_numerator);
  /*{
    a_quotBar.define(a_denomInv.box());
    T one = 1.0;
    forallInPlace(f_matVecProd,a_quotBar,a_denomInv,a_numerator,one);
    BoxData<T,C,MEM,1> corr(a_quotBar.box().grow(-1));
    corr.setToZero();
    // need to fix - this is actually a matrix multiplication.
    for (int dir = 0; dir < DIM; dir++)
      {
        BoxData<T,C,MEM,1> Dquotbar = Stencil<T>::Derivative(1,dir)(a_quotBar);
        BoxData<T,C,MEM,C> Ddenom = Stencil<T>::Derivative(1,dir)(a_denominator);
        BoxData<T,C,MEM,1> prod(Ddenom.box());
        T twelveth = -1.0/12.0;
        forallInPlace(f_matVecProd,prod,Ddenom,Dquotbar,twelveth);      
        corr += prod;
      }
    a_quot.define(corr.box());
    forallInPlace([] PROTO_LAMBDA(
                                  Var<T,C,MEM,1>& a_quot,
                                  const Var<T,C,MEM,1>& a_quotbar,
                                  const Var<T,C,MEM,C>& a_denominv,
                                  const Var<T,C,MEM,1>& a_corr)
                  {
                    Var<T,C,MEM,1> prod;
                    f_matVecProd(prod,a_denominv,a_corr,1.0);
                    for (int ind0 = 0; ind0 < C; ind0++)
                      {
                        a_quot(ind0) = a_quotbar(ind0)+prod(ind0);
                      }
                  },a_quot,a_quotBar,a_denomInv,corr);
                  }*/
  template<typename T,MemType MEM>
  inline void NTRowCS(
                      BoxData<T,DIM,MEM>& a_NTrow,
                      const BoxData<T,1,MEM>& a_adjD,
                      Box a_bx,
                      T a_h,
                      int a_block,
                      int a_centering);
  template<typename T,MemType MEM>
  inline void
  scaledDetA_A(
               BoxData<T,DIM,MEM,DIM>& a_scaledDetA_A,
               const BoxData<T,1,MEM>& a_radialScale,
               Box a_bx,
               T a_h,
               unsigned int a_block,
               unsigned int a_centering,
               bool a_scaled = true);
  template<typename T,MemType MEM>
  inline void
  detA(
       BoxData<T,DIM,MEM,DIM>& a_detA,
       Box a_bx,
       T a_h,
       unsigned int a_block,
       unsigned int a_centering);
  template<typename T,MemType MEM>
  inline void
  detA_A(
         BoxData<T,DIM,MEM,DIM>& a_detA_A,
         Box a_bx,
         T a_h,
         unsigned int a_block,
         unsigned int a_centering);  
  template <typename T, MemType MEM>
  inline void d1nCrossN(
                        BoxData<T,DIM,MEM>& a_d1n_x_n,
                        Box a_bx,
                        T a_h,
                        uint a_block);
  template <typename T, MemType MEM>
  inline void nCrossD2n(
                        BoxData<T,DIM,MEM>& a_n_x_d2n,
                        Box a_bx,
                        T a_h,
                        uint a_block);
  template <typename T, unsigned int C, MemType MEM>
  inline
  void cellMatrixQuotientMatAndInv(
                                   BoxData<T,C,MEM,1>& a_quot,
                                   BoxData<T,C,MEM,1>& a_quotBar,
                                   const BoxData<T,C,MEM,C>& a_denominator,
                                   const BoxData<T,C,MEM,C>& a_denomInv,
                                   const BoxData<T,C,MEM,1>& a_numerator);
  template<typename T,MemType MEM,unsigned int C, unsigned int D,unsigned int E>
  inline void 
  quadratureTemplate(
                     BoxData<T,C,MEM,D,E>& a_output,
                     Box a_bx,
                     T a_h,
                     unsigned int a_centering);
  template<typename T,MemType MEM,unsigned int C, unsigned int D,unsigned int E>
  inline void scaled2DTo3D_template(
                                    BoxData<T,C,MEM,D,E>& a_output,
                                    BoxData<T,1,MEM,1,1>& a_scale,
                                    Box a_bx,
                                    T a_h,
                                    unsigned int a_block,
                                    unsigned int a_centering,
                                    bool a_useScaling = false);
  template<typename T,MemType MEM,unsigned int C, unsigned int D,unsigned int E>
  inline void scaledDetA(
                         BoxData<T,C,MEM,D,E>& a_output,
                         const BoxData<T,1,MEM,1,1>& a_scale,
                         Box a_bx,
                         T a_h,
                         unsigned int a_block,
                         unsigned int a_centering,
                         bool a_useScaling = false);
  template<typename T,MemType MEM,unsigned int C, unsigned int D,unsigned int E>
  inline void 
  quadrature_detA(
                  BoxData<T,C,MEM,D,E>& a_output,
                  Box a_bx,
                  T a_h,
                  unsigned int a_centering);   
#include "implem/Proto_CubedSphereShellImplem.H"
} // end namespace Proto
