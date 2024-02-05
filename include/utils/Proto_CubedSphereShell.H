#pragma once
#ifndef _PROTO_CUBED_SPHERE_SHELL_H_
#define _PROTO_CUBED_SPHERE_SHELL_H_
#include "Proto.H"
#include "Proto_CubedSphereKernels.H"

using namespace Proto;
namespace Proto
{

namespace {
    ACCEL_DECORATION
    void f_cubedSphereMap(
            Array<double, DIM>& a_X,
            double& a_J,
            const Array<double, DIM>& a_X0,
            int a_block)
    {
        constexpr int r_dir = CUBED_SPHERE_SHELL_RADIAL_COORD;
        constexpr double r0 = CUBED_SPHERE_SHELL_R0;
        constexpr double r1 = CUBED_SPHERE_SHELL_R1;
        double r, xi, eta;
        switch (r_dir)
        {
            case 0:
                {
                    r   = r0 + (r1-r0)*a_X0[0];
                    xi  = -M_PI/4.0 + M_PI/2.0*a_X0[1];
                    eta = -M_PI/4.0 + M_PI/2.0*a_X0[2];
                    break;
                }
            case 2:
                {
                    r   = r0 + (r1-r0)*a_X0[2];
                    xi  = -M_PI/4.0 + M_PI/2.0*a_X0[0];
                    eta = -M_PI/4.0 + M_PI/2.0*a_X0[1];
                    break;
                }
        }
        double X0 = tan(xi);
        double Y0 = tan(eta);
        double d0 = sqrt(1+X0*X0+Y0*Y0);
        a_J = r*r*(1+X0*X0)*(1+Y0*Y0)/(d0*d0*d0);
        switch (a_block)
        {
            case 0:
                a_X[0] =  r*X0/d0;
                a_X[1] = -r*Y0/d0;
                a_X[2] = -r/d0;

                break;
            case 1:
                a_X[0] = r*X0/d0;
                a_X[1] = r*Y0/d0;
                a_X[2] = r/d0;
                break;
            case 2:
                a_X[0] =  r*X0/d0;
                a_X[1] = -r/d0;
                a_X[2] =  r*Y0/d0;
                break;
            case 3:
                a_X[0] = r/d0;
                a_X[1] = r*X0/d0;
                a_X[2] = r*Y0/d0;
                break;
            case 4:
                a_X[0] = -r*X0/d0;
                a_X[1] =  r/d0;
                a_X[2] =  r*Y0/d0;
                break;
            case 5:
                a_X[0] = -r/d0;
                a_X[1] = -r*X0/d0;
                a_X[2] =  r*Y0/d0;
                break;
            default:
                a_X[0] = 0;
                a_X[1] = 0;
                a_X[2] = 0;
                break;
        }
    }
    
    ACCEL_DECORATION
    void f_cartesianToPolar(
            Array<double, DIM>&       a_R,
            double&                   a_J,
            const Array<double, DIM>& a_X,
            int                       a_focalBlock)
    {
        constexpr int r_dir = CUBED_SPHERE_SHELL_RADIAL_COORD;
        
        double x0 = a_X[0];
        double y0 = a_X[1];
        double z0 = a_X[2];

        double x, y, z;
        switch (a_focalBlock)
        {
            case 0:
                {
                    x = -z0;
                    y = y0;
                    z = x0;
                    break;
                }
            case 1:
                {
                    x = z0;
                    y = y0;
                    z = -x0;
                    break;
                }
            case 2:
                {
                    x = -y0;
                    y = x0;
                    z = z0;
                    break;
                }
            case 3:
                {
                    x = x0;
                    y = y0;
                    z = z0;
                    break;
                }
            case 4:
                {
                    x = y0;
                    y = -x0;
                    z = z0;
                    break;
                }
            case 5:
                {
                    x = -x0;
                    y = -y0;
                    z = z0;
                    break;
                }
        }

        double R = sqrt(x*x + y*y + z*z);
        double r = sqrt(x*x + y*y);
        double phi = atan2(r,z);
        double theta = atan2(y,x);
        
        a_J = r*r*sin(phi);
        a_R[r_dir] = R;
        a_R[(r_dir + 1) % 3] = theta;
        a_R[(r_dir + 2) % 3] = phi;
    }
}

template<MemType MEM>
class MBMap_CubedSphereShell : public MBMapOp<MEM>
{
public:

    using MBMapOp<MEM>::MBMapOp;

    inline void apply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J,
            FluxBoxData<double, DIM, MEM>& a_NT)
    {
        auto& dx = this->dx();
        auto X0 = this->X(a_X.box(), dx);
        auto block = this->block();
        forallInPlace(
            [] PROTO_LAMBDA (
                Var<double, DIM, MEM>& a_X, Var<double>& a_J, Var<double, DIM, MEM>& a_X0,
                int a_block)
            {
                Array<double, DIM> X0;
                Array<double, DIM> X;
                double J;
                for (int dir = 0; dir < DIM; dir++) { X0[dir] = a_X0(dir); };
                f_cubedSphereMap(X, J, X0, a_block);
                for (int dir = 0; dir < DIM; dir++) { a_X(dir) = X[dir]; }
                a_J(0) = J;
            }, a_X, a_J, X0, block);
    }
    
    inline Array<double, DIM> apply(const Array<double, DIM>& a_Xi)
    {
        Array<double, DIM> X;
        double J;
        auto block = this->block();
        f_cubedSphereMap(X, J, a_Xi, block);
        return X;
    }
};

template<MemType MEM>
    class MBMap_CubedSphereShellPolar : public MBMapOp<MEM>
{
public:

    using MBMapOp<MEM>::MBMapOp;

    inline void apply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J,
            FluxBoxData<double, DIM, MEM>& a_NT)
    {
        Box b = a_X.box() & a_J.box();
        auto block = this->block();
        auto& dx = this->dx();
        auto X0 = this->X(b, dx);
        BoxData<double, DIM, MEM> Xc(b);
        BoxData<double, 1, MEM> Jc(b);
        
        // Mapped -> Cartesian
        forallInPlace(
            [] PROTO_LAMBDA (
                Var<double, DIM, MEM>& a_X, Var<double>& a_J, Var<double, DIM, MEM>& a_X0,
                int a_block)
            {
                Array<double, DIM> X0;
                Array<double, DIM> X;
                double J;
                for (int dir = 0; dir < DIM; dir++) { X0[dir] = a_X0(dir); };
                f_cubedSphereMap(X, J, X0, a_block);
                for (int dir = 0; dir < DIM; dir++) { a_X(dir) = X[dir]; }
                a_J(0) = J;
            }, Xc, Jc, X0, block);

        // Cartesian -> Spherical Polar
        int focalBlock = this->focalBlock();
        forallInPlace(
            [] PROTO_LAMBDA (
                Var<double, DIM, MEM>& a_X, Var<double>& a_J, Var<double, DIM, MEM>& a_X0,
                int a_focalBlock)
            {
                Array<double, DIM> X0;
                Array<double, DIM> X;
                double J;
                for (int dir = 0; dir < DIM; dir++) { X0[dir] = a_X0(dir); };
                f_cartesianToPolar(X, J, X0, a_focalBlock);
                for (int dir = 0; dir < DIM; dir++) { a_X(dir) = X[dir]; }
                a_J(0) = J;
            }, a_X, a_J, Xc, focalBlock);
        
        a_J *= Jc;
    }
    
    inline Array<double, DIM> apply(const Array<double, DIM>& a_Xi)
    {
        Array<double, DIM> Xc;
        Array<double, DIM> X;
        double J;
        auto block = this->block();
        auto focalBlock = this->focalBlock();
        f_cubedSphereMap(Xc, J, a_Xi, block);
        f_cartesianToPolar(X, J, Xc, focalBlock);
        return X;
    }
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
}; 
  template <typename T, MemType MEM>
  inline
  void cellMatrixQuotientMatAndInv(
                                   BoxData<T,DIM,MEM,1>& a_quot,
                                   BoxData<T,DIM,MEM,1>& a_quotBar,
                                   const BoxData<T,DIM,MEM,DIM>& a_denominator,
                                   const BoxData<T,DIM,MEM,DIM>& a_denomInv,
                                   const BoxData<T,DIM,MEM,1>& a_numerator);
  template <typename T, MemType MEM>
  inline
  void cellQuotientValAndInv(
                             BoxData<T,DIM,MEM,1>& a_quot,
                             BoxData<T,DIM,MEM,1>& a_quotBar,
                             const BoxData<T,DIM,MEM,DIM>& a_denominator,
                             const BoxData<T,DIM,MEM,DIM>& a_denomInv,
                             const BoxData<T,DIM,MEM,1>& a_numerator);
  
  template<typename T,MemType MEM>
  inline BoxData<T,1,MEM>
  detA(
       Box a_bx,
       T a_h,
       unsigned int a_centering);   
#include "implem/Proto_CubedSphereShellImplem.H"
} // end namespace Proto
#endif // end include guard
