#pragma once

#include "Proto.H"
#include "TestFunctions.H"

#define CUBED_SPHERE_SHELL_R0 1.0
#define CUBED_SPHERE_SHELL_R1 1.1
#define CUBED_SPHERE_SHELL_RADIAL_COORD 0

using namespace Proto;

template<MemType MEM>
class MBMap_CubeSphereShell : public MBMapOp<MEM>
{
    public:

    using MBMapOp<MEM>::MBMapOp;

    inline void apply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J,
            FluxBoxData<double, DIM, MEM>& a_NT)
    {
        auto& dx = this->dx();
        auto block = this->block();
#if 1
        auto C2C = Stencil<double>::CornersToCells(4);
        int r_dir = CUBED_SPHERE_SHELL_RADIAL_COORD;
        double r0 = CUBED_SPHERE_SHELL_R0;
        double r1 = CUBED_SPHERE_SHELL_R1;
        Box b = C2C.domain(a_X.box());
        BoxData<double, DIM, MEM> X0(b);
        BoxData<double, 1, MEM> J0(b);
        forallInPlace_p(f_cubeSphereMap, X0, J0, dx, r0, r1, block, r_dir);
        X0.copyTo(a_X);
        a_J |= C2C(J0);
#else
        BoxData<double, 1, MEM> radius(a_X.box());
        forallInPlace_p(
            [] PROTO_LAMBDA (Point& a_pt, Var<double, 1, MEM>& a_R, Array<double, DIM> a_dx)
            {
                a_R(0) = a_pt[0]*a_dx[0] + CUBED_SPHERE_SHELL_R0;
            }, radius, dx);
        Operator::cubedSphereGeometry(a_X, a_NT, a_J, radius, a_J.box().grow(Point::Ones()), dx[1], block);
#endif
    }
};

template<MemType MEM>
class MBMap_CubeSphereShellPolar : public MBMapOp<MEM>
{
    public:

    using MBMapOp<MEM>::MBMapOp;

    inline void apply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J,
            FluxBoxData<double, DIM, MEM>& a_NT)
    {
        auto& dx = this->dx();
        auto block = this->block();
#if 1
        int r_dir = CUBED_SPHERE_SHELL_RADIAL_COORD;
        double r0 = CUBED_SPHERE_SHELL_R0;
        double r1 = CUBED_SPHERE_SHELL_R1;
        BoxData<double, DIM, MEM> X0(a_X.box().grow(1));
        BoxData<double, 1, MEM> J0(X0.box());
        BoxData<double, DIM, MEM> X1(a_X.box().grow(1));
        BoxData<double, 1, MEM> J1(X1.box());
        forallInPlace_p(f_cubeSphereMap, X0, J0, dx, r0, r1, block, r_dir);
        int focalBlock = this->focalBlock();
        forallInPlace_p(f_cartesianToSphericalMap, X1, J1, X0, focalBlock, r_dir); 
        
        J0.copyTo(a_J);
        a_J *= J1;
        X1.copyTo(a_X);
#else
        
        BoxData<double, 1, MEM> radius(a_X.box());
        forallInPlace_p(
            [] PROTO_LAMBDA (Point& a_pt, Var<double, 1, MEM>& a_R, Array<double, DIM> a_dx)
            {
                a_R(0) = a_pt[0]*a_dx[0] + CUBED_SPHERE_SHELL_R0;
            }, radius, dx);
        Array<BoxData<double,DIM,MEM>, DIM> NT;
        Operator::cubedSphereGeometry(a_X, NT, a_J, radius, a_J.box().grow(Point::Ones()), dx[1], a_block);
#endif
    }
};

template<typename T, unsigned int C, MemType MEM, Centering CTR>
void interpBoundaries_CubeSphereShell(
        MBLevelBoxData<T, C, MEM, CTR>& a_data,
        unsigned int a_order = 4)
{
    std::vector<Point> footprint;
    Box K = Box::Kernel(2);
    for (auto pi : K)
    {
        // Assumes radial coordinate is the first
        if (pi.abs().sum() <= 2 && pi[CUBED_SPHERE_SHELL_RADIAL_COORD] == 0)
        {
            footprint.push_back(pi);
        }
    }
    MBLevelMap<MBMap_CubeSphereShellPolar, MEM> polarMaps[6];
    MBInterpOp op(a_data.ghost()[0], a_order);
    for (int bi = 0; bi < 6; bi++)
    {
        polarMaps[bi].define(a_data.layout(), a_data.ghost(), bi);
        op.define(polarMaps[bi], footprint, bi); 
    }
    op.apply(a_data, a_data);
}
