#pragma once

#include "Proto.H"
#include "Lambdas.H"

#define CUBE_SPHERE_SHELL_R0 1.0
#define CUBE_SPHERE_SHELL_R1 1.1
#define CUBE_SPHERE_SHELL_RADIAL_COORD 2

using namespace Proto;

template<MemType MEM>
class MBLevelMap_CubeSphereShell : public MBLevelMap<MEM>
{
    public:

    using MBLevelMap<MEM>::MBLevelMap;

    inline void apply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J,
            FluxBoxData<double, DIM, MEM>& a_NT,
            unsigned int a_block)
    {
        auto& layout = this->layout();
        auto& dx = this->dx(a_block);
#if 1
        auto C2C = Stencil<double>::CornersToCells(4);
        int r_dir = CUBE_SPHERE_SHELL_RADIAL_COORD;
        double r0 = CUBE_SPHERE_SHELL_R0;
        double r1 = CUBE_SPHERE_SHELL_R1;
        Box b = C2C.domain(a_X.box());
        BoxData<double, DIM, MEM> X0(b);
        BoxData<double, 1, MEM> J0(b);
        forallInPlace_p(f_cubeSphereMap, X0, J0, dx, r0, r1, a_block, r_dir);
        X0.copyTo(a_X);
        a_J |= C2C(J0);
#else
        
        BoxData<double, 1, MEM> radius(a_X.box());
        forallInPlace_p(
            [] PROTO_LAMBDA (Point& a_pt, Var<double, 1, MEM>& a_R, Array<double, DIM> a_dx)
            {
                a_R(0) = a_pt[0]*a_dx[0] + CUBE_SPHERE_SHELL_R0;
            }, radius, dx);
        Array<BoxData<double,DIM,MEM>, DIM> NT;
        Operator::cubedSphereGeometry(a_X, NT, a_J, radius, a_J.box().grow(Point::Ones()), dx[1], a_block);
#endif
    }
};
