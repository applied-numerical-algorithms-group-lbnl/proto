#pragma once

#include "Proto.H"
#include "Lambdas.H"

#define CUBE_SPHERE_SHELL_R0 1.0
#define CUBE_SPHERE_SHELL_R1 1.1

using namespace Proto;

template<MemType MEM>
class MBLevelMap_CubeSphereShell : public MBLevelMap<MEM>
{
    public:

    using MBLevelMap<MEM>::MBLevelMap;

    private:

    inline void apply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J,
            FluxBoxData<double, DIM, MEM>& a_NT,
            unsigned int a_block)
    {
        auto& layout = this->layout();
        auto& dx = this->dx(a_block);
#if 1
        double r0 = CUBE_SPHERE_SHELL_R0;
        double r1 = CUBE_SPHERE_SHELL_R1;
        forallInPlace_p(f_cubeSphereMap, a_X, a_J, dx, r0, r1, a_block);
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
