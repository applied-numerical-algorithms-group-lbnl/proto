#pragma once
#include "Proto.H"

#define CUBED_SPHERE_SHELL_R0 1.0
#define CUBED_SPHERE_SHELL_R1 1.1
#define CUBED_SPHERE_SHELL_RADIAL_COORD 0

namespace Proto
{

template<MemType MEM>
PROTO_KERNEL_START
void f_cubeSphereMapF(Point& a_pt, Var<double,DIM,MEM>& a_X, Var<double,1,MEM>& a_J,
        Array<double, DIM> a_dx, double a_r0, double a_r1, int a_block, int a_radialCoord)
{
    // subscript 0 are at corners
    Array<double, DIM> x  = a_pt;
    x += 0.5;
    x *= a_dx;
    Array<double, DIM> x0 = x - (0.5*a_dx);

    double r0, xi0, eta0, r, xi, eta;
    switch (a_radialCoord)
    {
        case 0:
        {
            r0   = a_r0 + (a_r1-a_r0)*x0[0];
            xi0  = -M_PI/4.0 + M_PI/2.0*x0[1];
            eta0 = -M_PI/4.0 + M_PI/2.0*x0[2];
            break;
        }
        case 2:
        {
            r0   = a_r0 + (a_r1-a_r0)*x0[2];
            xi0  = -M_PI/4.0 + M_PI/2.0*x0[0];
            eta0 = -M_PI/4.0 + M_PI/2.0*x0[1];
            break;
        }
    }
    
    double X0 = tan(xi0);
    double Y0 = tan(eta0);
    double d0 = sqrt(1+X0*X0+Y0*Y0);
    
    double X = tan(xi);
    double Y = tan(eta);
    double d = sqrt(1+X*X+Y*Y);

    a_J(0) = r0*r0*(1+X0*X0)*(1+Y0*Y0)/(d0*d0*d0);
    switch (a_block)
    {
        case 0:
            a_X(0) = r0*X0/d0;
            a_X(1) = -r0*Y0/d0;
            a_X(2) = -r0/d0;
            break;
        case 1:
            a_X(0) = r0*X0/d0;
            a_X(1) = r0*Y0/d0;
            a_X(2) = r0/d0;
            break;
        case 2:
            a_X(0) = r0*X0/d0;
            a_X(1) = -r0/d0;
            a_X(2) = r0*Y0/d0;
            break;
        case 3:
            a_X(0) = r0/d0;
            a_X(1) = r0*X0/d0;
            a_X(2) = r0*Y0/d0;
            break;
        case 4:
            a_X(0) = -r0*X0/d0;
            a_X(1) = r0/d0;
            a_X(2) = r0*Y0/d0;
            break;
        case 5:
            a_X(0) = -r0/d0;
            a_X(1) = -r0*X0/d0;
            a_X(2) = r0*Y0/d0;
            break;
        default:
            a_X(0) = 0;
            a_X(1) = 0;
            a_X(2) = 0;
            break;
    }
}
PROTO_KERNEL_END(f_cubeSphereMapF, f_cubeSphereMap);

template<MemType MEM>
PROTO_KERNEL_START
void f_cartesianToSphericalMapF(Point& a_pt, Var<double,DIM,MEM>& a_R, Var<double,1,MEM>& a_J,
        Var<double, DIM, MEM>& a_X, int a_block, int a_radialCoord)
{
    double x0 = a_X(0);
    double y0 = a_X(1);
    double z0 = a_X(2);

    double x, y, z;
    switch (a_block)
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
    // note that this J is at corners
    a_J(0) = r*r*sin(phi);
    a_R(a_radialCoord) = R;
    a_R((a_radialCoord + 1) % 3) = theta;
    a_R((a_radialCoord + 2) % 3) = phi;
}
PROTO_KERNEL_END(f_cartesianToSphericalMapF, f_cartesianToSphericalMap);

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
                [] PROTO_LAMBDA (
                    Point& a_pt, Var<double, 1, MEM>& a_R, Array<double, DIM> a_dx)
                {
                    a_R(0) = a_pt[0]*a_dx[0] + CUBED_SPHERE_SHELL_R0;
                }, radius, dx);
        Operator::cubedSphereGeometry(
                a_X, a_NT, a_J, radius, a_J.box().grow(Point::Ones()), dx[1], block);
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
                [] PROTO_LAMBDA (
                    Point& a_pt, Var<double, 1, MEM>& a_R, Array<double, DIM> a_dx)
                {
                a_R(0) = a_pt[0]*a_dx[0] + CUBED_SPHERE_SHELL_R0;
                }, radius, dx);
        Array<BoxData<double,DIM,MEM>, DIM> NT;
        Operator::cubedSphereGeometry(
                a_X, NT, a_J, radius, a_J.box().grow(Point::Ones()), dx[1], a_block);
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


class CubedSphereShell 
{
    public:

    static MBProblemDomain Domain(int a_domainSize, int a_thickness, int a_radialDir = 0)
    {
        MBProblemDomain domain(6);
        auto CCW = CoordPermutation::ccw(a_radialDir);
        auto CW = CoordPermutation::cw(a_radialDir);
        CoordPermutation RNorth, RSouth, RI;
        Point x = Point::X();
        Point y = Point::Y();
        Point z = Point::Z();
        for (int bi = 2; bi < 6; bi++)
        {
            int srcBlock = bi;
            int dstBlock = bi+1;
            if (dstBlock > 5) { dstBlock = 2; }
            switch (a_radialDir)
            {
                case 0:
                    domain.defineBoundary(srcBlock, dstBlock, y, RI);
                    domain.defineBoundary(srcBlock, 1, z, RNorth);
                    domain.defineBoundary(srcBlock, 0, -z, RSouth);
                    break;
                case 2:
                    domain.defineBoundary(srcBlock, dstBlock, x, RI);
                    domain.defineBoundary(srcBlock, 1, y, RNorth);
                    domain.defineBoundary(srcBlock, 0, -y, RSouth);
                    break;
                default:
                    MayDay<void>::Abort("NOT IMPLEMENTED");
            }
            RNorth = RNorth*CW;
            RSouth = RSouth*CCW;
        }
        domain.graph().fixRotations();
        Point domainSize = Point::Ones(a_domainSize);
        domainSize[a_radialDir] = a_thickness;
        for (int bi = 0; bi < 6; bi++)
        {
            domain.defineDomain(bi, domainSize);
        }
        return domain;
    }

    template<MemType MEM = MEMTYPE_DEFAULT>
    static MBLevelMap<MBMap_CubeSphereShell, MEM> Map(
            const MBDisjointBoxLayout& a_layout,
            const Array<Point, DIM+1> a_ghost)
    {
        MBLevelMap<MBMap_CubeSphereShell, HOST> map;
        map.define(a_layout, a_ghost);
        return map;
    }
    
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    static MBLevelMap<MBMap_CubeSphereShell, MEM> Map(
            const MBLevelBoxData<T,C,MEM,CTR>& a_data)
    {
        MBLevelMap<MBMap_CubeSphereShell, HOST> map;
        map.define(a_data.layout(), a_data.ghost());
        return map;
    }

    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    static void InterpBoundaries(
            MBLevelBoxData<T,C,MEM,CTR>& a_data,
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
};

} // end namespace Proto
