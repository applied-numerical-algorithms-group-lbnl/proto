#ifndef _PROTO_MBMAP_CUBED_SPHERE_SHELL_
#define _PROTO_MBMAP_CUBED_SPHERE_SHELL_

#include <cmath> //for atan2
#if 0
#define CUBED_SPHERE_SHELL_R0 1495978707000.0
#define CUBED_SPHERE_SHELL_R1 17951744484000.0
#define CUBED_SPHERE_SHELL_RADIAL_COORD 0
#else
#define CUBED_SPHERE_SHELL_R0 0.9
#define CUBED_SPHERE_SHELL_R1 1.0
#define CUBED_SPHERE_SHELL_RADIAL_COORD 2
#endif
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
    void f_cubedSphereInverseMap(
            Array<double, DIM>& a_XMap,
            const Array<double, DIM>& a_XCart,
            int a_block)
    {
        constexpr int r_dir = CUBED_SPHERE_SHELL_RADIAL_COORD;
        constexpr double r0 = CUBED_SPHERE_SHELL_R0;
        constexpr double r1 = CUBED_SPHERE_SHELL_R1;
        double xi = 0;
        double eta = 0;
        double x = a_XCart[0];
        double y = a_XCart[1];
        double z = a_XCart[2];
        double r = sqrt(x*x + y*y + z*z);
        switch (a_block)
        {
            case 0:
                xi = atan2(x,-z);
                eta = atan2(-y,-z); 
                break;
            case 1:
                xi = atan2(x,z);
                eta = atan2(y,z);
                break;
            case 2:
                xi = atan2(x,-y);
                eta = atan2(z,-y); 
                break;
            case 3:
                xi = atan2(y,x);
                eta = atan2(z,x);
                break;
            case 4:
                xi = atan2(-x,y);
                eta = atan2(z,y);
                break;
            case 5:
                xi = atan2(-y,-x);
                eta = atan2(z,-x);
                break;
            default:
                xi = 0;
                eta = 0;
                break;
        }
        r = (r - r0)/(r1-r0);
        xi = (xi + M_PI/4.0)*(2.0/M_PI);
        eta = (eta + M_PI/4.0)*(2.0/M_PI);
        switch (r_dir)
        {
            case 0:
                {
                    a_XMap[0] = r;
                    a_XMap[1] = xi;
                    a_XMap[2] = eta;
                    break;
                }
            case 2:
                {
                    a_XMap[0] = xi;
                    a_XMap[1] = eta;
                    a_XMap[2] = r;
                    break;
                }
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

    ACCEL_DECORATION
    void f_polarToCartesian(
            Array<double, DIM>&       a_X,
            const Array<double, DIM>& a_R,
            int                       a_focalBlock)
    {
        constexpr int r_dir = CUBED_SPHERE_SHELL_RADIAL_COORD;
        
        double R = a_R[r_dir];
        double theta = a_R[(r_dir + 1) % 3];
        double phi = a_R[(r_dir + 2) % 3];

        double x = R*sin(phi)*cos(theta);
        double y = R*sin(phi)*sin(theta);
        double z = R*cos(phi);

        double x0, y0, z0;
        switch (a_focalBlock)
        {
            case 0:
                {
                    z0 = -x;
                    y0 = y;
                    x0 = z;
                    break;
                }
            case 1:
                {
                    z0 = x;
                    y0 = y;
                    x0 = -z;
                    break;
                }
            case 2:
                {
                    y0 = -x;
                    x0 = y;
                    z0 = z;
                    break;
                }
            case 3:
                {
                    x0 = x;
                    y0 = y;
                    z0 = z;
                    break;
                }
            case 4:
                {
                    y0 = x;
                    x0 = -y;
                    z0 = z;
                    break;
                }
            case 5:
                {
                    x0 = -x;
                    y0 = -y;
                    z0 = z;
                    break;
                }
        }
        
        a_X[0] = x0;
        a_X[1] = y0;
        a_X[2] = z0;
    }
}

/// Cubed Sphere Map
/** Use this for the MAP template parameter for Cubed Sphere calculations
 */
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
        PROTO_ASSERT(a_X.box() == a_J.box(),
                "MBMap_CubedSphereShell::apply | Error: J and X should both have the same (node centered) box");
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
    
    inline void inverse(
            BoxData<double, DIM, MEM>& a_XMap,
            const BoxData<double, DIM, MEM>& a_XCart)
    {
        auto block = this->block();
        forallInPlace(
            [] PROTO_LAMBDA (
                Var<double, DIM, MEM>& a_X, Var<double, DIM, MEM>& a_X0,
                int a_block)
            {
                Array<double, DIM> X0;
                Array<double, DIM> X;
                for (int dir = 0; dir < DIM; dir++) { X0[dir] = a_X0(dir); };
                f_cubedSphereInverseMap(X, X0, a_block);
                for (int dir = 0; dir < DIM; dir++) { a_X(dir) = X[dir]; }
            }, a_XMap, a_XCart, block);
    }
    
    inline Array<double, DIM> inverse(const Array<double, DIM>& a_XCart)
    {
        Array<double, DIM> XMap;
        auto block = this->block();
        f_cubedSphereInverseMap(XMap, a_XCart, block);
        return XMap;
    }
};

/// Polar Cubed Sphere Map
/** Maps from the cubed sphere mapped space onto spherical polar coordinates centered on
 *  a specified block. This map is ONLY used for building the Cubed Sphere interpolation
 *  operators and is not recommended for public use.
 */
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
        PROTO_ASSERT(a_X.box() == a_J.box(),
                "MBMap_CubedSphereShell::apply | Error: J and X should both have the same (node centered) box");
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

    inline void inverse(
        BoxData<double, DIM, MEM>& a_XMap,
        const BoxData<double, DIM, MEM>& a_R)
    {
        auto block = this->block();
        auto focus = this->focalBlock();
        BoxData<double, DIM, MEM> XCart(a_R.box());
        forallInPlace(
            [] PROTO_LAMBDA (
                Var<double, DIM, MEM>& a_X, Var<double, DIM, MEM>& a_R,
                int a_focus)
            {
                Array<double, DIM> X;
                Array<double, DIM> R;
                for (int dir = 0; dir < DIM; dir++) { R[dir] = a_R(dir); };
                f_polarToCartesian(X, R, a_focus);
                for (int dir = 0; dir < DIM; dir++) { a_X(dir) = X[dir]; }
            }, XCart, a_R, focus);
        forallInPlace(
            [] PROTO_LAMBDA (
                Var<double, DIM, MEM>& a_X, Var<double, DIM, MEM>& a_XCart,
                int a_block)
            {
                Array<double, DIM> X0;
                Array<double, DIM> X;
                for (int dir = 0; dir < DIM; dir++) { X[dir] = a_XCart(dir); };
                f_cubedSphereInverseMap(X0, X, a_block);
                for (int dir = 0; dir < DIM; dir++) { a_X(dir) = X0[dir]; }
            }, a_XMap, XCart, block);
    }
    
    inline Array<double, DIM> inverse(const Array<double, DIM>& a_R)
    {
        Array<double, DIM> XMap;
        Array<double, DIM> X;
        auto focus = this->focalBlock();
        auto block = this->block();
        f_polarToCartesian(X, a_R, focus);
        f_cubedSphereInverseMap(XMap, X, block);
        return XMap;
    }
};

}
#endif //end include guard
