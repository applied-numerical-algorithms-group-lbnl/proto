#ifndef _PROTO_MBMAP_CUBED_SPHERE_SHELL_
#define _PROTO_MBMAP_CUBED_SPHERE_SHELL_

#define CUBED_SPHERE_SHELL_R0 0.1
#define CUBED_SPHERE_SHELL_R1 1.5
#define CUBED_SPHERE_SHELL_RADIAL_COORD 0

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
};

}
#endif //end include guard
