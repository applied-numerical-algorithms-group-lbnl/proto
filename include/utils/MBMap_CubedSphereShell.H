#ifndef _PROTO_MBMAP_CUBED_SPHERE_SHELL_
#define _PROTO_MBMAP_CUBED_SPHERE_SHELL_

#if 1
#define CUBED_SPHERE_SHELL_R0 1495978707000.0
#define CUBED_SPHERE_SHELL_R1 17951744484000.0
#define CUBED_SPHERE_SHELL_RADIAL_COORD 0
#else
#define CUBED_SPHERE_SHELL_R0 0.1
#define CUBED_SPHERE_SHELL_R1 2.0
#define CUBED_SPHERE_SHELL_RADIAL_COORD 0
#endif
#define CUBED_SPHERE_RADIAL_STRETCH 1.5
#include <cmath> //for atan2

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
        double r, xi, eta, C_rad, R_t;
        C_rad = CUBED_SPHERE_RADIAL_STRETCH;
        R_t = (r1 - r0) / (exp(C_rad) - 1.0);
        switch (r_dir)
        {
            case 0:
                {
                    // r   = r0 + (r1-r0)*a_X0[0];
				    r = r0 + R_t * (exp(C_rad * a_X0[0]) - 1.0);
                    xi  = -M_PI/4.0 + M_PI/2.0*a_X0[1];
                    eta = -M_PI/4.0 + M_PI/2.0*a_X0[2];
                    break;
                }
            case 2:
                {
                    // r   = r0 + (r1-r0)*a_X0[2];
				    r = r0 + R_t * (exp(C_rad * a_X0[2]) - 1.0);
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
                // xi = atan2(-x,z);
                // eta = atan2(y,z);
                xi = atan2(x,-z);
                eta = atan2(-y,-z); 
                break;
            case 1:
                xi = atan2(x,z);
                eta = atan2(y,z);
                break;
            case 2:
                // xi = atan2(-x,y);
                // eta = atan2(-z,y);
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
                // xi = atan2(y,x);
                // eta = atan2(-z,x);
                xi = atan2(-y,-x);
                eta = atan2(z,-x);
                break;
            default:
                xi = 0;
                eta = 0;
                break;
        }
        // r = (r - r0)/(r1-r0);
        double C_rad = CUBED_SPHERE_RADIAL_STRETCH;
        double R_t = (r1 - r0)/(exp(C_rad) - 1.0);
        r = (1.0/C_rad)*log(1.0+(r - r0)/R_t);
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
}

/// Cubed Sphere Map
/** Use this for the MAP template parameter for Cubed Sphere calculations
 */
template<MemType MEM>
class MBMap_CubedSphereShell : public MBMapOp<MEM>
{
public:
  using MBMapOp<MEM>::MBMapOp;
  Array<Array<uint,DIM>,6> m_permvals
        =  {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
  Array<Array<int,DIM>,6> m_signvals
        = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};
  Array<unsigned int,6> m_blockNormal = {2,2,1,0,1,0};
  
  inline Array<double,DIM> MBPointToNormal(MBPoint a_mbpt,
                                            double a_offset) const
  {
    Point pt = a_mbpt.point;
    int block = a_mbpt.block;
    auto perm = m_permvals[block];
    auto sign = m_signvals[block];
    Array<double,DIM> nlocal,normal;
    double dxtan = (this -> dx())[2];
    double angle1 = .5*M_PI*(pt[1]*dxtan + a_offset*dxtan - .5);
    double angle2 = .5*M_PI*(pt[2]*dxtan + a_offset*dxtan - .5);
    double tan1 = tan(angle1);
    double tan2 = tan(angle2);
    double normalization  = 1./sqrt(1.0 + tan1*tan1 + tan2*tan2);
    nlocal[0] = normalization;
    nlocal[1] = tan1*normalization;
    nlocal[2] = tan2*normalization;
    for (int dir = 0; dir < DIM; dir++)
    {
      normal[dir] = nlocal[perm[dir]]*sign[dir];
    }
    return normal;
  }
  inline MBPoint normalToMBPoint(Array<double,DIM> a_normal,
                                 MBProblemDomain a_domain) const
  {
    unsigned int retblock = 0;
    double offset = .5;
    Array<double,DIM> retNormal;
    
    for (unsigned int block = 0; block < 6; block++)
      {
        auto perm = m_permvals[block];
        auto sign = m_signvals[block];
        auto signBlockNormal = sign[m_blockNormal[retblock]];
        for (int dir = 0; dir < DIM; dir++)
          {
            retNormal[perm[dir]] = sign[dir]*a_normal[dir];
          }
        if ((abs(retNormal[1]) <= abs(retNormal[0])) &&
            (abs(retNormal[2]) <= abs(retNormal[0])) &&
            (a_normal[m_blockNormal[retblock]]*signBlockNormal > 0)
            )
          {
            break;
          }
        retblock++;
      }
     PROTO_ASSERT((retblock < 6),
     "MBMap_CubedSphereShell::normalToMBPoint | Error: normal not in any block");
    double dxi = (this->dx())[2];
    double angle1 = atan(retNormal[1]/retNormal[0]);
    double angle2 = atan(retNormal[2]/retNormal[0]);
    double xi1 = angle1*(2.0/M_PI) + .5;
    double xi2 = angle2*(2.0/M_PI) + .5;
    Point pt;   
    pt[0] = 0;
    pt[1] = xi1/dxi;
    pt[2] = xi2/dxi;
    Point highBlock = a_domain.getBlock(retblock).box().high();
    Point lowBlock = a_domain.getBlock(retblock).box().low();
    pt[1] = max(min(pt[1],highBlock[1]),lowBlock[1]);
    pt[2] = max(min(pt[2],highBlock[2]),lowBlock[2]);
    return MBPoint(pt,retblock);
  }
  inline MBPoint findValid(MBProblemDomain a_domain,
                                   MBPoint a_mbpt) const
  {
    auto block = a_mbpt.block;
    Box domainBox = a_domain.getBlock(block).box();
    Point pt = a_mbpt.point;
    //cout << "domainBox = " << domainBox << endl;
    //cout << "pt = " << pt << endl;
    if (domainBox.containsPoint(pt))
      {
        return a_mbpt;
      }
    else
      {
        auto normal = MBPointToNormal(a_mbpt,.5);
        //cout << "normal = " << normal << endl;
        auto retval = normalToMBPoint(normal,a_domain);
        //cout << "retval  = " << retval.point << " , " << retval.block << endl;
        retval.point[0] = a_mbpt.point[0];
        return retval;
      }
  }
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
  Array<Array<uint,DIM>,6> m_permvals
        =  {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
  Array<Array<int,DIM>,6> m_signvals
        = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};
  Array<unsigned int,6> m_blockNormal = {2,2,1,0,1,0};
  
  inline Array<double,DIM> MBPointToNormal(MBPoint a_mbpt,
                                            double a_offset) const
  {
    Point pt = a_mbpt.point;
    int block = a_mbpt.block;
    auto perm = m_permvals[block];
    auto sign = m_signvals[block];
    Array<double,DIM> nlocal,normal;
    double dxtan = (this -> dx())[2];
    double angle1 = .5*M_PI*(pt[1]*dxtan + a_offset*dxtan - .5);
    double angle2 = .5*M_PI*(pt[2]*dxtan + a_offset*dxtan - .5);
    double tan1 = tan(angle1);
    double tan2 = tan(angle2);
    double normalization  = 1./sqrt(1.0 + tan1*tan1 + tan2*tan2);
    nlocal[0] = normalization;
    nlocal[1] = tan1*normalization;
    nlocal[2] = tan2*normalization;
    for (int dir = 0; dir < DIM; dir++)
    {
      normal[dir] = nlocal[perm[dir]]*sign[dir];
    }
    return normal;
  }
  inline MBPoint normalToMBPoint(Array<double,DIM> a_normal,
                                 MBProblemDomain a_domain) const
  {
    unsigned int retblock = 0;
    double offset = .5;
    Array<double,DIM> retNormal;
    
    for (unsigned int block = 0; block < 6; block++)
      {
        auto perm = m_permvals[block];
        auto sign = m_signvals[block];
        auto signBlockNormal = sign[m_blockNormal[retblock]];
        for (int dir = 0; dir < DIM; dir++)
          {
            retNormal[perm[dir]] = sign[dir]*a_normal[dir];
          }
        if ((abs(retNormal[1]) <= abs(retNormal[0])) &&
            (abs(retNormal[2]) <= abs(retNormal[0])) &&
            (a_normal[m_blockNormal[retblock]]*signBlockNormal > 0)
            )
          {
            break;
          }
        retblock++;
      }
    
  
     PROTO_ASSERT((retblock < 6),
     "MBMap_CubedSphereShell::normalToMBPoint | Error: normal not in any block");
    double dxi = (this->dx())[2];
    double angle1 = atan(retNormal[1]/retNormal[0]);
    double angle2 = atan(retNormal[2]/retNormal[0]);
    double xi1 = angle1*(2.0/M_PI) + .5;
    double xi2 = angle2*(2.0/M_PI) + .5;
    Point pt;   
    pt[0] = 0;
    pt[1] = xi1/dxi;
    pt[2] = xi2/dxi;
    Point highBlock = a_domain.getBlock(retblock).box().high();
    Point lowBlock = a_domain.getBlock(retblock).box().low();
    pt[1] = max(min(pt[1],highBlock[1]),lowBlock[1]);
    pt[2] = max(min(pt[2],highBlock[2]),lowBlock[2]);
    return MBPoint(pt,retblock);
  }
  inline MBPoint findValid(MBProblemDomain a_domain,
                                   MBPoint a_mbpt) const override
  {
    auto block = a_mbpt.block;
    Box domainBox = a_domain.getBlock(block).box();
    Point pt = a_mbpt.point;
    //cout << "srcvalExtracted  = " << pt << " , " << block << endl;
    if (domainBox.containsPoint(pt))
      {
        return a_mbpt;
      }
    else
      {
        //cout << "srcval  = " << a_mbpt.point << " , " << a_mbpt.block << endl;
        auto normal = MBPointToNormal(a_mbpt,.5);
        //cout << "normal = " << normal << endl;
        auto retval = normalToMBPoint(normal,a_domain);
        //cout << "retval  = " << retval.point << " , " << retval.block << endl;
        retval.point[0] = a_mbpt.point[0];
        return retval;
      }
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
