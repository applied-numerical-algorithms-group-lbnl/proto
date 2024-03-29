#ifndef TEST_LAMBDAS
#define TEST_LAMBDAS
#include <cmath>
#include "Proto.H"

#define XPOINT_SIZE 3
#define XPOINT_RADIUS 2.0
#define CUBED_SPHERE_R0 1.0
#define CUBED_SPHERE_R1 1.1
#define SHEAR_SLOPE 1.0
using namespace Proto;

template<typename T, unsigned int C, MemType MEM>
PROTO_KERNEL_START
void f_MBConstF(Point& a_pt, Var<T,C,MEM>& a_data, unsigned int a_block, T a_value)
{
    for (int comp = 0; comp < C; comp++)
    {
        a_data(comp) = a_value;
    }
}
PROTO_KERNEL_END(f_MBConstF, f_MBConst);

template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_iotaCornerF(Point& a_pt, Var<T,DIM,MEM>& a_X, T a_dx)
{
    for (int dir = 0; dir < DIM; dir++)
    {
        a_X(dir) = a_pt[dir]*a_dx;
    }
}
PROTO_KERNEL_END(f_iotaCornerF, f_iotaCorner);

template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_iotaCenterF(Point& a_pt, Var<T,DIM,MEM>& a_X, T a_dx)
{
    for (int dir = 0; dir < DIM; dir++)
    {
        a_X(dir) = a_pt[dir]*a_dx + a_dx/2.0;
    }
}
PROTO_KERNEL_END(f_iotaCenterF, f_iotaCenter);

#ifdef PR_MMB

ACCEL_DECORATION
double MBPointID(unsigned int a_block, unsigned int a_comp, Point a_pt)
{
    Point x = a_pt;
#if DIM == 1
        return x[0]*1e-2;
#endif
#if DIM == 2
        return x[0]*1e-2 + x[1]*1e-4;
#endif
#if DIM > 2
        //return x[0]*1e-0 + x[1]*1e-1 + x[2]*1e-2;
        return x[1];
#endif
}

template<typename T, unsigned int C, MemType MEM>
PROTO_KERNEL_START
void f_MBPointIDF(Point& a_pt, Var<T,C,MEM>& a_data, unsigned int a_block)
{
    for (int comp = 0; comp < C; comp++)
    {
        //a_data(comp) = MBPointID(a_block, comp, a_pt);
        a_data(comp) = a_pt[comp];
    }
}
PROTO_KERNEL_END(f_MBPointIDF, f_MBPointID);

#ifdef DIM == 3
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

            //r   = a_r0 + (a_r1-a_r0)*x[0];
            //xi  = -M_PI/4.0 + M_PI/2.0*x[1];
            //eta = -M_PI/4.0 + M_PI/2.0*x[2];
            break;
        }
        case 2:
        {
            r0   = a_r0 + (a_r1-a_r0)*x0[2];
            xi0  = -M_PI/4.0 + M_PI/2.0*x0[0];
            eta0 = -M_PI/4.0 + M_PI/2.0*x0[1];

            //r   = a_r0 + (a_r1-a_r0)*x[2];
            //xi  = -M_PI/4.0 + M_PI/2.0*x[0];
            //eta = -M_PI/4.0 + M_PI/2.0*x[1];
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
    //a_R(0) = x;
    //a_R(1) = y;
    //a_R(2) = z;
}
PROTO_KERNEL_END(f_cartesianToSphericalMapF, f_cartesianToSphericalMap);

PROTO_KERNEL_START
void f_classicSphereMapF(Point& a_pt, Var<double,3>& a_X, Var<double,1>& a_J,
        Array<double, 3> a_dx, double a_r0, double a_r1)
{
    Array<double, DIM> x = a_pt;
    x += 0.5;
    x *= a_dx;
    Array<double, DIM> x0 = x - (0.5*a_dx);
    Array<double, DIM> x1 = x0 + a_dx;

    double r      = a_r0 + (a_r1-a_r0)*x[0];
    double phi    = M_PI/4.0 + M_PI/2.0*x[1];
    double theta  = 2.0*M_PI*x[2];
    
    double r0      = a_r0 + (a_r1-a_r0)*x0[0];
    double phi0    = M_PI/4.0 + M_PI/2.0*x0[1];
    double theta0  = 2.0*M_PI*x0[2];
    
    double r1      = a_r0 + (a_r1-a_r0)*x1[0];
    double phi1    = M_PI/4.0 + M_PI/2.0*x1[1];
    double theta1  = 2.0*M_PI*x1[2];
    
    a_J(0) = r*r*sin(phi);
    
    a_X(0) = r0*cos(theta0)*sin(phi0);
    a_X(1) = r0*sin(theta0)*sin(phi0);
    a_X(2) = r0*cos(phi0);

    //a_N(0,0) =  r0*sin(phi)*sin(theta);
    //a_N(1,0) = -r0*sin(phi)*cos(theta);
    //a_N(2,0) =  r0*cos(phi)*cos(theta);
    //a_N(0,1) = -r*cos(phi0)*cos(theta);
    //a_N(1,1) =  r*cos(phi0)*sin(theta);
    //a_N(2,1) =  0;
    //a_N(0,2) =  -r*sin(phi)*cos(theta0);
    //a_N(1,2) =  -r*sin(phi)*sin(theta0);
    //a_N(2,2) =  -r*sin(phi);
}
PROTO_KERNEL_END(f_classicSphereMapF, f_classicSphereMap);

MBProblemDomain buildIdentity(int a_domainSize)
{
    MBProblemDomain domain(1);
    domain.defineDomain(0, Point::Ones(a_domainSize));
    return domain;
}

MBProblemDomain buildXPoint(int a_domainSize, int a_numBlocks)
{
    MBProblemDomain domain(a_numBlocks);
    auto CCW = CoordPermutation::ccw();
    for (int ii = 0; ii < a_numBlocks; ii++)
    {
        domain.defineBoundary(ii, (ii+1) % a_numBlocks, 0, Side::Hi, CCW);
    }
    for (int bi = 0; bi < a_numBlocks; bi++)
    {
        domain.defineDomain(bi, Point::Ones(a_domainSize));
    }
    return domain;
}

MBProblemDomain buildShear(int a_domainSize)
{
    MBProblemDomain domain(4);
    auto I = CoordPermutation::identity();
    domain.defineBoundary(0,1,0,Side::Hi,I);
    domain.defineBoundary(1,2,1,Side::Hi,I);
    domain.defineBoundary(2,3,0,Side::Lo,I);
    domain.defineBoundary(3,0,1,Side::Lo,I);

    for (int ii = 0; ii < 4; ii++)
    {
        domain.defineDomain(ii, Point::Ones(a_domainSize));
    }
    return domain;
}
#endif

MBProblemDomain buildRing(Point a_domainSize)
{
    MBProblemDomain domain(4);
    CoordPermutation I;
    for (int bi = 0; bi < 4; bi++)
    {
        domain.defineBoundary(bi, (bi+1)%4, Point::X(), I);
    }
    domain.graph().fixRotations();
    for (int bi = 0; bi < 4; bi++)
    {
        domain.defineDomain(bi, a_domainSize);
    }
    return domain;
}

#if DIM > 2
MBProblemDomain buildCubeSphere(int a_domainSize)
{
    MBProblemDomain domain(6);
    auto R_theta = CoordPermutation::ccw(1);
    auto R_north = CoordPermutation::cw(0);
    auto R_south = CoordPermutation::ccw(0);
    CoordPermutation R_I;
    Point x = Point::Basis(0);
    Point y = Point::Basis(1);
    for (int bi = 2; bi < 6; bi++)
    {
        int srcBlock = bi;
        int dstBlock = bi+1;
        if (dstBlock > 5) { dstBlock = 2; }
        domain.defineBoundary(srcBlock, dstBlock, x, R_I);
        domain.defineBoundary(srcBlock, 1, y, R_north);
        domain.defineBoundary(srcBlock, 0, -y, R_south);
        R_north = R_north*R_theta;
        R_south = R_south*R_theta;
    }
    domain.graph().fixRotations();
    for (int bi = 0; bi < 6; bi++)
    {
        domain.defineDomain(bi, Point::Ones(a_domainSize));
    }
    return domain;
}
MBProblemDomain buildCubeSphereShell(int a_domainSize, int a_thickness, int a_radialDir = 0)
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
MBProblemDomain buildPolarShell(int a_domainSize, int a_thickness = 1)
{
    MBProblemDomain domain(4);
   
    CoordPermutation I;
    for (int bi = 0; bi < 4; bi++)
    {
        domain.defineBoundary(bi, (bi+1)%4, Point::Z(), I);
    }
    domain.graph().fixRotations();
    Point domainSize(a_thickness, a_domainSize, a_domainSize);
    for (int bi = 0; bi < 4; bi++)
    {
        domain.defineDomain(bi, domainSize);
    }
    return domain;
}
#endif


template<MemType MEM>
PROTO_KERNEL_START
void f_XPointMapF(Point& a_pt, Var<double,3,MEM>& a_data, unsigned int a_block, int numBlocks, int a_radius)
{
    double dtheta_mapped = 2*M_PI/numBlocks;
    double dtheta_unmapped = 2*M_PI/4;
    double theta_ratio = dtheta_mapped/dtheta_unmapped;
    Point ri = Point::Ones(a_radius) - a_pt;
    double r = max(ri[0], ri[1])*(1.0/a_radius);
    double theta = 3*M_PI/2.0 + atan2(ri[1]*1.0, ri[0]*1.0)*theta_ratio + dtheta_mapped*a_block;

    a_data(0) = r*cos(theta);
    a_data(1) = r*sin(theta);
    a_data(2) = 0;
}
PROTO_KERNEL_END(f_XPointMapF, f_XPointMap);

template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_JCubeSphere_(Point& a_pt, Var<T,1,MEM>& a_J, int a_block, Point a_blockSize)
{
    double r0 = CUBED_SPHERE_R0 - 0.1;
    double r1 = CUBED_SPHERE_R0 + 0.1;

    double L = -M_PI/4;
    double H = M_PI/4;
  
    double dX = 1.0/((double)a_blockSize[0]);
    double dY = 1.0/((double)a_blockSize[1]);
    double dZ = 1.0/((double)a_blockSize[2]);

    double X = a_pt[0]*dX + dX/2.0;
    double Y = a_pt[1]*dX + dY/2.0;
    double Z = a_pt[2]*dX + dZ/2.0;

    double xi = L + Y*(H - L);
    double eta = L + Z*(H - L);
    double R = r0 + X*(r1 - r0);

    double u = tan(xi);
    double v = tan(eta);
    double r = sqrt(1+u*u + v*v);

    double dR_dZ = (r1-r0);
    double du_dY = (H-L)/(cos(xi)*cos(xi));
    double dv_dZ = (H-L)/(cos(eta)*cos(eta));
    double drinv_dY = -u*(H-L)/(cos(eta)*cos(eta)*r*r*r);
    double drinv_dZ = -v*(H-L)/(cos(eta)*cos(eta)*r*r*r);

    double dx_dX =  dR_dZ*u/r;
    double dy_dX = -dR_dZ*v/r;
    double dz_dX =  dR_dZ/r;

    double dx_dY =  R*(drinv_dY + du_dY);
    double dy_dY = -R*v*(drinv_dY);
    double dz_dY =  R*(drinv_dY);
    
    double dx_dZ =  R*u*(drinv_dZ);
    double dy_dZ = -R*(drinv_dZ + dv_dZ);
    double dz_dZ =  R*(drinv_dZ);

    a_J(0) = dx_dX*dy_dY*dz_dZ + dx_dY*dy_dZ*dz_dX + dx_dZ*dy_dX*dz_dY
           - dx_dZ*dy_dY*dz_dX - dx_dX*dy_dZ*dz_dY - dx_dY*dy_dX*dz_dZ;

    a_J(0) *= -(dX*dY*dZ);
}
PROTO_KERNEL_END(f_JCubeSphere_, f_JCubeSphere);
#endif

ACCEL_DECORATION
int pointID(unsigned int a_comp, Point a_pt)
{
    Point x = a_pt + Point::Ones();
#if DIM == 1
        return (a_comp+1)*100 + x[0];
#endif
#if DIM == 2
        return (a_comp+1)*10000 + 100*x[0] + x[1];
#endif
#if DIM > 2
        return (a_comp+1)*100000 + 10000*x[0] + 100*x[1] + x[2];
#endif
}

template<typename T, unsigned int C, MemType MEM>
PROTO_KERNEL_START
void f_pointIDF(Point& a_pt, Var<T,C,MEM>& a_data)
{
    for (int comp = 0; comp < C; comp++)
    {
        a_data(comp) = pointID(comp, a_pt);
    }
}
PROTO_KERNEL_END(f_pointIDF, f_pointID);

PROTO_KERNEL_START
void halfPointSumF(Point p, Var<double> v)
{
  v(0) = 0;
  for (int ii = 0; ii < DIM; ii++)
  {
    v(0) += 0.5*p[ii];
  }
}
PROTO_KERNEL_END(halfPointSumF, halfPointSum)

PROTO_KERNEL_START
void pointSumF(Point p, Var<double> v)
{
  v(0) = 0;
  for (int ii = 0; ii < DIM; ii++)
  {
    v(0) += p[ii];
  }
}
PROTO_KERNEL_END(pointSumF, pointSum)

PROTO_KERNEL_START
void twicePointSumF(Point p, Var<double> v)
{
  v(0) = 0;
  for (int ii = 0; ii < DIM; ii++)
  {
    v(0) += 2.*p[ii];
  }
}
PROTO_KERNEL_END(twicePointSumF, twicePointSum)

PROTO_KERNEL_START
void cosxCosyFuncF(Point p, Var<double> v, double dx)
{
  double x = p[0]*dx;
  double y = p[1]*dx;
  v(0) = cos(x)*cos(y);
}
PROTO_KERNEL_END(cosxCosyFuncF,cosxCosyFunc)

PROTO_KERNEL_START
void cosxCosyPCosFuncF(Point p, Var<double> v, double dx)
{
  double x = p[0]*dx/2.0;
  double y = p[1]*dx/2.0;
  v(0) = cos(x)*cos(y) + cos(x);
}
PROTO_KERNEL_END(cosxCosyPCosFuncF,cosxCosyPCosFunc)

PROTO_KERNEL_START
void cosxFuncF(Point p, Var<double> v, double dx)
{
  double x = p[0]*dx/2.0;
  v(0) = cos(x);
}
PROTO_KERNEL_END(cosxFuncF,cosxFunc)

template<unsigned int C, MemType MEM>
PROTO_KERNEL_START
void f_polyMF(const Point& a_pt, Var<double, C, MEM>& a_data, unsigned int a_block, 
        Var<double, DIM, MEM>& a_X, Array<double, DIM> a_exp, Array<double, DIM> a_offset)
{
    double phi = 1;
    for (int dir = 0; dir < DIM; dir++)
    {
        if (a_exp[dir] != 0)
        {
            phi += pow(a_X(dir) - a_offset[dir], a_exp[dir]);
        }
    }
    for (int ii = 0; ii < C; ii++)
    {
        a_data(ii) = phi;
    }
}
PROTO_KERNEL_END(f_polyMF, f_polyM);

template<unsigned int C, MemType MEM>
PROTO_KERNEL_START
void f_phiMF(const Point& a_pt, Var<double, C, MEM>& a_data, unsigned int a_block, 
        Var<double, DIM, MEM>& a_X, Array<double, DIM> a_k, Array<double, DIM> a_offset)
{
    double phi = 1.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        double k = a_k[dir]*M_PI*2;
        phi *= sin(k*(a_X(dir) + a_offset[dir]));
    }
    for (int ii = 0; ii < C; ii++)
    {
        a_data(ii) = phi;
    }
}
PROTO_KERNEL_END(f_phiMF, f_phiM);

template<unsigned int C, MemType MEM>
PROTO_KERNEL_START
void f_LphiMF(const Point& a_pt, Var<double, C, MEM>& a_data, unsigned int a_block, 
        Var<double, DIM, MEM>& a_X, Array<double, DIM> a_k, Array<double, DIM> a_offset)
{
    double phi = 1.0;
    double k2 = 0.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        double k = a_k[dir]*M_PI*2;
        k2 -= k*k;
        phi *= sin(k*(a_X(dir) + a_offset[dir]));
    }
    for (int ii = 0; ii < C; ii++)
    {
        a_data(ii) = k2*phi;
    }
}
PROTO_KERNEL_END(f_LphiMF, f_LphiM);

template<unsigned int C, MemType MEM>
PROTO_KERNEL_START
void f_phiF(const Point& a_pt, Var<double, C, MEM>& a_data,
        Array<double, DIM> a_dx, Array<double, DIM> a_k, Array<double, DIM> a_offset)
{
    double x[DIM];
    double k[DIM];
    double phi = 1.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        x[dir] = a_pt[dir]*a_dx[dir] + a_dx[dir]/2.0 + a_offset[dir];
        k[dir] = (a_k[dir]*M_PI*2);
        phi *= sin(k[dir]*x[dir]);
    }
    
    for (int ii = 0; ii < C; ii++)
    {
        a_data(ii) = phi;
    }
}
PROTO_KERNEL_END(f_phiF, f_phi);

template<unsigned int C, MemType MEM>
PROTO_KERNEL_START
void f_phi_cornerF(const Point& a_pt, Var<double, C, MEM>& a_data, Array<double, DIM> a_dx,
        Array<double, DIM> a_k, Array<double, DIM> a_offset)
{
    double x[DIM];
    double k[DIM];
    double phi = 1.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        x[dir] = a_pt[dir]*a_dx[dir] + a_offset[dir];
        k[dir] = (a_k[dir]*M_PI*2);
        phi *= sin(k[dir]*x[dir]);
    }
    
    for (int ii = 0; ii < C; ii++)
    {
        a_data(ii) = phi;
    }
}
PROTO_KERNEL_END(f_phi_cornerF, f_phi_corner);

template<unsigned int C, MemType MEM>
PROTO_KERNEL_START
void f_phi_faceF(const Point& a_pt, Var<double, C, MEM>& a_data, double a_dx, Point a_k, Point a_offset, int a_dir)
{
    double x[DIM];
    double k[DIM];
    double phi = 1.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        x[dir] = (a_pt[dir] + a_offset[dir])*a_dx;
        if (dir != a_dir) { x[dir] += (a_dx/2.0);}
        k[dir] = (a_k[dir]*M_PI*2);
        phi *= sin(k[dir]*x[dir]);
    }
    
    for (int ii = 0; ii < C; ii++)
    {
        a_data(ii) = phi;
    }
}
PROTO_KERNEL_END(f_phi_faceF, f_phi_face);

template<unsigned int C, MemType MEM>
PROTO_KERNEL_START
void f_phi_avgF(const Point& a_pt, Var<double, C, MEM>& a_data,
        Array<double, DIM> a_dx, Array<double, DIM> a_k, Array<double, DIM> a_offset)
{
    double x0[DIM];
    double x1[DIM];
    double k[DIM];
    double kdx = 1.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        x0[dir] = a_pt[dir]*a_dx[dir] + a_offset[dir];
        x1[dir] = x0[dir] + a_dx[dir];
        k[dir] = (a_k[dir]*M_PI*2);
        kdx *= k[dir]*a_dx[dir];
    }
    

    for (int ii = 0; ii < C; ii++)
    {
        a_data(ii) = 1.0;
        for (int d = 0; d < DIM; d++)
        {
            a_data(ii) *= (-cos(k[d]*x1[d]) + cos(k[d]*x0[d]));
        }
        a_data(ii) *= 1.0/kdx;
    }
}
PROTO_KERNEL_END(f_phi_avgF, f_phi_avg);

template<unsigned int C, MemType MEM>
PROTO_KERNEL_START
void f_phi_face_avgF(const Point& a_pt, Var<double, C, MEM>& a_data, double a_dx, Point a_k, Point a_offset, int a_dir)
{
    double x0[DIM];
    double x1[DIM];
    double k[DIM];
    double kdx = 1.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        x0[dir] = (a_pt[dir] + a_offset[dir])*a_dx;
        x1[dir] = x0[dir] + a_dx;
        k[dir] = (a_k[dir]*M_PI*2);
        if (dir != a_dir)
        {
            kdx *= k[dir]*a_dx;
        }
    }
    
    for (int ii = 0; ii < C; ii++)
    {
        a_data(ii) = 1.0;
        for (int d = 0; d < DIM; d++)
        {
            if (d == a_dir) {continue;}
            a_data(ii)  *= (-cos(k[d]*x1[d]) + cos(k[d]*x0[d]));
        }
        a_data(ii) *= sin(k[a_dir]*x0[a_dir])/kdx;
    }
}
PROTO_KERNEL_END(f_phi_face_avgF, f_phi_face_avg);

template<unsigned int C, MemType MEM>
PROTO_KERNEL_START
void f_LphiF(const Point& a_pt, Var<double, C, MEM>& a_data,
        Array<double, DIM> a_dx, Array<double, DIM> a_k, Array<double, DIM> a_offset)
{
    double x[DIM];
    double k[DIM];
    double k2 = 0.0;
    double L0 = 1.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        //x[dir] = (a_pt[dir] + a_offset[dir])*a_dx + a_dx/2.0;
        //k[dir] = (a_k[dir]*M_PI*2);
        x[dir] = a_pt[dir]*a_dx[dir] + a_dx[dir]/2.0 + a_offset[dir];
        k[dir] = (a_k[dir]*M_PI*2);
        k2 += k[dir]*k[dir];
        L0 *= sin(k[dir]*x[dir]);
    }
    
    for (int cc = 0; cc < C; cc++)
    {
        a_data(cc) = -k2*L0;
    }
}
PROTO_KERNEL_END(f_LphiF, f_Lphi);

template<unsigned int C, MemType MEM>
PROTO_KERNEL_START
void f_Lphi_cornerF(const Point& a_pt, Var<double, C, MEM>& a_data, double a_dx, Point a_k, Point a_offset)
{
    double x[DIM];
    double k[DIM];
    double k2 = 0.0;
    double L0 = 1.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        x[dir] = (a_pt[dir] + a_offset[dir])*a_dx;
        k[dir] = (a_k[dir]*M_PI*2);
        k2 += k[dir]*k[dir];
        L0 *= sin(k[dir]*x[dir]);
    }
    
    for (int cc = 0; cc < C; cc++)
    {
        a_data(cc) = -k2*L0;
    }
}
PROTO_KERNEL_END(f_Lphi_cornerF, f_Lphi_corner);

template<unsigned int C, MemType MEM>
PROTO_KERNEL_START
void f_Lphi_avgF(const Point& a_pt, Var<double, C, MEM>& a_data, double a_dx, Point a_k, Point a_offset)
{
    double x0[DIM];
    double x1[DIM];
    double k[DIM];
    double kdx = 1.0;
    double k2 = 0.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        x0[dir] = (a_pt[dir] + a_offset[dir])*a_dx;
        x1[dir] = x0[dir] + a_dx;
        k[dir] = (a_k[dir]*M_PI*2);
        kdx *= k[dir]*a_dx;
        k2 += k[dir]*k[dir];
    }
    
    for (int ii = 0; ii < C; ii++)
    {
#if DIM == 2
        a_data(ii) =
            + cos(k[0]*x1[0])*cos(k[1]*x1[1])
            - cos(k[0]*x0[0])*cos(k[1]*x1[1])
            - cos(k[0]*x1[0])*cos(k[1]*x0[1])
            + cos(k[0]*x0[0])*cos(k[1]*x0[1]);
#elif DIM == 3
        a_data(ii) =
            + cos(k[0]*x1[0])*cos(k[1]*x1[1])*cos(k[2]*x1[2])
            - cos(k[0]*x0[0])*cos(k[1]*x1[1])*cos(k[2]*x1[2])
            - cos(k[0]*x1[0])*cos(k[1]*x0[1])*cos(k[2]*x1[2])
            + cos(k[0]*x0[0])*cos(k[1]*x0[1])*cos(k[2]*x1[2])
            - cos(k[0]*x1[0])*cos(k[1]*x1[1])*cos(k[2]*x0[2])
            + cos(k[0]*x0[0])*cos(k[1]*x1[1])*cos(k[2]*x0[2])
            + cos(k[0]*x1[0])*cos(k[1]*x0[1])*cos(k[2]*x0[2])
            - cos(k[0]*x0[0])*cos(k[1]*x0[1])*cos(k[2]*x0[2]);
#else
#endif
        a_data(ii) *= (-k2/kdx);
    }
}
PROTO_KERNEL_END(f_Lphi_avgF, f_Lphi_avg);

#ifdef PR_AMR
AMRGrid telescopingGrid(
        Point crseDomainSize,
        unsigned int numLevels,
        std::vector<Point> refRatios,
        std::vector<Point>   boxSizes,
        std::array<bool, DIM> periodicity)
{
    std::vector<DisjointBoxLayout> layouts;
    layouts.resize(numLevels);
    Box domainBox(crseDomainSize);
    ProblemDomain domain(crseDomainSize, periodicity);
    layouts[0].define(domain, domainBox, boxSizes[0]);
    for (int lvl = 1; lvl < numLevels; lvl++)
    {
        domain = domain.refine(refRatios[lvl-1]);
        domainBox = domainBox.grow(-domainBox.sizes()/4).refine(refRatios[lvl-1]);
        layouts[lvl].define(domain, domainBox, boxSizes[lvl]); 
    }
    return AMRGrid(layouts, refRatios, numLevels);
}

AMRGrid telescopingGrid(
        int domainSize,
        unsigned int numLevels,
        Point refRatio,
        Point boxSize)
{
    std::vector<Point> refRatios(numLevels-1, refRatio);
    std::vector<Point> boxSizes(numLevels, boxSize);
    Point crseDomainSize = Point::Ones(domainSize);
    std::array<bool, DIM> periodicity;
    periodicity.fill(true);
    return telescopingGrid(crseDomainSize, numLevels, refRatios, boxSizes, periodicity);
}
#endif //end if PR_AMR
#endif //end include guard
