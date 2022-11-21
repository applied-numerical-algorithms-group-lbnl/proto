#ifndef TEST_LAMBDAS
#define TEST_LAMBDAS
#include <cmath>
#include "Proto.H"

#define XPOINT_SIZE 5
#define XPOINT_RADIUS 2.0
#define CUBED_SPHERE_R0 1.0
#define CUBED_SPHERE_R1 2.0

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
        return (a_block+1)*1e2 + (a_comp+1)*1e-0 + x[0]*1e-2;
#endif
#if DIM == 2
        return (a_block+1)*1e2 + (a_comp+1)*1e-0 + x[0]*1e-2 + x[1]*1e-4;
#endif
#if DIM > 2
        //return (a_block+1)*1e2 + (a_comp+1)*1e-0 + x[0]*1e-2 + x[1]*1e-4 + x[2]*1e-6;
        return x[0]*1e-0 + x[1]*1e-1 + x[2]*1e-2;
#endif
}

template<typename T, unsigned int C, MemType MEM>
PROTO_KERNEL_START
void f_MBPointIDF(Point& a_pt, Var<T,C,MEM>& a_data, unsigned int a_block)
{
    for (int comp = 0; comp < C; comp++)
    {
        a_data(comp) = MBPointID(a_block, comp, a_pt);
    }
}
PROTO_KERNEL_END(f_MBPointIDF, f_MBPointID);

template<MemType MEM>
PROTO_KERNEL_START
void f_TestMapF(Point& a_pt, Var<double,DIM,MEM>& a_data, unsigned int a_block,
        double a_h, Point a_shift, Point a_k)
{
    double kx[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        kx[ii] = (a_shift[ii] + a_pt[ii])*(a_k[ii]*a_h);
    }
    for (int ii = 0; ii < DIM; ii++)
    {
        double prod = 1.0;
        for (int jj = 0; jj < DIM; jj++)
        {
            prod *= sin(2*M_PI*kx[jj]);
        }
        a_data(ii) = pow(-1,ii)*prod;
    }
}
PROTO_KERNEL_END(f_TestMapF, f_TestMap);

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
    for (int bi = 0; bi < 6; bi++)
    {
        domain.defineDomain(bi, Point::Ones(a_domainSize));
    }
    return domain;
}
#endif

MBProblemDomain buildIdentity(int a_domainSize)
{
    MBProblemDomain domain(1);
    domain.defineDomain(0, Point::Ones(a_domainSize));
    return domain;
}

MBProblemDomain buildXPoint(int a_domainSize)
{
    MBProblemDomain domain(XPOINT_SIZE);
    auto CCW = CoordPermutation::ccw();
    for (int ii = 0; ii < XPOINT_SIZE; ii++)
    {
        domain.defineBoundary(ii, (ii+1) % XPOINT_SIZE, 0, Side::Hi, CCW);
    }
    for (int bi = 0; bi < XPOINT_SIZE; bi++)
    {
        domain.defineDomain(bi, Point::Ones(a_domainSize));
    }
    return domain;
}


template<MemType MEM>
PROTO_KERNEL_START
void f_XPointMapF(Point& a_pt, Var<double,3,MEM>& a_data, unsigned int a_block, int numBlocks, int a_radius)
{
    double dtheta_mapped = 2*M_PI/numBlocks;
    double dtheta_unmapped = 2*M_PI/4;
    double theta_ratio = dtheta_mapped/dtheta_unmapped;
    Point ri = Point::Ones(a_radius) - a_pt;
    double r = max(ri[0], ri[1])*(1.0/a_radius);
    double theta = atan2(ri[1]*1.0, ri[0]*1.0)*theta_ratio + dtheta_mapped*a_block;

    a_data(0) = r*cos(theta);
    a_data(1) = r*sin(theta);
    a_data(2) = 0;
}
PROTO_KERNEL_END(f_XPointMapF, f_XPointMap);

PROTO_MAP_START(IdentityMap)
    x = X; y = Y; z = Z;
PROTO_MAP_END(IdentityMap)

PROTO_MAP_START(XPointMap)
    double numBlocks = XPOINT_SIZE;
    double RMax = XPOINT_RADIUS;
    
    double dtheta_mapped = 2*M_PI/numBlocks;
    double dtheta_unmapped = 2*M_PI/4;
    double theta_ratio = dtheta_mapped/dtheta_unmapped;
    
    double u = (1-X)*sqrt(1.0 - (1-Y)*(1-Y)/2.0);
    double v = (1-Y)*sqrt(1.0 - (1-X)*(1-X)/2.0);
    double r = sqrt(u*u + v*v)*RMax;

    double theta = atan2(v,u)*theta_ratio + dtheta_mapped*block;
   
    x = r*cos(theta);
    y = r*sin(theta);
    z = Z;
PROTO_MAP_END(XPointMap)

template<MemType MEM>
PROTO_KERNEL_START
void f_CubeSphereMapF(Point& a_pt, Var<double,3,MEM>& a_data, unsigned int a_block, Point a_blockSize, double a_r0, double a_r1)
{
    double L = -M_PI/4;
    double H = M_PI/4;
    
    double xi = L + (a_pt[0]*1.0)/(a_blockSize[0]*1.0)*(H - L);
    double eta = L + (a_pt[1]*1.0)/(a_blockSize[1]*1.0)*(H - L);
    double R = a_r0 + (a_pt[2]*1.0)/(a_blockSize[2]*1.0)*(a_r1 - a_r0);
    //double R = a_r0;
    double a = sqrt(3.0)/3.0*R;
    double x = a*tan(xi);
    double y = a*tan(eta);
    double r = sqrt(a*a + x*x + y*y);
    

    switch (a_block)
    {
        case 0:
            a_data(0) = R/r*x;
            a_data(1) = -R/r*y;
            a_data(2) = -R/r*a;
            break;
        case 1:
            a_data(0) = R/r*x;
            a_data(1) = R/r*y;
            a_data(2) = R/r*a;
            break;
        case 2:
            a_data(0) = R/r*x;
            a_data(1) = -R/r*a;
            a_data(2) = R/r*y;
            break;
        case 3:
            a_data(0) = R/r*a;
            a_data(1) = R/r*x;
            a_data(2) = R/r*y;
            break;
        case 4:
            a_data(0) = -R/r*x;
            a_data(1) = R/r*a;
            a_data(2) = R/r*y;
            break;
        case 5:
            a_data(0) = -R/r*a;
            a_data(1) = -R/r*x;
            a_data(2) = R/r*y;
            break;
        default:
            a_data(0) = 0;
            a_data(1) = 0;
            a_data(2) = 0;
            break;
    }
}
PROTO_KERNEL_END(f_CubeSphereMapF, f_CubeSphereMap);

PROTO_MAP_START(CubedSphereMap)
    double r0 = CUBED_SPHERE_R0;
    double r1 = CUBED_SPHERE_R1;

    double L = -M_PI/4;
    double H = M_PI/4;
    
    double xi = L + X*(H - L);
    double eta = L + Y*(H - L);
    double R = r0 + Z*(r1 - r0);
    double a = sqrt(3.0)/3.0*R;
    double u = a*tan(xi);
    double v = a*tan(eta);
    double r = sqrt(a*a + u*u + v*v);
    
    switch (block)
    {
        case 0:
            x =  R/r*u;
            y = -R/r*v;
            z = -R/r*a;
            break;
        case 1:
            x =  R/r*u;
            y =  R/r*v;
            z =  R/r*a;
            break;
        case 2:
            x =  R/r*u;
            y = -R/r*a;
            z =  R/r*v;
            break;
        case 3:
            x =  R/r*a;
            y =  R/r*u;
            z =  R/r*v;
            break;
        case 4:
            x = -R/r*u;
            y =  R/r*a;
            z =  R/r*v;
            break;
        case 5:
            x = -R/r*a;
            y = -R/r*u;
            z =  R/r*v;
            break;
        default:
            x = 0;
            y = 0;
            z = 0;
            break;
    }
PROTO_MAP_END(CubedSphereMap)

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
void f_phiMF(const Point& a_pt, Var<double, C, MEM>& a_data, Var<double, DIM, MEM>& a_X,
        Point a_k, Array<double, DIM> a_offset)
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
void f_phiF(const Point& a_pt, Var<double, C, MEM>& a_data,
        double a_dx, Point a_k, Point a_offset)
{
    double x[DIM];
    double k[DIM];
    double phi = 1.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        x[dir] = (a_pt[dir] + a_offset[dir])*a_dx + a_dx/2.0;
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
void f_phi_cornerF(const Point& a_pt, Var<double, C, MEM>& a_data, double a_dx, Point a_k, Point a_offset)
{
    double x[DIM];
    double k[DIM];
    double phi = 1.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        x[dir] = (a_pt[dir] + a_offset[dir])*a_dx;
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
void f_phi_avgF(const Point& a_pt, Var<double, C, MEM>& a_data, double a_dx, Point a_k, Point a_offset)
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
        kdx *= k[dir]*a_dx;
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
void f_LphiF(const Point& a_pt, Var<double, C, MEM>& a_data, double a_dx, Point a_k, Point a_offset)
{
    double x[DIM];
    double k[DIM];
    double k2 = 0.0;
    double L0 = 1.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        x[dir] = (a_pt[dir] + a_offset[dir])*a_dx + a_dx/2.0;
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
