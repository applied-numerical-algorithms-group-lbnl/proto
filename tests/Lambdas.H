#ifndef TEST_LAMBDAS
#define TEST_LAMBDAS
#include <cmath>
#include "Proto.H"

#ifdef PR_MMB
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

CUDA_DECORATION
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
        return (a_block+1)*1e2 + (a_comp+1)*1e-0 + x[0]*1e-2 + x[1]*1e-4 + x[2]*1e-6;
#endif
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
            a_data(1) = R/r*y;
            a_data(2) = -R/r*a;
            break;
        case 1:
            a_data(0) = R/r*x;
            a_data(1) = R/r*y;
            a_data(2) = R/r*a;
            break;
        case 2:
            a_data(2) = R/r*x;
            a_data(0) = R/r*y;
            a_data(1) = -R/r*a;
            break;
        case 3:
            a_data(2) = R/r*x;
            a_data(1) = R/r*y;
            a_data(0) = -R/r*a;
            break;
        case 4:
            a_data(2) = R/r*x;
            a_data(0) = R/r*y;
            a_data(1) = R/r*a;
            break;
        case 5:
            a_data(2) = R/r*x;
            a_data(1) = R/r*y;
            a_data(0) = R/r*a;
            break;
        default:
            a_data(0) = 0;
            a_data(1) = 0;
            a_data(2) = 0;
            break;
    }

}
PROTO_KERNEL_END(f_CubeSphereMapF, f_CubeSphereMap);
#endif

CUDA_DECORATION
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
void f_phiF(const Point& a_pt, Var<double, C, MEM>& a_data, double a_dx, double a_offset = 0.0)
{
    double x[DIM];
    
    double a = a_offset;
    for (int dir = 0; dir < DIM; dir++)
    {
        x[dir] = a_pt[dir]*a_dx + a + a_dx/2.0;
    }
    
    double k = M_PI*2;
    for (int ii = 0; ii < C; ii++)
    {
        a_data(ii) = sin(k*x[0])*sin(k*x[1]);
    }
}
PROTO_KERNEL_END(f_phiF, f_phi);

template<unsigned int C, MemType MEM>
PROTO_KERNEL_START
void f_phi_avgF(const Point& a_pt, Var<double, C, MEM>& a_data, double a_dx, double a_offset = 0.0)
{
    double x0[DIM];
    double x1[DIM];
    
    for (int dir = 0; dir < DIM; dir++)
    {
        x0[dir] = a_pt[dir]*a_dx + a_offset;
        x1[dir] = x0[dir] + a_dx;
    }
    
    double k = M_PI*2;
    for (int ii = 0; ii < C; ii++)
    {
        a_data(ii) = + cos(k*x1[0])*cos(k*x1[1])
            - cos(k*x0[0])*cos(k*x1[1])
            - cos(k*x1[0])*cos(k*x0[1])
            + cos(k*x0[0])*cos(k*x0[1]);
        a_data(ii) *= 1.0/(k*k*a_dx*a_dx);
    }
}
PROTO_KERNEL_END(f_phi_avgF, f_phi_avg);

template<unsigned int C, MemType MEM>
PROTO_KERNEL_START
void f_LphiF(const Point& a_pt, Var<double, C, MEM>& a_data, double a_dx, double a_offset = 0.0)
{
    double x[DIM];
    
    double a = a_offset;
    for (int dir = 0; dir < DIM; dir++)
    {
        x[dir] = a_pt[dir]*a_dx + a_dx/2.0 + a;
    }
    
    double k = M_PI*2;
    for (int cc = 0; cc < C; cc++)
    {
        a_data(cc) = -2.0*k*k*sin(k*x[0])*sin(k*x[1]);
    }
}
PROTO_KERNEL_END(f_LphiF, f_Lphi);

template<unsigned int C, MemType MEM>
PROTO_KERNEL_START
void f_Lphi_avgF(const Point& a_pt, Var<double, C, MEM>& a_data, double a_dx, double a_offset = 0.0)
{
    double x0[DIM];
    double x1[DIM];
    
    for (int dir = 0; dir < DIM; dir++)
    {
        x0[dir] = a_pt[dir]*a_dx + a_offset;
        x1[dir] = x0[dir] + a_dx;
    }
    
    double k = M_PI*2;
    for (int ii = 0; ii < C; ii++)
    {
        a_data(ii) = + cos(k*x1[0])*cos(k*x1[1])
            - cos(k*x0[0])*cos(k*x1[1])
            - cos(k*x1[0])*cos(k*x0[1])
            + cos(k*x0[0])*cos(k*x0[1]);
        a_data(ii) *= 1.0/(k*k*a_dx*a_dx);
        a_data(ii) *= -DIM*pow(k, 2.0);
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
#endif
#endif
