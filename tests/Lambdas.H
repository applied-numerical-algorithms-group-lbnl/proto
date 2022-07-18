#ifndef TEST_LAMBDAS
#define TEST_LAMBDAS
#include "Proto.H"

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

CUDA_DECORATION
int MBPointID(unsigned int a_block, unsigned int a_comp, Point a_pt)
{
    Point x = a_pt + Point::Ones();
#if DIM == 1
        return (a_block+1)*10000 + (a_comp+1)*100 + x[0];
#endif
#if DIM == 2
        return (a_block+1)*1000000 + (a_comp+1)*10000 + 100*x[0] + x[1];
#endif
#if DIM > 2
        return (a_block+1)*100000000 + (a_comp+1)*100000 + 10000*x[0] + 100*x[1] + x[2];
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

#endif
