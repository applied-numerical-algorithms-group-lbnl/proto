#include "Proto.H"

#define NUMCOMPS 1

using namespace Proto;

template<MemType MEM>
PROTO_KERNEL_START
void f_waveF(Point& a_pt, Var<double, NUMCOMPS, MEM>& a_data, double a_dx, int a_k = 1)
{
    double x[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        x[ii] = a_pt[ii]*a_dx + a_dx / 2.0;
    }
    
    double k = 2.0*M_PI*a_k;
    a_data(0) = sin(k*(x[0] + x[1]));
}
PROTO_KERNEL_END(f_waveF, f_wave);

template<MemType MEM>
PROTO_KERNEL_START
void f_LwaveF(Point& a_pt, Var<double, NUMCOMPS, MEM>& a_data, double a_dx, int a_k = 1)
{
    double x[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        x[ii] = a_pt[ii]*a_dx + a_dx / 2.0;
    }
    
    double k = 2.0*M_PI*a_k;
    a_data(0) = -2*k*k*sin(k*(x[0] + x[1]));
}
PROTO_KERNEL_END(f_LwaveF, f_Lwave);


// f(xi) = sum(xi)
template<MemType MEM>
PROTO_KERNEL_START
void rampF(const Point& a_pt, Var<double,NUMCOMPS, MEM>& a_data, double a_dx, double a_offset = 0.0)
{
    a_data(0) = 0.0;
    for (int ii = 0; ii < DIM; ii++)
    {
        a_data(0) += a_pt[ii]*a_dx + a_dx/2.0 - a_offset;
    }
}
PROTO_KERNEL_END(rampF, ramp);

template<MemType MEM>
PROTO_KERNEL_START
void pointIDF(const Point& a_pt, Var<int,NUMCOMPS, MEM>& a_data)
{
    Point x = a_pt + Point::Ones();
    for (int comp = 0; comp < NUMCOMPS; comp++)
    {
        a_data(comp) = (comp+1)*10 + x[0];
#if DIM > 1
        a_data(comp) = (comp+1)*100 + 10*x[0] + x[1];
#endif
#if DIM > 2
        a_data(comp) = (comp+1)*1000 + 100*x[0] + 10*x[1] + x[2];
#endif
    }
}
PROTO_KERNEL_END(pointIDF, pointID);

template<MemType MEM>
PROTO_KERNEL_START
void sinProdF(const Point& a_pt, Var<double, 1, MEM> a_data, double a_dx)
{
    double x0[DIM];
    double x1[DIM];
    
    double a = 0.125;
    for (int dir = 0; dir < DIM; dir++)
    {
        x0[dir] = a_pt[dir]*a_dx + a;
        x1[dir] = x0[dir] + a_dx;
    }
    
    double k = M_PI*2;
    a_data(0) = sin(k*x1[0])*sin(k*x1[1]);
}
PROTO_KERNEL_END(sinProdF, sinProd);

template<unsigned int C, MemType MEM>
PROTO_KERNEL_START
void sinProd_avgF(const Point& a_pt, Var<double, C, MEM> a_data, double a_dx, double a_offset = 0.0)
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
PROTO_KERNEL_END(sinProd_avgF, sinProd_avg);

template<unsigned int C, MemType MEM>
PROTO_KERNEL_START
void LsinProd_avgF(const Point& a_pt, Var<double, C, MEM> a_data, double a_dx, double a_offset = 0.0)
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
PROTO_KERNEL_END(LsinProd_avgF, LsinProd_avg);

template<MemType MEM>
PROTO_KERNEL_START
void sinx_avgF(const Point& a_pt, Var<double, 1, MEM> a_data, double a_dx, double a_offset = 0.0)
{
    double x0 = a_pt[0]*a_dx + a_offset;
    double x1 = x0 + a_dx;

    a_data(0) = (-cos(x1) + cos(x0))/a_dx;
}
PROTO_KERNEL_END(sinx_avgF, sinx_avg);
