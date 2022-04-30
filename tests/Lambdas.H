#ifndef TEST_LAMBDAS
#define TEST_LAMBDAS
#include "Proto.H"

#define COMPS 2

PROTO_KERNEL_START
void f_rampHostF(const Point& a_pt, Var<int,COMPS,HOST>& a_data)
{
    Point x = a_pt + Point::Ones();
    for (int comp = 0; comp < COMPS; comp++)
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
PROTO_KERNEL_END(f_rampHostF, f_rampHost);

PROTO_KERNEL_START
void
f_rampDeviF(const Point& a_pt, Var<int,COMPS,DEVICE>& a_data)
{
    Point x = a_pt + Point::Ones();
    for (int comp = 0; comp < COMPS; comp++)
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
PROTO_KERNEL_END(f_rampDeviF, f_rampDevi);

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

#endif
