#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>

#include "Proto.H"

#define PI 3.141592653589793
#define NUMCOMPS DIM+2

using namespace std;
using namespace Proto;

typedef Var<double,DIM> V;
typedef Var<double,NUMCOMPS> State;


//=================================================================================================
PROTO_KERNEL_START
void iotaFuncPointF(Point           & a_p,
               V               & a_X,
               double            a_h)
{
  for (int ii = 0; ii < DIM; ii++)
  {
    a_X(ii) = a_p[ii]*a_h + 0.5*a_h;
  }
}
PROTO_KERNEL_END(iotaFuncPointF,iotaFuncPoint)
//=================================================================================================
PROTO_KERNEL_START
void iotaFuncIntF(int               a_p[DIM],
                  V               & a_X,
                  double            a_h)
{
  for (int ii = 0; ii < DIM; ii++)
  {
    a_X(ii) = a_p[ii]*a_h + 0.5*a_h;
  }
}
PROTO_KERNEL_END(iotaFuncIntF,iotaFuncInt)
/***/
int main(int argc, char* argv[])
{

  int size1D = 16;
  Point lo = Point::Zeros();
  Point hi = Point::Ones(size1D - 1);
  Box dbx0(lo,hi);
  double s_dx = 1./size1D;
  int nGhost = 4;
  Box dbx = dbx0.grow(nGhost);
  Box dbx1 = dbx.grow(1);
  BoxData<double,NUMCOMPS> UBig(dbx1);
  BoxData<double,DIM> x(dbx1);

  unsigned long long int numflops = 0;
  printf(" calling int c-array iota function\n");
  forallInPlaceBaseOp_i(numflops, "int_version", iotaFuncInt,   dbx1, x, s_dx);

  printf(" calling point iota function\n");
  forallInPlaceBaseOp_p(numflops, "point_version",iotaFuncPoint, dbx1, x, s_dx);

  printf("done \n");

}
