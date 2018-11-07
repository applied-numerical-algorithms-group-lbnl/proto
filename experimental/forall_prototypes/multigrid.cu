#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>
#include "../../include/Proto.H"


using namespace std;
using namespace Proto;

typedef Var<double,   1> Scalar;


PROTO_KERNEL_START void initParabolaT(Point& p, Var<double>& data)
{
  data(0) = 0;
  for(int idir = 0; idir < DIM; idir ++)
  {
    data(0) += p[idir]*p[idir];
  }
}
PROTO_KERNEL_END(initParabolaT, initParabola);
/****************/
PROTO_KERNEL_START void setRHSF(Point&   a_p,  Var<double>& a_rhs)
{
  a_rhs(0) = 1.;
}
PROTO_KERNEL_END(setRHSF, setRHS) 

/****************/
void
multigridSolve()
{
  int nx = 16;
  Point lo = Point::Zeros();
  Point hi = Point::Ones(nx - 1);
  Bx domain(lo, hi);
  BoxData<double> rhs = forall_p<double>(setRHS, domain);
//  BoxData<double> rhs(domain);
//  forallInPlace_p(initParabola, domain, rhs);


  cout << "after setting rhs max  =  "<< rhs.max() << ", min = "<< rhs.min() << endl;

#ifdef PROTO_CUDA
  cudaError err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr, "cudaGetLastError() failed at %s:%i : %s\n",
            __FILE__, __LINE__, cudaGetErrorString(err));
  }
#endif

/**/
}
int main(int argc, char* argv[])
{

  multigridSolve();

}  
