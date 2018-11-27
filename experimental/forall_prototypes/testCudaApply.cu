//compile with 
// nvcc -DTHRUST_DEBUG -g -G -std=c++11 testCudaApply.cu
#define DIM 2
#define PROTO_CUDA 1
#include "../../include/Proto.H"
#include <cstdlib>
#include <cstdio>
#include <functional>
#include <iostream>

using namespace Proto;
using std::cout;
using std::endl;

PROTO_KERNEL_START void initParabolaT(Point& p, Var<double>& data)
{
  data(0) = 0;
  for(int idir = 0; idir < DIM; idir ++)
  {
    data(0) += p[idir]*p[idir];
  }
}
PROTO_KERNEL_END(initParabolaT, initParabola);

PROTO_KERNEL_START void initBogosityT(Point& p, Var<double>& data)
{
  data(0) = 123456789e10;
}
PROTO_KERNEL_END(initBogosityT, initBogosity);

int main(int argc, char** argv) 
{
  constexpr int nx = 4;
  Box domain = Box::Cube(nx);
  Box grrdom = domain.grow(1);
  BoxData<double> phi = forall_p<double>(initParabola, grrdom);
  BoxData<double> lap = forall_p<double>(initBogosity, domain);


  Stencil<double> lapsten = Stencil<double>::Laplacian();
  
  cout << "stencil given by:" << endl;
  lapsten.print();
  lapsten.cudaApplyBF(phi, lap, domain, true, 1.0);
  cout << "BF  laplacian (should be 2*DIM) max = " << lap.max() << ", min = " << lap.min() << endl;;


  lap.setVal(1.23456789e10);
  lapsten.cudaApply(phi, lap, domain, true, 1.0);
  cout << "new laplacian (should be 2*DIM) max = " << lap.max() << ", min = " << lap.min() << endl;;


  return 0;
}
