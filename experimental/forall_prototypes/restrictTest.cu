//compile with 
// nvcc -DTHRUST_DEBUG -g -G -std=c++11 restrictTest.cu
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

///this test is mainly to test if stencils work properly with a non-trivial srcRefRatio
int main(int argc, char** argv) 
{
  int refrat = 2;

  Stencil<double> aveSten;
  Box refbox = Box::Cube(refrat);
  int numpts = refbox.size();
  for(BoxIterator boxit = refbox.begin(); boxit != refbox.end(); ++boxit)
  {
    aveSten += (1.0/numpts)*Shift(*boxit);
  }
  aveSten.srcRatio() = Point::Ones(refrat);


  constexpr int nx = 8;
  Box fineDom = Box::Cube(nx);
  Box coarDom = fineDom.coarsen(2);
  BoxData<double> fineDat(fineDom);
  BoxData<double> coarDat(coarDom);
  fineDat.setVal(7.0);
  coarDat.setVal(1.23456789e10);
  
  aveSten.cudaApplyBF(fineDat, coarDat, coarDom, true, 1.0);
  cout << "BF  coar data (should be 7)  max = " << coarDat.max() << ", min = " << coarDat.min() << endl;;


  coarDat.setVal(1.23456789e10);
  aveSten.cudaApply(fineDat, coarDat, coarDom, true, 1.0);
  cout << "new coar data (should be 7)  max = " << coarDat.max() << ", min = " << coarDat.min() << endl;;


  return 0;
}
