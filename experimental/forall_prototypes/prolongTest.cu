//compile with 
// nvcc -DTHRUST_DEBUG -g -G -std=c++11 prolongTest.cu
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

///this test is mainly to test if stencils work properly with a non-trivial dstRefRatio
int main(int argc, char** argv) 
{
  int refrat = 2;

  Box refbox = Box::Cube(refrat);
  int numcolors = refbox.size();
  Stencil<double> proSten[numcolors];
  Point           pcolors[numcolors];
  int icolor = 0;
  for(BoxIterator boxit = refbox.begin(); boxit != refbox.end(); ++boxit)
  {
    proSten[icolor]  =  (1.0)*Shift(Point::Zeros());
    proSten[icolor].destRatio() = Point::Ones(refrat);
    pcolors[icolor] = *boxit;
    proSten[icolor].destShift() = pcolors[icolor];
    icolor++;
  }

  constexpr int nx = 8;
  Box fineDom = Box::Cube(nx);
  Box coarDom = fineDom.coarsen(2);
  BoxData<double> fineDat(fineDom);
  BoxData<double> coarDat(coarDom);
  coarDat.setVal(7.0);
  fineDat.setVal(1.23456789e10);

/*  for(int icolor = 0; icolor < numcolors;  icolor++)
  {
    proSten[icolor].cudaApplyBF(coarDat, fineDat, coarDom, true, 1.0);
  }
  cout << "BF  fine data (should be 7)  max = " << fineDat.max() << ", min = " << fineDat.min() << endl;;


  fineDat.setVal(1.23456789e10);*/
  for(int icolor = 0; icolor < numcolors;  icolor++)
  {
    proSten[icolor].cudaApply(coarDat, fineDat,  coarDom, true, 1.0);
  }
  cout << "new fine data (should be 7)  max = " << fineDat.max() << ", min = " << fineDat.min() << endl;;


  return 0;
}
