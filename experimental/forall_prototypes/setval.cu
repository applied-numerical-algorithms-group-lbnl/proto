#define DIM 2
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
using namespace Proto;
int main(int argc, char* argv[])
{

  int size1D  = 4;

  Point lo = Point::Zeros();
  Point hi = Point::Ones(size1D - 1);
  Box finedom(lo,hi);
  BoxData<int> finebd(finedom);
    
  finebd.setVal(-4);

  int maxfine = finebd.absMax();
  std::cout << "fine absMax value (should be  4)= ";
  std::cout << finebd.absMax() << std::endl;
  std::cout << "fine max    value (should be -4)= ";
  std::cout << finebd.max() << std::endl;
  std::cout << "fine min    value (should be -4)= ";
  std::cout << finebd.min() << std::endl;

  Box coardom = finedom.coarsen(2);
  BoxData<int, 1> coarbd(coardom);  
  coarbd.setVal(-2);

  std::cout << "coar absMax value (should be  2)= " ;
  std::cout << coarbd.absMax() << std::endl;
  std::cout << "coar max    value (should be -2)= " ;
  std::cout << coarbd.max() << std::endl;
  std::cout << "coar min    value (should be -2)= " ;
  std::cout << coarbd.min() << std::endl;

#ifndef PROTO_CUDA
  int* dataPtrCoar = coarbd.data(lo);
  int* dataPtrFine = finebd.data(lo);
  int loCoar = *dataPtrCoar;
  int loFine = *dataPtrFine;
  std::cout << "coar value at low (should be -2)= " << loCoar  << std::endl;
  std::cout << "fine value at low (should be -4)= " << loFine  << std::endl;
#endif

  BoxData<int> addedUp(finedom);
  addedUp.setVal(0.);
  addedUp += 4586;
  std::cout << "added           max, min (should be 4586)= " ;
  std::cout << addedUp.min() << " , " << addedUp.max() << std::endl;

  addedUp -= 4588;
  std::cout << "subracted    max, min value (should be -2)= " ;
  std::cout << addedUp.min() << " , " << addedUp.max() << std::endl;


  addedUp *= -7;
  std::cout << "multiplied        max, min  (should be 14)= " ;
  std::cout << addedUp.min() << " , " << addedUp.max() << std::endl;

  addedUp /= 2;
  std::cout << "divided            max, min  (should be 7)= " ;
  std::cout << addedUp.min() << " , " << addedUp.max() << std::endl;

  
  return 0;
    
}
