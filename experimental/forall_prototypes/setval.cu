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
  Bx finedom(lo,hi);
  BoxData<int> finebd(finedom);
    
  finebd.setVal(-4);

  int maxfine = finebd.absMax();
  std::cout << "fine absMax value (should be  4)= ";
  std::cout << finebd.absMax() << std::endl;
  std::cout << "fine max    value (should be -4)= ";
  std::cout << finebd.max() << std::endl;
  std::cout << "fine min    value (should be -4)= ";
  std::cout << finebd.min() << std::endl;

  Bx coardom = finedom.coarsen(2);
  BoxData<int, 1> coarbd(coardom);  
  coarbd.setVal(-2);

  std::cout << "coar absMax value (should be  2)= " ;
  std::cout << coarbd.absMax() << std::endl;
  std::cout << "coar max    value (should be -2)= " ;
  std::cout << coarbd.max() << std::endl;
  std::cout << "coar min    value (should be -2)= " ;
  std::cout << coarbd.min() << std::endl;

  return 0;
    
}
