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

  int size1D  = 64;

  Point lo = Point::Zeros();
  Point hi = Point::Ones(size1D - 1);
  Bx desmene(lo,hi);

  BoxData<double, 3> doublebd(desmene);
  BoxData<int   , 3>    intbd(desmene);
    
  doublebd.setVal(4.0);
  intbd.setVal(4);



  return 0;
    
}
