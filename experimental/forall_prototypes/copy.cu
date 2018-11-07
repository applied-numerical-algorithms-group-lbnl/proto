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
/****************/

int main(int argc, char* argv[])
{
  int nx = 8;
  Bx domain(Point::Zeros(), (nx-1)*Point::Ones());
  Bx grrdom = domain.grow(1);
  BoxData<double,1> phidom(domain);
  BoxData<double,1> phigrr(grrdom);
  phidom.setVal(-2.0);
  double phimax = phidom.max();
  double phimin = phidom.min();
  phigrr.setVal(-3.0);
  double grrmax = phigrr.max();
  double grrmin = phigrr.min();
  std::cout << "before copy phi (should be -2)  max = " << phimax << ", phi min = " << phimin << std::endl;
  std::cout << "before copy grr (should be -3)  max = " << grrmax << ", grr min = " << grrmin << std::endl;

//  phigrr.copyTo(phidom);
//
//  phimax = phidom.max();
//  phimin = phidom.min();
//  grrmax = phigrr.max();
//  grrmin = phigrr.min();
//
//  std::cout << "after copy phi (should be -3)  max = " << phimax << ", phi min = " << phimin << std::endl;
//  std::cout << "after copy grr (should be -3)  max = " << grrmax << ", grr min = " << grrmin << std::endl;
  

  return 0;
    
}
