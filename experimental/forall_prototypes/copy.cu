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
  Box domain(Point::Zeros(), (nx-1)*Point::Ones());
  Box grrdom = domain.grow(1);
  BoxData<double,1> psidom(domain);
  BoxData<double,1> phidom(domain);
  BoxData<double,1> phigrr(grrdom);
  phidom.setVal(-2.0);
  psidom.setVal(-4.0);
  double phimax = phidom.max();
  double phimin = phidom.min();
  phigrr.setVal(-3.0);
  double grrmax = phigrr.max();
  double grrmin = phigrr.min();
  double psimax = psidom.max();
  double psimin = psidom.min();
  std::cout << "before copy phi (should be -2)  max = " << phimax << ", phi min = " << phimin << std::endl;
  std::cout << "before copy psi (should be -4)  max = " << psimax << ", psi min = " << psimin << std::endl;
  std::cout << "before copy grr (should be -3)  max = " << grrmax << ", grr min = " << grrmin << std::endl;

  psidom.copyTo(phidom);
  phimax = phidom.max();
  phimin = phidom.min();


  std::cout << "after same size copy phi (should be -4)  max = " << phimax << ", phi min = " << phimin << std::endl;
  phigrr.copyTo(phidom);


  phimax = phidom.max();
  phimin = phidom.min();

  std::cout << "after diff size copy phi (should be -3)  max = " << phimax << ", phi min = " << phimin << std::endl;

  return 0;
    
}
