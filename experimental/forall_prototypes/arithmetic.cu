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
using std::cout;
using std::endl;
/****************/

void 
printStuff(BoxData<double, DIM> & a_data,
           std::string          a_location,
           double               a_expectedVal)
{
  cout << a_location << "(should be "  << a_expectedVal << "):" << endl; 
  for(int idir = 0; idir < DIM; idir++)
  {
    double phimax = a_data.max(idir);
    double phimin = a_data.min(idir);
    cout  << "max[" << idir << "]=" << phimax << ", min[" << idir<< "]=" << phimin << std::endl;
  }
}

int main(int argc, char* argv[])
{
  int nx = 8;
  Box domain(Point::Zeros(), (nx-1)*Point::Ones());

  BoxData<double,DIM> data1(domain);
  data1.setVal(2.);
  printStuff(data1, string("after setVal    "), 2.);

  data1 *= 4.;
  printStuff(data1, string("after scalarMult"), 8.);
             
  data1 /= 2.;
  printStuff(data1, string("after scalarDiv "), 4.);

  data1 += 4.;
  printStuff(data1, string("after scalarAdd "), 8.);

  data1 -= 2.;
  printStuff(data1, string("after scalarSub "), 6.);
             

  BoxData<double,DIM> data2(domain);
  data2.setVal(2.);

  data1 /= data2;                        
  printStuff(data1, string("after bdDiv     "), 3.);

  data1 *= data2;
  printStuff(data1, string("after bdMult    "), 6.);
                                         
  data1 += data2;                        
  printStuff(data1, string("after bdAdd     "), 8.);
                                         
                                         
  data1 -= data2;                        
  printStuff(data1, string("after bdSub     "), 6.);

  return 0;
    
}
