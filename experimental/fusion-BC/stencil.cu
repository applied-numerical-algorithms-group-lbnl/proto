#include <iostream>
#include <Proto.H>

#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "Proto.H"
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"

using namespace std;
using namespace Proto;

typedef Var<double,1> State;
typedef Var<double,1> V;


PROTO_KERNEL_START
unsigned int InitializeStateF(State& a_U, State& a_U2)
{
    a_U(0) = 1;
    a_U2(0) = 0;
    return 0;
}
PROTO_KERNEL_END(InitializeStateF, InitializeState)

void WriteData( BoxData<double, 1>&a_state, int it)
{
    char basename[1024];
    sprintf(basename,"euler.%06d",it);

    const char* varnames[1];
    varnames[0] = "martin";
    double origin[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        origin[ii] = 0.0;
    }
    WriteBoxData(basename,a_state,varnames,origin,1);
};

int main()
{
  std::cout << " This code works only on GPU " << std::endl;
  std::cout << " Dim = 2" << std::endl;
  unsigned int size1D = 8;
  unsigned int size2D= size1D*size1D;

  Stencil<double> ret2 = ((double)(2))*Shift::Zeros();
  Stencil<double> ret3 = ((double)(3))*Shift::Zeros();
  Stencil<double> ret4 = ((double)(4))*Shift::Zeros();

  Box b = Box::Cube(size1D);
  Box bminus= b.grow(-1);   
  Box bminus2= b.grow(-2);   

  std::cout << " b Low : " << b.low() <<  " b High : " << b.high() << std::endl;
  std::cout << " bminus Low : " << bminus.low() <<  " bminus High : " << bminus.high() << std::endl;


  BoxData<double,1> myBoxDatain(b);
  BoxData<double,1> myBoxDataout(b);
  BoxData<double,1> myBoxDataoutfused(b);

  std::vector<Stencil<double>> sten;
  std::vector<Box> bx;
  
  for(int i=0;i<4;i++) sten.push_back(ret2);
  for(int i=0;i<4;i++) sten.push_back(ret3);

  forallInPlace(InitializeState, b, myBoxDatain, myBoxDataout);


  Box FxInf = b.faceBox(0,Side::Lo); 
  Box FxSup = b.edge(Point::Basis(1),1); 
  Box FyInf = b.faceBox(1,Side::Lo); 
  Box FySup = b.edge(Point::Basis(0),1); 

  bx.push_back(FxInf);
  bx.push_back(FxSup);
  bx.push_back(FyInf);
  bx.push_back(FySup);

  Box FxInfminus = bminus.faceBox(0,Side::Lo); 
  Box FxSupminus = bminus.edge(Point::Basis(1),1); 
  Box FyInfminus = bminus.faceBox(1,Side::Lo); 
  Box FySupminus = bminus.edge(Point::Basis(0),1); 

  bx.push_back(FxInfminus);
  bx.push_back(FxSupminus);
  bx.push_back(FyInfminus);
  bx.push_back(FySupminus);

  for( int it = 0 ; it < sten.size() ; it++)
    sten[it].apply(myBoxDatain, myBoxDataout, bx[it], true, 1.0);
  
  ret4.apply(myBoxDatain, myBoxDataout, bminus2, true, 1.0);

  WriteBoxData(myBoxDataout, 1);

  std::cout << " fused stencil " << std::endl;
  
  forallInPlace(InitializeState, b, myBoxDatain, myBoxDataoutfused);
  FusedStencil<double> myTest;
  
  auto& ptr = myTest.getStencils();

  for(int i=0;i<4;i++) ptr.push_back(&ret2);
  for(int i=0;i<4;i++) ptr.push_back(&ret3);
 
  myTest.copyInfo();

  myTest.cudaApplyFused(myBoxDatain, myBoxDataoutfused, bx.data(), true, 1.0); 
  ret4.apply(myBoxDatain, myBoxDataoutfused, bminus2, true, 1.0);
  WriteBoxData(myBoxDataoutfused,2);


  return 0;
}


