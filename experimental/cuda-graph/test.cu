
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
#include "cuda-graph.H"

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
  std::cout << " Dim = " << DIM << std::endl;
  unsigned int size1D = 32;
  unsigned int l      = 1;
  unsigned int iter   = 1;
 
  std::vector<Stencil<double>*> sten;
  std::vector<Box> bx;
  std::vector<Stencil<double>> ret;
  ret.resize(l);

  for(int it = 0; it < l ; it++)
  {
    ret[it] = ((double)(it))*Shift::Zeros();
    for(int i=0;i<2*DIM;i++) sten.push_back(&(ret[it]));
  }
  std::cout << " number of Stencils: " << sten.size() << std::endl; 

  Box base = Box::Cube(size1D);

  std::vector<Box> b;

  for(int it = 0; it < l ; it++)
  {
    Box tmp = base.grow(-it);
    b.push_back(tmp);
  }
  std::cout << " base Low : " << base.low() <<  " base High : " << base.high() << std::endl;

  BoxData<double,1> myBoxDatain(base);
  BoxData<double,1> myBoxDataout(base);
  BoxData<double,1> myBoxDataoutfused(base);

  
  forallInPlace(InitializeState, base, myBoxDatain, myBoxDataout);

  bx.resize(l*2*DIM);

  for(int it=0 ; it<l ; it++)
  {
    bx[it*4]   = b[it].faceBox(0,Side::Lo);
    bx[it*4+1] = b[it].edge(Point::Basis(1),1);
    bx[it*4+2] = b[it].faceBox(1,Side::Lo);
    bx[it*4+3] = b[it].edge(Point::Basis(0),1);
  }

  std::cout << " number of boxes: " << bx.size() << std::endl; 
  std::cout << " fused stencil " << std::endl;
  
  forallInPlace(InitializeState, base, myBoxDatain, myBoxDataoutfused);
 
  cudaEvent_t start, stop;
 
 {
 float milliseconds = 0;
  FusedStencil<double> myTest1;
  
  {
    auto& tmp = myTest1.getStencils();
    for(int it = 0 ; it<sten.size() ; it++)
      tmp.push_back(sten[it]);  
  }
  myTest1.copyInfo();

  std::cout << " agressive Merge " << std::endl;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start);
  for(int i=0; i<iter ; i++)
    myTest1.cudaApplyFused(myBoxDatain, myBoxDataout, bx.data(), true, 1.0); 
  cudaEventRecord(stop);

  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&milliseconds, start, stop);
  std::cout << " time: " << milliseconds << " ms" << std::endl;
  }
  //WriteBoxData(myBoxDataout,3);
 
  {
 float milliseconds = 0;
  FusedStencilGraph<double> myTest;
  {
    auto& tmp = myTest.getStencils();
    for(int it = 0 ; it<sten.size() ; it++)
      tmp.push_back(sten[it]);  
  }
  myTest.copyInfo();

  std::cout << " cuda graph " << std::endl;
 
  cudaEventRecord(start);
  for(int i=0; i<iter ; i++)
    myTest.cudaGraph(myBoxDatain, myBoxDataoutfused, bx.data(), true, 1.0); 
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&milliseconds, start, stop);
  std::cout << " time: " << milliseconds << " ms" << std::endl;
  }
 // WriteBoxData(myBoxDataoutfused,2);

  std::cout << " classic " << std::endl;

  {
 float milliseconds = 0;
  cudaEventRecord(start);
  for(int i=0; i<iter ; i++)
    for(int it = 0; it < sten.size() ; it++)
      sten[it]->apply(myBoxDatain, myBoxDataoutfused, bx[it], true, 1.0); 
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&milliseconds, start, stop);
  std::cout << " time: " << milliseconds << " ms" << std::endl;
  }
  return 0;
}


