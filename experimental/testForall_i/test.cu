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
unsigned int InitoneF(State& a_U)
{
    a_U(0) = 1;
    return 0;
}
PROTO_KERNEL_END(InitoneF, Initone)


PROTO_KERNEL_START
unsigned int InittwoF(State& a_U)
{
    a_U(0) = 2;
    return 0;
}
PROTO_KERNEL_END(InittwoF, Inittwo)

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

bool checkAnswer(double *ptr, unsigned int size1D)
{
  //edge = 1
  for(int i = 0; i<size1D ; i++)
	for(int j = 0 ; j<size1D ; j++)
		if( (i==0 || i == size1D-1) && (j==0 || j==size1D-1) )
			if(ptr[i+j*size1D]!=1)
				return false;
  //inside = 2
  for(int i = 1; i<size1D-1 ; i++)
	for(int j = 1 ; j<size1D-1 ; j++)
		if(ptr[i+j*size1D]!=2)
			return false;

  return true;
}

int main()
{
  cudaSetDevice(1);
  std::cout << " This code works only on GPU " << std::endl;
  std::cout << " Dim = 2" << std::endl;
  unsigned int size1D = 16;
  unsigned int size2D= size1D*size1D;

  Box b = Box::Cube(size1D);
  Box bminus= b.grow(-1);   
  Box bminus2= b.grow(-2);   

  BoxData<double,1> myBoxDatain(b);

  std::cout << " fill from " << b.low() << " to " << b.high() << " with 1 "<<std::endl;
  forallInPlace(Initone, b , myBoxDatain);
  std::cout << " fill from " << bminus.low() << " to " << bminus.high() << " with 2 "<<std::endl;
  forallInPlace(Inittwo, bminus , myBoxDatain);


  double *h_ptr = new double[size2D];
  double *d_ptr = myBoxDatain.dataPtr();
  unsigned int sizeBox = myBoxDatain.box().size();
  assert(size2D == sizeBox);
  unsigned int nBytes = sizeBox * sizeof(double);

  cudaMemcpy(h_ptr, d_ptr, nBytes, cudaMemcpyDeviceToHost);

  cudaDeviceSynchronize();

  bool check = checkAnswer(h_ptr, size1D);

  if(check)
	std::cout << " The result is correct " << std::endl;
  else 
	std::cout << " The result is wrong " << std::endl;

#ifdef ZASA
  WriteBoxData(myBoxDatain, 1);
#endif
  return 0;
}


