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
#include "Proto_Timer.H"

using namespace std;
using namespace Proto;

typedef Var<double,1> State;
typedef Var<double,1> V;


PROTO_KERNEL_START
unsigned int InitF(State& a_U, double a_val)
{
    a_U(0) = a_val;
    return 0;
}
PROTO_KERNEL_END(InitF, Init)

PROTO_KERNEL_START
unsigned int Init_pV2F(Point p, State& a_U)
{
    a_U(0) = p[0]+p[1]*10;
    return 0;
}
PROTO_KERNEL_END(Init_pV2F, Init_pV2)


PROTO_KERNEL_START
unsigned int Init_pF(Point p, State& a_U, double a_val)
{
    a_U(0) = a_val;
    return 0;
}
PROTO_KERNEL_END(Init_pF, Init_p)

PROTO_KERNEL_START
unsigned int Init_iV2F(int p[DIM], State& a_U)
{
    a_U(0) = p[0]+p[1]*10;
    return 0;
}
PROTO_KERNEL_END(Init_iV2F, Init_iV2)

PROTO_KERNEL_START
unsigned int Init_iF(int p[DIM], State& a_U, double a_val)
{
    a_U(0) = a_val;
    return 0;
}
PROTO_KERNEL_END(Init_iF, Init_i)


void print(double *ptr, unsigned int size1D)
{
  //edge = 1
  for(int i = 0; i<size1D ; i++)
  {
	for(int j = 0 ; j<size1D ; j++)
		std::cout << ptr[i+j*size1D] << " ";
        std::cout << std::endl;
  }				
  std::cout << std::endl;
}

bool checkAnswer(double *ptr, unsigned int size1D)
{
  //edge = 1
  for(int i = 0; i<size1D ; i++)
	for(int j = 0 ; j<size1D ; j++)
		if( (i==0 || i == size1D-1) && (j==0 || j==size1D-1) )
			if(ptr[i+j*size1D]!=1)
			{
				std::cout << " error [" << i << "," << j << "] =" << ptr[i+j*size1D] << " != 1 " <<std::endl;
				return false;
			}
  //inside = 2
  for(int i = 1; i<size1D-1 ; i++)
	for(int j = 1 ; j<size1D-1 ; j++)
		if(ptr[i+j*size1D]!=2)
		{
			std::cout << " error [" << i << "," << j << "] =" << ptr[i+j*size1D] << " != 2 " <<std::endl;
			return false;
		}
  return true;
}


bool checkAnswer_p(double *ptr, unsigned int size1D)
{
  //edge = 1
  for(int i = 0; i<size1D ; i++)
	for(int j = 0 ; j<size1D ; j++)
		if( (i==0 || i == size1D-1) && (j==0 || j==size1D-1) )
			if(ptr[i+j*size1D]!=1)
			{
				std::cout << " error [" << i << "," << j << "] =" << ptr[i+j*size1D] << " != 1 " <<std::endl;
				return false;
			}
  //inside = 2
  for(int i = 1; i<size1D-1 ; i++)
	for(int j = 1 ; j<size1D-1 ; j++)
		if(ptr[i+j*size1D]!=i+j*10)
		{
			std::cout << " error [" << i << "," << j << "] =" << ptr[i+j*size1D] << " != 2 " <<std::endl;
			return false;
		}
  return true;
}

bool run_test_forall()
{
  unsigned int size1D = 16;
  unsigned int size2D= size1D*size1D;

  Box b = Box::Cube(size1D);
  Box bminus= b.grow(-1);   

  BoxData<double,1> myBoxDatain(b);
  double a = 1;
  forallInPlace(Init, b, myBoxDatain, a);
  a=2;
  forallInPlace(Init, bminus, myBoxDatain, a);

  double *h_ptr = new double[size2D];
  double *d_ptr=myBoxDatain.dataPtr();
  unsigned int sizeBox = myBoxDatain.box().size();
  assert(size2D == sizeBox);
  unsigned int nBytes = sizeBox * sizeof(double);

  protoMemcpy(h_ptr, d_ptr, nBytes, protoMemcpyDeviceToHost);
  protoDeviceSynchronize();

  bool check = checkAnswer(h_ptr, size1D);
//  print(h_ptr,size1D);

  assert(check);
  return check;
}

bool run_test_forall_p()
{
  unsigned int size1D = 16;
  unsigned int size2D= size1D*size1D;
  double *h_ptr = new double[size2D];
  double *d_ptr;

  Box b = Box::Cube(size1D);
  Box bminus= b.grow(-1);   

  BoxData<double,1> myboxdataforall_p(b);
  unsigned int sizeBox = myboxdataforall_p.box().size();
  assert(size2D == sizeBox);
  unsigned int nBytes = sizeBox * sizeof(double);

  double val=1;
  forallInPlace_p(Init_p, b, myboxdataforall_p, val);
  forallInPlace_p(Init_pV2, bminus, myboxdataforall_p);

  d_ptr = myboxdataforall_p.dataPtr();
  protoMemcpy(h_ptr, d_ptr, nBytes, protoMemcpyDeviceToHost);

  protoDeviceSynchronize();

  bool check = checkAnswer_p(h_ptr, size1D);
 
  assert(check); 
  return check;
}


bool run_test_forall_i()
{
  unsigned int size1D = 16;
  unsigned int size2D= size1D*size1D;
  double *h_ptr = new double[size2D];

  Box b = Box::Cube(size1D);
  Box bminus= b.grow(-1);   

  BoxData<double,1> myboxdataforall_i(b);
  unsigned int sizeBox = myboxdataforall_i.box().size();
  assert(size2D == sizeBox);
  unsigned int nBytes = sizeBox * sizeof(double);

  double val = 1;
  forallInPlace_i(Init_i, b, myboxdataforall_i, val);
  forallInPlace_i(Init_iV2, bminus, myboxdataforall_i);

  double * d_ptr = myboxdataforall_i.dataPtr();
  protoMemcpy(h_ptr, d_ptr, nBytes, protoMemcpyDeviceToHost);

  protoDeviceSynchronize();

  bool check = checkAnswer_p(h_ptr, size1D);

  assert(check);
  return check;
}


