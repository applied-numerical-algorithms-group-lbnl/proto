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
#include <iostream>

#include "Proto.H"
#include "Proto_Timer.H"

using namespace std;
using namespace Proto;

typedef Var<double,1> State;
typedef Var<double,1> V;


PROTO_KERNEL_START
unsigned int test_forall_initF(State& a_U, double a_val)
{
    a_U(0) = a_val;
    return 0;
}
PROTO_KERNEL_END(test_forall_initF, test_forall_init)

PROTO_KERNEL_START
unsigned int test_forall_init_pV2F(Point p, State& a_U)
{
    a_U(0) = p[0]+p[1]*10;
#if DIM == 3
    a_U(0)+= p[2]*100;
#endif
    return 0;
}
PROTO_KERNEL_END(test_forall_init_pV2F, test_forall_init_pV2)


PROTO_KERNEL_START
unsigned int test_forall_init_pF(Point p, State& a_U, double a_val)
{
    a_U(0) = a_val;
    return 0;
}
PROTO_KERNEL_END(test_forall_init_pF, test_forall_init_p)

PROTO_KERNEL_START
unsigned int test_forall_init_iV2F(int p[DIM], State& a_U)
{
    a_U(0) = p[0]+p[1]*10;
#if DIM == 3
    a_U(0)+= p[2]*100;
#endif
    return 0;
}
PROTO_KERNEL_END(test_forall_init_iV2F, test_forall_init_iV2)

PROTO_KERNEL_START
unsigned int test_forall_init_iF(int p[DIM], State& a_U, double a_val)
{
    a_U(0) = a_val;
    return 0;
}
PROTO_KERNEL_END(test_forall_init_iF, test_forall_init_i)


void test_forall_print(double *ptr, unsigned int size1D)
{
  //edge = 1
  for(size_t i = 0; i<size1D ; i++)
  {
	for(size_t j = 0 ; j<size1D ; j++)
		std::cout << ptr[i+j*size1D] << " ";
        std::cout << std::endl;
  }				
  std::cout << std::endl;
}

bool test_forall_check_answer(double *ptr, unsigned int size1D)
{
#if DIM == 3
  //edge = 1
  for(size_t i = 0; i<size1D ; i++)
	for(size_t j = 0 ; j<size1D ; j++)
		for(size_t k = 0 ; k<size1D ; k++)
		if( (i==0 || i == size1D-1) && (j==0 || j==size1D-1) && (k==1||k==size1D-1))
			if(ptr[i+j*size1D+k*size1D*size1D]!=1)
			{
				std::cout << " error [" << i << "," << j << "," << k << "] =" << ptr[i+j*size1D+k*size1D*size1D] << " != 1 " <<std::endl;
				return false;
			}
  //inside = 2
  for(size_t i = 1; i<size1D-1 ; i++)
	for(size_t j = 1 ; j<size1D-1 ; j++)
		for(size_t k = 1 ; k<size1D-1 ; k++)
		if(ptr[i+j*size1D+k*size1D*size1D]!=2)
		{
			std::cout << " error [" << i << "," << j << "," << k << "] =" << ptr[i+j*size1D+k*size1D*size1D] << " != 2 " <<std::endl;
			return false;
		}
#else
  //edge = 1
  for(size_t i = 0; i<size1D ; i++)
	for(size_t j = 0 ; j<size1D ; j++)
		if( (i==0 || i == size1D-1) && (j==0 || j==size1D-1) )
			if(ptr[i+j*size1D]!=1)
			{
				std::cout << " error [" << i << "," << j << "] =" << ptr[i+j*size1D] << " != 1 " <<std::endl;
				return false;
			}
  //inside = 2
  for(size_t i = 1; i<size1D-1 ; i++)
	for(size_t j = 1 ; j<size1D-1 ; j++)
		if(ptr[i+j*size1D]!=2)
		{
			std::cout << " error [" << i << "," << j << "] =" << ptr[i+j*size1D] << " != 2 " <<std::endl;
			return false;
		}
#endif
  return true;
}


bool test_forall_check_answer_p(double *ptr, unsigned int size1D)
{
#if DIM == 3
  //edge = 1
  for(size_t i = 0; i<size1D ; i++)
	for(size_t j = 0 ; j<size1D ; j++)
		for(size_t k = 0 ; k<size1D ; k++)
		if( (i==0 || i == size1D-1) && (j==0 || j==size1D-1) && (k==0 || k==size1D-1))
			if(ptr[i+j*size1D+k*size1D*size1D]!=1)
			{
				std::cout << " error [" << i << "," << j << "," << k << "] =" << ptr[i+j*size1D+k*size1D*size1D] << " != 1 " <<std::endl;
				return false;
			}
  //inside = 2
  for(size_t i = 1; i<size1D-1 ; i++)
	for(size_t j = 1 ; j<size1D-1 ; j++)
		for(size_t k = 1 ; k<size1D-1 ; k++)
		if(ptr[i+j*size1D+k*size1D*size1D]!=i+j*10+k*100)
		{
			std::cout << " error [" << i << "," << j << ","<< k << "] =" << ptr[i+j*size1D+k*size1D*size1D] << " != 2 " <<std::endl;
			return false;
		}
#else
  //edge = 1
  for(size_t i = 0; i<size1D ; i++)
	for(size_t j = 0 ; j<size1D ; j++)
		if( (i==0 || i == size1D-1) && (j==0 || j==size1D-1) )
			if(ptr[i+j*size1D]!=1)
			{
				std::cout << " error [" << i << "," << j << "] =" << ptr[i+j*size1D] << " != 1 " <<std::endl;
				return false;
			}
  //inside = 2
  for(size_t i = 1; i<size1D-1 ; i++)
	for(size_t j = 1 ; j<size1D-1 ; j++)
		if(ptr[i+j*size1D]!=i+j*10)
		{
			std::cout << " error [" << i << "," << j << "] =" << ptr[i+j*size1D] << " != 2 " <<std::endl;
			return false;
		}
#endif
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
  forallInPlace(test_forall_init, b, myBoxDatain, a);
  a=2;
  forallInPlace(test_forall_init, bminus, myBoxDatain, a);

  double *d_ptr=myBoxDatain.data();
  unsigned int sizeBox = myBoxDatain.box().size();
#if DIM == 3
  assert(size2D*size1D == sizeBox);
  double *h_ptr = new double[size2D*size1D];
#else
  assert(size2D == sizeBox);
  double *h_ptr = new double[size2D];
#endif
  unsigned int nBytes = sizeBox * sizeof(double);

  protoMemcpy(MEMTYPE_DEFAULT,h_ptr, d_ptr, nBytes, protoMemcpyDeviceToHost);
  protoDeviceSynchronize(MEMTYPE_DEFAULT);

  bool check = test_forall_check_answer(h_ptr, size1D);
  if(!check) test_forall_print(h_ptr,size1D);

  assert(check);
  return check;
}

bool run_test_forall_p()
{
  unsigned int size1D = 16;
  unsigned int size2D= size1D*size1D;
  double *d_ptr;

  Box b = Box::Cube(size1D);
  Box bminus= b.grow(-1);   

  BoxData<double,1> myboxdataforall_p(b);
  unsigned int sizeBox = myboxdataforall_p.box().size();
#if DIM == 3
  assert(size2D*size1D == sizeBox);
  double *h_ptr = new double[size2D*size1D];
#else
  assert(size2D == sizeBox);
  double *h_ptr = new double[size2D];
#endif
  unsigned int nBytes = sizeBox * sizeof(double);

  double val=1;
  forallInPlace_p(test_forall_init_p, b, myboxdataforall_p, val);
  forallInPlace_p(test_forall_init_pV2, bminus, myboxdataforall_p);

  d_ptr = myboxdataforall_p.data();

  protoMemcpy(MEMTYPE_DEFAULT,h_ptr, d_ptr, nBytes, protoMemcpyDeviceToHost);
  protoDeviceSynchronize(MEMTYPE_DEFAULT);

  bool check = test_forall_check_answer_p(h_ptr, size1D);
  if(!check) test_forall_print(h_ptr,size1D);
 
  assert(check); 
  return check;
}


bool run_test_forall_i()
{
  unsigned int size1D = 16;
  unsigned int size2D= size1D*size1D;

  Box b = Box::Cube(size1D);
  Box bminus= b.grow(-1);   

  BoxData<double,1> myboxdataforall_i(b);
  unsigned int sizeBox = myboxdataforall_i.box().size();
#if DIM == 3
  assert(size2D*size1D == sizeBox);
  double *h_ptr = new double[size2D*size1D];
#else
  assert(size2D == sizeBox);
  double *h_ptr = new double[size2D];
#endif
  unsigned int nBytes = sizeBox * sizeof(double);

  double val = 1;
  forallInPlace_i(test_forall_init_i, b, myboxdataforall_i, val);
  forallInPlace_i(test_forall_init_iV2, bminus, myboxdataforall_i);

  double * d_ptr = myboxdataforall_i.data();
  protoMemcpy(MEMTYPE_DEFAULT,h_ptr, d_ptr, nBytes, protoMemcpyDeviceToHost);
  protoDeviceSynchronize(MEMTYPE_DEFAULT);

  bool check = test_forall_check_answer_p(h_ptr, size1D);
  if(!check) test_forall_print(h_ptr,size1D);

  assert(check);
  return check;
}


