#pragma once

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

#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"

using namespace std;
using namespace Proto;

typedef Var<double,1> State;
typedef Var<double,1> V;


PROTO_KERNEL_START
unsigned int test_stencil_initialize_stateF(State& a_U)
{
    a_U(0) = 0;
    return 0;
}
PROTO_KERNEL_END(test_stencil_initialize_stateF, test_stencil_initialize_state)

PROTO_KERNEL_START
unsigned int test_stencil_valueF(State& a_U, double a_val)
{
    a_U(0) = a_val;
    return 0;
}
PROTO_KERNEL_END(test_stencil_valueF, test_stencil_value)

PROTO_KERNEL_START
unsigned int test_stencil_escalierF(Point p, State& a_U)
{
    a_U(0) = p[0]+p[1]*10;
#if DIM == 3
    a_U(0)+= p[2]*100;
#endif
    return 0;
}
PROTO_KERNEL_END(test_stencil_escalierF, test_stencil_escalier)

bool test_stencil_test_value(double *a_ptr, unsigned int a_size, double a_val)
{
  for(int i = 0 ; i < a_size ; i++)
    if(a_ptr[i] != a_val)
      return false;

  return true;
}


void test_stencil_print_mesh(double *a_ptr, unsigned int a_size)
{
#if DIM == 3
  for(int k = 0 ; k < a_size ; k++)
  {
#endif
  for(int j = 0 ; j < a_size ; j++)
  {
    for(int i = 0 ; i < a_size ; i++)
    {
#if DIM == 3
      unsigned int idx = i + (j + k * a_size) * a_size;
#else    
      unsigned int idx = i + j * a_size;
#endif
      std::cout << a_ptr[idx] << " ";
    }
    std::cout << std::endl;
  }
#if DIM == 3
    std::cout << std::endl;
    std::cout << std::endl;
  }
#endif
}

bool test_stencil_one_point_and_full()
{
  unsigned int size1D = 16;
  Box b = Box::Cube(size1D);
  BoxData<double,1> meshin(b);
  BoxData<double,1> meshout(b);

  double val = 2;
  double val_sten = 2;
  forallInPlace(test_stencil_initialize_state, b, meshout);
  forallInPlace(test_stencil_value, b, meshin, val);

  Stencil<double> sten = ((double)(val_sten))*Shift::Zeros();

  sten.apply(meshin, meshout, b, true, 1.0);


  double *d_ptr= meshout.data();
  unsigned int sizeBox = meshout.box().size();
  double *h_ptr = new double[sizeBox];
  unsigned int nBytes = sizeBox * sizeof(double);

  protoMemcpy(MEMTYPE_DEFAULT,h_ptr, d_ptr, nBytes, protoMemcpyDeviceToHost);
  protoDeviceSynchronize(MEMTYPE_DEFAULT);

  // result = 2*a_val 
  bool check = test_stencil_test_value(h_ptr,sizeBox, val*val_sten);
  if(!check) test_stencil_print_mesh(h_ptr,size1D);

  return check;
}

template<typename F>
bool test_stencil_solution(double *a_ptr, Box all, Box bx, F solution)
{
#if DIM == 3
  for(int k = bx.low()[2] ; k < bx.high()[2] ; k++)
#endif
  for(int j = bx.low()[1] ; j < bx.high()[1] ; j++)
    for(int i = bx.low()[0] ; i < bx.high()[0] ; i++)
#if DIM == 3
      if(!solution(a_ptr, all, i, j, k))
#else    
      if(!solution(a_ptr, all, i, j))
#endif
        return false;

  return true;
}


template<typename Func>
bool test_stencil_boxdata_box(Stencil<double> &a_sten, BoxData<double,1> &a_in, BoxData<double,1> &a_out, Box &a_bx, Func &solution)
{
  a_sten.apply(a_in, a_out, a_bx, true, 1.0);
  double *d_ptr= a_out.data();
  unsigned int sizeBox = a_out.box().size();
  double *h_ptr = new double[sizeBox];
  unsigned int nBytes = sizeBox * sizeof(double);

  protoMemcpy(MEMTYPE_DEFAULT,h_ptr, d_ptr, nBytes, protoMemcpyDeviceToHost);
  protoDeviceSynchronize(MEMTYPE_DEFAULT);
  bool check = test_stencil_solution(h_ptr,a_out.box(), a_bx, solution);
  if(!check) test_stencil_print_mesh(h_ptr,a_out.box().size(0));

  return check;
}

bool test_stencil_one_point_and_full_v2()
{
  unsigned int size1D = 16;
  Box b = Box::Cube(size1D);
  BoxData<double,1> meshin(b);
  BoxData<double,1> meshout(b);

  double val = 3;
  double val_sten = 2;
  forallInPlace(test_stencil_initialize_state, b, meshout);
  forallInPlace(test_stencil_value, b, meshin, val);

  Stencil<double> sten = ((double)(val_sten))*Shift::Zeros();

#if DIM == 3
  auto sol = [val,val_sten](double* ptr, Box all, int i, int j, int k)->bool 
	{ 
		if(ptr[i+(j+k*all.size(1))*all.size(0)] == val * val_sten) return true;
		else return false;
	};
#else
auto sol = [val,val_sten](double* ptr, Box all, int i, int j)->bool 
      {
	      if(ptr[i+j*all.size(0)] == val * val_sten) return true;
	      else return false;
      };
#endif
  return test_stencil_boxdata_box(sten, meshin, meshout, b, sol);
}


bool test_stencil_one_point_and_sub_box()
{
  unsigned int size1D = 16;
  Box b = Box::Cube(size1D);
  BoxData<double,1> meshin(b);
  BoxData<double,1> meshout(b);

  double val = 3;
  double val_sten = 2;
  forallInPlace(test_stencil_initialize_state, b, meshout);
  forallInPlace(test_stencil_value, b, meshin, val);

  Box subbox = b.grow(-2);

  Stencil<double> sten = ((double)(val_sten))*Shift::Zeros();

#if DIM == 3
  auto sol = [val,val_sten](double* ptr, Box all, int i, int j, int k)->bool 
	{ 
		if(ptr[i+(j+k*all.size(1))*all.size(0)] == val * val_sten) return true;
		else return false;
	};
#else
auto sol = [val,val_sten](double* ptr, Box all, int i, int j)->bool 
      {
	      if(ptr[i+j*all.size(0)] == val * val_sten) return true;
	      else return false;
      };
#endif
  return test_stencil_boxdata_box(sten, meshin, meshout, subbox, sol);
}


bool test_stencil_two_point_and_sub_box()
{
  unsigned int size1D = 16;
  Box b = Box::Cube(size1D);
  BoxData<double,1> meshin(b);
  BoxData<double,1> meshout(b);

  double val = 3;
  double val_sten = 2;
  forallInPlace(test_stencil_initialize_state, b, meshout);
  forallInPlace_p(test_stencil_escalier, b, meshin);

  Box subbox = b.grow(-3);
  Stencil<double> sten = ((double)(val_sten))*Shift::Zeros();
  sten += ((double)(1.5*val_sten))*Shift::Basis(0, 1);

#if DIM == 3
  auto sol = [val,val_sten](double* ptr, Box all, int i, int j, int k)->bool 
	{ 
		if(ptr[i+(j+k*all.size(1))*all.size(0)] == val_sten*((i+10*j+100*k)+1.5*(i+1+10*(j)+100*(k)))) return true;
		else return false;
	};
#else
auto sol = [val,val_sten](double* ptr, Box all, int i, int j)->bool 
      {
		if(ptr[i+j*all.size(0)] == val_sten*((i+10*j)+1.5*(i+1+10*j))) return true;
	      	else return false;
      };
#endif
  return test_stencil_boxdata_box(sten, meshin, meshout, subbox, sol);
}


bool test_stencil_laplacian_constant_and_sub_box()
{
  unsigned int size1D = 16;
  Box b = Box::Cube(size1D);
  BoxData<double,1> meshin(b);
  BoxData<double,1> meshout(b);

  double val = 3.213414;
  double val_sten = 2;
  forallInPlace(test_stencil_value, b, meshout,val);
  forallInPlace(test_stencil_value, b, meshin, val);

  Box subbox = b.grow(-2);

  Stencil<double> sten = Stencil<double>::Laplacian();

#if DIM == 3
  auto sol = [val,val_sten](double* ptr, Box all, int i, int j, int k)->bool 
	{ 
		if(ptr[i+(j+k*all.size(1))*all.size(0)] == 0) return true;
		else return false;
	};
#else
auto sol = [val,val_sten](double* ptr, Box all, int i, int j)->bool 
      {
	      if(ptr[i+j*all.size(0)] == 0) return true;
	      else return false;
      };
#endif
  return test_stencil_boxdata_box(sten, meshin, meshout, subbox, sol);
}


bool test_stencil_laplacian_escalier_and_sub_box()
{
  unsigned int size1D = 16;
  Box b = Box::Cube(size1D);
  BoxData<double,1> meshin(b);
  BoxData<double,1> meshout(b);

  double val = 3.213414;
  double val_sten = 2;
  //forallInPlace(test_stencil_value, b, meshout,val);
  forallInPlace_p(test_stencil_escalier, b, meshin);
  forallInPlace_p(test_stencil_escalier, b, meshout);

  Box subbox = b.grow(-2);

  Stencil<double> sten = Stencil<double>::Laplacian();

#if DIM == 3
  auto sol = [val,val_sten](double* ptr, Box all, int i, int j, int k)->bool 
	{ 
		if(ptr[i+(j+k*all.size(1))*all.size(0)] == 0) return true;
		else return false;
	};
#else
auto sol = [val,val_sten](double* ptr, Box all, int i, int j)->bool 
      {
	      if(ptr[i+j*all.size(0)] == 0) return true;
	      else return false;
      };
#endif
  return test_stencil_boxdata_box(sten, meshin, meshout, subbox, sol);
}
