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
#include <iostream>

#include "Proto.H"
#include "Proto_Timer.H"

using namespace Proto;

bool test_boxdata_operators_check_value(double* a_ptr, double a, unsigned int a_size)
{
  for(int id = 0 ; id < a_size ; id++)
    if(a_ptr[id] != a)
    {
//      std::cout << " a_ptr[" << id << "]= " << a_ptr[id] << " != " << a << std::endl;
      return false;
    }

  return true;
}

bool test_boxdata_operators_check_value_box(double* a_ptr, double a, Proto::Box b, Proto::Box all)
{
  for(int id = 0 ; id < b.size() ; id++)
  {
    unsigned int idx = all.index(b[id]);
    if( a_ptr[idx] != a)
    {
//      std::cout << " a_ptr[" << idx << "]= " << a_ptr[idx] << " != " << a << " (id ----> " << id << ")" << std::endl;
      return false;
    }
  }
  return true;
}

bool test_boxdata_operators_set_value(double a_init)
{
  unsigned size1D=16;
  Box b = Box::Cube(size1D);
  BoxData<double,1> myBoxData(b);
  myBoxData.setVal(a_init);
  
  double *host = new double[myBoxData.size()];
  protoMemcpy(MEMTYPE_DEFAULT,host,myBoxData.data(),myBoxData.size()*sizeof(double),protoMemcpyDeviceToHost);
  bool check = test_boxdata_operators_check_value(host,a_init,myBoxData.size());
  return check;
}


bool test_boxdata_operators_copy(double a_val, double a_varBoxSize)
{
  unsigned size1D=16;
  Box b = Box::Cube(size1D);
  Box bis = b.grow(a_varBoxSize);
  BoxData<double,1> to(bis);
  BoxData<double,1> from(b);
  to.setVal(a_varBoxSize);
  from.setVal(a_val);

  from.copyTo(to);  

  double *host = new double[to.size()];

  protoMemcpy(MEMTYPE_DEFAULT,host,to.data(),to.size()*sizeof(double),protoMemcpyDeviceToHost);
  bool check    = test_boxdata_operators_check_value(host,a_val,to.size()) ;

  bool checkbis = test_boxdata_operators_check_value_box(host,a_val,bis,to.box());

  assert(check);
  assert(checkbis);

  return check && checkbis;
}


bool test_boxdata_operators_copy_to_box(double a_val, double a_varBoxSize)
{
  unsigned size1D=16;
  Box b = Box::Cube(size1D);
  Box bis = b.grow(a_varBoxSize);
  BoxData<double,1> to(b);
  BoxData<double,1> from(b);
  to.setVal(a_varBoxSize);
  from.setVal(a_val);

  from.copyTo(to,bis,Proto::Point(0));  

  double *host = new double[to.size()];
  protoMemcpy(MEMTYPE_DEFAULT,host,to.data(),to.size()*sizeof(double),protoMemcpyDeviceToHost);
  bool check    = !test_boxdata_operators_check_value(host,a_val,to.size());
  bool checkbis = test_boxdata_operators_check_value_box(host,a_val,bis,to.box());

  if(!checkbis) 
  {
	  for(int i = 0 ; i < size1D*size1D; i++)
		  std::cout << host[i] << " ";
	  std::cout << std::endl;
  }

  assert(check);
  assert(checkbis);
  return check && checkbis;
}

bool test_boxdata_operators_copy_full()
{
  return test_boxdata_operators_copy(56.233452,0);
}

bool test_boxdata_operators_copy_to_smaller_box_data()
{
  return test_boxdata_operators_copy(56.233452,-1);
}

bool test_boxdata_operators_copy_to_smaller_box()
{
  return test_boxdata_operators_copy_to_box(56.233452,-1);
}

bool test_boxdata_operators_set_value_zero()
{
  return test_boxdata_operators_set_value(0);
}

bool test_boxdata_operators_set_value_minus_three()
{
  return test_boxdata_operators_set_value(-3);
}

bool test_boxdata_operators_set_value_two()
{
  return test_boxdata_operators_set_value(2);
}

double computeSolution(double a, double b,  Proto::BoxDataOp op)
{
  if(op == Proto::BoxDataOp::Add)
    return a+b;
  if(op == Proto::BoxDataOp::Subtract)
    return a-b;
  if(op == Proto::BoxDataOp::Multiply)
    return a*b;
  if(op == Proto::BoxDataOp::Divide)
  {
    assert(b!=0);
    return a/b;
  }
  std::abort();
}

bool test_boxdata_operators_op(double a_init, double a_val, Proto::BoxDataOp op)
{
  unsigned size1D=16;
  assert(size1D > 0);
  Box b = Box::Cube(size1D);
  BoxData<double,1> myBoxData(b);
  BoxData<double,1> myBoxData2(b);
  myBoxData.setVal(a_init);
  myBoxData2.setVal(a_val);

  switch(op)
  {
    case BoxDataOp::Add:
      myBoxData += myBoxData2;
      break;
    case BoxDataOp::Subtract:
      myBoxData -= myBoxData2;
      break;
    case BoxDataOp::Multiply:
      myBoxData *= myBoxData2;
      break;
    case BoxDataOp::Divide:
      myBoxData /= myBoxData2;
      break;
    default:
      PR_error("invalid operator given to test_boxdata_operators_op");
  }

  double res = computeSolution(a_init,a_val,op);  
  double *host = new double[myBoxData.size()];
  protoMemcpy(MEMTYPE_DEFAULT,host,myBoxData.data(),myBoxData.size()*sizeof(double),protoMemcpyDeviceToHost);
  bool check = test_boxdata_operators_check_value(host,res,myBoxData.size());
  return check;
}

bool test_boxdata_operators_op_scalar(double a_init, double a_val, Proto::BoxDataOp op)
{
  unsigned size1D=16;
  assert(size1D > 0);
  Box b = Box::Cube(size1D);
  BoxData<double,1> myBoxData(b);
  myBoxData.setVal(a_init);

  switch(op)
  {
    case BoxDataOp::Add:
      myBoxData += a_val;
      break;
    case BoxDataOp::Subtract:
      myBoxData -= a_val;
      break;
    case BoxDataOp::Multiply:
      myBoxData *= a_val;
      break;
    case BoxDataOp::Divide:
      myBoxData /= a_val;
      break;
    default:
      PR_error("invalid operator given to test_boxdata_operators_op_scalar");
  }

  double res = computeSolution(a_init,a_val,op);  
  double *host = new double[myBoxData.size()];
  protoMemcpy(MEMTYPE_DEFAULT,host,myBoxData.data(),myBoxData.size()*sizeof(double),protoMemcpyDeviceToHost);
  bool check = test_boxdata_operators_check_value(host,res,myBoxData.size());
  return check;
}

bool test_boxdata_operators_add_one_and_two()
{
  return test_boxdata_operators_op(1,2,Proto::BoxDataOp::Add);
}

bool test_boxdata_operators_subtract_one_and_two()
{
  return test_boxdata_operators_op(1,2,Proto::BoxDataOp::Subtract);
}

bool test_boxdata_operators_multiply_one_and_two()
{
  return test_boxdata_operators_op(1,2,Proto::BoxDataOp::Multiply);
}

bool test_boxdata_operators_divide_one_and_two()
{
  return test_boxdata_operators_op(1,2,Proto::BoxDataOp::Divide);
}

bool test_boxdata_operators_add_neg()
{
  return test_boxdata_operators_op(-1.342,-65.33,Proto::BoxDataOp::Add);
}

bool test_boxdata_operators_subtract_neg()
{
  return test_boxdata_operators_op(-1.342,-65.33,Proto::BoxDataOp::Subtract);
}

bool test_boxdata_operators_multiply_neg()
{
  return test_boxdata_operators_op(-1.342,-65.33,Proto::BoxDataOp::Multiply);
}

bool test_boxdata_operators_divide_neg()
{
  return test_boxdata_operators_op(-1.342,-65.33,Proto::BoxDataOp::Divide);
}

bool test_boxdata_operators_scalar_add_one_and_two()
{
  return test_boxdata_operators_op_scalar(1,2,Proto::BoxDataOp::Add);
}

bool test_boxdata_operators_scalar_subtract_one_and_two()
{
  return test_boxdata_operators_op_scalar(1,2,Proto::BoxDataOp::Subtract);
}

bool test_boxdata_operators_scalar_multiply_one_and_two()
{
  return test_boxdata_operators_op_scalar(1,2,Proto::BoxDataOp::Multiply);
}

bool test_boxdata_operators_scalar_divide_one_and_two()
{
  return test_boxdata_operators_op_scalar(1,2,Proto::BoxDataOp::Divide);
}

bool test_boxdata_operators_scalar_add_neg()
{
  return test_boxdata_operators_op_scalar(-1.342,-65.33,Proto::BoxDataOp::Add);
}

bool test_boxdata_operators_scalar_subtract_neg()
{
  return test_boxdata_operators_op_scalar(-1.342,-65.33,Proto::BoxDataOp::Subtract);
}

bool test_boxdata_operators_scalar_multiply_neg()
{
  return test_boxdata_operators_op_scalar(-1.342,-65.33,Proto::BoxDataOp::Multiply);
}

bool test_boxdata_operators_scalar_divide_neg()
{
  return test_boxdata_operators_op_scalar(-1.342,-65.33,Proto::BoxDataOp::Divide);
}
