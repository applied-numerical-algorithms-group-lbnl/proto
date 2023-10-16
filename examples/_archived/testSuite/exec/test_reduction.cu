#pragma once
#include<Proto.H>

void test_reduction_linear(double* ptr, double a_coef, double a_pt, unsigned int a_size)
{
  for(int i = 0 ; i < a_size ; i++)
  {
    ptr[i] = i * a_coef + a_pt;
  }
}

bool test_reduction_min_linear_init_val(double a_val, double a_pt)
{
  unsigned int size = 32;
  double * data = new double[size];
 
  test_reduction_linear(data, a_val, a_pt, size);

  double * device;

  protoMallocGPU(device,sizeof(double)*size); 
  protoMemcpyGPU(device, data, sizeof(double)*size, protoMemcpyHostToDevice);

  Reduction<double,Operation::Min> red;
  red.reset();

  red.reduce(device, size);
 
  double result = red.fetch();
 
  bool check = false;

  double good;
  if(a_val < 0) 
    good = ((size-1) * a_val + a_pt);
  else
    good = a_pt;

  check = result == good;

  if(!check)
    std::cout << " correct result: " << good << " | what we get " << result << std::endl; 
 
  return check;
}

bool test_reduction_max_linear_init_val(double a_val, double a_pt)
{
  unsigned int size = 32;
  double * data = new double[size];
 
  test_reduction_linear(data, a_val, a_pt, size);

  double * device;
  protoMallocGPU(device,sizeof(double)*size); 
  protoMemcpyGPU(device, data, sizeof(double)*size, protoMemcpyHostToDevice);

  Reduction<double,Operation::Max> red;
  red.reset();

  red.reduce(device, size);
 
  double result = red.fetch();
 
  bool check = false;

  if(a_val > 0) 
    check = result == ((size-1) * a_val + a_pt);
  else
    check = result == a_pt;

  if(!check)
    std::cout << result << std::endl; 
 
  return check;
}


bool test_reduction_abs_linear_init_val(double a_val, double a_pt)
{
  unsigned int size = 32;
  double * data = new double[size];
 
  test_reduction_linear(data, a_val, a_pt, size);

  double * device;
  protoMallocGPU(device,sizeof(double)*size); 
  protoMemcpyGPU(device, data, sizeof(double)*size, protoMemcpyHostToDevice);

  Reduction<double,Operation::Abs> red;
  red.reset();

  red.reduce(device, size);
 
  double result = red.fetch();
 
  bool check = false;

  double good = std::max( std::abs(a_pt),std::abs((size-1) * a_val+a_pt) );

  check = result == good;

  if(!check)
    std::cout << result << std::endl; 
 
  return check;
}

bool test_reduction_min_linear_init_1()
{
  return test_reduction_min_linear_init_val(1,-1.01);
}

bool test_reduction_min_linear_init_minus_2()
{
  return test_reduction_min_linear_init_val(-2,1.01);
}

bool test_reduction_max_linear_init_1()
{
  return test_reduction_max_linear_init_val(1,-1.01);
}

bool test_reduction_max_linear_init_minus_2()
{
  return test_reduction_max_linear_init_val(-2,1.01);
}


bool test_reduction_abs_linear_init_1()
{
  return test_reduction_abs_linear_init_val(1,-12.01);
}

bool test_reduction_abs_linear_init_minus_2()
{
  return test_reduction_abs_linear_init_val(-2,12.01);
}

bool test_reduction_reset_min()
{
  unsigned int size = 32;
  double * data = new double[size];
  Reduction<double,Operation::Min> red;
  red.reset();
 
  double result = red.fetch();
 
  bool check = result == std::numeric_limits<double>::max();
  return check;
}


bool test_reduction_reset_max()
{
  unsigned int size = 32;
  double * data = new double[size];
  Reduction<double,Operation::Max> red;
  red.reset();
 
  double result = red.fetch();
 
  bool check = result == std::numeric_limits<double>::min();
  return check;
}

bool test_reduction_reset_abs()
{
  unsigned int size = 32;
  double * data = new double[size];
  Reduction<double,Operation::Abs> red;
  red.reset();
 
  double result = red.fetch();
 
  bool check = result == double(0);
  return check;
}
