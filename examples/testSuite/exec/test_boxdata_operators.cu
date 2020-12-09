#pragma once

bool test_boxdata_operators_check_value(double* a_ptr, double a, unsigned int a_size)
{
  for(int id = 0 ; id < a_size ; id++)
    if(a_ptr[id] != a)
    {
      std::cout << " a_ptr[" << id << "]= " << a_ptr[id] << " != " << a << std::endl;
      return false;
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
  protoMemcpy(host,myBoxData.data(),myBoxData.size()*sizeof(double),protoMemcpyDeviceToHost);
  bool check = test_boxdata_operators_check_value(host,a_init,myBoxData.size());
  return check;
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
  }

  double res = computeSolution(a_init,a_val,op);  
  double *host = new double[myBoxData.size()];
  protoMemcpy(host,myBoxData.data(),myBoxData.size()*sizeof(double),protoMemcpyDeviceToHost);
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
