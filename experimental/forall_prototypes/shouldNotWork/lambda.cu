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


//compile with 
//nvcc -g -G -DTHRUST_DEBUG -DPROTO_MEM_CHECK --expt-extended-lambda -DPROTO_CUDA=1 -std=c++11 -DDIM=2 lambda.cu

using namespace std;
using namespace Proto;


PROTO_KERNEL_START void xeqylambdaF(Point& p, Var<int, 1>& err, 
                                    Var<double, DIM>& xv,
                                    Var<double, DIM>& yv)
{  
  err(0) = 0;
  err(0) = 0;
  for (int ii = 0; ii < DIM; ii++)
  {
    if(xv(ii) != yv(ii))
    {
      err(0) = 1;
    }
  }
}
PROTO_KERNEL_END(xeqylambdaF,xeqylambda)
/****************/
void
multigridSolve()
{
  int nx = 16;
  Point lo = Point::Zeros();
  Point hi = Point::Ones(nx - 1);
  Bx B(lo, hi);
  BoxData<double, DIM> X(B);
  BoxData<double, DIM> Y(B);
  X.setVal(4.);
  Y.setVal(7.);

  //first try without the lambda function
  BoxData<int> errf = forall_p<int>(xeqylambda, B, X, Y);
  
  cout << "after standard forall err (should be 1) max   =  "<< errf.max() << ", min = "<< errf.min() << endl;

  //now try with a lambda
  BoxData<int> errg = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                      Var<double, DIM> xv,
                                                      Var<double, DIM> yv) 
                                    {  
                                      err(0) = 0;
                                      for (int ii = 0; ii < DIM; ii++)
                                      {
                                        if(xv(ii) != yv(ii))
                                        {
                                          err(0) = 1;
                                        }
                                      }
                                    }, B, X, Y);


  cout << "after lamda forall err (should be 1) max   =  "<< errg.max() << ", min = "<< errg.min() << endl;
/**/
}
int main(int argc, char* argv[])
{

  multigridSolve();

}  
