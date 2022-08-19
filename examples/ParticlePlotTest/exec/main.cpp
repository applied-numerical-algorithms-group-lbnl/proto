#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <string>
#include <memory>
#include <array>
#include "Proto.H"
#include "writers.h"

using namespace std;

int main(int argc, char* argv[])
{
  unsigned int N;
  double size;
  cout << "input log_2(number of grid points), domain size" << endl;
  cin >> N ;
  cout << "N = " << N << endl;
  Box bx(Point::Zeros(),Point::Ones(N-1));
  double h = 1.0/N;
  vector<array<double,DIM+1> > particles;
  for (auto bxint = bx.begin();bxint.ok();++bxint)
    {
      array<double,DIM+1> temp;
      temp[0] = (*bxint)[0]*h;
      temp[1] = (*bxint)[1]*h;
      temp[2] = 1.0;
      particles.push_back(temp);
    }
  string testname="test";
  const char *const varnames[] = {"x_1", "x_2", "strength", "var3", "var4",
                                  "var5", "var6",
                                  "var7", "var8",
                                  "var9"};
  PWrite<DIM+1>(particles,varnames,testname,0);
}
