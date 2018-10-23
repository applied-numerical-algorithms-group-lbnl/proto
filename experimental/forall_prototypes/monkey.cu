
#define DIM 2
#define PROTO_CUDA 1
#include "../../include/Proto.H"
#include <cstdlib>
#include <cstdio>
#include <functional>
#include <iostream>
using namespace Proto;
using std::cout;
using std::endl;

void func_p(Point a_p)
{
  cout << a_p << endl;
}
int
int_ref(BoxData<int>& a_data, const Point& a_p){return a_data(a_p);}

int main(int argc, char** argv) 
{
  printf("hello\n");
  constexpr int nx = 4;
  Bx domain = Bx::Cube(nx);
  BoxData<int> data(domain);
  int val = 0;
  for(BxIterator bit = domain.begin(); bit != domain.end(); ++bit)
  {
    data(*bit) = val;
    ++val;
  }


  for(BxIterator bit = domain.begin(); bit != domain.end(); ++bit)
  {
    cout << int_ref(data, *bit) << endl;
  }

  printf("goodbye\n");
  return 0;
}
