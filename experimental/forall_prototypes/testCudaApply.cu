
#define DIM 3
#define PROTO_CUDA 1
#include "../../include/Proto.H"
#include <cstdlib>
#include <cstdio>
#include <functional>
#include <iostream>

using namespace Proto;

PROTO_KERNEL_START void initParabolaT(Point& p, Var<double>& data)
{
  data(0) = 0;
  for(int idir = 0; idir < DIM; idir ++)
  {
    data(0) += p[idir]*p[idir];
  }
}
PROTO_KERNEL_END(initParabolaT, initParabola);

PROTO_KERNEL_START void initBogosityT(Point& p, Var<double>& data)
{
  data(0) = 123456789e10;
}
PROTO_KERNEL_END(initBogosityT, initBogosity);

int main(int argc, char** argv) 
{
  constexpr int nx = 16;
  Bx domain = Bx::Cube(nx);
  Bx grrdom = domain.grow(1);
  BoxData<double, 1> phi(grrdom);
  BoxData<double, 1> lap(domain);
  cudaForall_p(initParabola,  grrdom,  phi);
  cudaForall_p(initBogosity,  domain,  lap);
//  phi.cudaSyncHost();
//  lap.cudaSyncHost();

//  printf("phi init :\n");
//  phi.printData();
//  printf("lap init :\n");
//  lap.printData();
  


  Stencil<double> lapsten = Stencil<double>::Laplacian(2);
//  Stencil<double> lapsten(Shift(Point::Zeros()), 1.0);
  
  printf("stencil given by\n");
  lapsten.print();
//  lapsten.apply(phi, lap, domain, true, 1.0);
//  printf("phi after host apply :\n");
//  phi.printData();
//  printf("lap after host apply :\n");
//  lap.printData();
//  phi.cudaSyncDevice();
//  lap.cudaSyncDevice();

  lapsten.cudaApply(phi, lap, domain, true, 1.0);
  cudaDeviceSynchronize();
  phi.cudaSyncHost();
  lap.cudaSyncHost();
  
//  printf("phi after device apply :\n");
//  phi.printData();
//  printf("lap after device apply :\n");
//  lap.printData();
  using std::cout;
  using std::endl;
  double tol = 1.0e-12;
  for(BxIterator bit = domain.begin(); bit != domain.end(); ++bit)
  {
    double lapval = lap(*bit, 0);
    double corval = 2*DIM;
    if(std::abs(lapval - corval) > tol)
    {
      cout << "Test Failed: laplacian has wrong value of "<< lapval << " at point " << *bit << endl;
      return -3;
    }
  }


  lapsten.apply(phi, lap, domain, true, 1.0);
  cudaDeviceSynchronize();
  phi.cudaSyncHost();
  lap.cudaSyncHost();
  
//  printf("phi after device apply :\n");
//  phi.printData();
//  printf("lap after device apply :\n");
//  lap.printData();


  for(BxIterator bit = domain.begin(); bit != domain.end(); ++bit)
  {
    double lapval = lap(*bit, 0);
    double corval = 2*DIM;
    if(std::abs(lapval - corval) > tol)
    {
      cout << "Test Failed 2: laplacian has wrong value of "<< lapval << " at point " << *bit << endl;
      return -5;
    }
  }
  cout << "Test Passed!" << endl;
  return 0;
}
