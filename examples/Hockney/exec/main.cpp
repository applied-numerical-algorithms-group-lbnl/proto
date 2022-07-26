#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <string>
#include <memory>
#include <array>
#include "Proto.H"
#include "Hockney.H"
#include "FFTMD.H"

using namespace std;

int main(int argc, char* argv[])
{
  unsigned int M;
  unsigned int N;
  double size;
  cout << "input log_2(number of grid points), domain size" << endl;
  cin >> M >> size;
  cout << "M = " << M << ", size = " << size << endl;
  N = Power(2,M);
  double h = size/N;
  Hockney hockney(h,M);
  Box bx(Point::Zeros(),Point::Ones(N-1));
  BoxData<double> charge(bx);
  forallInPlace_p(
                  [N] PROTO_LAMBDA
                  (Point& a_pt, Var<double,1>& a_charge, double& a_h, double& a_sz) 
                   {
                     double x = a_pt[0]*a_h - a_sz/2;
                     double y = a_pt[1]*a_h - a_sz/2;
                     double radius = sqrt(x*x + y*y); 
                     if (radius < .25)
                     //if (a_pt == Point::Ones(N/2))
                       {
                         a_charge(0) = pow(cos(M_PI*radius/.25/2),6);
                         //a_charge(0) = 1.0/a_h/a_h;
                       }                     
                     else
                       {
                         a_charge(0) = 0.;
                       }
                   }
                  ,charge,h,size);

  BoxData<double> sv(bx);
  charge.copyTo(sv);
  HDF5Handler h5;
  std::string name = "charge"+to_string(N)+"_"+to_string(size);
  h5.writePatch(h,charge,name);
  hockney.convolve(charge);
  BoxData<double> LPhi = Stencil<double>::Laplacian()(charge,bx.grow(-2),1.0/(h*h));
  LPhi -= sv;
  /*forallInPlace_p(
                  [ ] PROTO_LAMBDA
                  (Point& a_pt, Var<double,1>& a_field, double& a_h, double& a_sz) 
                   {
                     double x = a_pt[0]*a_h - a_sz/2;
                     double y = a_pt[1]*a_h - a_sz/2;
                     double radius = sqrt(x*x + y*y); 
                     if (radius!=0)
                       {
                         a_field(0) /= (.5/M_PI)*log(radius);
                       }
                     else
                       {
                         a_field(0) = 0;                    
                       }
                  }
                  ,charge,h,size);*/
  std::string nameField = "field"+to_string(N)+"_"+to_string(size);
  h5.writePatch(h,charge,nameField);
  std::string nameLPhi = "LPhi"+to_string(N)+"_"+to_string(size);
  h5.writePatch(h,LPhi,nameLPhi);
}
