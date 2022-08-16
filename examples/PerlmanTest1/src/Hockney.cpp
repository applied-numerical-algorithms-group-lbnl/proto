#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include "PowerItoI.H"
#include "Hockney.H"

using namespace std;
using namespace Proto;
Hockney::Hockney(const double& a_h,int a_M)
{ 
  shared_ptr<FFT1DW> p_fftw1d = shared_ptr<FFT1DW>(new FFT1DW(a_M+1));
  shared_ptr<FFT1D> p_fft = dynamic_pointer_cast<FFT1D>(p_fftw1d);
  m_fftmd.define(p_fft);
  m_h = a_h;
  m_M = a_M;
  m_N = Power(2,a_M);
  Point high = Point::Ones(m_N-1);
  Point low = Point::Ones(-m_N);
  Box ddomain(low,high);
  m_kernelHat =BoxData<complex<double> >(ddomain);
   forallInPlace_p(
                   [ ] PROTO_LAMBDA
                   (Point& a_pt,Var<std::complex<double >,1>& a_kernel,double& a_h) 
                   {
                     if (a_pt != Point::Zeros())
                       {
                         double distsq = a_h*a_h*(a_pt[0]*a_pt[0] + a_pt[1]*a_pt[1]);
                         double realKernel = (1.0/M_PI)*log(distsq)*a_h*a_h;
                         a_kernel(0) = complex<double>(realKernel,0.);
                       }
                     else
                       {
                         a_kernel(0) = complex<double>(0.,0.);
                       }
                   }
                   ,m_kernelHat,m_h);
   m_fftmd.forwardCCcen(m_kernelHat);
};
void Hockney::convolve(BoxData<double>& a_rhs)
{
  Box rhsDomain = a_rhs.box();
  Point low = rhsDomain.low();
  Point high = rhsDomain.high();
 
  assert(low == Point::Zeros());
  assert(high == Point::Ones()*(m_N-1));

  low = high*(-1) + Point::Ones(-1);
  Box ddomain(low,high);
  BoxData<complex<double> > rhsDouble(ddomain);
  complex<double> cxzero(0.,0.);
  rhsDouble.setVal(cxzero);
  double scale = 1./pow(m_N*4.,DIM);
   forallInPlace(
                 [ ] PROTO_LAMBDA
                 (Var<std::complex<double >,1>& a_rhsdouble,
                  Var<double,1>& a_rhs)
                 {
                   a_rhsdouble(0) = complex<double>(a_rhs(0),0.);
                 }
                 ,rhsDouble,a_rhs);
  BoxData<double > realOut(ddomain);
  m_fftmd.forwardCCcen(rhsDouble);
  forallInPlace_p(
                 [ ] PROTO_LAMBDA
                 (Point& a_p,
                  Var<std::complex<double >,1>& a_rhsdouble,
                  Var<std::complex<double >,1>& a_kernelhat)
                 {
                   a_rhsdouble(0) *= a_kernelhat(0);
                 }
                 ,rhsDouble,m_kernelHat);
  m_fftmd.inverseCCcen(rhsDouble);
  a_rhs.setVal(0.);
  forallInPlace(
                [scale] PROTO_LAMBDA
                (Var<double,1>& a_rhs,
                 Var<std::complex<double>,1>& a_rhsdouble, double& a_h)
                {
                  a_rhs(0) = a_rhsdouble(0).real()*scale;
                }
                ,a_rhs,rhsDouble,m_h);
}


