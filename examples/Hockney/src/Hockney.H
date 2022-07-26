#ifndef _HOCKNEY_
#define _HOCKNEY_
#include "Proto.H"
#include "FFT1DW.H"
#include "FFT1D.H"
#include "FFTMD.H"
#include <array>
#include <memory>
using namespace std;
using namespace Proto;
class Hockney
{
public:
  Hockney();
  Hockney(const double& a_h,int a_M);
  void convolve(BoxData<double>& a_rhs);
  ~Hockney(){};
protected:
  double m_h;
  int m_M,m_N;
  FFTMD m_fftmd;
  bool m_isDefined;
  BoxData<complex<double> > m_kernelHat;
};
#endif
  
