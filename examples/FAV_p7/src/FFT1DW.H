#ifndef _FFT1DW_H_
#define _FFT1DW_H_
#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
#include "PowerItoI.H"
#include "FFT1D.H"
#include "fftw3.h"

using namespace std;

class FFT1DW:public FFT1D
{
public:
  FFT1DW(unsigned int a_M);
  virtual ~FFT1DW();
  virtual void forwardFFTCC(vector<complex<double> > & a_fHat, 
                          const vector<complex<double> >& f) const;
  virtual void inverseFFTCC(vector<complex<double> > & a_f, 
                          const vector<complex<double> > & a_fHat) const;
protected:
  mutable vector<complex<double> > m_in;
  vector<complex<double> > m_out;
  fftw_plan m_forward, m_inverse;
};
#endif
