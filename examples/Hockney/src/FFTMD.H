#ifndef _FFTMD_H_
#define _FFTMD_H_
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <memory>
#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_BoxData.H"
#include "Proto_WriteBoxData.H"
#include "FFT1D.H"
using namespace std;
using namespace Proto;
class FFTMD
{
public:
  FFTMD(){};
  FFTMD(shared_ptr<FFT1D > a_fft1dPtr);
  void define(shared_ptr<FFT1D > a_fft1dPtr);
  void forwardCC(BoxData<complex<double> > & a_f) const;
  void inverseCC(BoxData<complex<double> > & a_fHat) const;
  void forwardCCcen(BoxData<complex<double> > & a_f) const;
  void inverseCCcen(BoxData<complex<double> > & a_fHat) const;
  const int& getN() const;
  const int& getM() const;
private:
  int m_N;
  int m_M;
  shared_ptr<FFT1D> m_fft1dPtr;
};
#endif
