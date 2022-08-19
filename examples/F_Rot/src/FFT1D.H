#ifndef _FFT1D_H_
#define _FFT1D_H_
#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
#include "PowerItoI.H"
using namespace std;
class FFT1D
{
public:
  // Interface class for complex-to-complex power-of-two FFT on the unit interval. 
  FFT1D(){};
  // Constructor. argument a_M specifies number of points is N= 2^{a_M}
  FFT1D(unsigned int a_M){m_M = a_M; m_N = Power(2,m_M);}
  virtual ~FFT1D() { }
  // Forward FFT: a_fHat[k] = \sum_j=0^{N-1} a_f[j] z^{j k}, z = e^{-2 \pi \iota /m_N}
  virtual void forwardFFTCC(vector<complex<double> > & a_fHat, 
                            const vector<complex<double> >& f) const = 0;
  // inverse FFT: a_f[j] = \sum_{k=0}^{N-1} a_fHat[k] z^{j k}, z = e^{2 \pi \iota /m_N}
  virtual void inverseFFTCC(vector<complex<double> > & a_f, 
                          const vector<complex<double> > & a_fHat) const = 0;
  // Access functions.
  const unsigned int& getN(){return m_N;} 
  const unsigned int& getM(){return m_M;}
protected:
  unsigned int m_M, m_N;
};
#endif 
