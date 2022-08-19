#include <cstdio>
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "Proto.H"
#include "FFT1D.H"
#include "FFTMD.H"
using namespace std;
using namespace Proto;

FFTMD::FFTMD(shared_ptr<FFT1D> a_fft1dPtr)
{
 
  m_fft1dPtr = a_fft1dPtr;
  m_M = m_fft1dPtr->getM();
  m_N = m_fft1dPtr->getN();
}; 
void FFTMD::define(shared_ptr<FFT1D> a_fft1dPtr)
{
 
  m_fft1dPtr = a_fft1dPtr;
  m_M = m_fft1dPtr->getM();
  m_N = m_fft1dPtr->getN();
};
void FFTMD::forwardCC(BoxData<complex<double> > & a_f) const
{ 
  vector<complex<double> > f1d(m_N);
  vector<complex<double> > fHat1d(m_N);
  PR_TIMERS("fftmdForward");
  PR_TIMER("fft1d",t1);
  int low[DIM],high[DIM];
  for (int dir = 0;dir < DIM; dir++)
    {
      for (int dir2 = 0;dir2 < DIM;dir2++)
        {
          low[dir2]= 0;
          high[dir2] = m_N-1;
        }
      high[dir]=0;
      Point lo(low),hi(high);
      cout << DIM << " " << m_N << " " << hi[0] << " " << hi[1] << " " << hi[2] << endl;
      Box base(lo,hi);
      int offset = Power(m_N+1,dir);
      Point edir = Point::Basis(dir);
      for (auto it=base.begin();it.ok();++it)
        {
          Point pt = *it;
          for (int l = 0 ; l < m_N;l++)
            {
              f1d[l] = a_f(pt + edir*l);
            }
          PR_START(t1);
          m_fft1dPtr->forwardFFTCC(fHat1d,f1d);
          PR_STOP(t1);
          for (int l = 0 ; l < m_N;l++)
            {
              a_f(pt + edir*l) = fHat1d[l];
            }
          //cout << "In FFTMD: " << ztot 
          //   << " , " << a_f[pt] << " , dir = " << dir << endl;
          //cout << counter << endl;
        }
    }        
};
void FFTMD::inverseCC(BoxData<complex<double> > & a_fHat) const
{int low[DIM],high[DIM];
  PR_TIMERS("fftmdInverse");
  PR_TIMER("fft1drev",t3);
  vector<complex<double> > f1d(m_N);
  vector<complex<double> > fHat1d(m_N);
  for (int dir = 0;dir < DIM ; dir++)
    {
      for (int dir2 = 0;dir2 < DIM;dir2++)
        {
          low[dir2]= 0;
          high[dir2] = m_N-1;
        }
      high[dir]=0;
      int offset = Power(m_N+1,dir);
      Box base(low,high);
      Point edir = Point::Basis(dir);
      for (auto it=base.begin();it.ok();++it)
        {
          Point pt = *it;
          //complex<double>* foo = &a_fHat.getPointer()[base.getIndex(pt)];
          for (int l = 0 ; l < m_N;l++)
            {
              fHat1d[l] =a_fHat(pt + edir*l);
            }
          PR_START(t3);
          m_fft1dPtr->inverseFFTCC(f1d,fHat1d);
          PR_STOP(t3);
          for (int l = 0 ; l < m_N;l++)
            {
              a_fHat(pt + edir*l) = f1d[l];
            }
        }
    }
};
const int& FFTMD::getN() const
{
  return m_N;
  
};
const int& FFTMD::getM() const
{
  return m_M;
};
void FFTMD::forwardCCcen(BoxData<complex<double> > & a_f) const
{
  int low[DIM],high[DIM];
  vector<complex<double> > f1d(m_N);
  vector<complex<double> > fHat1d(m_N);
  PR_TIMERS("fftmdForward");
  PR_TIMER("fft1d",t1);
  for (int dir = 0;dir < DIM; dir++)
    {
      for (int dir2 = 0;dir2 < DIM;dir2++)
        {
          low[dir2]= -m_N/2;
          high[dir2] = m_N/2-1;
        }
      high[dir]=0;
      low[dir]=0;
      Box base(low,high);
      Point edir = Point::Basis(dir);
      for (auto it = base.begin();it.ok();++it)
        {
          Point pt = *it;
          for (int l = 0 ; l < m_N/2;l++)
            {
              f1d[l] = a_f(pt + edir*l);
              f1d[m_N - l-1] = a_f(pt-edir*(l+1)); 
            }
          PR_START(t1);
          m_fft1dPtr->forwardFFTCC(fHat1d,f1d);
          PR_STOP(t1);
          
          for (int l = 0 ; l < m_N/2;l++)
            {
              a_f(pt + edir*l) = fHat1d[l];
              a_f(pt-edir*(l+1)) = fHat1d[m_N-l-1];
            }
        }
      
    }     
};
void FFTMD::inverseCCcen(BoxData<complex<double> > & a_fHat) const
{int low[DIM],high[DIM];
  PR_TIMERS("fftmdInverse");
  PR_TIMER("fft1drev",t3);
  vector<complex<double> > f1d(m_N);
  vector<complex<double> > fHat1d(m_N);
  for (int dir = 0;dir < DIM ; dir++)
    {
      for (int dir2 = 0;dir2 < DIM;dir2++)
        {
          low[dir2]= -m_N/2;
          high[dir2] = m_N/2-1;
        }
      high[dir]=0;
      low[dir]=0;
      int offset = Power(m_N+1,dir);
      Box base(low,high);
      Point edir = Point::Basis(dir);
      for (auto it = base.begin();it.ok();++it)
        {
          Point pt = *it;
          for (int l = 0 ; l < m_N/2;l++)
            {
              fHat1d[l] = a_fHat(pt + edir*l);
              fHat1d[m_N - l-1] = a_fHat(pt-edir*(l+1)); 
            }
          PR_START(t3);
          m_fft1dPtr->inverseFFTCC(f1d,fHat1d);
          PR_STOP(t3);
          for (int l = 0 ; l < m_N/2;l++)
            {
              a_fHat(pt + edir*l) = f1d[l];
              a_fHat(pt-edir*(l+1)) = f1d[m_N - l-1];
            }
        }
    }
};
