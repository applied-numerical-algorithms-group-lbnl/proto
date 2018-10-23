#pragma once

#include <cmath>
#include <utility>
//#include <dvec.h>
#include <cstdio>
#include <iostream>

// example using the Intel Vector class.  just not available most places.
#ifdef __AVX512BW__
struct Real: public F64vec8
{
  Real() {;}
  Real(const double& val):F64vec8(val,val,val,val,val,val,val,val){;}
  Real(F64vec8&& val):F64vec8(std::move(val)){;}
};
#define VLENGTH 8

#elif defined __AVX__
#include <immintrin.h>

typedef double Real __attribute__ (( vector_size(8*4) ));

__m256i mask[4]={
  {0,0,0,0},
  {~0,0,0,0},
  {~0,~0,0,0},
  {~0,~0,~0,0}};

#define VLENGTH 4
constexpr Real one{1,1,1,1};
constexpr Real zero{0,0,0,0};
// done as multiply to get around bug in Intel C++ compiler FIXME
inline void set(Real& rr, double val){rr=one*val;}
inline Real sqrt(const Real& s)
{
  return _mm256_sqrt_pd(s);
}
inline Real maskLoad(const Real& src, int count)
{
  return _mm256_maskload_pd ((double const *) &src, mask[count]);
}
inline void maskStore(Real& dest, const Real& src, int count)
{
  _mm256_maskstore_pd ((double *)&dest, mask[count], src);
}

#else
#error This code needs to have some VLENGTH defined
 #include <stophere>
#endif


std::ostream& operator<<(std::ostream& os, const Real& r)
{
  const double* d = (const double*)&r;
  for(int i=0; i<VLENGTH; ++i) os<<d[i]<<" ";
  return os;
}

template<unsigned char D> inline void pincr(Real** ptr)
{
  ++ptr[D-1];
  pincr<D-1>(ptr);
}
template<> inline void pincr<1>(Real** ptr){ ++(ptr[0]);}

template <unsigned int C>
struct Var
{
  Real* m_ptr[C];
  inline Real& operator()(int i){ return *(m_ptr[i]);}
  inline const Real& operator()(int i) const { return *(m_ptr[i]);}
  inline Var& operator++(){ pincr<C>(m_ptr); return *this;}
  const double& operator()(int comp, int i) const {return ((const double*)m_ptr[comp])[i];}
};

template <unsigned int C>
std::ostream& operator<<(std::ostream& os, const Var<C>& v)
{
  for(int i=0; i<C; ++i) os<<*(v.m_ptr[i]);
  return os;
}
