#pragma once

#include "PR_Real.h"
#include "forall.h"

double s_gamma = 1.4;
double s_invdx = 4.0;

#define DIM 3
#define NUMCOMPS 5



static inline void init(Var<NUMCOMPS>& a_v)
{ 
  a_v(0) = one*0.05;
  a_v(1) = one; 
  a_v(2) = one*1.5;
  a_v(3) = zero;
  a_v(4) = one*450.0;
}

static inline
void consToPrim(Var<NUMCOMPS>& a_W, 
                const Var<NUMCOMPS>& a_U)
{
    Real rho = a_U(0);
    Real v2 = zero;

    a_W(0) = rho;

    for (int i = 1; i <= DIM; i++)
    {
        Real v;
        v = a_U(i) / rho;

        a_W(i) = v;
        v2 += v*v;
    }
    Real E = (a_U(NUMCOMPS-1) -  rho * v2* 0.5) * (s_gamma - 1.0);
    a_W(NUMCOMPS-1) = E;
}

int dir=0; // cheating a variable into the function through global
Real v_gamma = one*s_gamma/(s_gamma-1.0);
static inline void getFlux(Var<NUMCOMPS>& a_F, const Var<NUMCOMPS>& a_W)
{
    Real F0 = a_W(dir+1)*a_W(0);
    Real W2 = zero;

    a_F(0) = F0;

    for (int d = 1; d <= DIM; d++)
    {
        Real Wd = a_W(d);

        a_F(d) = Wd*F0;
        W2 += Wd*Wd;
    }

    a_F(dir+1) += a_W(NUMCOMPS-1);
    a_F(NUMCOMPS-1) = v_gamma * a_W(dir+1) * a_W(NUMCOMPS-1) + F0 * W2 * 0.5;
}

Real vinvdx = one*(-s_invdx);
static inline void scaleInvDx(Var<NUMCOMPS>& v)
{
  for(int i=0; i<NUMCOMPS; ++i)
    {
      v(i)*=vinvdx;
    }
}

