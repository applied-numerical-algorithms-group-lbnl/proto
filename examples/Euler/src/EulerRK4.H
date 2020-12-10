#ifndef _EULERRK4_H_
#define _EULERRK4_H_

#include "EulerOp.H"
#include "Proto.H"
#include "Proto_WriteBoxData.H"
#include <iostream>
using namespace Proto;

class EulerDX;
class EulerState
{
public:
  EulerState(){};
  ~EulerState(){};
    EulerState(const Box& a_box,
               const double a_dx,
               const double a_gamma);
    double m_dx;
    double m_gamma;
  BoxData<double,NUMCOMPS> m_U;
  Box m_dbx0;
  int m_stage = 0;
  Box m_bdry[2*DIM];
  Point m_shift[2*DIM];
  Reduction<double> m_Rxn;
  void increment(const EulerDX& a_DX);
};

class EulerDX
{
public:
  EulerDX(){};
  ~EulerDX(){};
  BoxData<double,NUMCOMPS> m_DU;
  Box m_box;
  void increment(double& a_weight,const EulerDX& a_incr);
  void init(EulerState& a_State);
  void operator*=(const double& a_weight);
};

class EulerRK4Op
{
public:
    EulerRK4Op(){};
    ~EulerRK4Op(){};
    static int s_count;
    void operator()(
                  EulerDX& a_DX,
                  double a_time,
                  double a_dt,
                  EulerState& a_State);
};
#endif
