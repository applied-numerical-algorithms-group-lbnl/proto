#ifndef _ADVECTION_RK4_H_
#define _ADVECTION_RK4_H_

#include "Proto.H"
#include <iostream>
using namespace Proto;

PROTO_KERNEL_START
void evaluatePhiFace_p_temp(Point& a_p,
                            Var<double>& phi,
                            const double& time,
                            const double& vel,
                            const double& dx)
{
  //phi(0)=1.0;
  double x=(a_p[0]-0.5)*dx;
  double R=std::abs(x-vel*time-0.5);
  double R0=0.15;
  double pi_div_2=1.57079632679;
  if(R<=R0)
    phi(0)=pow(cos(pi_div_2*(R/R0)),8);
  else
    phi(0)=0.0;
}
PROTO_KERNEL_END(evaluatePhiFace_p_temp,evaluatePhiFace_p)

class AdvectionDX;
class AdvectionState
{
public:
  AdvectionState(const double& domain_length,
                 const int& n_cells,
                 const double& vel);

  void increment(const AdvectionDX& incr);
  static void setBoundaryConditions(BoxData<double>& state_ext);

  double m_L;
  int m_N;
  double m_dx;
  double m_vel;
  BoxData<double,1> m_phi;
};

class AdvectionDX
{
public:
  AdvectionDX();
  ~AdvectionDX();

  void init(AdvectionState& state);
  void increment(const double& weight, const AdvectionDX& incr);
  void operator*=(const double& weight);

  BoxData<double> m_dF;
};

class AdvectionOp
{
public:
  AdvectionOp();
  ~AdvectionOp();
  void operator()(AdvectionDX& k, double time, double& dt, AdvectionState& state); 
};

#endif
