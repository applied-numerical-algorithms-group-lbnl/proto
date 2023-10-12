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
  //phi(0)=sin(2*M_PI*(x-vel*time));
  double R=std::abs(x-vel*time-0.5);
  double R0=0.15;
  double pi_div_2=1.57079632679;
  if(R<=R0)
    phi(0)=pow(cos(pi_div_2*(R/R0)),8);
  else
    phi(0)=0.0;
}
PROTO_KERNEL_END(evaluatePhiFace_p_temp,evaluatePhiFace_p)

PROTO_KERNEL_START
void evaluateExactDiv_p_temp(Point& a_p,
                             Var<double>& phi,
                             const double& time,
                             const double& vel,
                             const double& dx)
{
  //phi(0)=1.0;
  double xp=(a_p[0]+0.5)*dx;
  double xm=(a_p[0]-0.5)*dx;
  phi(0)=vel*(sin(2*M_PI*(xm-vel*time))-sin(2*M_PI*(xp-vel*time)))/dx;
  //double R=std::abs(x-vel*time-0.5);
  //double R0=0.15;
  //double pi_div_2=1.57079632679;
  //if(R<=R0)
  //  phi(0)=pow(cos(pi_div_2*(R/R0)),8);
  //else
  //  phi(0)=0.0;
}
PROTO_KERNEL_END(evaluateExactDiv_p_temp,evaluateExactDiv_p)

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

class AdvectionRK4
{
public:
  AdvectionRK4();
  ~AdvectionRK4();
  void advance(double time, double& dt, AdvectionState& state);
  void advance(BoxData<double>& k_new, BoxData<double>& k_prev, double time, double& dt, AdvectionState& state);
  void operator()(AdvectionDX& k, double time, double& dt, AdvectionState& state); 
};

class AdvectionOp
{
public:
    AdvectionOp();
    ~AdvectionOp();
    static void RK4Step(BoxData<double>& dF, BoxData<double>& curr_phi, double time, double dt, double dx, double vel);
    static void ComputeFlux(BoxData<double>& flux, BoxData<double>& phi, double time, double& dt, double vel);
    static void ComputeDivergenceFlux(BoxData<double>& dF, BoxData<double>& flux, double dx, double vel);
};

#endif
