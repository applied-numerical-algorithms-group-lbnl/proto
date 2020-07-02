#ifndef _ADVECTION_RK4_H_
#define _ADVECTION_RK4_H_

#include "Proto.H"
#include <iostream>
using namespace Proto;

class AdvectionDX;
class AdvectionState
{
public:
  AdvectionState(const double& domain_length,
                 const int& n_cells,
                 const double& vel,
                 const std::string& init_case);

  void increment(const AdvectionDX& incr);

  double m_L;
  int m_N;
  double m_dx;
  double m_vel;
  BoxData<double,1> m_phi;  

private:
  //void InitializePhi(const std::string& init_case);
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
