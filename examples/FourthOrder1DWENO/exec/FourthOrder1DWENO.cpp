#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "Proto.H"
#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"
#include "AdvectionRK4.H"

using namespace std;
using namespace Proto;

PROTO_KERNEL_START
void evaluatePhi_p_temp(Point& a_p,
                        Var<double>& phi,
                        const double& time,
                        const double& vel,
                        const double& dx)
{
  double x=a_p[0]*dx;
  double R=std::abs(x-vel*time-0.5);
  double R0=0.15;
  double pi_div_2=1.57079632679;
  if(R<=R0)
    phi(0)=pow(cos(pi_div_2*(R/R0)),8);
  else
    phi(0)=0.0;
}
PROTO_KERNEL_END(evaluatePhi_p_temp,evaluatePhi_p)

int main(int argc, char* argv[])
{
  int Ncells=64;
  double vel=1.0;
  double time=0.0;
  AdvectionState state(1.0,Ncells,vel);
  forallInPlace_p(evaluatePhi_p,state.m_phi,time,state.m_vel,state.m_dx);
  std::cout << "Max initial phi: " << state.m_phi.absMax() << std::endl;
  RK4<AdvectionState,AdvectionOp,AdvectionDX> rk4_timestepper;

  double dt=0.5*vel/Ncells;
  int maxStep=10000;
  double tStop=1;
  for(int k=0; k<maxStep && time<tStop; k++)
    {
      rk4_timestepper.advance(time,dt,state);
      time+=dt;
      std::cout << "Time: , max initial phi: " << time << " , " << state.m_phi.absMax() << std::endl;
    }
  BoxData<double> exact_solution(state.m_phi.box());
  exact_solution.setVal(0.0);
  forallInPlace_p(evaluatePhi_p,exact_solution,time,state.m_vel,state.m_dx);
  BoxData<double> error(state.m_phi.box());
  state.m_phi.copyTo(error);
  error-=exact_solution;
  std::cout << "Abs max computed solution: " << state.m_phi.absMax() << std::endl;
  std::cout << "Max error: " << error.absMax() << std::endl;


  return 0;
}
