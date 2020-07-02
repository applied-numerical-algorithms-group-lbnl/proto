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
void initialPhi_p_temp(Point& a_p,
                       Var<double>& phi)
{
  double R=std::abs(a_p[0]-0.5);
  double R0=0.15;
  double pi_div_2=1.57079632679;
  if(R<=R0)
    phi(0)=pow(cos(pi_div_2*(R/R0)),8);
  else
    phi(0)=0.0;
}
PROTO_KERNEL_END(initialPhi_p_temp,initialPhi_p)

int main(int argc, char* argv[])
{
  int Ncells=64;
  double vel=1.0;
  AdvectionState state(1.0,Ncells,vel);
  forallInPlace_p(initialPhi_p,state.m_phi);
  RK4<AdvectionState,AdvectionOp,AdvectionDX> rk4_timestepper;
  double time=0.0;
  double dt=0.5*vel/Ncells;
  int maxStep=10000;
  double tStop=1;
  for(int k=0; k<maxStep && time<=tStop; k++)
    {
      rk4_timestepper.advance(time,dt,state);
      time+=dt;
    }


  return 0;
}
