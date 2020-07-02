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

int main(int argc, char* argv[])
{
  int Ncells=64;
  double vel=1.0;
  AdvectionState state(1.0,Ncells,vel,"");
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
