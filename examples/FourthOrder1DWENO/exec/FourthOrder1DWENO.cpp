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
  AdvectionState state(1.0,Ncells,1.0,"");
  RK4<AdvectionState,AdvectionOp,AdvectionDX> rk4_timestepper;
  double time=0.0;
  double dt=0.1;
  rk4_timestepper.advance(time,dt,state);


  return 0;
}
