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

using namespace std;
using namespace Proto;

/**
 * @brief Define and fill phi_cell with cell-averages of phi=1.5*x^3+2*x^2+3 on the domain
 * [0,1].
 * @param phi_cell defined on the box [-Nghost,Ncell+Nghost-1] and filled with cell-averages of phi.
 * @param Ncells number of cells in the discretization of [0,1]
 * @param Nghost number of ghost cells
 * @returns defines and fills phi_cell
 */
void InitializePhiCell(BoxData<double>& phi_cell,
                       const int Ncells,
                       const int Nghost)
{
  Box box(Point({0}),Point({Ncells}));
  phi_cell.define(box.grow(Nghost));
  //TODO: set phi=1.5*x^3+2*x^2+3
  phi_cell.setVal(1.0);
}

void InitializeVelocity(BoxData<double>& vel,
                        const int Ncells)
{
  Box box(Point({0}),Point({Ncells+1}));
  vel.define(box);
  vel.setVal(1.0);
}

PROTO_KERNEL_START
void smoothnessFactors_temp(Var<double>& wl,
                            Var<double>& wr,
                            const Var<double>& cl1,
                            const Var<double>& cl2,
                            const Var<double>& cr1,
                            const Var<double>& cr2,
                            const double& eps)
{
  //std::cout << cl1(0) << ", " << cl2(0) << std::endl;
  //bl(0)=cl1(0)*cl2(0);
  //std::cout << bl(0) << std::endl;
  //br(0)=cr1(0)*cr2(0);
  double bl=(4.0/3.0)*cl1(0)*cl1(0)+(1.0/2.0)*cl1(0)*cl2(0)+(1.0/4.0)*cl2(0)*cl2(0);
  double br=(4.0/3.0)*cr1(0)*cr1(0)-(1.0/2.0)*cr1(0)*cr2(0)+(1.0/4.0)*cr2(0)*cr2(0);
  double temp1=(eps+bl)*(eps+bl)+(eps+br)*(eps+br);
  double temp2=(eps+bl)*(eps+bl)/temp1;
  double temp3=(eps+br)*(eps+br)/temp2;
  double al=temp2*(0.75+temp2*temp2-1.5*temp2);
  double ar=temp3*(0.75+temp3*temp3-1.5*temp3);
  wl(0)=al/(ar+al);
  wr(0)=ar/(ar+al);
}
PROTO_KERNEL_END(smoothnessFactors_temp,smoothnessFactors)

PROTO_KERNEL_START
void computePhiFaceAve_temp(Var<double>& phi_face,
                            const Var<double>& vel_face,
                            const Var<double>& wl,
                            const Var<double>& wr,
                            const Var<double>& fl,
                            const Var<double>& fr)
{
  double max_w=std::max(wl(0),wr(0));
  double min_w=std::min(wl(0),wr(0));
  if(vel_face(0)>0)
    phi_face(0)=max_w*fl(0)+min_w*fr(0);
  else
    phi_face(0)=max_w*fr(0)+min_w*fl(0);
}
PROTO_KERNEL_END(computePhiFaceAve_temp,computePhiFaceAve)
                            

/***/
int main(int argc, char* argv[])
{
  //cos^6(2*pi*x)
  //by-hand periodic boundary conditions
  //add time-stepping (RK4) --> use template RK4 pattern in Proto
  //retime with Milo's changes
  int Ncells=64;
  int Nghost=2;
  BoxData<double> phi_cell;
  InitializePhiCell(phi_cell,Ncells,Nghost);

  BoxData<double> vel;
  InitializeVelocity(vel,Ncells);

  Stencil<double> S_c1=1.0*Shift::Zeros()-2.0*Shift::Basis(0,-1)+1.0*Shift::Basis(0,-2);
  Stencil<double> S_c2=1.0*Shift::Zeros()-1.0*Shift::Basis(0,-2);
  BoxData<double> cl1=S_c1(phi_cell);
  BoxData<double> cl2=S_c2(phi_cell);
  BoxData<double> cr1=alias(cl1,Point::Ones(-1));
  BoxData<double> cr2=alias(cl2,Point::Ones(-1));
  std::cout << "Phi cell domain: " << phi_cell.box() << std::endl;
  std::cout << "cl1 domain: " << cl1.box() << std::endl;
  std::cout << "cl2 domain: " << cl2.box() << std::endl;
  std::cout << "cr1 domain: " << cr1.box() << std::endl;
  std::cout << "cr2 domain: " << cr2.box() << std::endl;

  BoxData<double> wl(phi_cell.box().grow(-Nghost));
  wl.setVal(0.0);
  BoxData<double> wr(phi_cell.box().grow(-Nghost));
  wr.setVal(0.0);
  double eps=1e-5;
  //cl1.setVal(1.0);
  //cl2.setVal(2.0);
  //cr1.setVal(1.5);
  //cr2.setVal(3.0);
  forallInPlace(smoothnessFactors,wl,wr,cl1,cl2,cr1,cr2,eps);
  //std::cout << b1(Point({0})) << b1(Point({3}))<<std::endl;

  Stencil<double> S_fl=(1.0/6.0)*(5.0*Shift::Basis(0,-1)+2.0*Shift::Zeros()-1.0*Shift::Basis(0,-2));
  Stencil<double> S_fr=(1.0/6.0)*(2.0*Shift::Basis(0,-1)+5.0*Shift::Zeros()-1.0*Shift::Basis(0,1));
  BoxData<double> fl=S_fl(phi_cell);
  BoxData<double> fr=S_fr(phi_cell);

  BoxData<double> phi_face(vel.box());
  forallInPlace(computePhiFaceAve,phi_face,vel,wl,wr,fl,fr);

  return 0;
}
