#include "AdvectionRK4.H"

AdvectionState::AdvectionState(const double& domain_length,
                               const int& n_cells,
                               const double& vel):
  m_L(domain_length),
  m_N(n_cells),
  m_vel(vel)
{
  //TODO: add assertions for m_N<=0
  m_dx=m_L/m_N;
  Box box(Point({0}),Point({m_N-1}));
  m_phi.define(box);
  m_phi.setVal(0.0);
}

void AdvectionState::setBoundaryConditions(BoxData<double>& state_ext)
{
  //std::cout << "state_ext max before: " << state_ext.absMax() << std::endl;
  //std::cout << "Before" << std::endl;
  //state_ext.print();
  Box box_ext=state_ext.box();
  Box box_valid=box_ext.grow(-2);
  Point shift=Point::Basis(0)*box_valid.size();
  Box inter=box_ext.shift(shift)&box_valid;
  //std::cout << "Intersection box: " << inter << std::endl;
  //std::cout << "Shifted intersection box: " << inter.shift(-1.0*shift) << std::endl;
  state_ext.copyTo(state_ext,box_ext.shift(shift)&box_valid,-1.0*shift);
  state_ext.copyTo(state_ext,box_ext.shift(-1.0*shift)&box_valid,shift);
  //std::cout << "After" << std::endl;
  //state_ext.print();
  //std::cout << "state_ext max after: " << state_ext.absMax() << std::endl; 
}

void AdvectionState::increment(const AdvectionDX& incr)
{
  m_phi+=incr.m_dF;
}

AdvectionDX::AdvectionDX()
{}

AdvectionDX::~AdvectionDX()
{}

void AdvectionDX::init(AdvectionState& state)
{
  m_dF.define(state.m_phi.box());
  m_dF.setVal(0.0);
}

void AdvectionDX::increment(const double& weight, const AdvectionDX& incr)
{
  BoxData<double> temp(incr.m_dF.box());
  (incr.m_dF).copyTo(temp);
  temp*=weight;
  m_dF+=temp;
}

void AdvectionDX::operator*=(const double& weight)
{
  m_dF*=weight;
}

AdvectionOp::AdvectionOp()
{
}

AdvectionOp::~AdvectionOp()
{
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
                            const double& vel,
                            const Var<double>& wl,
                            const Var<double>& wr,
                            const Var<double>& fl,
                            const Var<double>& fr)
{
  //double max_w=std::max(wl(0),wr(0));
  //double min_w=std::min(wl(0),wr(0));
  //if(vel>0)
  //  phi_face(0)=max_w*fl(0)+min_w*fr(0);
  //else
  //  phi_face(0)=max_w*fr(0)+min_w*fl(0);
  if(vel>0)
    phi_face(0)=fl(0);
  else
    phi_face(0)=fr(0);
}
PROTO_KERNEL_END(computePhiFaceAve_temp,computePhiFaceAve)

void AdvectionOp::advance(double time, double& dt, AdvectionState& state)
{
  BoxData<double> k0(state.m_phi.box());
  k0.setVal(0.0);
  BoxData<double> k1(state.m_phi.box());
  advance(k1,k0,time,dt,state);
  state.m_phi+=k1;

  /*
  BoxData<double> k1(state.m_phi.box());
  advance(k1,k0,time,dt,state);
  k1*=0.5;
  BoxData<double> k2(state.m_phi.box());
  double new_time=time+0.5*dt;
  advance(k2,k1,new_time,dt,state);
  k2*=0.5;
  BoxData<double> k3(state.m_phi.box());
  advance(k3,k2,new_time,dt,state);
  BoxData<double> k4(state.m_phi.box());
  new_time=time+dt;
  advance(k4,k3,new_time,dt,state);
  state.m_phi+=(1.0/6.0)*k1;
  state.m_phi+=(1.0/3.0)*k2;
  state.m_phi+=(1.0/3.0)*k3;
  state.m_phi+=(1.0/6.0)*k4;
  */
  /*
  BoxData<double> phi_face;
  BoxData<double> dF(state.m_phi.box());
  dF.setVal(0.0);
  ComputeFaceAve(phi_face,dF,time,dt,state);
  //BoxData<double> phi_face(Box(Point({0}),Point({state.m_phi.box().size()})));
  //std::cout << phi_face.box() << std::endl;
  //phi_face.setVal(0.0);
  //forallInPlace_p(evaluatePhiFace_p,phi_face,time,state.m_vel,state.m_dx);
  BoxData<double> exact_div(state.m_phi.box());
  exact_div.setVal(0.0);
  double dx_inv=1;
  dx_inv/=state.m_dx;
  double c=state.m_vel*dx_inv;
  //std::cout << "c: " << c << std::endl;
  Stencil<double> S_div=c*Shift::Zeros()-c*Shift::Basis(0,1);
  exact_div+=S_div(phi_face);
  //forallInPlace_p(evaluateExactDiv_p,exact_div,time,state.m_vel,state.m_dx);
  exact_div*=dt;
  state.m_phi+=exact_div;
  */
}

void AdvectionOp::advance(BoxData<double>& k_new, BoxData<double>& k_prev, double time, double& dt, AdvectionState& state)
{
  BoxData<double> phi_face;
  ComputeFaceAve(phi_face,k_prev,time,dt,state);
  double dx_inv=1;
  dx_inv/=state.m_dx;
  double c=state.m_vel*dx_inv;
  //std::cout << "c: " << c << std::endl;
  Stencil<double> S_div=c*Shift::Zeros()-c*Shift::Basis(0,1);
  k_new.setVal(0.0);
  k_new+=S_div(phi_face);
  k_new*=dt;
}

void AdvectionOp::ComputeFaceAve(BoxData<double>& phi_face, BoxData<double>& dF, double time, double& dt, AdvectionState& state)
{
  //k contains the previous intermediate step weighed by the current step weight and dt.
  //The current state at which we compute the flux is state+k
  BoxData<double> curr_state(state.m_phi.box().grow(2));
  dF.copyTo(curr_state);
  //curr_state*=dt;
  curr_state+=state.m_phi;
  //Enforce periodic boundary conditions
  AdvectionState::setBoundaryConditions(curr_state);

  //Compute Flux
  //TODO (debugging):
  //1) Simplify PhiFaceAve to not use WENO
  //2) Run out to time=0.125
  //3) Change example to phi0=1 (checks periodic boundary conditions)
  //4) Flip sign of velocity
  //Should see O(h^3) accuracy
  Stencil<double> S_c1=1.0*Shift::Zeros()-2.0*Shift::Basis(0,-1)+1.0*Shift::Basis(0,-2);
  Stencil<double> S_c2=1.0*Shift::Zeros()-1.0*Shift::Basis(0,-2);
  BoxData<double> cl1=S_c1(curr_state);
  BoxData<double> cl2=S_c2(curr_state);
  BoxData<double> cr1=alias(cl1,Point::Ones(-1));
  BoxData<double> cr2=alias(cl2,Point::Ones(-1));
  //std::cout << "Phi cell domain: " << curr_state.box() << std::endl;
  //std::cout << "cl1 domain: " << cl1.box() << std::endl;
  //std::cout << "cl2 domain: " << cl2.box() << std::endl;
  //std::cout << "cr1 domain: " << cr1.box() << std::endl;
  //std::cout << "cr2 domain: " << cr2.box() << std::endl;

  Box wbox=cl1.box()&cr1.box();
  BoxData<double> wl(wbox);
  wl.setVal(0.0);
  BoxData<double> wr(wbox);
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
  BoxData<double> fl=S_fl(curr_state);
  BoxData<double> fr=S_fr(curr_state);

  phi_face=forall<double>(computePhiFaceAve,state.m_vel,wl,wr,fl,fr);
}

void AdvectionOp::operator()(AdvectionDX& k, double time, double& dt, AdvectionState& state)
{
  //k contains the previous intermediate step weighed by the current step weight and dt.
  //The current state at which we compute the flux is state+k
  BoxData<double> curr_state(state.m_phi.box().grow(2));
  (k.m_dF).copyTo(curr_state);
  //curr_state*=dt;
  curr_state+=state.m_phi;
  //Enforce periodic boundary conditions
  AdvectionState::setBoundaryConditions(curr_state);

  //Compute Flux
  //TODO (debugging):
  //1) Simplify PhiFaceAve to not use WENO
  //2) Run out to time=0.125
  //3) Change example to phi0=1 (checks periodic boundary conditions)
  //4) Flip sign of velocity
  //Should see O(h^3) accuracy
  Stencil<double> S_c1=1.0*Shift::Zeros()-2.0*Shift::Basis(0,-1)+1.0*Shift::Basis(0,-2);
  Stencil<double> S_c2=1.0*Shift::Zeros()-1.0*Shift::Basis(0,-2);
  BoxData<double> cl1=S_c1(curr_state);
  BoxData<double> cl2=S_c2(curr_state);
  BoxData<double> cr1=alias(cl1,Point::Ones(-1));
  BoxData<double> cr2=alias(cl2,Point::Ones(-1));
  //std::cout << "Phi cell domain: " << curr_state.box() << std::endl;
  //std::cout << "cl1 domain: " << cl1.box() << std::endl;
  //std::cout << "cl2 domain: " << cl2.box() << std::endl;
  //std::cout << "cr1 domain: " << cr1.box() << std::endl;
  //std::cout << "cr2 domain: " << cr2.box() << std::endl;

  Box wbox=cl1.box()&cr1.box();
  BoxData<double> wl(wbox);
  wl.setVal(0.0);
  BoxData<double> wr(wbox);
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
  BoxData<double> fl=S_fl(curr_state);
  BoxData<double> fr=S_fr(curr_state);

  BoxData<double> phi_face=forall<double>(computePhiFaceAve,state.m_vel,wl,wr,fl,fr);
  //std::cout << "phi_face box: " << phi_face.box() << std::endl;
  phi_face.setVal(0.0);
  forallInPlace_p(evaluatePhiFace_p,phi_face,time,state.m_vel,state.m_dx);
  //phi_face.print();

  //Compute divergence
  (k.m_dF).setVal(0.0);
  double dx_inv=1;
  dx_inv/=state.m_dx;
  double c=state.m_vel*dx_inv;
  //std::cout << "c: " << c << std::endl;
  Stencil<double> S_div=c*Shift::Zeros()-c*Shift::Basis(0,1);
  (k.m_dF)+=S_div(phi_face);

  BoxData<double> exact_div((k.m_dF).box());
  forallInPlace_p(evaluateExactDiv_p,exact_div,time,state.m_vel,state.m_dx);
  BoxData<double> error(exact_div.box());
  error.setVal(0.0);
  //std::cout << "exact div box: " << exact_div.box() << std::endl;
  exact_div.copyTo(error);
  error-=(k.m_dF);
  //error.print();

  (k.m_dF)*=dt;
  //(k.m_dF).print();
  //std::cout << "flux domain: " << k.m_dF.box() << std::endl;
}
