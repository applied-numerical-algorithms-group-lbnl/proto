#include "AdvectionRK4.H"

AdvectionState::AdvectionState(const double& domain_length,
                               const int& n_cells,
                               const double& vel,
                               const std::string& init_case):
  m_L(domain_length),
  m_N(n_cells),
  m_vel(vel)
{
  //TODO: add assertions for m_N<=0
  m_dx=m_L/m_N;
  Box box(Point({0}),Point({m_N-1}));
  m_phi.define(box.grow(2));
  m_phi.setVal(0.0);
  //InitializePhi(init_case);
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
  double max_w=std::max(wl(0),wr(0));
  double min_w=std::min(wl(0),wr(0));
  if(vel>0)
    phi_face(0)=max_w*fl(0)+min_w*fr(0);
  else
    phi_face(0)=max_w*fr(0)+min_w*fl(0);
}
PROTO_KERNEL_END(computePhiFaceAve_temp,computePhiFaceAve)

void AdvectionOp::operator()(AdvectionDX& k, double time, double& dt, AdvectionState& state)
{
  //k contains the previous intermediate step weighed by the current step weight.
  //The current state at which we compute the flux is state+dt*k 
  BoxData<double> curr_state(state.m_phi.box());
  (k.m_dF).copyTo(curr_state);
  curr_state*=dt;
  curr_state+=state.m_phi;
  
  //Compute Flux
  Stencil<double> S_c1=1.0*Shift::Zeros()-2.0*Shift::Basis(0,-1)+1.0*Shift::Basis(0,-2);
  Stencil<double> S_c2=1.0*Shift::Zeros()-1.0*Shift::Basis(0,-2);
  BoxData<double> cl1=S_c1(curr_state);
  BoxData<double> cl2=S_c2(curr_state);
  BoxData<double> cr1=alias(cl1,Point::Ones(-1));
  BoxData<double> cr2=alias(cl2,Point::Ones(-1));
  std::cout << "Phi cell domain: " << curr_state.box() << std::endl;
  std::cout << "cl1 domain: " << cl1.box() << std::endl;
  std::cout << "cl2 domain: " << cl2.box() << std::endl;
  std::cout << "cr1 domain: " << cr1.box() << std::endl;
  std::cout << "cr2 domain: " << cr2.box() << std::endl;

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

  //Compute divergence
  (k.m_dF).setVal(0.0);
  double dx_inv=1;
  dx_inv/=state.m_dx;
  double c=state.m_vel*dx_inv;
  Stencil<double> S_div=c*Shift::Zeros()-c*Shift::Basis(0,1);
  (k.m_dF)+=S_div(phi_face);  
}
