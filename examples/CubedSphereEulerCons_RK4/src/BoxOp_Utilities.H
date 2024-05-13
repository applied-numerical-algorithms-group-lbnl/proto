template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::radialMetrics(
                                             BoxData<T, 1, MEM>& a_radius,
                                             BoxData<T,DIM,MEM>& a_Diagr,
                                             BoxData<T,DIM,MEM>& a_adjDr,
                                             BoxData<T, 1, MEM>& a_dVolr,
                                             Box a_bx) const
{
  int r_dir = CUBED_SPHERE_SHELL_RADIAL_COORD;
  double r0 = CUBED_SPHERE_SHELL_R0;
  double r1 = CUBED_SPHERE_SHELL_R1;
  double C_rad = CUBED_SPHERE_RADIAL_STRETCH;
  T dxi0 = this->dx()[r_dir];
  double dr = (r1-r0)*dxi0;
  forallInPlace_p([]PROTO_LAMBDA(Point a_pt,
                                 Var<double>& a_rad,
                                 Var<double,DIM>& a_diagr,
                                 Var<double,DIM>& a_adjdiagr,
                                 Var<double>& a_dVolr,
                                 double a_r0,
                                 double a_r1,
                                 double a_C_rad,
                                 double a_dr,
                                 double a_dxi0,
                                 int a_rdir)
                  {
                    // a_rad(0) = a_r0 + a_dr*a_pt[a_rdir];
                    // T rlow = a_rad(0);
                    // T rhigh = a_rad(0)+a_dr;

                    T R_t = (a_r1 - a_r0) / (exp(a_C_rad) - 1.0);
                    T etalow = a_pt[a_rdir]*a_dxi0;
				            T rlow = a_r0 + R_t * (exp(a_C_rad * etalow) - 1.0);
                    T etahigh = (a_pt[a_rdir] + 1.0)*a_dxi0;
                    T rhigh = a_r0 + R_t * (exp(a_C_rad * etahigh) - 1.0);

                    // a_rad(0) = rlow;
                    a_rad(0) = (rlow+rhigh)/2.0;
                    a_dVolr(0) = (rhigh*rhigh*rhigh - rlow*rlow*rlow)/3.0/a_dxi0;
                    a_diagr(0) = (rhigh - rlow)/a_dxi0;
                    a_diagr(1) = (rhigh + rlow)/2.0;
                    a_diagr(2) = (rhigh + rlow)/2.0;
                    a_adjdiagr(0) = rlow*rlow;
                    a_adjdiagr(1) = .5*(rhigh*rhigh - rlow*rlow)/a_dxi0;
                    a_adjdiagr(2) = .5*(rhigh*rhigh - rlow*rlow)/a_dxi0;
                  },a_bx,a_radius,a_Diagr,a_adjDr,a_dVolr,r0,r1,C_rad,dr,dxi0,r_dir);
}
template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline BoxData<T,NUMCOMPS,MEM> 
BoxOp_EulerCubedSphere<T,MAP,MEM>::checkRadial(State& a_inputs,
                                               T a_h,
                                               unsigned int a_block) const
{
  BoxData<T,NUMCOMPS,MEM> outputs(a_inputs.box());
  outputs.setToZero();
  Box bx = a_inputs.box();
  BoxData<T,1,MEM> deta = CubedSphereShell::detA<T,MEM>(bx.grow(2),a_h,0); 
  BoxData<T,1,MEM> detaInv = forall<T,1,MEM>
    ([ ] PROTO_LAMBDA
     (Var<T,1,MEM>& a_output, const Var<T,1,MEM>& a_deta)
     {
       a_output(0) = 1.0/a_deta(0);
     }
     ,deta);
  BoxData<T,1,MEM> rho = slice(outputs,iRHO);
  BoxData<T,1,MEM> rhoIn = slice(a_inputs,iRHO);
  rho += rhoIn*detaInv;
  BoxData<T,1,MEM> rhoE = slice(outputs,iE);
  BoxData<T,1,MEM> rhoEIn = slice(a_inputs,iE);
  rhoE += rhoEIn*detaInv;
  
  BoxData<T,DIM,MEM,DIM> deta_amat = CubedSphereShell::detA_A<T,MEM>(bx.grow(2),a_h,a_block,0); 
  BoxData<T,DIM,MEM,DIM> invcoef(bx);
  invcoef.setToZero();
  forallInPlace(f_matinv3by3,invcoef,deta_amat);

  BoxData<T,DIM,MEM> Mom = slice<T,NUMCOMPS,DIM,MEM>(outputs,iMOMX);
  BoxData<T,DIM,MEM> MomIn = slice<T,NUMCOMPS,DIM,MEM>(a_inputs,iMOMX);
  BoxData<T,DIM,MEM> tmpmom(outputs.box());
  T one = 1.0;  
  forallInPlace(f_matVecProd,tmpmom,invcoef,MomIn,one);
  tmpmom.copyTo(Mom);

  BoxData<T,DIM,MEM> B = slice<T,NUMCOMPS,DIM,MEM>(outputs,iBX);
  BoxData<T,DIM,MEM> BIn = slice<T,NUMCOMPS,DIM,MEM>(a_inputs,iBX);
  BoxData<T,DIM,MEM> tmpB(outputs.box());
  forallInPlace(f_matVecProd,tmpB,invcoef,BIn,one);
  tmpB.copyTo(B);

  return outputs;
}
template<
  typename T,
  template<MemType> typename MAP,
  MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::primToCons(
                                              State& a_JU,
                                              const State& a_WPoint,
                                              const BoxData<T>& a_dVolr,
                                              T a_gamma,
                                              T a_h,
                                              uint a_block) const
{
  PR_TIMERS("primToCons");
  Array<Array<uint,DIM>,6>
    permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
  Array<Array<int,DIM>,6>
    sign = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};
  Stencil<T> convolveStencil = (1.0*Shift(Point::Zeros())
                                + (1.0/24.0)*Stencil<T>::Laplacian());
  Stencil<T> deconvolveStencil = (1.0*Shift(Point::Zeros())
                                  - (1.0/24.0)*Stencil<T>::Laplacian());
  Stencil<T> convolveStencilRadial = (1.0*Shift(Point::Zeros())
                                      + (1.0/24.0)*Stencil<T>::Derivative(2,0,2));
  Stencil<T> deconvolveStencilRadial = (1.0*Shift(Point::Zeros())
                                        - (1.0/24.0)*Stencil<T>::Derivative(2,0,2));
  
  // Primitive variables at points -> <dVolr*(rho,vel_sph,rho*E)>

  BoxData<T,1,MEM> dVolPoint = deconvolveStencilRadial(a_dVolr);
  BoxData<T,1,MEM> ke = kineticEnergyW(a_WPoint,a_block,-1);
  BoxData<T,1,MEM> me = magneticEnergyW(a_WPoint,a_block,-1);
  BoxData<T,NUMCOMPS,MEM> UPoint =
    forall<T,NUMCOMPS,MEM>([ ] PROTO_LAMBDA
                           (Var<T,NUMCOMPS,MEM>& a_upoint,
                            const Var<T,NUMCOMPS>& a_wpoint,
                            const Var<T,1,MEM>& a_dvolpoint,
                            const Var<T,1,MEM>& a_ke,
                            const Var<T,1,MEM>& a_me,
                            T a_gamma)
                           {
                             a_upoint(iRHO) = a_wpoint(iRHO)*a_dvolpoint(0);
                             a_upoint(iMOMX) = a_wpoint(iVX)*a_upoint(iRHO);
                             a_upoint(iMOMY) = a_wpoint(iVY)*a_upoint(iRHO);
                             a_upoint(iMOMZ) = a_wpoint(iVZ)*a_upoint(iRHO);
                             a_upoint(iE) = (a_wpoint(iP)/(a_gamma - 1.0) + a_ke(0) + a_me(0))
                             *a_dvolpoint(0);
                             a_upoint(iBX) = a_wpoint(iBX)*a_dvolpoint(0);
                             a_upoint(iBY) = a_wpoint(iBY)*a_dvolpoint(0);
                             a_upoint(iBZ) = a_wpoint(iBZ)*a_dvolpoint(0);
                           },a_WPoint,dVolPoint,ke,me,a_gamma);
 
  BoxData<T,NUMCOMPS,MEM> U = convolveStencil(UPoint);
  
  // spherical velocities -> Cartesian velocities:
  // multiply <dvolr*rho u>, <dVolr*B> by det(A)*A using fourth-order product rule.
  
  Box bxOut = this->layout()[this->index()];
  a_JU.define(bxOut);
  Box bx = bxOut.grow(1);
  BoxData<T,DIM,MEM> JUMom = slice<T,NUMCOMPS,DIM,MEM>(a_JU,iMOMX);
  BoxData<T,DIM,MEM> UMom = slice<T,NUMCOMPS,DIM,MEM>(U,iMOMX);
  BoxData<T,DIM,MEM> JUB = slice<T,NUMCOMPS,DIM,MEM>(a_JU,iBX);
  BoxData<T,DIM,MEM> UB = slice<T,NUMCOMPS,DIM,MEM>(U,iBX);
  BoxData<T,DIM,MEM,DIM> amat = CubedSphereShell::detA_A<T,MEM>(bx,a_h,a_block,0);
  
  T one = 1.0;
  unsigned int block = this->block();
  Box blockBox = this->layout().blockDomainBox(block);  
  T twelveth = 1.0/12.0;
  forallInPlace(f_matVecProd,JUMom,amat,UMom,one);
  forallInPlace(f_matVecProd,JUB,amat,UB,one);
  for (int dir = 1; dir < DIM; dir++)
    {
      BoxData<T,DIM,MEM,DIM> damat = Stencil<T>::Derivative(1,dir,2)(amat);
      BoxData<T,DIM,MEM> dumom =
        CubedSphereShell::deriv3case(UMom,bxOut,blockBox,dir);
      BoxData<T,DIM,MEM> produ = forall<T,DIM,MEM>(f_matVecProd,damat,dumom,twelveth);
      JUMom += produ;
      BoxData<T,DIM,MEM> dub =
        CubedSphereShell::deriv3case(UB,bxOut,blockBox,dir);
      BoxData<T,DIM,MEM> prodb = forall<T,DIM,MEM>(f_matVecProd,damat,dub,twelveth);
      JUB += prodb;
    }
  // I'm not sure these copies are necessary.
  JUMom.copyTo(UMom); 
  JUB.copyTo(UB);

  // multiply <dvolr*rho>, <dVolr*rho*E> by det(A) using fourth-order product rule.

  BoxData<T,1,MEM> deta = CubedSphereShell::detA<T,MEM>(bx,a_h,0);
  {
    BoxData<T,1,MEM> JURho = slice(a_JU,iRHO);
    BoxData<T,1,MEM> URho = slice(U,iRHO);
    BoxData<T,1,MEM> tmpRho = deta*URho;
    tmpRho.copyTo(JURho);
    BoxData<T,1,MEM> JUE = slice(a_JU,iE);
    BoxData<T,1,MEM> UE = slice(U,iE);
    BoxData<T,1,MEM> tmpE = UE * deta;
    tmpE.copyTo(JUE);   
    for (int dir = 1; dir < DIM; dir++)
      {        
        BoxData<T,1,MEM> ddeta = Stencil<T>::Derivative(1,dir,2)(deta,twelveth);
        BoxData<T,1,MEM> drho =
          CubedSphereShell::deriv3case(URho,bxOut,blockBox,dir);
        JURho += ddeta*drho;
        BoxData<T,1,MEM> de =
          CubedSphereShell::deriv3case(UE,bxOut,blockBox,dir);
        JUE += ddeta*de;
      }
  }
}
template<
  typename T,
  template<MemType> typename MAP,
  MemType MEM>
BoxData<T,1,MEM> 
BoxOp_EulerCubedSphere<T,MAP,MEM>::kineticEnergyW(
                                              const State& a_WPoint, 
                                              uint a_block,
                                              uint a_centering) const
{
  PR_TIMERS("kineticEnergyW");
  Array<Array<uint,DIM>,6>
     permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
  Array<Array<int,DIM>,6>
    sign = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};  
  // Compute fluxes at face centers.
  
  // Compute Cartesian UPoint.
  T dxiPerp = this->dx()[1];
  Box bx = a_WPoint.box();
  Point high = bx.high();
  Point low = bx.low();
  high[0] = low[0];
  Box bx0(low,high);
  BoxData<T,DIM,MEM,DIM> Amat(bx);
  T half = .5;
  T offseta = half;
  T offsetb = half;
  if (a_centering == 1)
    {
      offseta = 0.0;
    }
  if (a_centering == 2)
    {
      offsetb = 0.0;
    }
  forallInPlace_p(f_Amatrix,bx0,Amat,permute[a_block],sign[a_block],
                  dxiPerp,offseta,offsetb);
  spreadSlice(Amat);
  BoxData<T,1,MEM> ke =
    forall<T,1,MEM>([ ] PROTO_LAMBDA
                    (                
                     Var<T,1,MEM> a_ke,
                     const Var<T,DIM,MEM,DIM>& a_amat,
                     const Var<T,NUMCOMPS>& a_wpoint)
                    {
                      vecType velSph;
                      vecType velCart;
                  
                      velSph[0] = a_wpoint(iVX);       
                      velSph[1] = a_wpoint(iVY);
                      velSph[2] = a_wpoint(iVZ);
                      velCart[0] = 0.;       
                      velCart[1] = 0.;
                      velCart[2] = 0.;
                      T half = .5;
                      T one = 1.0;
                      a_ke(0) = 0.;
                      for (int dir0 = 0; dir0 < DIM; dir0++)
                        {
                          for (int dir1 = 0; dir1 < DIM; dir1++)
                            {
                              velCart[dir0] += a_amat(dir0,dir1)*velSph[dir1];
                            }
                          a_ke(0) += velCart[dir0]*velCart[dir0]/2;
                        }
                      a_ke(0) *= a_wpoint(iRHO);
                    },Amat,a_WPoint);
  return ke;
}

template<
  typename T,
  template<MemType> typename MAP,
  MemType MEM>
BoxData<T,1,MEM> 
BoxOp_EulerCubedSphere<T,MAP,MEM>::magneticEnergyW(
                                              const State& a_WPoint, 
                                              uint a_block,
                                              uint a_centering) const
{
  PR_TIMERS("magneticEnergyW");
  Array<Array<uint,DIM>,6>
     permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
  Array<Array<int,DIM>,6>
    sign = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};  
  // Compute fluxes at face centers.
  
  // Compute Cartesian UPoint.
  T dxiPerp = this->dx()[1];
  Box bx = a_WPoint.box();
  Point high = bx.high();
  Point low = bx.low();
  high[0] = low[0];
  Box bx0(low,high);
  BoxData<T,DIM,MEM,DIM> Amat(bx);
  T half = .5;
  T offseta = half;
  T offsetb = half;
  if (a_centering == 1)
    {
      offseta = 0.0;
    }
  if (a_centering == 2)
    {
      offsetb = 0.0;
    }
  forallInPlace_p(f_Amatrix,bx0,Amat,permute[a_block],sign[a_block],
                  dxiPerp,offseta,offsetb);
  spreadSlice(Amat);
  BoxData<T,1,MEM> me =
    forall<T,1,MEM>([ ] PROTO_LAMBDA
                    (                
                     Var<T,1,MEM> a_me,
                     const Var<T,DIM,MEM,DIM>& a_amat,
                     const Var<T,NUMCOMPS>& a_wpoint)
                    {
                      vecType bSph;
                      vecType bCart;
                  
                      bSph[0] = a_wpoint(iBX);       
                      bSph[1] = a_wpoint(iBY);
                      bSph[2] = a_wpoint(iBZ);
                      bCart[0] = 0.;       
                      bCart[1] = 0.;
                      bCart[2] = 0.;
                      T half = .5;
                      T one = 1.0;
                      a_me(0) = 0.;
                      for (int dir0 = 0; dir0 < DIM; dir0++)
                        {
                          for (int dir1 = 0; dir1 < DIM; dir1++)
                            {
                              bCart[dir0] += a_amat(dir0,dir1)*bSph[dir1];
                            }
                          a_me(0) += bCart[dir0]*bCart[dir0]/8.0/M_PI;
                        }
                    },Amat,a_WPoint);
  return me;
}

template<
  typename T,
  template<MemType> typename MAP,
  MemType MEM>
BoxData<T,1,MEM> 
BoxOp_EulerCubedSphere<T,MAP,MEM>::vdotbW(
                                        const State& a_WPoint, 
                                        uint a_block,
                                        uint a_centering) const
{
  PR_TIMERS("vdotbW");
  Array<Array<uint,DIM>,6>
     permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
  Array<Array<int,DIM>,6>
    sign = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};  
  // Compute fluxes at face centers.
  
  // Compute Cartesian UPoint.
  T dxiPerp = this->dx()[1];
  Box bx = a_WPoint.box();
  Point high = bx.high();
  Point low = bx.low();
  high[0] = low[0];
  Box bx0(low,high);
  BoxData<T,DIM,MEM,DIM> Amat(bx);
  T half = .5;
  T offseta = half;
  T offsetb = half;
  if (a_centering == 1)
    {
      offseta = 0.0;
    }
  if (a_centering == 2)
    {
      offsetb = 0.0;
    }
  forallInPlace_p(f_Amatrix,bx0,Amat,permute[a_block],sign[a_block],
                  dxiPerp,offseta,offsetb);
  spreadSlice(Amat);
  BoxData<T,1,MEM> vb =
    forall<T,1,MEM>([ ] PROTO_LAMBDA
                    (                
                     Var<T,1,MEM> a_vb,
                     const Var<T,DIM,MEM,DIM>& a_amat,
                     const Var<T,NUMCOMPS>& a_wpoint)
                    {
                      vecType velSph;
                      vecType velCart;
                      vecType bSph;
                      vecType bCart;

                      velSph[0] = a_wpoint(iVX);       
                      velSph[1] = a_wpoint(iVY);
                      velSph[2] = a_wpoint(iVZ);
                      velCart[0] = 0.;       
                      velCart[1] = 0.;
                      velCart[2] = 0.;
                  
                      bSph[0] = a_wpoint(iBX);       
                      bSph[1] = a_wpoint(iBY);
                      bSph[2] = a_wpoint(iBZ);
                      bCart[0] = 0.;       
                      bCart[1] = 0.;
                      bCart[2] = 0.;
                      T half = .5;
                      T one = 1.0;
                      a_vb(0) = 0.;
                      for (int dir0 = 0; dir0 < DIM; dir0++)
                        {
                          for (int dir1 = 0; dir1 < DIM; dir1++)
                            {
                              velCart[dir0] += a_amat(dir0,dir1)*velSph[dir1];
                              bCart[dir0] += a_amat(dir0,dir1)*bSph[dir1];
                            }
                          a_vb(0) += velCart[dir0]*bCart[dir0];
                        }
                    },Amat,a_WPoint);
  return vb;
}

template<
  typename T,
  template<MemType> typename MAP,
  MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::primToSph(
                                             State& a_USph,
                                             const State& a_WPoint,
                                             const BoxData<T,1,MEM>& a_dVolr, 
                                             uint a_block) const
{
  PR_TIMERS("primToSph");
  Stencil<T> convolveStencil = (1.0*Shift(Point::Zeros())
                                + (1.0/24.0)*Stencil<T>::Laplacian());
  Stencil<T> deconvolveStencil = (1.0*Shift(Point::Zeros())
                                  - (1.0/24.0)*Stencil<T>::Laplacian());
  Stencil<T> convolveStencilRadial = (1.0*Shift(Point::Zeros())
                                      + (1.0/24.0)*Stencil<T>::Derivative(2,0,2));
  Stencil<T> deconvolveStencilRadial = (1.0*Shift(Point::Zeros())
                                        - (1.0/24.0)*Stencil<T>::Derivative(2,0,2));
  // Primitive variables -> V_r*(rho,vel_sph,rho*E).

  BoxData<T,1,MEM> dVolPoint = deconvolveStencilRadial(a_dVolr);
  BoxData<T,1,MEM> ke = kineticEnergyW(a_WPoint,a_block,-1);
  BoxData<T,1,MEM> me = magneticEnergyW(a_WPoint,a_block,-1);
  BoxData<T,NUMCOMPS,MEM> UPoint =
    forall<T,NUMCOMPS,MEM>([ ] PROTO_LAMBDA
                           (Var<T,NUMCOMPS,MEM>& a_upoint,
                            const Var<T,NUMCOMPS>& a_wpoint,
                            const Var<T,1,MEM>& a_dvolpoint,
                            const Var<T,1,MEM>& a_ke,
                            T a_gamma)
                           {
                             a_upoint(iRHO) = a_wpoint(iRHO)*a_dvolpoint(0);
                             a_upoint(iMOMX) = a_wpoint(iVX)*a_upoint(iRHO);
                             a_upoint(iMOMY) = a_wpoint(iVY)*a_upoint(iRHO);
                             a_upoint(iMOMZ) = a_wpoint(iVZ)*a_upoint(iRHO);
                             a_upoint(iE) = (a_wpoint(iP)/(a_gamma - 1.0) + a_ke(0))
                             *a_dvolpoint(0);
                             a_upoint(iBX) = a_wpoint(iBX)*a_dvolpoint(0);
                             a_upoint(iBY) = a_wpoint(iBY)*a_dvolpoint(0);
                             a_upoint(iBZ) = a_wpoint(iBZ)*a_dvolpoint(0);
                           },a_WPoint,dVolPoint,ke,m_gamma);
 
  a_USph |= convolveStencil(UPoint);
  cout << "UPoint, USph boxes = "  << UPoint.box() << " , " << a_USph.box()  << endl;
}

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline BoxData<T,NUMCOMPS,MEM> 
BoxOp_EulerCubedSphere<T,MAP,MEM>::consToSph(State& a_input,
                                          Box a_bx,
                                          unsigned int a_block) const
{
  T dxiPerp = this->dx()[1];
  BoxData<T,NUMCOMPS,MEM> output(a_bx);
  output.setToZero();
  Box bx = a_bx;
  Box bxMetrics = a_bx.grow(Point::Ones() - Point::Basis(0));
  BoxData<T,1,MEM> deta = CubedSphereShell::detA<T,MEM>(bxMetrics,dxiPerp,0); 
  BoxData<T,1,MEM> detaInv = forall<T,1,MEM>
    ([ ] PROTO_LAMBDA
     (Var<T,1,MEM>& a_output, const Var<T,1,MEM>& a_deta)
     {
       a_output(0) = 1.0/a_deta(0);
     }
     ,deta);
  BoxData<T,1,MEM> rho = slice(output,iRHO);
  BoxData<T,1,MEM> rhoIn = slice(a_input,iRHO);
  
  rho += rhoIn*detaInv;
  BoxData<T,1,MEM> rhoE = slice(output,iE);
  BoxData<T,1,MEM> rhoEIn = slice(a_input,iE);
  rhoE += rhoEIn*detaInv;
  
  BoxData<T,DIM,MEM,DIM> deta_amat =
    CubedSphereShell::detA_A<T,MEM>(bxMetrics,dxiPerp,a_block,0); 
  BoxData<T,DIM,MEM,DIM> invcoef(bx);
  invcoef.setToZero();
  forallInPlace(f_matinv3by3,invcoef,deta_amat);
  BoxData<T,DIM,MEM> Mom = slice<T,NUMCOMPS,DIM,MEM>(output,iMOMX);
  BoxData<T,DIM,MEM> MomBar;
  BoxData<T,DIM,MEM> MomIn = slice<T,NUMCOMPS,DIM,MEM>(a_input,iMOMX);
  BoxData<T,DIM,MEM> tmpmom(output.box());
  T one = 1.0;  
  // forallInPlace(f_matVecProd,tmpmom,invcoef,MomIn,one);
  cellMatrixQuotientMatAndInv(Mom,MomBar,deta_amat,invcoef,MomIn);
  // tmpmom.copyTo(Mom);

  BoxData<T,DIM,MEM> B = slice<T,NUMCOMPS,DIM,MEM>(output,iBX);
  BoxData<T,DIM,MEM> BBar;
  BoxData<T,DIM,MEM> BIn = slice<T,NUMCOMPS,DIM,MEM>(a_input,iBX);
  BoxData<T,DIM,MEM> tmpB(output.box());
  // forallInPlace(f_matVecProd,tmpB,invcoef,BIn,one);
  cellMatrixQuotientMatAndInv(B,BBar,deta_amat,invcoef,BIn);
  // tmpB.copyTo(B);

  return output;
}
template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline BoxData<T,NUMCOMPS,MEM> 
BoxOp_EulerCubedSphere<T,MAP,MEM>::sphToCart(State& a_input,
                                         Box a_bx,
                                         unsigned int a_block) const
{
  T dxiPerp = this->dx()[1];
  BoxData<T,NUMCOMPS,MEM> output(a_bx);
  output.setToZero();
  Box bx = a_bx;
  BoxData<T,1,MEM> deta = CubedSphereShell::detA<T,MEM>(bx,dxiPerp,0); 
 
  BoxData<T,1,MEM> rho = slice(output,iRHO);
  BoxData<T,1,MEM> rhoIn = slice(a_input,iRHO);
  rho += rhoIn*deta;
  BoxData<T,1,MEM> rhoE = slice(output,iE);
  BoxData<T,1,MEM> rhoEIn = slice(a_input,iE);
  rhoE += rhoEIn*deta;
  
  BoxData<T,DIM,MEM,DIM> deta_amat = CubedSphereShell::detA_A<T,MEM>(bx,dxiPerp,a_block,0); 
  BoxData<T,DIM,MEM> Mom = slice<T,NUMCOMPS,DIM,MEM>(output,iMOMX);
  BoxData<T,DIM,MEM> MomIn = slice<T,NUMCOMPS,DIM,MEM>(a_input,iMOMX);
  BoxData<T,DIM,MEM> tmpmom(output.box());
  T one = 1.0;  
  forallInPlace(f_matVecProd,tmpmom,deta_amat,MomIn,one);
  tmpmom.copyTo(Mom);

  BoxData<T,DIM,MEM> B = slice<T,NUMCOMPS,DIM,MEM>(output,iBX);
  BoxData<T,DIM,MEM> BIn = slice<T,NUMCOMPS,DIM,MEM>(a_input,iBX);
  BoxData<T,DIM,MEM> tmpB(output.box());
  forallInPlace(f_matVecProd,tmpB,deta_amat,BIn,one);
  tmpB.copyTo(B);

  return output;
}

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::cartToSph(
                                             State& a_USph,
                                             const State& a_JU,
                                             const BoxData<T,1,MEM>& a_dVolr,
                                             uint a_block) const
{
  // Get spherical conserved variables.
  PR_TIMERS("cartToSph");
  T dxiPerp = this->dx()[1];
  Box bx = a_JU.box();
  auto domainBox = this->layout().domain().getBlock(a_block).box();
  
  Stencil<T> convolveStencil = (1.0*Shift(Point::Zeros())
                                + (1.0/24.0)*Stencil<T>::Laplacian());
  Stencil<T> deconvolveStencil = (1.0*Shift(Point::Zeros())
                                  - (1.0/24.0)*Stencil<T>::Laplacian());
  Stencil<T> convolveStencilRadial = (1.0*Shift(Point::Zeros())
                                      + (1.0/24.0)*Stencil<T>::Derivative(2,0,2));
  Stencil<T> deconvolveStencilRadial = (1.0*Shift(Point::Zeros())
                                        - (1.0/24.0)*Stencil<T>::Derivative(2,0,2));
  Box bxMetric = a_USph.box().grow(Point::Ones() - Point::Basis(0));
  BoxData<T,1,MEM> deta = CubedSphereShell::detA<T,MEM>(bxMetric,dxiPerp,0); 
  BoxData<T,1,MEM> detaInv = forall<T,1,MEM>
    ([ ] PROTO_LAMBDA
     (Var<T,1,MEM>& a_output, const Var<T,1,MEM>& a_deta)
     {
       a_output(0) = 1.0/a_deta(0);
     }
     ,deta);
  
  BoxData<T,DIM,MEM,DIM> deta_amat =
    CubedSphereShell::detA_A<T,MEM>(bxMetric,dxiPerp,a_block,0); 
  BoxData<T,DIM,MEM,DIM> invcoef(bx);
  invcoef.setToZero();
  forallInPlace(f_matinv3by3,invcoef,deta_amat);
 
  // Set up storage for radial variables.
  Box bx4 = bx.grow(-1);
  BoxData<T,NUMCOMPS,MEM> U(bx4);
  BoxData<T,NUMCOMPS,MEM> UBar(bx);
  // JU -> V_r*(rho,rho*vel_cart,rho*E) = JU / det(A);
  {
    PR_TIMERS("ConsToRadialCons");
    BoxData<T,1,MEM> tmp(bx); 
    BoxData<T,1,MEM> JRho = slice(a_JU,iRHO);
    BoxData<T,1,MEM> rho = slice(a_USph,iRHO);
    //tmp = JRho*detaInv;
    //tmp.copyTo(rho);
    cellQuotientValAndInv(rho,deta,detaInv,JRho,domainBox);
    BoxData<T,1,MEM> JRhoE = slice(a_JU,iE);
    BoxData<T,1,MEM> rhoE = slice(a_USph,iE);
    //tmp = JRhoE*detaInv;
    //tmp.copyTo(rhoE);
    cellQuotientValAndInv(rhoE,deta,detaInv,JRhoE,domainBox);
    // V_r*rho*v_cart -> V_r*rho*v_sphere = (det(A) A)^{-1} V_r*rho*v_cart;
    BoxData<T,DIM,MEM> JUMom = slice<T,NUMCOMPS,DIM,MEM>(a_JU,iMOMX);
    BoxData<T,DIM,MEM> USphMom = slice<T,NUMCOMPS,DIM,MEM>(a_USph,iMOMX);    
    cellMatrixQuotientMatAndInv(USphMom,deta_amat,invcoef,JUMom,domainBox);

    BoxData<T,DIM,MEM> JUB = slice<T,NUMCOMPS,DIM,MEM>(a_JU,iBX);
    BoxData<T,DIM,MEM> USphB = slice<T,NUMCOMPS,DIM,MEM>(a_USph,iBX);    
    cellMatrixQuotientMatAndInv(USphB,deta_amat,invcoef,JUB,domainBox);
  }
}
template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::sphToPrim(
                                          State& a_W,
                                          State& a_WBar,
                                          const State& a_USph,
                                          uint a_block) const
{
  // Get radial conserved variables.
  PR_TIMERS("USphToPrim");
  
  Stencil<T> convolveStencil = (1.0*Shift(Point::Zeros())
                                + (1.0/24.0)*Stencil<T>::Laplacian());
  Stencil<T> deconvolveStencil = (1.0*Shift(Point::Zeros())
                                  - (1.0/24.0)*Stencil<T>::Laplacian());
         
  BoxData<T,NUMCOMPS,MEM> UPoint = deconvolveStencil(a_USph);
  a_WBar = forall<T,NUMCOMPS,MEM>(f_sphToPrim,a_USph);
  a_W = forall<T,NUMCOMPS,MEM>(f_sphToPrim,UPoint);
  // Compute pressure.
  BoxData<T,1,MEM> keBar = kineticEnergyW(a_WBar,a_block,-1);
  BoxData<T,1,MEM> kePoint = kineticEnergyW(a_W,a_block,-1);
  BoxData<T,1,MEM> meBar = magneticEnergyW(a_WBar,a_block,-1);
  BoxData<T,1,MEM> mePoint = magneticEnergyW(a_W,a_block,-1);
    forallInPlace([ ] PROTO_LAMBDA
                (Var<T,NUMCOMPS>& a_w,
                 const Var<T,1,MEM>& a_ke,
                 const Var<T,1,MEM>& a_me,
                                  T a_gamma)
                {
                  a_w(iP) = (a_w(iP) - a_ke(0) - a_me(0))*(a_gamma - 1.0);
                },a_W,kePoint,mePoint,m_gamma);
  forallInPlace([ ] PROTO_LAMBDA
                (Var<T,NUMCOMPS>& a_w,
                 const Var<T,1,MEM>& a_ke,
                 const Var<T,1,MEM>& a_me,
                 T a_gamma)
                {
                  a_w(iP) = (a_w(iP) - a_ke(0) - a_me(0))*(a_gamma - 1.0);
                },a_WBar,keBar,meBar,m_gamma);
  T coef = 1.0/24.0;
  a_W += Stencil<T>::Laplacian()(a_WBar,coef);
}

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::sphToPrim_2O(
                                          State& a_W,
                                          const State& a_USph,
                                          uint a_block) const
{
  // Get radial conserved variables.
  PR_TIMERS("USphToPrim_2O");
   
  a_W = forall<T,NUMCOMPS,MEM>(f_sphToPrim,a_USph);
  // Compute pressure.
  BoxData<T,1,MEM> kePoint = kineticEnergyW(a_W,a_block,-1);
  BoxData<T,1,MEM> mePoint = magneticEnergyW(a_W,a_block,-1);
    forallInPlace([ ] PROTO_LAMBDA
                (Var<T,NUMCOMPS>& a_w,
                 const Var<T,1,MEM>& a_ke,
                 const Var<T,1,MEM>& a_me,
                                  T a_gamma)
                {
                  a_w(iP) = (a_w(iP) - a_ke(0) - a_me(0))*(a_gamma - 1.0);
                },a_W,kePoint,mePoint,m_gamma);
}