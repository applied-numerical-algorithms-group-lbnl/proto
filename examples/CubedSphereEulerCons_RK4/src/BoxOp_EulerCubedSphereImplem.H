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
  // T dxi0 = this->dx()[1];
  T dxi0 = this->dx()[r_dir];
  double dr = (r1-r0)*dxi0;
  forallInPlace_p([]PROTO_LAMBDA(Point a_pt,
                                 Var<double>& a_rad,
                                 Var<double,DIM>& a_diagr,
                                 Var<double,DIM>& a_adjdiagr,
                                 Var<double>& a_dVolr,
                                 double a_r0,
                                 double a_r1,
                                 double a_dr,
                                 double a_dxi0,
                                 int a_rdir)
                  {
                    a_rad(0) = a_r0 + a_dr*a_pt[a_rdir];
                    T rlow = a_rad(0);
                    T rhigh = a_rad(0)+a_dr;
                    // a_rad(0) = rlow;
                    a_rad(0) = (rlow+rhigh)/2.0;
                    a_dVolr(0) = (rhigh*rhigh*rhigh - rlow*rlow*rlow)/3.0/a_dxi0;
                    a_diagr(0) = (rhigh - rlow)/a_dxi0;
                    a_diagr(1) = (rhigh + rlow)/2.0;
                    a_diagr(2) = (rhigh + rlow)/2.0;
                    a_adjdiagr(0) = rlow*rlow;
                    a_adjdiagr(1) = .5*(rhigh*rhigh - rlow*rlow)/a_dxi0;
                    a_adjdiagr(2) = .5*(rhigh*rhigh - rlow*rlow)/a_dxi0;
                  },a_bx,a_radius,a_Diagr,a_adjDr,a_dVolr,r0,r1,dr,dxi0,r_dir);
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

// ====================================================================================
// New factorization of operator implems start here.
// ====================================================================================
template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::operator()
  (
   State&                          a_Rhs,
   Array<State, DIM>&              a_fluxes,
   const State&                    a_USph,
   int                             a_block,
   T                               a_scale) const
{
  PR_TIME("BoxOp_Euler::operator()");
  HDF5Handler h5;
  Box bxRad = a_USph.box().grow(Point::Basis(0,2));
  BoxData<T,DIM,MEM> Dr(bxRad);
  BoxData<T,1,MEM> radius(bxRad);
  BoxData<T,DIM,MEM> adjDr(bxRad);
  BoxData<T,1,MEM> dVolr(bxRad);
  radialMetrics(radius,Dr,adjDr,dVolr,bxRad);
  auto domainBox = this->layout().domain().getBlock(a_block).box();
  // Primitive Variables. 
  // Transform to primitive variables, spherical velocities.
  Box bx = a_USph.box();
  State W;
  State WBar;
  //cout << "inputbox = " << bx << endl;
  
  sphToPrim(W,WBar,a_USph,a_block);
  // Compute Fluxes, rhs.
  a_Rhs.setVal(0.0);
  Box outputbox = bx;
  // Loop over directions.
  T one = 1.0;
  T half = 0.5;
  T dxiR = this->dx()[0];
  T dxiPerp = this->dx()[1];  
  for (int dir = 0; dir < DIM; dir++)
    {
      PR_TIME("BoxOp_Euler::operator()::directional flux calculation");        
      // State WFace = Stencil<T>::CellToFace(dir)(W);
      // State WFaceBar = Stencil<T>::CellToFace(dir,Side::Lo,2)(WBar);      

      State WFace_L = Stencil<T>::CellToFace(dir)(W);
      State WFaceBar_L = Stencil<T>::CellToFace(dir,Side::Lo,2)(WBar);    

      State WFace_R = Stencil<T>::CellToFace(dir)(W);
      State WFaceBar_R = Stencil<T>::CellToFace(dir,Side::Lo,2)(WBar);   
 
      bool limiter = true;
      if (limiter){
        Stencil<double> idOP = 1.0*Shift(Point::Zeros());
        State WFace_L_temp = idOP(WFace_L);
        State WFaceBar_L_temp = idOP(WFaceBar_L);    
        State WFace_R_temp = idOP(WFace_R);
        State WFaceBar_R_temp = idOP(WFaceBar_R); 
        Limiters_4O(WFace_L,WFace_R,WFace_L_temp,WFace_R_temp,W,W,dir,dxiR, dxiPerp, dxiPerp);
        Limiters_4O(WFaceBar_L,WFaceBar_R,WFaceBar_L_temp,WFaceBar_R_temp,WBar,WBar,dir,dxiR, dxiPerp, dxiPerp);
        // Limiters_minmod(WFace_L,WFace_R,W,this->dx(),dir);
        // Limiters_minmod(WFaceBar_L,WFaceBar_R,WBar,this->dx(),dir);
      }

      // Fluxes on dir face.
      // Average of advective fluxes in spherical variables at the
      // face center, given point values of W.
      BoxData<T,NUMCOMPS> physFlux_L,physFlux_R;
      BoxData<T,1,MEM> adjDrDir = slice(adjDr,dir);  
      primToFlux(physFlux_L,WFace_L,WFaceBar_L,adjDrDir,a_block,dir);
      primToFlux(physFlux_R,WFace_R,WFaceBar_R,adjDrDir,a_block,dir);

      //Required mapping variables.
      Box bxface = physFlux_L.box();
      BoxData<T,DIM,MEM,DIM> deta_a = CubedSphereShell::detA_A<T,MEM>(bxface,dxiPerp,a_block,dir);
      BoxData<T,1,MEM> deta = CubedSphereShell::detA<T,MEM>(bxface,dxiPerp,dir); 
      Box bx = W.box();
      Box bxdir = bx.extrude(Point::Basis(dir));
      BoxData<T,DIM,MEM> ntrow = CubedSphereShell::adjA<T,MEM>(bx,dxiPerp,a_block,dir);
      Array<Array<uint,DIM>,6> permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
      Array<Array<int,DIM>,6> sign = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}}; 
      Point high = bx.high();
      Point low = bx.low();
      high[0] = low[0];
      Box bx0(low,high);
      BoxData<T,DIM,MEM,DIM> A_matrix(bx);
      T offseta = half;
      T offsetb = half;
      if (dir == 1) offseta = 0.0;
      if (dir == 2) offsetb = 0.0;
      forallInPlace_p(f_Amatrix,bx0,A_matrix,permute[a_block],sign[a_block],
                      dxiPerp,offseta,offsetb);
      spreadSlice(A_matrix);

      State fluxdir_L(bxface),fluxdir_R(bxface),fluxdir(bxface),fluxdir_Rus(bxface);
      primToMappedFlux(fluxdir_L,physFlux_L,WFace_L,WFaceBar_L,adjDrDir,deta_a,deta,ntrow,a_block,dir);
      primToMappedFlux(fluxdir_R,physFlux_R,WFace_R,WFaceBar_R,adjDrDir,deta_a,deta,ntrow,a_block,dir);
      primToMappedRusanovFlux(fluxdir_Rus,WFace_L,WFaceBar_L,WFace_R,WFaceBar_R,adjDrDir,deta_a,deta,A_matrix,ntrow,a_block,dir);

      fluxdir.setVal(0.0);
      fluxdir += fluxdir_R;
      fluxdir += fluxdir_L;
      fluxdir -= fluxdir_Rus;
      fluxdir *= 0.5;

      // Save fluxes, update RHS.
      fluxdir.copyTo(a_fluxes[dir]);
      BoxData<T,NUMCOMPS,MEM> rhsIncr;
      T dxInv = 1.0/(this->dx()[dir]);
      rhsIncr = Stencil<T>::FluxDivergence(dir)(fluxdir,dxInv);
      a_Rhs += rhsIncr;
      outputbox &= rhsIncr.box();
      
    }
  if (a_block == BLOCKTEST)
    {
      HDF5Handler h5;
      BoxData<T,NUMCOMPS,MEM> rhsRadial = checkRadial(a_Rhs,dxiPerp,a_block);
      h5.writePatch(dxiPerp,rhsRadial,"rhsRadial");
    }
  a_Rhs *= (a_scale);
  if (procID() == 0)cout << "output box = " << outputbox << endl;
}

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::primToMappedFlux(
                                          State& a_fluxdir,
                                          const State& a_physFlux,
                                          const State& a_WFace,
                                          const State& a_WFaceBar,
                                          const BoxData<T,1,MEM>& a_adjDrDir,
                                          const BoxData<T,DIM,MEM,DIM>& a_deta_a,
                                          const BoxData<T,1,MEM>& a_deta,
                                          const BoxData<T,DIM,MEM>& a_ntrow,
                                          uint a_block,
                                          uint a_dir) const
{
  // Multiply momentum and B fluxes by det(A) A.
  BoxData<T,DIM,MEM> fluxdirMom = slice<T,NUMCOMPS,DIM,MEM>(a_fluxdir,iMOMX);
  BoxData<T,DIM,MEM> physFluxMom = slice<T,NUMCOMPS,DIM,MEM>(a_physFlux,iMOMX);
  BoxData<T,DIM,MEM> tmp =
    Operator::_faceMatrixProductAB(a_deta_a,physFluxMom,a_deta_a,physFluxMom,a_dir);
  tmp.copyTo(fluxdirMom);

  BoxData<T,DIM,MEM> fluxdirB = slice<T,NUMCOMPS,DIM,MEM>(a_fluxdir,iBX);
  BoxData<T,DIM,MEM> physFluxB = slice<T,NUMCOMPS,DIM,MEM>(a_physFlux,iBX);
  BoxData<T,DIM,MEM> tmp2 =
    Operator::_faceMatrixProductAB(a_deta_a,physFluxB,a_deta_a,physFluxB,a_dir);
  tmp2.copyTo(fluxdirB);

  // Scale rho, rho E fluxes by det(A),
  forallInPlace([ ] PROTO_LAMBDA
                (Var<T,NUMCOMPS,MEM>& a_advflux,
                  const Var<T,NUMCOMPS,MEM>& a_fluxin,
                  const Var<T,1,MEM>& a_deta,
                  const Var<T,1,MEM>& a_adjddir)
                {                    
                  a_advflux(iRHO) = a_deta(0)*a_fluxin(iRHO);
                  a_advflux(iE) =
                    a_deta(0)*a_fluxin(iE);  
                },a_fluxdir,a_physFlux,a_deta,a_adjDrDir);
    
  for (int comp = 0; comp < NUMCOMPS;comp++)
    {
      BoxData<T,1,MEM> advfluxcomp = slice(a_fluxdir,comp);
      BoxData<T,1,MEM> advfluxprod = faceProduct(advfluxcomp,a_adjDrDir,a_dir);
      advfluxprod.copyTo(advfluxcomp);
    }
    
  // Momentum fluxes - pressure forces.      
  {
    BoxData<T,1,MEM> pslice = slice(a_WFace,iP);
    BoxData<T,DIM,MEM> tmp2(a_ntrow.box());
    a_ntrow.copyTo(tmp2);
    BoxData<T,1,MEM> tmp3 = faceProduct(a_adjDrDir,pslice,a_dir);
    
    for (int dirvec = 0; dirvec < DIM; dirvec++)
      {
        BoxData<T> ntslice = slice(tmp2,dirvec);          
        ntslice *= tmp3;
      }
    fluxdirMom += tmp2;
  }
}

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::primToMappedRusanovFlux(
                                  State& a_fluxdir,
                                  const State& a_WFace_L,
                                  const State& a_WFaceBar_L,
                                  const State& a_WFace_R,
                                  const State& a_WFaceBar_R,
                                  const BoxData<T,1,MEM>& a_adjDdir,
                                  const BoxData<T,DIM,MEM,DIM>& a_deta_a,
                                  const BoxData<T,1,MEM>& a_deta,
                                  const BoxData<T,DIM,MEM,DIM>& a_A,
                                  const BoxData<T,DIM,MEM>& a_ntrow,
                                  uint a_block,
                                  uint a_dir) const
{

  // Compute fluxes at face centers.
  PR_TIMERS("primToMappedRusanovFlux");
  Stencil<T> deconvolveStencilRadial = (1.0*Shift(Point::Zeros())
                                    - (1.0/24.0)*Stencil<T>::Derivative(2,0,2));
  // Deconvolve W to obtain point values.
  BoxData<T,NUMCOMPS,MEM> WFacePoint_L = _deconvolveFace(a_WFace_L,a_WFaceBar_L,a_dir);  
  BoxData<T,1,MEM> kePoint_L = kineticEnergyW(WFacePoint_L,a_block,a_dir);
  BoxData<T,1,MEM> mePoint_L = magneticEnergyW(WFacePoint_L,a_block,a_dir);

  BoxData<T,1,MEM> adjDdirPoint;
  if (a_dir != 0)
  {
    adjDdirPoint = deconvolveStencilRadial(a_adjDdir);
  }
  else
  {
    adjDdirPoint.define(a_adjDdir.box());
    a_adjDdir.copyTo(adjDdirPoint);
  }
 
  State q_L = forall<T,NUMCOMPS>
    ([ ] PROTO_LAMBDA
     (
      Var<T,NUMCOMPS>& a_u,
      const Var<T,NUMCOMPS>& a_wface,      
      const Var<T,1>& a_ke,
      const Var<T,1>& a_me,
      const Var<T>& a_adjddir,
      const T a_gamma)
     {
        a_u(iRHO ) = a_adjddir(0)*a_wface(iRHO);
        a_u(iMOMX) = a_adjddir(0)*a_wface(iRHO)*a_wface(iVX);
        a_u(iMOMY) = a_adjddir(0)*a_wface(iRHO)*a_wface(iVY);
        a_u(iMOMZ) = a_adjddir(0)*a_wface(iRHO)*a_wface(iVZ);
        a_u(iE   ) = a_adjddir(0)*(a_ke(0) + a_me(0) + a_wface(iP)/(a_gamma - 1.0));
        a_u(iBX  ) = a_adjddir(0)*a_wface(iBX);
        a_u(iBY  ) = a_adjddir(0)*a_wface(iBY);
        a_u(iBZ  ) = a_adjddir(0)*a_wface(iBZ);
     },WFacePoint_L,kePoint_L,mePoint_L,adjDdirPoint,m_gamma);

  // Compute Average conserved variables from pointwise values.
  State q_avg_L = (1.0*Shift(Point::Zeros()) +
                (1.0/24.0)*Stencil<T>::LaplacianFace(a_dir))(q_L);


  // Deconvolve W to obtain point values.
  BoxData<T,NUMCOMPS,MEM> WFacePoint_R = _deconvolveFace(a_WFace_R,a_WFaceBar_R,a_dir);  
  BoxData<T,1,MEM> kePoint_R = kineticEnergyW(WFacePoint_R,a_block,a_dir);
  BoxData<T,1,MEM> mePoint_R = magneticEnergyW(WFacePoint_R,a_block,a_dir);
 
  State q_R = forall<T,NUMCOMPS>
    ([ ] PROTO_LAMBDA
     (
      Var<T,NUMCOMPS>& a_u,
      const Var<T,NUMCOMPS>& a_wface,      
      const Var<T,1>& a_ke,
      const Var<T,1>& a_me,
      const Var<T>& a_adjddir,
      const T a_gamma)
     {
        a_u(iRHO ) = a_adjddir(0)*a_wface(iRHO);
        a_u(iMOMX) = a_adjddir(0)*a_wface(iRHO)*a_wface(iVX);
        a_u(iMOMY) = a_adjddir(0)*a_wface(iRHO)*a_wface(iVY);
        a_u(iMOMZ) = a_adjddir(0)*a_wface(iRHO)*a_wface(iVZ);
        a_u(iE   ) = a_adjddir(0)*(a_ke(0) + a_me(0) + a_wface(iP)/(a_gamma - 1.0));
        a_u(iBX  ) = a_adjddir(0)*a_wface(iBX);
        a_u(iBY  ) = a_adjddir(0)*a_wface(iBY);
        a_u(iBZ  ) = a_adjddir(0)*a_wface(iBZ);
     },WFacePoint_R,kePoint_R,mePoint_R,adjDdirPoint,m_gamma);

  // Compute Average conserved variables from pointwise values.
  State q_avg_R = (1.0*Shift(Point::Zeros()) +
                (1.0/24.0)*Stencil<T>::LaplacianFace(a_dir))(q_R);
  
  // Difference of conserved variables.
  q_avg_R -= q_avg_L;

  T one = 1.0;
  BoxData<T,DIM,MEM> fluxdirMom = slice<T,NUMCOMPS,DIM,MEM>(a_fluxdir,iMOMX);
  BoxData<T,DIM,MEM> Con_diff_Mom = slice<T,NUMCOMPS,DIM,MEM>(q_avg_R,iMOMX);
  forallInPlace(f_matVecProd,fluxdirMom,a_deta_a,Con_diff_Mom,one);

  BoxData<T,DIM,MEM> fluxdirB = slice<T,NUMCOMPS,DIM,MEM>(a_fluxdir,iBX);
  BoxData<T,DIM,MEM> Con_diff_B = slice<T,NUMCOMPS,DIM,MEM>(q_avg_R,iBX);
  forallInPlace(f_matVecProd,fluxdirB,a_deta_a,Con_diff_B,one);
    
  forallInPlace([ ] PROTO_LAMBDA
                (Var<T,NUMCOMPS,MEM>& a_Rusflux,
                  const Var<T,NUMCOMPS,MEM>& a_fluxin,
                  const Var<T,1,MEM>& a_deta)
                {                    
                  a_Rusflux(iRHO) = a_deta(0)*a_fluxin(iRHO);
                  a_Rusflux(iE  ) = a_deta(0)*a_fluxin(iE);  
                },a_fluxdir,q_avg_R,a_deta);

  // Find (|w_d| + c_fast*||A||_d) and multiply with a_fluxdir.
  forallInPlace([ ] PROTO_LAMBDA
               (Var<T,NUMCOMPS,MEM>& a_Rusflux,
                const Var<T,NUMCOMPS>& a_wface_L,
                const Var<T,NUMCOMPS>& a_wface_R,
                const Var<T,DIM,MEM,DIM>& a_Amatrix,
                const T a_gamma,
                const int a_d)
                {        
                  double rho = 0.5*(a_wface_L(iRHO)+a_wface_R(iRHO));
                  double w_d = 0.5*(a_wface_L(iVX+a_d)+a_wface_R(iVX+a_d));
                  double B_x = 0.5*(a_wface_L(iBX)+a_wface_R(iBX));
                  double B_y = 0.5*(a_wface_L(iBY)+a_wface_R(iBY));
                  double B_z = 0.5*(a_wface_L(iBZ)+a_wface_R(iBZ));
                  double B_mag = sqrt(B_x*B_x + B_y*B_y + B_z*B_z);
                  double p = 0.5*(a_wface_L(iP)+a_wface_R(iP));
                  if (p < 0.0) p = 0.0;
                  double af = sqrt(a_gamma*p/rho + B_mag*B_mag/4.0/M_PI/rho);  
                  double A_d_mag = sqrt(a_Amatrix(0,a_d)*a_Amatrix(0,a_d) + a_Amatrix(1,a_d)*a_Amatrix(1,a_d) + a_Amatrix(2,a_d)*a_Amatrix(2,a_d));
                  
                  for (int comp = 0; comp < NUMCOMPS;comp++)
                  {
                    a_Rusflux(comp) *= (abs(w_d) + af*A_d_mag);
                  }
                },a_fluxdir,a_WFace_L,a_WFace_R,a_A,m_gamma,a_dir);
}

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void BoxOp_EulerCubedSphere<T,MAP,MEM>::
primToFlux(
           State& a_physFlux,
           const State& a_WFace,
           const State& a_WFaceBar,
           const BoxData<T,1,MEM>& a_adjDdir,
           uint a_block,
           uint a_centering) const
{
  // Compute fluxes at face centers.
  PR_TIMERS("primToFlux");
  T dxiPerp = this->dx()[1];
  Stencil<T> deconvolveStencilRadial = (1.0*Shift(Point::Zeros())
                                    - (1.0/24.0)*Stencil<T>::Derivative(2,0,2));
  // Deconvolve W to obtain point values.
  BoxData<T,NUMCOMPS,MEM> WFacePoint = _deconvolveFace(a_WFace,a_WFaceBar,a_centering);  
  BoxData<T,1,MEM> kePoint = kineticEnergyW(WFacePoint,a_block,a_centering);
  BoxData<T,1,MEM> mePoint = magneticEnergyW(WFacePoint,a_block,a_centering);
  BoxData<T,1,MEM> vbPoint = vdotbW(WFacePoint,a_block,a_centering);
 
  BoxData<T,1,MEM> adjDdirPoint;
  if (a_centering != 0)
    {
      adjDdirPoint = deconvolveStencilRadial(a_adjDdir);
    }
  else
    {
      adjDdirPoint.define(a_adjDdir.box());
      a_adjDdir.copyTo(adjDdirPoint);
    }
  State fluxPoint = forall<T,NUMCOMPS>
    ([ ] PROTO_LAMBDA
     (
      Var<T,NUMCOMPS>& a_flux,
      const Var<T,NUMCOMPS>& a_wface,      
      const Var<T,1>& a_ke,
      const Var<T,1>& a_me,
      const Var<T,1>& a_vb,
      const Var<T>& a_adjddir,
      const T a_gamma,
      uint a_block,
      uint a_centering)
     {
       T velnorm = a_wface(a_centering+iVX);
       T veltan1 = a_wface((a_centering+1)%DIM + iVX);
       T veltan2 = a_wface((a_centering+2)%DIM + iVX);
       T bnorm = a_wface(a_centering+iBX);
       T btan1 = a_wface((a_centering+1)%DIM + iBX);
       T btan2 = a_wface((a_centering+2)%DIM + iBX);
       a_flux(iRHO) = a_wface(iRHO)*velnorm;
       a_flux(a_centering+iMOMX) = a_flux(iRHO)*velnorm - bnorm*bnorm/4.0/M_PI;
       a_flux((a_centering+1)%DIM + iMOMX) = a_flux(iRHO)*veltan1 - bnorm*btan1/4.0/M_PI;
       a_flux((a_centering+2)%DIM + iMOMX) = a_flux(iRHO)*veltan2 - bnorm*btan2/4.0/M_PI;
       a_flux(iE) = (2.0*a_me(0) + a_ke(0) + a_wface(iP)*a_gamma/(a_gamma - 1.0))*velnorm - bnorm*a_vb(0)/4.0/M_PI;
       a_flux(iBX) = velnorm*a_wface(iBX) - a_wface(iVX)*bnorm;
       a_flux(iBY) = velnorm*a_wface(iBY) - a_wface(iVY)*bnorm;
       a_flux(iBZ) = velnorm*a_wface(iBZ) - a_wface(iVZ)*bnorm;
     },WFacePoint,kePoint,mePoint,vbPoint,adjDdirPoint,m_gamma,a_block,a_centering);
   if ((a_block == BLOCKTEST) && a_centering == 0)
    {
      HDF5Handler h5;
      h5.writePatch(this->dx(),kePoint,"kePoint" + to_string(a_centering));
      h5.writePatch(this->dx(),mePoint,"mePoint" + to_string(a_centering));
      h5.writePatch(this->dx(),WFacePoint,"WPointPF" + to_string(a_centering));          
    }
  // Compute Average fluxes from pointwise values.
  a_physFlux = (1.0*Shift(Point::Zeros()) +
                (1.0/24.0)*Stencil<T>::LaplacianFace(a_centering))(fluxPoint);
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

// ====================================================================================
// New factorization of operator implems ends here.
// ====================================================================================


// ====================================================================================
// Linear Viscocity implems start here.
// ====================================================================================

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::LinearVisc(
                        State& a_Rhs,
                        const State& a_USph,
                        int a_block,
                        T a_scale) const
{
  a_Rhs.setVal(0.0);
  PR_TIME("BoxOp_Euler::LinearVisc()");
  HDF5Handler h5;
  Box bxRad = a_USph.box().grow(Point::Basis(0,2));
  BoxData<T,DIM,MEM> Dr(bxRad);
  BoxData<T,1,MEM> radius(bxRad);
  BoxData<T,DIM,MEM> adjDr(bxRad);
  BoxData<T,1,MEM> dVolr(bxRad);
  radialMetrics(radius,Dr,adjDr,dVolr,bxRad);
  auto domainBox = this->layout().domain().getBlock(a_block).box();
  // Primitive Variables. 
  // Transform to primitive variables, spherical velocities.
  Box bx = a_USph.box();
  State W;
  State WBar;
  //cout << "inputbox = " << bx << endl;
  
  sphToPrim(W,WBar,a_USph,a_block);
  // Compute Fluxes, rhs.
  
  Box outputbox = bx;
  // Loop over directions.
  T one = 1.0;
  T half = 0.5;
  T dxiPerp = this->dx()[1];  
  for (int dir = 0; dir < DIM; dir++)
    {
      PR_TIME("BoxOp_Euler::operator()::directional flux calculation");        
      State WFace = Stencil<T>::CellToFace(dir)(W);
      State WFaceBar = Stencil<T>::CellToFace(dir,Side::Lo,2)(WBar);      
    
      BoxData<T,1,MEM> adjDrDir = slice(adjDr,dir);  

      Stencil<double> SPlus = 1.0*Shift(Point::Zeros()); // This SPlus, SMinus, and IdOp should give D5 at lower face.
      Stencil<double> SMinus = 1.0*Shift(-2*Point::Basis(dir));
      Stencil<double> IdOp = 1.0*Shift(-Point::Basis(dir));
      Stencil<double> D1 = SPlus - IdOp;
      Stencil<double> D2 = SPlus - 2.0*IdOp + SMinus;
      Stencil<double> D5 = (-1.0) * D1 * D2  * D2;
      State physFlux = D5(a_USph, 1.0/64);

      //Required mapping variables.
      Box bxface = physFlux.box();
      BoxData<T,DIM,MEM,DIM> deta_a = CubedSphereShell::detA_A<T,MEM>(bxface,dxiPerp,a_block,dir);
      BoxData<T,1,MEM> deta = CubedSphereShell::detA<T,MEM>(bxface,dxiPerp,dir); 
      Box bx = W.box();
      Box bxdir = bx.extrude(Point::Basis(dir));
      BoxData<T,DIM,MEM> ntrow = CubedSphereShell::adjA<T,MEM>(bx,dxiPerp,a_block,dir);
      Array<Array<uint,DIM>,6> permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
      Array<Array<int,DIM>,6> sign = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}}; 
      Point high = bx.high();
      Point low = bx.low();
      high[0] = low[0];
      Box bx0(low,high);
      BoxData<T,DIM,MEM,DIM> A_matrix(bx);
      T offseta = half;
      T offsetb = half;
      if (dir == 1) offseta = 0.0;
      if (dir == 2) offsetb = 0.0;
      forallInPlace_p(f_Amatrix,bx0,A_matrix,permute[a_block],sign[a_block],
                      dxiPerp,offseta,offsetb);
      spreadSlice(A_matrix);

      State fluxdir(bxface);
      primToMappedViscFlux(fluxdir,physFlux,WFace,WFaceBar,adjDrDir,deta_a,deta,A_matrix,ntrow,a_block,dir);

      // update RHS.
      BoxData<T,NUMCOMPS,MEM> rhsIncr;
      T dxInv = 1.0/(this->dx()[dir]);
      rhsIncr = Stencil<T>::FluxDivergence(dir)(fluxdir,dxInv);

      a_Rhs += rhsIncr;
      outputbox &= rhsIncr.box();   
    }
  if (a_block == BLOCKTEST)
    {
      HDF5Handler h5;
      BoxData<T,NUMCOMPS,MEM> rhsRadial = checkRadial(a_Rhs,dxiPerp,a_block);
      h5.writePatch(dxiPerp,rhsRadial,"rhsViscRadial");
    }
  a_Rhs *= (a_scale);

}

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::primToMappedViscFlux(
                                  State& a_fluxdir,
                                  const State& a_physFlux,
                                  const State& a_WFace,
                                  const State& a_WFaceBar,
                                  const BoxData<T,1,MEM>& a_adjDdir,
                                  const BoxData<T,DIM,MEM,DIM>& a_deta_a,
                                  const BoxData<T,1,MEM>& a_deta,
                                  const BoxData<T,DIM,MEM,DIM>& a_A,
                                  const BoxData<T,DIM,MEM>& a_ntrow,
                                  uint a_block,
                                  uint a_dir) const
{
  // Multiply momentum and B fluxes by det(A) A.
  BoxData<T,DIM,MEM> fluxdirMom = slice<T,NUMCOMPS,DIM,MEM>(a_fluxdir,iMOMX);
  BoxData<T,DIM,MEM> physFluxMom = slice<T,NUMCOMPS,DIM,MEM>(a_physFlux,iMOMX);
  BoxData<T,DIM,MEM> tmp =
    Operator::_faceMatrixProductAB(a_deta_a,physFluxMom,a_deta_a,physFluxMom,a_dir);
  tmp.copyTo(fluxdirMom);

  BoxData<T,DIM,MEM> fluxdirB = slice<T,NUMCOMPS,DIM,MEM>(a_fluxdir,iBX);
  BoxData<T,DIM,MEM> physFluxB = slice<T,NUMCOMPS,DIM,MEM>(a_physFlux,iBX);
  BoxData<T,DIM,MEM> tmp2 =
    Operator::_faceMatrixProductAB(a_deta_a,physFluxB,a_deta_a,physFluxB,a_dir);
  tmp2.copyTo(fluxdirB);

  // Scale rho, rho E fluxes by det(A),
  forallInPlace([ ] PROTO_LAMBDA
                (Var<T,NUMCOMPS,MEM>& a_fluxout,
                  const Var<T,NUMCOMPS,MEM>& a_fluxin,
                  const Var<T,1,MEM>& a_deta)
                {                    
                  a_fluxout(iRHO) = a_deta(0)*a_fluxin(iRHO);
                  a_fluxout(iE) = a_deta(0)*a_fluxin(iE);  
                },a_fluxdir,a_physFlux,a_deta);
    
  for (int comp = 0; comp < NUMCOMPS;comp++)
    {
      BoxData<T,1,MEM> advfluxcomp = slice(a_fluxdir,comp);
      BoxData<T,1,MEM> advfluxprod = faceProduct(advfluxcomp,a_adjDdir,a_dir);
      advfluxprod.copyTo(advfluxcomp);
    }

  // Find (|w_d| + c_fast*Adj(A)_{d,.}) and multiply with a_fluxdir.
  forallInPlace([ ] PROTO_LAMBDA
               (Var<T,NUMCOMPS,MEM>& a_fluxdir,
                const Var<T,NUMCOMPS>& a_wface,
                const Var<T,DIM,MEM>& a_ntrow,
                const Var<T,DIM,MEM,DIM>& a_Amatrix,
                const T a_gamma,
                const int a_d)
                {        
                  double rho = a_wface(iRHO);
                  double w_d = a_wface(iVX+a_d);
                  double B_x = a_wface(iBX);
                  double B_y = a_wface(iBY);
                  double B_z = a_wface(iBZ);
                  double B_mag = sqrt(B_x*B_x + B_y*B_y + B_z*B_z);
                  double p = a_wface(iP);
                  if (p < 0.0) p = 0.0;
                  double af = sqrt(a_gamma*p/rho + B_mag*B_mag/4.0/M_PI/rho); 
                  double AdjA = sqrt(a_ntrow(0)*a_ntrow(0) + a_ntrow(1)*a_ntrow(1) + a_ntrow(2)*a_ntrow(2)); 
                  double A_d_mag = sqrt(a_Amatrix(0,a_d)*a_Amatrix(0,a_d) + a_Amatrix(1,a_d)*a_Amatrix(1,a_d) + a_Amatrix(2,a_d)*a_Amatrix(2,a_d));
                  
   
                  for (int comp = 0; comp < NUMCOMPS;comp++)
                  {
                    a_fluxdir(comp) *= (abs(w_d) + af*A_d_mag);
                  }
                },a_fluxdir,a_WFace,a_ntrow,a_A,m_gamma,a_dir);
}

// ====================================================================================
// Linear Viscocity implems end here.
// ====================================================================================


// ====================================================================================
// Limiter implems start here.
// ====================================================================================

template<typename T,MemType MEM>
PROTO_KERNEL_START
void del2_W_c_calcF(Var<T,NUMCOMPS,MEM>& del2_W_c,
                    const Var<T,NUMCOMPS,MEM>& a_W_ave,
                    const Var<T,NUMCOMPS,MEM>& a_W_ave_behind,
                    const Var<T,NUMCOMPS,MEM>& a_W_ave_ahead)
{
  for (int i=0; i< NUMCOMPS; i++) {
    del2_W_c(i) = a_W_ave_behind(i) - 2.0*a_W_ave(i) + a_W_ave_ahead(i);
  }
}
PROTO_KERNEL_END(del2_W_c_calcF, del2_W_c_calc)

template<typename T,MemType MEM>
PROTO_KERNEL_START
void del3_W_calcF(Var<T,NUMCOMPS,MEM>& a_del3_W,
                  const Var<T,NUMCOMPS,MEM>& a_del2_W_c,
                  const Var<T,NUMCOMPS,MEM>& a_del2_W_c_behind)
{
  for (int i=0; i< NUMCOMPS; i++) {
    a_del3_W(i) = a_del2_W_c(i) - a_del2_W_c_behind(i);
  }
}
PROTO_KERNEL_END(del3_W_calcF, del3_W_calc)

const bool printlim = false;

template<typename T,MemType MEM>
PROTO_KERNEL_START
void limiter_calcF(
  const Point& a_pt,
  Var<T,NUMCOMPS,MEM>& a_W_ave_L_ahead_limited,
  Var<T,NUMCOMPS,MEM>& a_W_ave_R_limited,
  const Var<T,NUMCOMPS,MEM>& a_W_ave_L_ahead,
  const Var<T,NUMCOMPS,MEM>& a_W_ave_R,
  const Var<T,NUMCOMPS,MEM>& a_W_ave,
  const Var<T,NUMCOMPS,MEM>& a_W_ave_ahead,
  const Var<T,NUMCOMPS,MEM>& a_W_ave_ahead2,
  const Var<T,NUMCOMPS,MEM>& a_W_ave_behind,
  const Var<T,NUMCOMPS,MEM>& a_W_ave_behind2,
  const Var<T,NUMCOMPS,MEM>& a_del2_W_c,
  const Var<T,NUMCOMPS,MEM>& a_del2_W_c_ahead,
  const Var<T,NUMCOMPS,MEM>& a_del2_W_c_behind,
  const Var<T,NUMCOMPS,MEM>& a_del3_W,
  const Var<T,NUMCOMPS,MEM>& a_del3_W_L,
  const Var<T,NUMCOMPS,MEM>& a_del3_W_R,
  const Var<T,NUMCOMPS,MEM>& a_del3_W_behind,
  const Var<T,NUMCOMPS,MEM>& a_del3_W_ahead,
  const Var<T,NUMCOMPS,MEM>& a_del3_W_ahead2,
  int a_dir,
  const double a_dx,
  const double a_dy,
  const double a_dz)
{
  double rhs = 1.0 - 1.0e-12;
  double rhs_test, lhs_test, rhs_test2, lhs_test2,rhs_test_a, lhs_test_a, rhs_rho, lhs_rho;
  double a_del_W_f_m, a_del_W_f_p, a_del2_W_f, a_del3_W_C, a_del2_W_lim, a_rho_i, a_del3_W_min, a_del3_W_max;
  for (int i=0; i< NUMCOMPS; i++)
  {
    a_del_W_f_m = a_W_ave(i) - a_W_ave_R(i);
    a_del_W_f_p = a_W_ave_L_ahead(i) - a_W_ave(i);
    a_del2_W_f = 6.0*(a_W_ave_R(i) - 2.0*a_W_ave(i) + a_W_ave_L_ahead(i));
    a_del3_W_C = (a_del2_W_c_ahead(i) - a_del2_W_c_behind(i))/2.0;
    if ((a_del2_W_c_behind(i) >= 0.0) && (a_del2_W_c(i) >= 0.0) && (a_del2_W_c_ahead(i) >= 0.0) && (a_del2_W_f >= 0.0)) {
      a_del2_W_lim = std::min({std::abs(a_del2_W_f), 1.25*std::abs(a_del2_W_c_behind(i)),
                                1.25*std::abs(a_del2_W_c(i)), 1.25*std::abs(a_del2_W_c_ahead(i))});
    } else if ((a_del2_W_c_behind(i) < 0.0) && (a_del2_W_c(i) < 0.0) && (a_del2_W_c_ahead(i) < 0.0) && (a_del2_W_f < 0.0)) {
      a_del2_W_lim = -1.0*std::min({std::abs(a_del2_W_f), 1.25*std::abs(a_del2_W_c_behind(i)),
                                    1.25*std::abs(a_del2_W_c(i)), 1.25*std::abs(a_del2_W_c_ahead(i))});
    } else {
      a_del2_W_lim = 0.0;
    }
    rhs_rho = 1.0e-12*std::max({std::abs(a_W_ave_behind2(i)), std::abs(a_W_ave_behind(i)), std::abs(a_W_ave(i)), std::abs(a_W_ave_ahead(i)), std::abs(a_W_ave_ahead2(i))});
    lhs_rho = std::abs(a_del2_W_f);
    if (lhs_rho <= rhs_rho) {
      a_rho_i = 0.0;
    } else {
      a_rho_i = a_del2_W_lim/a_del2_W_f;
    }
    a_del3_W_min = std::min({a_del3_W_behind(i), a_del3_W(i), a_del3_W_ahead(i), a_del3_W_ahead2(i)});
    a_del3_W_max = std::max({a_del3_W_behind(i), a_del3_W(i), a_del3_W_ahead(i), a_del3_W_ahead2(i)});

    if (printlim)
    {
      if (i==1 && a_dir == 0){
        if (a_pt[0]  >= 16 && a_pt[0]  <= 25 && a_pt [1] == 45 ) {
          cout << setw(30) << setprecision(18) << a_pt[0]
                << setw(30) << setprecision(18) << a_pt[1]
              << setw(30) << setprecision(18) << a_W_ave(i)
                << setw(30) << setprecision(18) << a_W_ave_L_ahead(i)
                << setw(30) << setprecision(18) << a_W_ave_R(i)
                << setw(30) << setprecision(18) << a_del_W_f_p
                << setw(30) << setprecision(18) << a_del_W_f_m
                << setw(30) << setprecision(18) << a_del2_W_c(i)
                << setw(30) << setprecision(18) << a_del2_W_f
                << setw(30) << setprecision(18) << a_del3_W(i)
                << setw(30) << setprecision(18) << a_del2_W_lim
                << setw(30) << setprecision(18) << a_rho_i
                << setw(30) << setprecision(18) << a_del3_W_min
                << setw(30) << setprecision(18) << a_del3_W_max
                << setw(30) << setprecision(18) << a_del3_W_L(i)
                << setw(30) << setprecision(18) << a_del3_W_R(i)
                << setw(30) << setprecision(18) << a_del3_W_C
              << endl;
        }
      }
    }

    if ((a_del_W_f_m * a_del_W_f_p <= 0.0) || ((a_W_ave(i) - a_W_ave_behind2(i))*(a_W_ave_ahead2(i) - a_W_ave(i)) <= 0.0)) {
      if (a_rho_i < rhs) {
        rhs_test = a_del3_W_max - a_del3_W_min;
        lhs_test = 0.1*std::max({std::abs(a_del3_W_min), std::abs(a_del3_W_max)});
        lhs_test_a = std::min({std::abs(a_del3_W_L(i)),std::abs(a_del3_W_R(i))})*1.1;
        rhs_test_a = std::abs(a_del3_W_C);
        if (lhs_test <= rhs_test){
        //if (!(lhs_test_a > rhs_test_a && a_del3_W_L(i)*a_del3_W_R(i) > 0.)) {
          rhs_test2 = 2.0*std::abs(a_del_W_f_m);
          lhs_test2 = std::abs(a_del_W_f_p);
          if (a_del_W_f_m * a_del_W_f_p < 0.0) {
            a_W_ave_L_ahead_limited(i) = a_W_ave(i) + a_rho_i*a_del_W_f_p;
          } else if (lhs_test2 >= rhs_test2) {
            a_W_ave_L_ahead_limited(i) = a_W_ave(i) + 2.0*(1.0-a_rho_i)*a_del_W_f_m + a_rho_i*a_del_W_f_p;
          }
          rhs_test2 = 2.0*std::abs(a_del_W_f_p);
          lhs_test2 = std::abs(a_del_W_f_m);
          if (a_del_W_f_m * a_del_W_f_p < 0.0) {
            a_W_ave_R_limited(i) = a_W_ave(i) - a_rho_i*a_del_W_f_m;
          } else if (lhs_test2 >= rhs_test2) {
            a_W_ave_R_limited(i) = a_W_ave(i) - 2.0*(1.0-a_rho_i)*a_del_W_f_p - a_rho_i*a_del_W_f_m;
          }

        }
      }
    } else {
      rhs_test2 = 2.0*std::abs(a_del_W_f_m);
      lhs_test2 = std::abs(a_del_W_f_p);
      if (lhs_test2 >= rhs_test2) {
        a_W_ave_L_ahead_limited(i) = a_W_ave(i) + 2.0*a_del_W_f_m;
      }
      rhs_test2 = 2.0*std::abs(a_del_W_f_p);
      lhs_test2 = std::abs(a_del_W_f_m);
      if (lhs_test2 >= rhs_test2) {
        a_W_ave_R_limited(i) = a_W_ave(i) - 2.0*a_del_W_f_p;
      }
    }
  }
}
PROTO_KERNEL_END(limiter_calcF, limiter_calc)


template<typename T,MemType MEM>
PROTO_KERNEL_START
void eta_tilde_d_calcF(Var<double,1>& a_eta_tilde_d,
                        const Var<T,NUMCOMPS,MEM>& a_W_ave,
                        const Var<T,NUMCOMPS,MEM>& a_W_ave_ahead,
                        const Var<T,NUMCOMPS,MEM>& a_W_ave_ahead2,
                        const Var<T,NUMCOMPS,MEM>& a_W_ave_behind,
                        const Var<T,NUMCOMPS,MEM>& a_W_ave_behind2,
                        int a_d)
{
  double z1=0.85, z0=0.75, delta = 0.33, si;
  double p, p_ahead, p_ahead2, p_behind, p_behind2;

  p = a_W_ave(iP) + (a_W_ave(iBX)*a_W_ave(iBX)+a_W_ave(iBY)*a_W_ave(iBY)+a_W_ave(iBZ)*a_W_ave(iBZ))/8.0/M_PI;
  p_ahead = a_W_ave_ahead(iP) + (a_W_ave_ahead(iBX)*a_W_ave_ahead(iBX)+a_W_ave_ahead(iBY)*a_W_ave_ahead(iBY)+a_W_ave_ahead(iBZ)*a_W_ave_ahead(iBZ))/8.0/M_PI;
  p_ahead2 = a_W_ave_ahead2(iP) + (a_W_ave_ahead2(iBX)*a_W_ave_ahead2(iBX)+a_W_ave_ahead2(iBY)*a_W_ave_ahead2(iBY)+a_W_ave_ahead2(iBZ)*a_W_ave_ahead2(iBZ))/8.0/M_PI;
  p_behind = a_W_ave_behind(iP) + (a_W_ave_behind(iBX)*a_W_ave_behind(iBX)+a_W_ave_behind(iBY)*a_W_ave_behind(iBY)+a_W_ave_behind(iBZ)*a_W_ave_behind(iBZ))/8.0/M_PI;
  p_behind2 = a_W_ave_behind2(iP) + (a_W_ave_behind2(iBX)*a_W_ave_behind2(iBX)+a_W_ave_behind2(iBY)*a_W_ave_behind2(iBY)+a_W_ave_behind2(iBZ)*a_W_ave_behind2(iBZ))/8.0/M_PI;

  double arg = std::abs(p_ahead-p_behind)/std::abs(p_ahead2-p_behind2);
  if (arg > z1) {si = 0.0;}
  if (arg < z0) {si = 1.0;}
  if ((arg<=z1) && (arg>=z0)) {si = 1.0 - (arg-z0)/(z1-z0);}

  double lhs1 = a_W_ave_behind(1+a_d)-a_W_ave_ahead(1+a_d);
  double lhs2 = std::abs(p_ahead-p_behind)/std::min({p_ahead,p_behind});
  if ((lhs1 > 0.0) && (lhs2 > delta)) {
    a_eta_tilde_d(0) = si;
  } else {
    a_eta_tilde_d(0) = 0.0;
  }
}
PROTO_KERNEL_END(eta_tilde_d_calcF, eta_tilde_d_calc)

template<typename T,MemType MEM>
PROTO_KERNEL_START
void eta_d_calcF(Var<T,1,MEM>& a_eta_d,
                  const Var<T,1,MEM>& a_eta_tilde_d,
                  const Var<T,1,MEM>& a_eta_tilde_d_ahead,
                  const Var<T,1,MEM>& a_eta_tilde_d_behind)
{
  a_eta_d(0) = std::min({a_eta_tilde_d(0),a_eta_tilde_d_ahead(0),a_eta_tilde_d_behind(0)});
}
PROTO_KERNEL_END(eta_d_calcF, eta_d_calc)

template<typename T,MemType MEM>
PROTO_KERNEL_START
void eta_calcF(Var<T,1,MEM>& a_eta,
                const Var<T,1,MEM>& a_eta_d,
                const Var<T,1,MEM>& a_eta_d_old)
{
  a_eta(0) = std::min({a_eta_d(0),a_eta_d_old(0)});
}
PROTO_KERNEL_END(eta_calcF, eta_calc)

template<typename T,MemType MEM>
PROTO_KERNEL_START
void Flat_calcF(Var<T,NUMCOMPS,MEM>& a_flattened,
                const Var<T,NUMCOMPS,MEM>& a_not_flattened,
                const Var<T,NUMCOMPS,MEM>& a_W_ave,
                const Var<double,1>& a_eta,
                int a_dir)
{
  for (int i=0; i< NUMCOMPS; i++) {
    a_flattened(i) = (1.0-a_eta(0))*a_not_flattened(i) + a_eta(0)*a_W_ave(i);       //This works
  }
}
PROTO_KERNEL_END(Flat_calcF, Flat_calc)

template<typename T,MemType MEM>
PROTO_KERNEL_START
void compare_calcF(const Point& a_pt,
        Var<T,NUMCOMPS,MEM>& a_W_diff,
                const Var<T,NUMCOMPS,MEM>& a_W_ave_R_limited,
                const Var<T,NUMCOMPS,MEM>& a_W_ave)
{
  for (int i=0; i< NUMCOMPS; i++) {
    if (a_W_ave_R_limited(i) != a_W_ave(i)) cout << a_pt[0] << a_pt[1] << endl;
  }
}
PROTO_KERNEL_END(compare_calcF, compare_calc)


template<
  typename T,
  template<MemType> typename MAP,
  MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::Limiters_4O(BoxData<double,NUMCOMPS>& a_W_ave_L_lim_flat,
                  BoxData<double,NUMCOMPS>& a_W_ave_R_lim_flat,
                  BoxData<double,NUMCOMPS>& a_W_ave_L,
                  BoxData<double,NUMCOMPS>& a_W_ave_R,
                  BoxData<double,NUMCOMPS>& a_W_ave,
                  BoxData<double,NUMCOMPS>& a_W_bar,
                  const int a_d,
                  const double a_dx,
                  const double a_dy,
                  const double a_dz) const
{
  static Stencil<double> m_copy;
  m_copy = 1.0*Shift(Point::Zeros());
  //Limiter Starts here
  State W_ave_L_ahead = alias(a_W_ave_L,Point::Basis(a_d)*(-1));
  //State W_ave_L_ahead_limited = alias(W_ave_L_ahead);  // This works, but won't be able to compare diffrences as if W_ave_L_ahead_limited changes, so does W_ave_L_ahead
  State W_ave_L_ahead_limited = m_copy(W_ave_L_ahead);
  //State W_ave_R_limited = alias(a_W_ave_R);
  State W_ave_R_limited = m_copy(a_W_ave_R);
  int limiter_apply = 1;
  int slope_flattening_apply = 1;
  if (limiter_apply == 1) {
    State W_ave_ahead = alias(a_W_ave,Point::Basis(a_d)*(-1));
    State W_ave_ahead2 = alias(a_W_ave,Point::Basis(a_d)*(-2));
    State W_ave_behind = alias(a_W_ave,Point::Basis(a_d)*(1));
    State W_ave_behind2 = alias(a_W_ave,Point::Basis(a_d)*(2));
    State del2_W_c  = forall<double,NUMCOMPS,MEM>( del2_W_c_calc, a_W_ave, W_ave_behind, W_ave_ahead);
    State del2_W_c_ahead = alias(del2_W_c,Point::Basis(a_d)*(-1));
    State del2_W_c_behind = alias(del2_W_c,Point::Basis(a_d)*(1));
    State del3_W = forall<double,NUMCOMPS,MEM>( del3_W_calc, del2_W_c, del2_W_c_behind);
    State del3_W_L = forall<double,NUMCOMPS,MEM>( del3_W_calc, del2_W_c, del2_W_c_behind);
    State del3_W_R = forall<double,NUMCOMPS,MEM>( del3_W_calc, del2_W_c_ahead, del2_W_c);
    State del3_W_behind = alias(del3_W,Point::Basis(a_d)*(1));
    State del3_W_ahead = alias(del3_W,Point::Basis(a_d)*(-1));
    State del3_W_ahead2 = alias(del3_W,Point::Basis(a_d)*(-2));
    forallInPlace_p( limiter_calc,W_ave_L_ahead_limited, W_ave_R_limited, W_ave_L_ahead, a_W_ave_R,
                      a_W_ave,W_ave_ahead, W_ave_ahead2, W_ave_behind, W_ave_behind2, del2_W_c, del2_W_c_ahead, del2_W_c_behind, del3_W,
                      del3_W_L, del3_W_R, del3_W_behind, del3_W_ahead, del3_W_ahead2, a_d, a_dx, a_dy, a_dz);
  }
  Box dbx0 = W_ave_L_ahead_limited.box();
  Box dbx1 = W_ave_L_ahead.box();
  State W_diff(dbx0);
  //forallInPlace_p(compare_calc, W_diff, W_ave_R_limited, a_W_ave_R);
  // Slope flattening starts here
  if (slope_flattening_apply == 1) {
    Scalar eta;
    Scalar eta_old;
    for (int d = 0; d < DIM; d++)
    {
      State W_bar_ahead = alias(a_W_bar,Point::Basis(d)*(-1));
      State W_bar_ahead2 = alias(a_W_bar,Point::Basis(d)*(-2));
      State W_bar_behind = alias(a_W_bar,Point::Basis(d)*(1));
      State W_bar_behind2 = alias(a_W_bar,Point::Basis(d)*(2));
      Scalar eta_tilde_d = forall<double,1,MEM>( eta_tilde_d_calc, a_W_bar, W_bar_ahead, W_bar_ahead2,
                                            W_bar_behind, W_bar_behind2, d);
      Scalar eta_tilde_d_ahead = alias(eta_tilde_d,Point::Basis(d)*(-1));
      Scalar eta_tilde_d_behind = alias(eta_tilde_d,Point::Basis(d)*(1));
      Scalar eta_d = forall<double,1,MEM>( eta_d_calc, eta_tilde_d, eta_tilde_d_ahead, eta_tilde_d_behind);
      if (d>0) {
        eta = forall<double,1,MEM>(eta_calc, eta_d, eta_old);
      } else {
        eta = alias(eta_d);
      }
      eta_old = alias(eta);
    }
    State W_ave_L_ahead_lim_flat = forall<double,NUMCOMPS,MEM>(Flat_calc, W_ave_L_ahead_limited, a_W_ave, eta, a_d);
    a_W_ave_R_lim_flat = forall<double,NUMCOMPS,MEM>(Flat_calc, W_ave_R_limited, a_W_ave, eta, a_d);
    a_W_ave_L_lim_flat = alias(W_ave_L_ahead_lim_flat,Point::Basis(a_d)*(1));
    //Slope flattening ends here
  } else {
    State W_ave_L_ahead_lim_flat = alias(W_ave_L_ahead_limited);
    a_W_ave_R_lim_flat = alias(W_ave_R_limited);
    a_W_ave_L_lim_flat = alias(W_ave_L_ahead_lim_flat,Point::Basis(a_d)*(1));
  }
}

template<typename T,MemType MEM>
PROTO_KERNEL_START
void minmod_calcF(const Point& a_pt,
                  Var<T,NUMCOMPS,MEM>& W_low_ahead_limited,
                  Var<T,NUMCOMPS,MEM>& W_high_limited,
                  const Var<T,NUMCOMPS,MEM>& a_W,
                  const Var<T,NUMCOMPS,MEM>& a_dWL,
                  const Var<T,NUMCOMPS,MEM>& a_dWR,
                  Var<T,DIM,MEM>& a_dx_sph,
                  Var<T,DIM,MEM>& a_dx_sph_behind,
                  Var<T,DIM,MEM>& a_dx_sph_ahead,
                  int a_d)
{
  double c1     = 1.0/(a_dx_sph(a_d) + a_dx_sph_ahead(a_d));
  double c2     = 1.0/(a_dx_sph(a_d) + a_dx_sph_behind(a_d));

  for (int i = 0; i< NUMCOMPS; i++){ 
    double WP     = a_W(i);
    double dWR    = a_dWR(i);
    double dWL    = a_dWL(i);

    // SLOPE  = dx*minmod( dWR*c1, dWL*c2 )

    // double signA = /abs(dWR*c1);
    double minmod;
    double A = dWR*c1;
    double B = dWL*c2;
    if (A*B < 0){
      minmod = 0;
    } else if (abs(A) < abs(B)){
      minmod = A;
    } else {
      minmod = B;
    }

    // double minmod = signA*max(0.0,min(abs(dWR*c1),signA*dWL*c2));
    double SLOPE = a_dx_sph(a_d)*minmod;

    W_low_ahead_limited(i) = WP + SLOPE;
    W_high_limited(i) = WP - SLOPE;

    // Reduce order if density or p become negative at faces
    // if (i==0 || i==4){
    if (i==4){
      if (W_low_ahead_limited(i) < 0) W_low_ahead_limited(i) = 0.;//WP;
      if (W_high_limited(i) < 0) W_high_limited(i) = 0.;//WP;
    }
  }
}
PROTO_KERNEL_END(minmod_calcF, minmod_calc)




template<
  typename T,
  template<MemType> typename MAP,
  MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::Limiters_minmod(BoxData<double,NUMCOMPS>& a_W_low_lim,
                  BoxData<double,NUMCOMPS>& a_W_high_lim,
                  BoxData<double,NUMCOMPS>& a_W,
                  BoxData<double,DIM>& a_dx_sph,
                  const int a_d) const
{
  static Stencil<double> m_difference;
  m_difference = -1.0*Shift(Point::Zeros()) + 1.0*Shift(Point::Basis(a_d)); 
  static Stencil<double> m_copy;
  m_copy = 1.0*Shift(Point::Zeros());

  State W_low_ahead = alias(a_W,Point::Basis(a_d)*(-1));
  State W_low_ahead_limited = m_copy(W_low_ahead);
  State W_high_limited = m_copy(a_W);
  int limiter_apply = 1;
  if (limiter_apply == 1) {
    State dWR = m_difference(a_W);
    State dWL = alias(dWR,Point::Basis(a_d)*(1));
    BoxData<double,DIM> a_dx_sph_behind = alias(a_dx_sph,Point::Basis(a_d)*(1));
    BoxData<double,DIM> a_dx_sph_ahead = alias(a_dx_sph,Point::Basis(a_d)*(-1));
    forallInPlace_p(minmod_calc,W_low_ahead_limited, W_high_limited, a_W, dWL,
                      dWR, a_dx_sph, a_dx_sph_behind, a_dx_sph_ahead, a_d);
  }
  a_W_high_lim = alias(W_high_limited);
  a_W_low_lim = alias(W_low_ahead_limited,Point::Basis(a_d)*(1));		
}


// ====================================================================================
// Limiter implems end here.
// ====================================================================================

#if 0
template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
static inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::
preStage(
         MBLevelBoxData<T,NUMCOMPS,MEM>& a_UStage)
{
   CubedSphereShell::consToSphInterpEuler(a_UStage,a_iop,dVolrLev, dx, 4);
}
template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
static inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::
postStage(
          State& a_RhsStage,
          MBBoundaryRegister<T,2*NUMCOMPS,MEM,PR_CELL>& a_blockReg,
          const MBLevelBoxData<T,NUMCOMPS,MEM>& a_dVolrLev)
{  
}
template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
static inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::radialGhost(State& a_U,
                                                BoxData<T,1,MEM>& a_radius,
                                                T a_dxi0,
                                                Box a_bxGhosted,
                                                Box a_blockBox
                                                )
{
  if (a_blockBox.high()[0] < a_bxGhosted.high()[0])
    {
      Point low = a_bxGhosted.low();
      low[0] = a_blockBox.high()[0] + 1;
      Box bdryBoxHigh(low, a_bxGhosted.high());
      forallInPlace_p(f_radialBCs, bdryBoxHigh, a_USph, a_radius, a_dxi0, m_gamma, m_eps, thickness);
    }
  if (a_blockBox.low()[0] > a_bxGhosted.low()[0])
    {
      Point high = a_bxGhosted.high();
      high[0] = a_blockBox.low()[0] + 1;
      Box bdryBoxLow(a_bxGhosted.low(), high);
      forallInPlace_p(f_radialBCs, bdryBoxLow, a_USph, radius, a_dxi0, m_gamma, m_eps, thickness);
    }
}
#endif
