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
                    a_rad(0) = rlow;
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
  tmpmom.copyTo(B);

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
  
  // spherical velocities -> Cartesian velocities.
  a_JU.define(U.box());
  Box bx = U.box();
  BoxData<T,DIM,MEM> JUMom = slice<T,NUMCOMPS,DIM,MEM>(a_JU,iMOMX);
  BoxData<T,DIM,MEM> UMom = slice<T,NUMCOMPS,DIM,MEM>(U,iMOMX);
  BoxData<T,DIM,MEM> JUB = slice<T,NUMCOMPS,DIM,MEM>(a_JU,iBX);
  BoxData<T,DIM,MEM> UB = slice<T,NUMCOMPS,DIM,MEM>(U,iBX);
  BoxData<T,DIM,MEM,DIM> amat = CubedSphereShell::detA_A<T,MEM>(bx,a_h,a_block,0);
  
  T one = 1.0;
  forallInPlace(f_matVecProd,JUMom,amat,UMom,one);
  JUMom.copyTo(UMom); //Are these copies necessary?
  forallInPlace(f_matVecProd,JUB,amat,UB,one);
  JUB.copyTo(UB);

  // scale by det(A)

  BoxData<T,1,MEM> deta = CubedSphereShell::detA<T,MEM>(bx,a_h,0);
  {
    BoxData<T,1,MEM> JUcomp = slice(a_JU,iRHO);
    BoxData<T,1,MEM> Ucomp = slice(U,iRHO);
    BoxData<T,1,MEM> tmp;
    tmp = Ucomp * deta;
    tmp.copyTo(JUcomp);
  }
  {
    BoxData<T,1,MEM> JUcomp = slice(a_JU,iE);
    BoxData<T,1,MEM> Ucomp = slice(U,iE);
    BoxData<T,1,MEM> tmp;
    tmp = Ucomp * deta;
    tmp.copyTo(JUcomp);
  }
}
#if 0
template<
  typename T,
  template<MemType> typename MAP,
  MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::primToFlux(
                                          State& a_physFlux,
                                          const State& a_WFace,
                                          const State& a_WFaceBar,
                                          const BoxData<T,1,MEM>& a_adjDdir,
                                          uint a_block,
                                          uint a_centering) const
{
  // Compute fluxes at face centers.
  PR_TIMERS("primToFlux");
  Stencil<T> deconvolveStencilRadial = (1.0*Shift(Point::Zeros())
                                    - (1.0/24.0)*Stencil<T>::Derivative(2,0,2));
  // Deconvolve W to obtain point values.
  BoxData<T,NUMCOMPS,MEM> WFacePoint = _deconvolveFace(a_WFace,a_WFaceBar,a_centering);  
  BoxData<T,1,MEM> kePoint = kineticEnergyW(WFacePoint,a_block,a_centering);
 
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
      const Var<T>& a_adjddir,
      const T a_gamma,
      uint a_block,
      uint a_centering)
     {
       T velnorm = a_wface(a_centering+iVX);
       a_flux(iRHO) = a_wface(iRHO)*velnorm*a_adjddir(0);
       T veltan1 = a_wface((a_centering+1)%DIM + iVX);
       T veltan2 = a_wface((a_centering+2)%DIM + iVX);
       a_flux(a_centering+iMOMX) = a_flux(iRHO)*velnorm;
       a_flux((a_centering+1)%DIM + iMOMX) = a_flux(iRHO)*veltan1;
       a_flux((a_centering+2)%DIM + iMOMX) = a_flux(iRHO)*veltan2;
       a_flux(iE) =
       (a_ke(0) + a_wface(iP)*a_gamma/(a_gamma - 1.0))*a_adjddir(0)*velnorm;
     },WFacePoint,kePoint,adjDdirPoint,m_gamma,a_block,a_centering);
   if ((a_block == BLOCKTEST) && a_centering == 0)
    {
      HDF5Handler h5;
      h5.writePatch(this->dx(),kePoint,"kePoint" + to_string(a_centering));
      h5.writePatch(this->dx(),WFacePoint,"WPointPF" + to_string(a_centering));          
    }
  // Compute Average fluxes from pointwise values.
  a_physFlux = (1.0*Shift(Point::Zeros()) +
                (1.0/24.0)*Stencil<T>::LaplacianFace(a_centering))(fluxPoint);
}
#endif
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
  T half = .5;
  Box bx = a_WPoint.box();
  Point high = bx.high();
  Point low = bx.low();
  high[0] = low[0];
  Box bx0(low,high);
  BoxData<T,DIM,MEM,DIM> Amat(bx);
  forallInPlace_p(f_Amatrix,bx0,Amat,permute[a_block],sign[a_block],dxiPerp,half,half);
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
  T half = .5;
  Box bx = a_WPoint.box();
  Point high = bx.high();
  Point low = bx.low();
  high[0] = low[0];
  Box bx0(low,high);
  BoxData<T,DIM,MEM,DIM> Amat(bx);
  forallInPlace_p(f_Amatrix,bx0,Amat,permute[a_block],sign[a_block],dxiPerp,half,half);
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
#if 0
template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::consToPrim(
                                          State& a_W,
                                          State& a_WBar,
                                          const State& a_JU,
                                          const BoxData<T,1,MEM>& a_dVolr,
                                          T a_gamma,
                                          T a_h,
                                          uint a_block) const
{
  // Get radial conserved variables.
  PR_TIMERS("ConsToPrim");
  a_WBar.define(a_JU.box());
  
  a_WBar.setToZero();
  
  Stencil<T> convolveStencil = (1.0*Shift(Point::Zeros())
                                + (1.0/24.0)*Stencil<T>::Laplacian());
  Stencil<T> deconvolveStencil = (1.0*Shift(Point::Zeros())
                                  - (1.0/24.0)*Stencil<T>::Laplacian());
  Stencil<T> convolveStencilRadial = (1.0*Shift(Point::Zeros())
                                      + (1.0/24.0)*Stencil<T>::Derivative(2,0,2));
  Stencil<T> deconvolveStencilRadial = (1.0*Shift(Point::Zeros())
                                        - (1.0/24.0)*Stencil<T>::Derivative(2,0,2));
  Box bx = a_JU.box();
  BoxData<T,1,MEM> deta = CubedSphereShell::detA<T,MEM>(bx.grow(2),a_h,0); 
  BoxData<T,1,MEM> detaInv = forall<T,1,MEM>
    ([ ] PROTO_LAMBDA
     (Var<T,1,MEM>& a_output, const Var<T,1,MEM>& a_deta)
     {
       a_output(0) = 1.0/a_deta(0);
     }
     ,deta);
  
  BoxData<T,DIM,MEM,DIM> deta_amat =
    CubedSphereShell::detA_A<T,MEM>(bx.grow(2),a_h,a_block,0); 
  BoxData<T,DIM,MEM,DIM> invcoef(bx);
  invcoef.setToZero();
  forallInPlace(f_matinv3by3,invcoef,deta_amat);
 
  // Set up storage for radial variables.
  Box bx4 = bx.grow(-1);
  a_W.define(bx4);
  a_W.setToZero();
  BoxData<T,NUMCOMPS,MEM> U(a_JU.box());
  // JU -> V_r*(rho,rho*vel_cart,rho*E) = JU / det(A);
  {
    PR_TIMERS("ConsToRadialCons");
    BoxData<T,1,MEM> tmp(a_JU.box()); 
    BoxData<T,1,MEM> JRho = slice(a_JU,iRHO);
    BoxData<T,1,MEM> rho = slice(U,iRHO);
    tmp = JRho*detaInv;
    tmp.copyTo(rho);
    BoxData<T,1,MEM> JRhoE = slice(a_JU,iE);
    BoxData<T,1,MEM> rhoE = slice(U,iE);
    tmp = JRhoE*detaInv;
    tmp.copyTo(rhoE);
    // V_r*rho*v_cart -> V_r*rho*v_sphere = (detA A)^{-1} V_r*rho*v_cart;
    BoxData<T,DIM,MEM> JUMom = slice<T,NUMCOMPS,DIM,MEM>(a_JU,iMOMX);
    BoxData<T,DIM,MEM> UMom = slice<T,NUMCOMPS,DIM,MEM>(U,iMOMX);
    BoxData<T,DIM,MEM> tmpmom(JUMom.box());
    T one = 1.0;  
    forallInPlace(f_matVecProd,tmpmom,invcoef,JUMom,one);
    tmpmom.copyTo(UMom);
    // Deconvolve-pointvise-convolve to go from volr rho Ur to W.
    // divide by dVolr and deconvolve.
    BoxData<T> tmp2(U.box());
    
    for (int comp = 0; comp < NUMCOMPS;comp++)
      {      
        BoxData<T> Ucomp = slice(U,comp);
        tmp2 = cellQuotient(Ucomp,a_dVolr);
        tmp2.copyTo(Ucomp);
    }
  }
  {
    PR_TIMERS("Cons Spherical to W, Wbar");
    BoxData<T,NUMCOMPS,MEM> UPoint = deconvolveStencil(U);    
    // Pointwise UrPoint -> WPoint.
    BoxData<T,NUMCOMPS,MEM> WPoint =
      forall<T,NUMCOMPS,MEM>([ ] PROTO_LAMBDA
                             (Var<T,NUMCOMPS,MEM>& a_wpoint,
                              Var<T,NUMCOMPS,MEM>& a_upoint)                           
                             {
                               Var<T,NUMCOMPS,MEM> wpoint;
                               a_wpoint(iRHO) = a_upoint(iRHO);
                               a_wpoint(iVX) = a_upoint(iMOMX)/a_upoint(iRHO);
                               a_wpoint(iVY) = a_upoint(iMOMY)/a_upoint(iRHO);
                               a_wpoint(iVZ) = a_upoint(iMOMZ)/a_upoint(iRHO);
                               a_wpoint(iP) = a_upoint(iE);                             
                             },UPoint);
    // Compute pressure.
    BoxData<T,1,MEM> kePoint = kineticEnergyW(WPoint,a_h,a_block,-1); 
    forallInPlace([ ] PROTO_LAMBDA
                  (Var<T,NUMCOMPS>& a_wpoint,
                   const Var<T,1,MEM>& a_ke,
                   T a_gamma)
                  {
                    a_wpoint(iP) = (a_wpoint(iP) - a_ke(0))*(a_gamma - 1.0);
                  },WPoint,kePoint,a_gamma);
    a_W = convolveStencil(WPoint);
    WPoint.copyTo(a_WBar);
  }
}
template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::primToFlux(
                                          State& a_physFlux,
                                          const State& a_WFace,
                                          const State& a_WFaceBar,
                                          const BoxData<T,1,MEM>& a_adjDdir,
                                          T a_h,
                                          T a_gamma,
                                          uint a_block,
                                          uint a_centering) const
{
  // Compute fluxes at face centers.
  PR_TIMERS("primToFlux");
  Stencil<T> deconvolveStencilRadial = (1.0*Shift(Point::Zeros())
                                    - (1.0/24.0)*Stencil<T>::Derivative(2,0,2));
  // Deconvolve W to obtain point values.
  BoxData<T,NUMCOMPS,MEM> WFacePoint = _deconvolveFace(a_WFace,a_WFaceBar,a_centering);  
  BoxData<T,1,MEM> kePoint = kineticEnergyW(WFacePoint,a_h,a_block,a_centering);
 
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
      const Var<T>& a_adjddir,
      const T a_gamma,
      uint a_block,
      uint a_centering)
     {
       T velnorm = a_wface(a_centering+iVX);
       a_flux(iRHO) = a_wface(iRHO)*velnorm*a_adjddir(0);
       T veltan1 = a_wface((a_centering+1)%DIM + iVX);
       T veltan2 = a_wface((a_centering+2)%DIM + iVX);
       a_flux(a_centering+iVX) = a_flux(iRHO)*velnorm;
       a_flux((a_centering+1)%DIM + iVX) = a_flux(iRHO)*veltan1;
       a_flux((a_centering+2)%DIM + iVX) = a_flux(iRHO)*veltan2;
       a_flux(iE) =
       (a_ke(0) + a_wface(iP)*a_gamma/(a_gamma - 1.0))*a_adjddir(0)*velnorm;
     },WFacePoint,kePoint,adjDdirPoint,a_gamma,a_block,a_centering);
   if ((a_block == BLOCKTEST) && a_centering == 0)
    {
      HDF5Handler h5;
      h5.writePatch(a_h,kePoint,"kePoint" + to_string(a_centering));
      h5.writePatch(a_h,WFacePoint,"WPointPF" + to_string(a_centering));          
    }
  // Compute Average fluxes from pointwise values.
  a_physFlux = (1.0*Shift(Point::Zeros()) +
                (1.0/24.0)*Stencil<T>::LaplacianFace(a_centering))(fluxPoint);
}

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::operator()
  (
   State&                          a_Rhs,
   Array<State, DIM>&              a_fluxes,
   const State&                    a_JU,
   T                               a_scale) const
{
    auto block = this->block();
    int r_dir = CUBED_SPHERE_SHELL_RADIAL_COORD;
    auto domainBox = this->layout().domain().getBlock(block).box();
    int thickness = domainBox.sizes()[r_dir];
    auto B0 = this->box();
    auto B1 = B0.grow(6);
    
    BoxData<double> radius(B1);
    BoxData<double,DIM,HOST> Dr(B1);
    BoxData<double,DIM,HOST> adjDr(B1);
    BoxData<double,1,HOST> dVolr(B1);

    radialMetrics(radius, Dr, adjDr, dVolr, B1, thickness);
    
    double dx1 = this->dx()[(r_dir + 1) % DIM];
    double dx2 = this->dx()[(r_dir + 2) % DIM];
    PROTO_ASSERT(dx1 == dx2,
            "BoxOp_EulerCubedSphere | Error: Anisotropic grid spacing detected");

    double dx = dx1;

    this->operator()(a_Rhs, a_fluxes, a_JU, Dr, adjDr, dVolr, dx, block, a_scale);
}

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::operator()
  (
   State&                          a_Rhs,
   Array<State, DIM>&              a_fluxes,
   const State&                    a_JU,
   const BoxData<T,DIM,MEM>&       a_Dr,
   const BoxData<T,DIM,MEM>&       a_adjDr,
   const BoxData<T,1,MEM>&         a_dVolr,
   T                               a_h,
   int                             a_block,
   T                               a_scale) const
{
  PR_TIME("BoxOp_Euler::operator()");
  HDF5Handler h5;
  // Primitive Variables.

  // Transform to primitive variables, spherical velocities.
  Box bx = a_JU.box();
  State W;
  State WBar;    
  consToPrim(W,WBar,a_JU,a_dVolr,m_gamma,a_h,a_block);
  // Compute Fluxes, rhs.
  a_Rhs.setVal(0.0);
  // Loop over directions.
    
  for (int dir = 0; dir < DIM; dir++)
    {
      PR_TIME("BoxOp_Euler::operator()::directional flux calculation");        
      State WFace = Stencil<T>::CellToFace(dir)(W);
      State WFaceBar = Stencil<T>::CellToFace(dir,Side::Lo,2)(WBar);       
        
      // Fluxes on dir face.

      // Physical fluxes.
     
      // Average of advective fluxes in spherical variables at the
      // face center, given point values of W.
      BoxData<T,NUMCOMPS> physFlux;
      BoxData<T,1,MEM> adjDrDir = slice(a_adjDr,dir);  
      primToFlux(physFlux,WFace,WFaceBar,adjDrDir,a_h,m_gamma,a_block,dir);
      Box bxface = physFlux.box();
      State fluxdir(bxface);
        
      // Advective fluxes for Cartesian variables.
      // Scalar fluxes.
      // Geometric quantities.
      // Multiply momentum fluxes by det(A) A.
      // Advective fluxes for Cartesian momentum.
      T one = 1.0;
      BoxData<T,DIM,MEM,DIM> deta_a =
        CubedSphereShell::detA_A<T,MEM>(bxface,a_h,a_block,dir);
      BoxData<T,DIM,MEM> fluxdirMom = slice<T,NUMCOMPS,DIM,MEM>(fluxdir,iMOMX);
      BoxData<T,DIM,MEM> physFluxMom = slice<T,NUMCOMPS,DIM,MEM>(physFlux,iMOMX);
      BoxData<T,DIM,MEM> tmp =
        Operator::_faceMatrixProductAB(deta_a,physFluxMom,deta_a,physFluxMom,dir);
      tmp.copyTo(fluxdirMom);
      // Scale rho, rho E fluxes by det(A),
        
      BoxData<T,1,MEM> deta = CubedSphereShell::detA<T,MEM>(bxface,a_h,dir);   
      forallInPlace([ ] PROTO_LAMBDA
                    (Var<T,NUMCOMPS,MEM>& a_advflux,
                     const Var<T,NUMCOMPS,MEM>& a_fluxin,
                     const Var<T,1,MEM>& a_deta,
                     const Var<T,1,MEM>& a_adjddir)
                    {                    
                      a_advflux(iRHO) = a_deta(0)*a_fluxin(iRHO);
                      a_advflux(iE) = a_deta(0)*a_fluxin(iE);  
                    },fluxdir,physFlux,deta,adjDrDir);
        
      for (int comp = 0; comp < NUMCOMPS;comp++)
        {
          BoxData<T,1,MEM> advfluxcomp = slice(fluxdir,comp);
          BoxData<T,1,MEM> advfluxprod = faceProduct(advfluxcomp,adjDrDir,dir);
          advfluxprod.copyTo(advfluxcomp);
        }
        
      // Momentum fluxes - pressure forces.
      {
        BoxData<T,1,MEM> pslice = slice(WFace,iP);
        Box bx = W.box();
        Box bxdir = bx.extrude(Point::Basis(dir));
        BoxData<T,DIM,MEM> ntrow = CubedSphereShell::adjA<T,MEM>(bx,a_h,a_block,dir);
        BoxData<T,DIM,MEM> tmp2(ntrow.box());
        BoxData<T,1,MEM> tmp3 = faceProduct(adjDrDir,pslice,dir);
        for (int dirvec = 0; dirvec < DIM; dirvec++)
          {
            BoxData<T> ntslice = slice(ntrow,dirvec);          
            ntslice *= tmp3;
          }
       
        fluxdirMom += ntrow;
        
        // Save fluxes, update RHS.
        fluxdir.copyTo(a_fluxes[dir]);
        BoxData<T,NUMCOMPS,MEM> rhsIncr;
        T dx = this->dx()[dir];
        rhsIncr = Stencil<T>::FluxDivergence(dir)(fluxdir,1.0/dx);
        a_Rhs += rhsIncr;         
      }
    }
  if (a_block == BLOCKTEST)
    {
      HDF5Handler h5;
      BoxData<T,NUMCOMPS,MEM> rhsRadial = checkRadial(a_Rhs,a_h,a_block);
      h5.writePatch(a_h,rhsRadial,"rhsRadial");
    }
  a_Rhs *= (a_scale);
}
#endif
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
  T dxiPerp = this->dx()[1];  
  for (int dir = 0; dir < DIM; dir++)
    {
      PR_TIME("BoxOp_Euler::operator()::directional flux calculation");        
      State WFace = Stencil<T>::CellToFace(dir)(W);
      State WFaceBar = Stencil<T>::CellToFace(dir,Side::Lo,2)(WBar);       
        
      // Fluxes on dir face.

      // Physical fluxes.
     
      // Average of advective fluxes in spherical variables at the
      // face center, given point values of W.
      BoxData<T,NUMCOMPS> physFlux;
      BoxData<T,1,MEM> adjDrDir = slice(adjDr,dir);  
      primToFlux(physFlux,WFace,WFaceBar,adjDrDir,a_block,dir);
      Box bxface = physFlux.box();
      h5.writePatch(dxiPerp,physFlux,
                    "physFlux:dir:"+ to_string(dir) + "Block"+to_string(a_block));
      State fluxdir(bxface);
        

      BoxData<T,DIM,MEM,DIM> deta_a =
        CubedSphereShell::detA_A<T,MEM>(bxface,dxiPerp,a_block,dir);
      BoxData<T,DIM,MEM> fluxdirMom = slice<T,NUMCOMPS,DIM,MEM>(fluxdir,iMOMX);
      BoxData<T,DIM,MEM> physFluxMom = slice<T,NUMCOMPS,DIM,MEM>(physFlux,iMOMX);
      BoxData<T,DIM,MEM> tmp =
        Operator::_faceMatrixProductAB(deta_a,physFluxMom,deta_a,physFluxMom,dir);
      tmp.copyTo(fluxdirMom);

      BoxData<T,DIM,MEM> fluxdirB = slice<T,NUMCOMPS,DIM,MEM>(fluxdir,iBX);
      BoxData<T,DIM,MEM> physFluxB = slice<T,NUMCOMPS,DIM,MEM>(physFlux,iBX);
      BoxData<T,DIM,MEM> tmp2 =
        Operator::_faceMatrixProductAB(deta_a,physFluxB,deta_a,physFluxB,dir);
      tmp2.copyTo(fluxdirB);

      // Scale rho, rho E fluxes by det(A),
        
      BoxData<T,1,MEM> deta = CubedSphereShell::detA<T,MEM>(bxface,dxiPerp,dir);   
      forallInPlace([ ] PROTO_LAMBDA
                    (Var<T,NUMCOMPS,MEM>& a_advflux,
                     const Var<T,NUMCOMPS,MEM>& a_fluxin,
                     const Var<T,1,MEM>& a_deta,
                     const Var<T,1,MEM>& a_adjddir)
                    {                    
                      a_advflux(iRHO) = a_deta(0)*a_fluxin(iRHO);
                      a_advflux(iE) =
                        a_deta(0)*a_fluxin(iE);  
                    },fluxdir,physFlux,deta,adjDrDir);
        
      for (int comp = 0; comp < NUMCOMPS;comp++)
        {
          BoxData<T,1,MEM> advfluxcomp = slice(fluxdir,comp);
          BoxData<T,1,MEM> advfluxprod = faceProduct(advfluxcomp,adjDrDir,dir);
          advfluxprod.copyTo(advfluxcomp);
        }
        
      // Momentum fluxes - pressure forces.      
      {
        BoxData<T,1,MEM> pslice = slice(WFace,iP);
        Box bx = W.box();
        Box bxdir = bx.extrude(Point::Basis(dir));
        BoxData<T,DIM,MEM> ntrow = CubedSphereShell::adjA<T,MEM>(bx,dxiPerp,a_block,dir);
        BoxData<T,DIM,MEM> tmp2(ntrow.box());
        ntrow.copyTo(tmp2);
        BoxData<T,1,MEM> tmp3 = faceProduct(adjDrDir,pslice,dir);
        for (int dirvec = 0; dirvec < DIM; dirvec++)
          {
            // BoxData<T> ntslice = slice(ntrow,dirvec); // This will edit a_ntrow too. We don't want that when we implement Rusanov.         
            BoxData<T> ntslice = slice(tmp2,dirvec);          
            ntslice *= tmp3;
          }
       
        // fluxdirMom += ntrow;
        fluxdirMom += tmp2;
        
        // Save fluxes, update RHS.
        fluxdir.copyTo(a_fluxes[dir]);
        BoxData<T,NUMCOMPS,MEM> rhsIncr;
        T dxInv = 1.0/(this->dx()[dir]);
        rhsIncr = Stencil<T>::FluxDivergence(dir)(fluxdir,dxInv);
        a_Rhs += rhsIncr;
        outputbox &= rhsIncr.box();
      }
    }
  if (a_block == BLOCKTEST)
    {
      HDF5Handler h5;
      BoxData<T,NUMCOMPS,MEM> rhsRadial = checkRadial(a_Rhs,dxiPerp,a_block);
      h5.writePatch(dxiPerp,rhsRadial,"rhsRadial");
    }
  a_Rhs *= (a_scale);
  cout << "output box = " << outputbox << endl;
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
