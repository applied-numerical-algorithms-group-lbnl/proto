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
   T                               a_dt,
   T                               a_scale) const
{
  PR_TIMERS("BoxOp_Euler::operator()");
  int a_block = this->block();
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
  // a_Rhs.setVal(0.0);
  T dxiPerp = this->dx()[1];
#if 0
  if (procID() == 0)
    {
      HDF5Handler h5;
      //Box bx0 = a_Rhs.box();
      //BoxData<T,NUMCOMPS,MEM> dummy(bx0),dummyface(bxface);
      h5.writePatch(dxiPerp,a_USp,"USph_"+to_string(procID()));
      h5.writePatch(dxiPerp,W,"W_"+to_string(procID()));
      h5.writePatch(dxiPerp,WBar,"WBar_"+to_string(procID()));
    }
#endif
  Box outputbox = bx;
  // Loop over directions.
  T one = 1.0;
  T half = 0.5;
  T dxiR = this->dx()[0];
  
  BoxData<T,1,MEM> deta0(bx);
  BoxData<T,DIM,MEM> ntrow0(bx);
  BoxData<T,DIM,MEM,DIM> A_matrix0(bx);
  BoxData<T,1,MEM>speedMax(this->box());
  speedMax.setVal(0.0);
  
  for (int dir = 0; dir < DIM; dir++)
    {
      PR_TIME("BoxOp_Euler::operator()::directional flux calculation");        
      State WFace_L = Stencil<T>::CellToFace(dir)(W);   
      State WFace_R = Stencil<T>::CellToFace(dir)(W);
      
      //State WFace_L = Stencil<T>::CellToFace(dir,Side::Lo,5)(W);      
      //State WFace_R = ((1.0*Shift(Point::Basis(dir,-1)))*Stencil<T>::CellToFace(dir,Side::Hi,5))(W);
      State WFaceBar_L = Stencil<T>::CellToFace(dir,Side::Lo,2)(WBar);    
      State WFaceBar_R = Stencil<T>::CellToFace(dir,Side::Lo,2)(WBar);   

      //Required mapping variables.
      Box bxface = WFace_L.box().grow(-2);     
      BoxData<T,1,MEM> adjDrDir = slice(adjDr,dir);  
      BoxData<T,DIM,MEM,DIM> deta_a = CubedSphereShell::detA_A<T,MEM>(bxface,dxiPerp,a_block,dir);
      BoxData<T,1,MEM> deta = CubedSphereShell::detA<T,MEM>(bxface,dxiPerp,dir); 
      if (dir ==0) deta.copyTo(deta0);
      Box bx = W.box();
      Box bxdir = bx.extrude(Point::Basis(dir));
      BoxData<T,DIM,MEM> ntrow = CubedSphereShell::adjA<T,MEM>(bx,dxiPerp,a_block,dir);
      if (dir == 0) ntrow.copyTo(ntrow0);
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
      if (dir == 0) A_matrix.copyTo(A_matrix0);
      
      int limiter = ParseInputs::get_limiter();
      if (limiter != 1 && limiter != 2)
      {
        limiter = 1;
      }

      if (limiter == 1){
        PR_TIMERS("BoxOp_Euler::operator()::limiters");
        Stencil<double> idOP = 1.0*Shift(Point::Zeros());
        State WFace_L_temp = idOP(WFace_L);
        State WFaceBar_L_temp = idOP(WFaceBar_L);    
        State WFace_R_temp = idOP(WFace_R);
        State WFaceBar_R_temp = idOP(WFaceBar_R); 
        Limiters_4O(WFace_L,WFace_R,WFace_L_temp,WFace_R_temp,W,W,dir,dxiR, dxiPerp, dxiPerp);
        Limiters_4O(WFaceBar_L,WFaceBar_R,WFaceBar_L_temp,WFaceBar_R_temp,WBar,WBar,dir,dxiR, dxiPerp, dxiPerp);
      }

      if (limiter == 2){
        BoxData<T,DIM> dx(bxface);
        forallInPlace([ ] PROTO_LAMBDA
        ( Var<T,DIM,MEM>& a_dx_d,
          const Var<T,DIM,MEM>& a_adjA,
          const Var<T,1,MEM>& a_adjDrDir,
          const Var<T,1,MEM>& a_detA,
          const Var<T,1,MEM>& a_dVolr,
          const Var<T,DIM,MEM>& a_Dr,
          const T a_dxiPerp,
          const T a_dxiR,
          const int a_d)
        { 
          if (a_d > 0){
            T adjA_mag = sqrt(a_adjA(0)*a_adjA(0) + a_adjA(1)*a_adjA(1) + a_adjA(2)*a_adjA(2));
            a_dx_d(a_d) = a_dVolr(0) * a_detA(0) * a_dxiPerp / a_adjDrDir(0) / adjA_mag;
          } 
          if (a_d == 0){
            a_dx_d(a_d) = a_Dr(0) * a_dxiR;
          }                   
        },dx,ntrow,adjDrDir,deta,dVolr,Dr,dxiPerp,dxiR,dir);
        Limiters_minmod(WFace_L,WFace_R,W,dx,dir);
        Limiters_minmod(WFaceBar_L,WFaceBar_R,WBar,dx,dir);
      }

      // Fluxes on dir face.
      // Average of advective fluxes in spherical variables at the
      // face center, given point values of W.
      
      BoxData<T,NUMCOMPS> physFlux_L,physFlux_R;
      primToFlux(physFlux_L,WFace_L,WFaceBar_L,adjDrDir,a_block,dir);
      primToFlux(physFlux_R,WFace_R,WFaceBar_R,adjDrDir,a_block,dir);
#if 0
      if (procID() == 0)
        {
          HDF5Handler h5;
          //Box bx0 = a_Rhs.box();
          //BoxData<T,NUMCOMPS,MEM> dummy(bx0),dummyface(bxface);
          h5.writePatch(dxiPerp,physFlux_L,"physFlux_L_"+to_string(dir) + "_"+to_string(procID()));
          h5.writePatch(dxiPerp,physFlux_R,",physFlux_R_"+to_string(dir) + "_"+to_string(procID()));
        }
#endif
      Box bxface0 = a_fluxes[dir].box();
      State fluxdir_L(bxface),fluxdir_R(bxface),fluxdir(bxface0),fluxdir_Rus(bxface);
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
      T dxInv = -a_scale/(this->dx()[dir]);
      //DEBUG
#ifndef PR_DISABLE_INIT_CHECKS
      // bool infOrNanBlock =fluxdir.containsInfOrNAN();
      // if (infOrNanBlock)
      // {
      //   h5.writePatch(1.0,fluxdir,"MHD_IO_InfOrNAN_"+to_string(dir) + "_" +to_string(a_block));
      //   cout << "Inf or Nan at Point " << fluxdir.findInfOrNAN()
      //        << ", block = " << a_block
      //        << ", dir = " << dir << endl;
      // }
#endif
      rhsIncr = Stencil<T>::FluxDivergence(dir)(fluxdir,dxInv);
      a_Rhs += rhsIncr;
      outputbox &= rhsIncr.box();
#if 0
      {
        HDF5Handler h5;
        Box bx0 = a_Rhs.box();
        bxface = bx0.extrude(dir);
        BoxData<T,NUMCOMPS,MEM> dummy(bx0),dummyface(bxface);
        a_fluxes[dir].copyTo(dummyface);
        rhsIncr.copyTo(dummy);
        h5.writePatch(dxiPerp,dummyface,"fluxes_"+to_string(dir)+"_"+to_string(procID()));
        h5.writePatch(dxiPerp,dummy,"rhsIncr_"+to_string(dir)+"_"+to_string(procID()));
      }
#endif
    }
  // a_dt_next = 1.0/speedMax.absMax();
  // if (a_block == BLOCKTEST)
  //   {
  //     HDF5Handler h5;
      
  //     BoxData<T,NUMCOMPS,MEM> rhsRadial = checkRadial(a_Rhs,dxiPerp,a_block);
  //     h5.writePatch(dxiPerp,rhsRadial,"rhsRadial");
  //   }
  // if (procID() == 0)cout << "output box = " << outputbox << endl;
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
  PR_TIMERS("primToMappedFlux");
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
  // Compute Average fluxes from pointwise values.
  a_physFlux = (1.0*Shift(Point::Zeros()) +
                (1.0/24.0)*Stencil<T>::LaplacianFace(a_centering))(fluxPoint);
}


// ====================================================================================
// New factorization of operator implems ends here.
// ====================================================================================


// ====================================================================================
// Pre and Post implems start here.
// ====================================================================================
template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline T BoxOp_EulerCubedSphere<T,MAP,MEM>::
dtCFL(MBLevelBoxData<T,NUMCOMPS,MEM>& a_JU,
      MBInterpOp& a_iop,
      MBLevelBoxData<T,1,MEM>& a_dVolr)
{
  auto layout = a_JU.layout();
  auto map = CubedSphereShell::Map(layout, ghost());
  auto eulerOp = CubedSphereShell::Operator<BoxOp_EulerCubedSphere, double, HOST>(map);
  MBLevelBoxData<T,NUMCOMPS,MEM> JUTemp(layout,a_JU.ghost());
  a_JU.copyTo(JUTemp);
  Reduction<T,Operation::Max,MEM> dtinv;
  dtinv.reset();
  CubedSphereShell::consToSphInterpEuler(JUTemp,a_iop,a_dVolr,4);
  for (auto dit : layout)
    {
      int block = layout.block(dit);
      BoxData<T,NUMCOMPS,MEM> W(layout[dit]);
      //BoxData<T,NUMCOMPS,MEM> WBar(JUTemp[dit].box());
      BoxData<T,NUMCOMPS,MEM>& sph_i = JUTemp[dit];
      BoxData<T,NUMCOMPS,MEM> sphTemp(layout[dit]);
      sph_i.copyTo(sphTemp);
      eulerOp[dit].sphToPrim_2O(W,sphTemp,block);
      eulerOp[dit].dtInv(dtinv,W);
    }
  return 1.0/dtinv.fetch(); 
}
template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::
StageLevel(MBLevelBoxData<T,NUMCOMPS,MEM>& a_k,
           MBLevelBoxData<T,NUMCOMPS,MEM>& a_UStage,
           MBBoundaryRegister<T,NUMREG,MEM>& a_blockreg,
           Array<Reduction<T,Operation::Sum,MEM>,NUMCOMPS>& a_consSum,
           T a_time, 
           T a_dt,
           T a_weight,
           int a_kstage)
{
  PR_TIMERS("BoxOp_Euler::StageLevel");
  auto layout = a_UStage.layout();
  auto map = CubedSphereShell::Map(layout, ghost());
  auto eulerOp = CubedSphereShell::Operator<BoxOp_EulerCubedSphere, double, HOST>(map);
  T dt_min = 1.0e20;
  T dt_next_local = 1.e20;
  T wgtxdt = a_weight*a_dt;
  for (auto dit : a_UStage.layout())
  {
    PROTO_ASSERT(a_k[dit].box() == a_k.layout()[dit],"k_i isn't exactly the size of a rhs");
    auto &k_i = a_k[dit];
    auto &USph_i = a_UStage[dit];
    Array<BoxData<double, NUMCOMPS>, DIM> fluxes;
    fluxes[0].define(k_i.box().extrude(0));
    fluxes[1].define(k_i.box().extrude(1));
    fluxes[2].define(k_i.box().extrude(2));    
    fluxes[0].setToZero();
    fluxes[1].setToZero();
    fluxes[2].setToZero();
    
    eulerOp[dit](k_i, fluxes, USph_i, 1.0);
    auto dx = eulerOp[dit].dx();
    unsigned int block = layout.block(dit);
    Box blockBox = layout.getBlock(block).domain().box();
    //auto mbdata = a_blockreg.bounds(dit);
    a_blockreg.increment(fluxes,dit,wgtxdt);
    CubedSphereShell::radialCons<T,NUMCOMPS,MEM>(fluxes[0],a_consSum,dx[0],blockBox,wgtxdt);
  }
}
template<
  typename T,
  template<MemType> typename MAP,
  MemType MEM>
Array<T,NUMCOMPS> 
BoxOp_EulerCubedSphere<T,MAP,MEM>::
  applyOp(MBLevelBoxData<T,NUMCOMPS,MEM>& a_JU,
          MBInterpOp& a_iop,
          MBLevelBoxData<T,1,MEM>& a_dVolrLev) 
{
  T weight = 1.0;
  auto layout = a_JU.layout();
  MBBoundaryRegister<T,BoxOp_EulerCubedSphere<T,MAP,MEM>::numState(),MEM>
    blockreg(layout,1,Point::Zeros());
  blockreg.clear();
  Array<Reduction<T,Operation::Sum,MEM>,BoxOp_EulerCubedSphere<T,MAP,MEM>::numState()> consSum, consRhs;
  MBLevelBoxData<T,NUMCOMPS,MEM> JUTemp(layout,a_JU.ghost());
  MBLevelBoxData<T,NUMCOMPS,MEM> rhs(layout,Point::Zeros());
  T time = 0.;
  T dt = 1.0;
  int kstage = 0;
  a_JU.copyTo(JUTemp);
  kstage = 0;
  int opswitchiv = ParseInputs::get_operator_switchiv();
  bool flagi = (opswitchiv%4) != 1;
  BoxOp_EulerCubedSphere<T,MAP,MEM>::
    preStageLevel(JUTemp, rhs, a_dVolrLev, a_iop, blockreg, consSum,
                  time,dt, weight, 0);
  if (procID() == 0) cout << opswitchiv << ": flagi = " << flagi << endl;  
  if (flagi) 
  BoxOp_EulerCubedSphere<T,MAP,MEM>::
    StageLevel(rhs, JUTemp, blockreg, consSum,
               time, dt, weight, kstage);
  blockreg.reflux(rhs);
  
  for (int comp = 0; comp < NUMCOMPS; comp++)
    {
      consRhs[comp].reset();
      for (auto dit : layout)
        {
          auto rhscomp = slice(rhs[dit],comp);
          rhscomp.reduce(consRhs[comp]);        
        }
    }
  Array<T,NUMCOMPS> outcons;
  for (int comp = 0; comp < NUMCOMPS; comp++)
    {
      outcons[comp] = 0.;
      outcons[comp] += consRhs[comp].fetch();
      outcons[comp] -= consSum[comp].fetch();
    }
  rhs.copyTo(a_JU);
  return outcons;
}
template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::dtInv(Reduction<T,Operation::Max,MEM>& a_dtInv,
                                         const State& a_W) const
{
  Box bx= this->box();  
  T dxiR = this->dx()[0];
  T dxiPerp = this->dx()[1];
  int block = this->block();
  Box bxRadius = bx;
  BoxData<T,DIM,MEM> Dr(bxRadius);
  BoxData<T,1,MEM> radius(bxRadius);
  BoxData<T,DIM,MEM> adjDr(bxRadius);
  BoxData<T,1,MEM> dVolr(bxRadius);
  BoxData<T,1,MEM> speedMax(bxRadius);
  speedMax.setToZero();
  radialMetrics(radius,Dr,adjDr,dVolr,bxRadius);
  T half = .5;
  for (int dir = 0 ; dir < DIM; dir++)
    {
      Box bxface = bx;
      BoxData<T,1,MEM> adjDrDir = slice(adjDr,dir); 
      BoxData<T,1,MEM> deta =
        CubedSphereShell::detA<T,MEM>(bxface,dxiPerp,dir);
      BoxData<T,DIM,MEM> ntrow = CubedSphereShell::adjA<T,MEM>(bxface,dxiPerp,block,dir);
      Array<Array<uint,DIM>,6> permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
      Array<Array<int,DIM>,6> sign = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}}; 
      Point high = bxface.high();
      Point low = bxface.low();
      high[0] = low[0];
      Box bx0(low,high);
      BoxData<T,DIM,MEM,DIM> A_matrix(bx);
      T offseta = half;
      T offsetb = half;
      if (dir == 1) offseta = 0.0;
      if (dir == 2) offsetb = 0.0;
      forallInPlace_p(f_Amatrix,bx0,A_matrix,permute[block],sign[block],
                      dxiPerp,offseta,offsetb);
      spreadSlice(A_matrix);
      forallInPlace(f_speedmax,
                    speedMax,a_W,A_matrix,
                    ntrow,adjDrDir,deta,
                    dVolr,Dr,dxiPerp,
                    dxiR,dir,m_gamma);
    }
  speedMax.reduce(a_dtInv);
}
#if 0
template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
// Use case: increment boundary registers. 
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::
PostStagePatch(
               std::vector<MBBoundaryData<T,NUMCOMPS,MEM> >& a_mbdata,
               Array<BoxData<T,NUMCOMPS,MEM>,DIM>& a_fluxes, 
               T a_time,
               T a_dt,
               T a_weight)
{
  // Increment boundary registers.
  PR_TIMERS("BoxOp_Euler::postStagePatch");
  abort(); // Deprecated. Can resurrect if we need it.
}
#endif
// Impose flux matching, apply artificial viscosity, divergence cleaning.
template<
  typename T,
  template<MemType> typename MAP,
  MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::
postStepLevel(MBLevelBoxData<T,NUMCOMPS,MEM>& a_rhsTotal,             
              MBBoundaryRegister<T,NUMCOMPS,MEM,PR_CELL>& a_blockReg)
{
  abort();
}
#if 0
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

// ====================================================================================
// Pre and Post implems end here.
// ====================================================================================


