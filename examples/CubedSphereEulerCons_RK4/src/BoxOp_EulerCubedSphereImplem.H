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
   T&                              a_dt_next,
   const State&                    a_USph,
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
  a_Rhs.setVal(0.0);
  Box outputbox = bx;
  // Loop over directions.
  T one = 1.0;
  T half = 0.5;
  T dxiR = this->dx()[0];
  T dxiPerp = this->dx()[1];  
  BoxData<T,1,MEM> deta0(bx);
  BoxData<T,DIM,MEM> ntrow0(bx);
  BoxData<T,DIM,MEM,DIM> A_matrix0(bx);
  BoxData<T,1,MEM> dt_min(this->box());
  dt_min.setVal(1.0e10);
  for (int dir = 0; dir < DIM; dir++)
    {
      PR_TIME("BoxOp_Euler::operator()::directional flux calculation");        
      // State WFace = Stencil<T>::CellToFace(dir)(W);
      // State WFaceBar = Stencil<T>::CellToFace(dir,Side::Lo,2)(WBar);      

      State WFace_L = Stencil<T>::CellToFace(dir)(W);
      State WFaceBar_L = Stencil<T>::CellToFace(dir,Side::Lo,2)(WBar);    

      State WFace_R = Stencil<T>::CellToFace(dir)(W);
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
 
      int limiter = 1;

      if (limiter == 1){
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
      T dxInv = -1.0/(this->dx()[dir]);
      rhsIncr = Stencil<T>::FluxDivergence(dir)(fluxdir,dxInv);
      a_Rhs += rhsIncr;
      outputbox &= rhsIncr.box();

      // Compute dt_min.
      forallInPlace([ ] PROTO_LAMBDA
        ( Var<T,1,MEM>& a_dt_min,
          const Var<T,NUMCOMPS,MEM>& a_W,
          const Var<T,DIM,MEM,DIM>& a_amat,
          const Var<T,DIM,MEM>& a_adjA,
          const Var<T,1,MEM>& a_adjDrDir,
          const Var<T,1,MEM>& a_detA,
          const Var<T,1,MEM>& a_dVolr,
          const Var<T,DIM,MEM>& a_Dr,
          const T a_dxiPerp,
          const T a_dxiR,
          const int a_d,
          const T a_gamma)
        { 
          double dx_d, dt_d;
          if (a_d > 0){
            T adjA_mag = sqrt(a_adjA(0)*a_adjA(0) + a_adjA(1)*a_adjA(1) + a_adjA(2)*a_adjA(2));
            dx_d = a_dVolr(0) * a_detA(0) * a_dxiPerp / a_adjDrDir(0) / adjA_mag;
          } 
          if (a_d == 0){
            dx_d = a_Dr(0) * a_dxiR;
          }   

          vecType velSph;
          vecType velCart;
          vecType bSph;
          vecType bCart;

          velSph[0] = a_W(iVX);       
          velSph[1] = a_W(iVY);
          velSph[2] = a_W(iVZ);
          velCart[0] = 0.;       
          velCart[1] = 0.;
          velCart[2] = 0.;
      
          bSph[0] = a_W(iBX);       
          bSph[1] = a_W(iBY);
          bSph[2] = a_W(iBZ);
          bCart[0] = 0.;       
          bCart[1] = 0.;
          bCart[2] = 0.;
          T half = .5;
          T one = 1.0;
          T v2 = 0.;
          T b2 = 0.;
          for (int dir0 = 0; dir0 < DIM; dir0++)
          {
            for (int dir1 = 0; dir1 < DIM; dir1++)
            {
              velCart[dir0] += a_amat(dir0,dir1)*velSph[dir1];
              bCart[dir0] += a_amat(dir0,dir1)*bSph[dir1];
            }
            v2 += velCart[dir0]*velCart[dir0];
            b2 += bCart[dir0]*bCart[dir0];
          }
          T rho = a_W(iRHO);
          T B_mag = sqrt(b2);
          T p = a_W(iP);
          if (p < 0.0) p = 0.0;
          T af = sqrt(a_gamma*p/rho + b2/4.0/M_PI/rho);  
          a_dt_min(0) = min(a_dt_min(0),dx_d/(af + sqrt(v2)));
        },dt_min,W,A_matrix0,ntrow0,adjDrDir,deta0,dVolr,Dr,dxiPerp,dxiR,dir,m_gamma);
      
    }
  a_dt_next = dt_min.min();
  if (a_block == BLOCKTEST)
    {
      HDF5Handler h5;
      BoxData<T,NUMCOMPS,MEM> rhsRadial = checkRadial(a_Rhs,dxiPerp,a_block);
      h5.writePatch(dxiPerp,rhsRadial,"rhsRadial");
    }
  a_Rhs *= (a_scale);
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
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::
StageLevel(MBLevelBoxData<T,NUMCOMPS,MEM>& a_rhs, 
           T& a_dt_next,
           MBLevelBoxData<T,NUMCOMPS,MEM>& a_UStage,
           MBBoundaryRegister<T,NUMREG,MEM>& a_blockreg,
           T a_time, 
           T a_dt,
           T a_weight,
           int a_kstage)
{
  auto map = CubedSphereShell::Map(a_UStage.layout(), ghost());
  auto eulerOp = CubedSphereShell::Operator<BoxOp_EulerCubedSphere, double, HOST>(map);
  T dt_min = 1.0e10;
  T wgtxdt = a_weight*a_dt;
  for (auto dit : a_UStage.layout())
  {
    auto &rhs_i = a_rhs[dit];
    auto &USph_i = a_UStage[dit];
    Array<BoxData<double, NUMCOMPS>, DIM> fluxes;
    fluxes[0].define(rhs_i.box().extrude(0));
    fluxes[1].define(rhs_i.box().extrude(1));
    fluxes[2].define(rhs_i.box().extrude(2));    
    fluxes[0].setToZero();
    fluxes[1].setToZero();
    fluxes[2].setToZero();
    
    eulerOp[dit](rhs_i, fluxes, a_dt_next, USph_i, a_dt);
    //cout << procID()  << ", " << fluxes[1].absMax(0) << endl;
    //cout << procID()  << ", " << fluxes[2].absMax(0) << endl;
    auto mbdata = a_blockreg.bounds(dit);
    //a_blockreg.increment(fluxes,dit,wgtxdt);
    if (a_dt_next < dt_min && a_kstage == 0) dt_min = a_dt_next;
  }
  if (a_kstage == 0) a_dt_next = dt_min;

}
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
  abort();
}

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


