template<
  typename T,
  template<MemType> typename MAP,
  MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::
viscFlux(
         State& a_Rhs,
         Array<BoxData<double, NUMCOMPS>, DIM>& a_fluxes,
         const State& a_USph,
         T a_scale,
         T a_dt) const
{
  int block = this->block();
  State rhs(a_Rhs.box());
  rhs.setVal(0.0);
  PR_TIME("BoxOp_Euler::viscosity terms()");
  HDF5Handler h5;
  Box bxRad = a_USph.box().grow(Point::Basis(0,2));
  BoxData<T,DIM,MEM> Dr(bxRad);
  BoxData<T,1,MEM> radius(bxRad);
  BoxData<T,DIM,MEM> adjDr(bxRad);
  BoxData<T,1,MEM> dVolr(bxRad);
  radialMetrics(radius,Dr,adjDr,dVolr,bxRad);
  auto domainBox = this->layout().domain().getBlock(block).box();
  auto dxi = this->dx();
  // Transform to primitive variables, spherical velocities.
  Box bx = a_Rhs.box();
  State W;
  State WBar; 
  sphToPrim(W,WBar,a_USph,block);
  
  
  // Compute Fluxes, rhs.  
  Box outputbox = bx;
  // Loop over directions.
  T one = 1.0;
  T half = 0.5;
  T dxiPerp = this->dx()[1];
  
  // Setup for nonlinear viscosity.
  // compute and store volumetric flow rates.
  int opswitchiv = ParseInputs::get_operator_switchiv();
  int opswitchln = ParseInputs::get_operator_switchln();
  bool viscl = (opswitchln % 4 != 1);
  bool viscnl = (opswitchln % 4 != 2);
  //if (procID() == 0) cout << opswitchln << ": viscl = " << viscl << endl;
  //if (procID() == 0) cout << opswitchln << ": viscnl = " << viscnl << endl;
  
  for (int dir = 0; dir < DIM; dir++)
    {
      PR_TIME("BoxOp_Euler::viscous operator");
      // high order linear viscosity.
      State WFace = Stencil<T>::CellToFace(dir)(W);
      BoxData<T,1,MEM> adjDrDir = slice(adjDr,dir);  
      Stencil<double> D2 = 1.0*Shift(Point::Basis(dir)) + (-2.0)*Shift(Point::Zeros())
        + 1.0*Shift(Point::Basis(dir,-1));
      Stencil<double> D1 = (1.0)*Shift(Point::Zeros()) + (-1.0)*Shift(Point::Basis(dir,-1));
      Stencil<double> D5 =  D1 * D2 * D2;
      Box bxface = a_Rhs.box().extrude(Point::Basis(dir));
      State physFlux(bxface);
      physFlux.setToZero();
      State fluxdir(bxface);
      fluxdir.setToZero();
      
      //Required mapping variables.
      BoxData<T,DIM,MEM,DIM> deta_a = CubedSphereShell::detA_A<T,MEM>(bxface,dxiPerp,block,dir);
      BoxData<T,1,MEM> deta = CubedSphereShell::detA<T,MEM>(bxface,dxiPerp,dir); 
      Box bx = W.box();
      Box bxdir = bx.extrude(Point::Basis(dir));
      BoxData<T,DIM,MEM> ntrow = CubedSphereShell::adjA<T,MEM>(bx,dxiPerp,block,dir);
#if 0
      // Debugging diagnostics, #if'ed out for now.
      if (dir == 0)
        {
          auto vr = slice(W,iVX);
          auto adjDrslice = slice(adjDr,0);
          BoxData<T,1,MEM> adjdrcomp =
            ((0.5)*Shift(Point::Basis(0)) + (0.5)*Shift(Point::Zeros()))(adjDrslice);
          auto temp = vr*adjdrcomp;
          //temp *= deta;
          //BoxData<T,1,MEM> dvdr = Stencil<T>::Derivative(1,0,2)(temp,a_Rhs.box(),1./dxi[0]);
          BoxData<T,1,MEM> dvdr =
            ((-1.0)*Shift(Point::Basis(0,-1)) + 1.0*Shift(Point::Basis(0)))
            (temp,bxface,0.5/dxi[0]);
          // cout << procID() << ", dAvrdxi = " << dvdr.absMax() << endl;
          forallInPlace([ ] PROTO_LAMBDA(
                                         Var<T,1,MEM>& a_numer,
                                         const Var<T,1,MEM>& a_denom1,
                                         const Var<T,1,MEM>& a_denom2)
                        {
                          a_numer(0) = a_numer(0)/(a_denom1(0));
                        },dvdr,dVolr,deta);
  
          // cout << procID() << ", dAvrdVolr = " << dvdr.absMax() << endl;
        }
#endif // endif for #if 0.
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
      forallInPlace_p(f_Amatrix,bx0,A_matrix,permute[block],sign[block],
                      dxiPerp,offseta,offsetb);
      spreadSlice(A_matrix);
      BoxData<T,NUMCOMPS,MEM> linearViscFlux;
      BoxData<T,NUMCOMPS,MEM> NLViscFlux;
      BoxData<T,NUMCOMPS,MEM> rhsIncr(bxface);
      if (viscl)
        {
          // cout << "viscl = " << viscl << endl;
          State fluxtemp(bxface);
          fluxtemp.setToZero();
          linearViscFlux = D5(a_USph, 1.0/128);
          primToMappedViscFlux(fluxtemp,linearViscFlux,WFace,adjDrDir,
                               deta_a,deta,A_matrix,ntrow,block,dir);
          BoxData<T,1,MEM> linvisccoef(bxface);
          LinearViscCoef(linvisccoef,WFace,deta,A_matrix,ntrow,block,dir);
          // cout << procID() << ": max lin visc coef = " << linvisccoef.absMax() << endl;
          for (int comp = 0 ; comp < NUMCOMPS; comp++)
            {
              auto fluxtempcomp = slice(fluxtemp,comp);
              auto fluxcomp = slice(fluxdir,comp);
              fluxtempcomp *= linvisccoef;
              fluxcomp += fluxtempcomp;
            }
        }
      if (viscnl) 
        {
          // cout << "viscnl = " << viscnl << endl;
          State fluxtemp(bxface);
          fluxtemp.setToZero();
          NLViscFlux = D1(a_USph);
          primToMappedViscFlux(fluxtemp,NLViscFlux,WFace,adjDrDir,
                               deta_a,deta,A_matrix,ntrow,block,dir);
          BoxData<T,1,MEM> nlvisccoef(bxface);
          NonLinearViscCoef(nlvisccoef,W,adjDr,A_matrix,ntrow,dVolr,deta,dxi,block,dir);
          // if (procID() == 0) cout << "dir =" << dir << ", max nlvisc coef = " << nlvisccoef.absMax() << endl;
          for (int comp = 0 ; comp < NUMCOMPS; comp++)
            {
              auto fluxtempcomp = slice(fluxtemp,comp);
              auto fluxcomp = slice(fluxdir,comp);
              fluxtempcomp *= nlvisccoef;
              fluxcomp += fluxtempcomp;
            }
        }
      // set fluxes, update RHS.
     
      fluxdir.copyTo(a_fluxes[dir]);
      int zerowidth = 1;
      // DEBUG?? Switch off artificial viscosity on radial boundaries.
      if (dir == 0)
        {
          auto layout = this->layout();
          unsigned int block = this->block();
          Box blockBox = layout.getBlock(block).domain().box();
          Box highBoundary0 = blockBox.face(0,Side::Hi,zerowidth).shift(0,1);
          Box lowBoundary0 = blockBox.face(0,Side::Lo,zerowidth);
          if (!((highBoundary0 & a_fluxes[0].box()).empty()))
            {
              a_fluxes[0].setVal(0.,a_fluxes[0].box().face(0,Side::Hi,zerowidth));
            }
          if (!((lowBoundary0 & a_fluxes[0].box()).empty()))
            {
              a_fluxes[0].setVal(0.,a_fluxes[0].box().face(0,Side::Lo,zerowidth));
            }
        }
    
      T dxInv = 1.0/dxi[dir];
      rhsIncr = Stencil<T>::FluxDivergence(dir)(a_fluxes[dir],dxInv);
      rhs += rhsIncr;
      
      outputbox &= rhsIncr.box();
#if 0
      // << procID() << ", " << outputbox << ", " << rhsIncr.absMax(0) <<  endl;
      if (procID() == 0)
        {
           T c_MP = 67262192e-24;
           Array<T,NUMCOMPS> scalevars =
             {c_MP,1.e5*c_MP,1.e5*c_MP,1.e5*c_MP,1.e12,1.e6,1.e6,1.e6};
          auto J = dVolr*deta;
          rhsIncr *= a_dt;  
          for (int comp; comp < NUMCOMPS; comp++)
            {
              auto rhsincrcomp = slice(rhsIncr,comp);
              forallInPlace([ ] PROTO_LAMBDA(Var<T,1,MEM>& a_rhs,
                                             const Var<T,1,MEM>& a_j,
                                             T a_scalevar)
                            {
                              a_rhs(0) /= (a_j(0)*a_scalevar);
                            },rhsincrcomp,J,scalevars[comp]);
              
            }
          if (procID() == 0)
            {
              cout << "dt = " << a_dt << ", rhsViscmax = " << rhsIncr.absMax(0) << endl;
            }
          HDF5Handler h5;
          h5.writePatch(dxiPerp,rhsIncr,"rhsVisc"+to_string(dir));
        }
#endif
    }
  rhs *= (a_scale);
  a_Rhs += rhs;
}
template<
  typename T,
  template<MemType> typename MAP,
  MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::
primToMappedViscFlux(
                     State& a_fluxdir,
                     const State& a_physFlux,
                     const State& a_WFace,
                     const BoxData<T,1,MEM>& a_adjDdir,
                     const BoxData<T,DIM,MEM,DIM>& a_deta_a,
                     const BoxData<T,1,MEM>& a_deta,
                     const BoxData<T,DIM,MEM,DIM>& a_A,
                     const BoxData<T,DIM,MEM>& a_ntrow,
                     uint a_block,
                     uint a_dir) const
{
  T one = 1.0;
  // Multiply momentum and B fluxes by det(A) A.
  BoxData<T,DIM,MEM> fluxdirMom = slice<T,NUMCOMPS,DIM,MEM>(a_fluxdir,iMOMX);
  BoxData<T,DIM,MEM> physFluxMom = slice<T,NUMCOMPS,DIM,MEM>(a_physFlux,iMOMX);
  BoxData<T,DIM,MEM> tmp = forall<T,DIM,MEM>(f_matVecProd,a_deta_a,physFluxMom,one);
  //Operator::_faceMatrixProductAB(a_deta_a,physFluxMom,a_deta_a,physFluxMom,a_dir);
  //CHANGE to simple product.
  tmp.copyTo(fluxdirMom);

  BoxData<T,DIM,MEM> fluxdirB = slice<T,NUMCOMPS,DIM,MEM>(a_fluxdir,iBX);
  BoxData<T,DIM,MEM> physFluxB = slice<T,NUMCOMPS,DIM,MEM>(a_physFlux,iBX);
  BoxData<T,DIM,MEM> tmp2 = forall<T,DIM,MEM>(f_matVecProd,a_deta_a,physFluxB,one);
  // Operator::_faceMatrixProductAB(a_deta_a,physFluxB,a_deta_a,physFluxB,a_dir);
  // CHANGE to simple product.
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
      BoxData<T,1,MEM> advfluxprod = advfluxcomp*a_adjDdir;
      //BoxData<T,1,MEM> advfluxprod = faceProduct(advfluxcomp,a_adjDdir,a_dir);
      advfluxprod.copyTo(advfluxcomp);
    }        
}
template<
  typename T,
  template<MemType> typename MAP,
  MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::LinearViscCoef(
                                                  BoxData<T,1,MEM>& a_LNViscCoef,
                                                  const State& a_WFace,
                                                  const BoxData<T,1,MEM>& a_deta,
                                                  const BoxData<T,DIM,MEM,DIM>& a_A,
                                                  const BoxData<T,DIM,MEM>& a_ntrow,
                                                  uint a_block,
                                                  uint a_dir) const
{
  // Find (|w_d| + c_fast*Adj(A)_{d,.}) and multiply with a_fluxdir.
  forallInPlace([ ] PROTO_LAMBDA
                (Var<T,1,MEM>& a_linfluxcoef,
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
                  
                  a_linfluxcoef(0) = (abs(w_d) + af*A_d_mag);
                },a_LNViscCoef,a_WFace,a_ntrow,a_A,m_gamma,a_dir);
}
template<
  typename T,
  template<MemType> typename MAP,
  MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::
NonLinearViscCoef(
                  BoxData<T,1,MEM>& a_NLViscCoef,
                  const State& a_W,
                  const BoxData<T,DIM,MEM>& a_adjDr,
                  const BoxData<T,DIM,MEM,DIM>& a_A,
                  const BoxData<T,DIM,MEM>& a_ntrow,
                  const BoxData<T,1,MEM>& a_dVolr,
                  const BoxData<T,1,MEM>& a_deta,
                  Array<T,DIM> a_dxi,
                  uint a_block,
                  uint a_dir) const
{
  // Nonlinear artificial viscosity.
  PR_TIMERS("BoxOp_EulerCubedSphere::NL Viscosity coefficient");
  HDF5Handler h5;
  // Compute divergence of the velocity.
  BoxData<T,DIM,MEM> volflow(a_W.box());
  BoxData<T,DIM,MEM> vel = slice<T,NUMCOMPS,DIM,MEM>(a_W,iVX);
  Array<BoxData<T,1,MEM>,DIM> dvoldn;
  BoxData<T,1,MEM> deta0 = CubedSphereShell::detA<T,MEM>(a_W.box(),a_dxi[1],0);
  for (int dir = 0; dir < DIM; dir++)
    {
      BoxData<T,1,MEM> veldir = slice(vel,dir);
      BoxData<T,1,MEM> volflowdir(veldir.box());
      veldir.copyTo(volflowdir);
      auto adjDrdir = slice(a_adjDr,dir);
      BoxData<T,1,MEM> adjDrdircen(a_adjDr.box());
      if (dir == 0)
        {
          adjDrdircen = (0.5*Shift(Point::Zeros()) + 0.5*Shift(Point::Basis(0)))(adjDrdir);
        }
      else
        {
         adjDrdir.copyTo(adjDrdircen);
        }
      T dxidir = a_dxi[dir];
      T scale = 1.0/dxidir;
      volflowdir *= adjDrdircen;
      volflowdir *= deta0;
      BoxData<T,1,MEM> temp = (1.0*Shift(Point::Zeros()) - 1.0*Shift(Point::Basis(0,-1)))
        (volflowdir,scale);
      dvoldn[dir].define(temp.box());
      temp.copyTo(dvoldn[dir]);
    }
  int dirnorm = a_dir;
  int dirtan1 = (a_dir + 1)%DIM;
  int dirtan2 = (a_dir + 2)%DIM;
  Stencil<T> faceToCorner1 = 
    (.5*Shift(Point::Zeros()) + .5*Shift(Point::Basis(dirnorm,-1)))*
    (.5*Shift(Point::Zeros()) + .5*Shift(Point::Basis(dirtan1,1)));
  Stencil<T> faceToCorner2 = 
    (.5*Shift(Point::Zeros()) + .5*Shift(Point::Basis(dirnorm,-1)))*
    (.5*Shift(Point::Zeros()) + .5*Shift(Point::Basis(dirtan2,1)));
  BoxData<T,1,MEM> dvft1dtan1 = faceToCorner1(dvoldn[dirtan1]);
  BoxData<T,1,MEM> dvft1dtan2 = faceToCorner2(dvoldn[dirtan2]);
  BoxData<T,1,MEM> divu =
    forall<T,1,MEM>([] PROTO_LAMBDA(
                                    Var<T,1,MEM>&            a_divu,
                                    Var<T,1,MEM>&            a_dvolfndn,
                                    Var<T,1,MEM>&            a_dvolft1dt1,
                                    Var<T,1,MEM>&            a_dvolft2dt2,
                                    const Var<T,1,MEM>&      a_deta0,
                                    const Var<T,1,MEM>&      a_dvolr)
                    {
                      //a_dvolfndn(0) /= a_deta0(0)*a_dvolr(0);
                      //a_dvolft1dt1(0) /= a_deta0(0)*a_dvolr(0);
                      //a_dvolft2dt2(0) /= a_deta0(0)*a_dvolr(0);
                      a_divu(0) =
                      (a_dvolfndn(0) + a_dvolft1dt1(0) + a_dvolft2dt2(0))
                      /(a_deta0(0)*a_dvolr(0));
                    },dvoldn[a_dir],dvft1dtan1,dvft1dtan2,deta0,a_dVolr);
  
  int nrad = int (1.0/a_dxi[0]);
  #if 0
  if (procID() == 0)
    {
      //h5.writePatch(a_dxi[0],divu,
      //             "divu_"+to_string(a_dir) + to_string(nrad) +
      //             "patch_"+to_string(procID()));
      Box bxface = a_NLViscCoef.box();
      auto adjDrdir = slice(a_adjDr,a_dir);
      BoxData<T,1,MEM> dummy(bxface),dummy0(bxface),dummy1(bxface),
        dummy2(bxface),dummy3(bxface),dummy4(bxface),dummy5(bxface);
      divu.copyTo(dummy);
      dvoldn[a_dir].copyTo(dummy0);
      //dvoldn[dirtan1].copyTo(dummy1);
      //dvoldn[dirtan2].copyTo(dummy2);
      dvft1dtan1.copyTo(dummy1);
      dvft1dtan2.copyTo(dummy2);
      a_dVolr.copyTo(dummy3);
      deta0.copyTo(dummy4);
      adjDrdir.copyTo(dummy5);
      
      // cout << procID() << ",max divu = " << dummy.absMax() <<
      //", dxi[dir] = " << a_dxi[a_dir] <<
        // ", dvolf1dtan1 = " << dummy1.absMax() <<
        //", dvolf1dtan2 = " << dummy2.absMax() <<
        // ", a_dVolr = " << dummy3.absMax() <<
        //", deta0 = " << dummy4.absMax() <<
        //", adjDrdir = " << dummy5.absMax() <<
        //endl;
    }
#endif
  {
    PR_TIMERS("BoxOp_EulerCubedSphere::NL Viscosity coefficient - final");
    divu *= a_dxi[a_dir];
    T alpha = .3;
    T beta = .3;
    BoxData<T,1,MEM> csq = forall<T,1,MEM>(f_csq,a_W,a_A,m_gamma);

    //cout << procID() << ", lambda max = "
    //    << divu.absMax() << ", csq max = " << csq.absMax() << endl;
    BoxData<T,1,MEM> csqshiftr = alias(csq,Point::Basis(a_dir));
    BoxData<T,1,MEM> csqshiftr2 = alias(csq,Point::Basis(a_dir,2));
    BoxData<T,1,MEM> csqshiftl = alias(csq,Point::Basis(a_dir,-1));
    forallInPlace([] PROTO_LAMBDA(
                                  Var<T,1,MEM>&            a_coef,
                                  const Var<T,1,MEM>&      a_divu,
                                  const Var<T,1,MEM>&      a_csql,
                                  const Var<T,1,MEM>&      a_csqr,
                                  const Var<T,1,MEM>&      a_csql2,
                                  const Var<T,1,MEM>&      a_csqr2,
                                  T                        a_dx,
                                  T                        a_alpha,
                                  T                        a_beta)
                  {
                    T lambda = abs(a_divu(0));
                    a_coef(0) = a_alpha*lambda;
                    min(1.,lambda*lambda/
                        (min(a_csql(0),min(a_csqr(0),
                                           min(a_csql2(0),a_csqr2(0))))
                         *a_beta));
                  },a_NLViscCoef,divu,
                  csq,csqshiftr,csqshiftl,csqshiftr2,
                  a_dxi[a_dir],alpha,beta);
#if 0    
    if (procID() == 0)
    {
      h5.writePatch(a_dxi[0],a_NLViscCoef,
                   "divu_"+to_string(a_dir) + to_string(nrad) +
                   "patch_"+to_string(procID()));
    }
#endif
  }
}
