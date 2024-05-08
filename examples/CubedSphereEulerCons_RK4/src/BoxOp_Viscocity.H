template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::LinearVisc(
                        State& a_Rhs,
                        Array<BoxData<double, NUMCOMPS>, DIM>& a_fluxes,
                        const State& a_USph,
                        T a_scale) const
{
  int a_block = this->block();
  State rhs(a_Rhs.box());
  rhs.setVal(0.0);
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
      Stencil<double> D2 = 1.0*Shift(Point::Basis(dir)) + (-2.0)*Shift(Point::Zeros())
       + 1.0*Shift(-Point::Basis(dir));
      Stencil<double> D1 = 1.0*Shift(Point::Zeros()) + (-1.0)*Shift(-Point::Basis(dir));
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
      T dxInv = -1.0/(this->dx()[dir]);
      fluxdir.copyTo(a_fluxes[dir]);
      rhsIncr = Stencil<T>::FluxDivergence(dir)(fluxdir,dxInv);

      rhs += rhsIncr;
      outputbox &= rhsIncr.box();   
    }
  if (a_block == BLOCKTEST)
    {
      HDF5Handler h5;
      BoxData<T,NUMCOMPS,MEM> rhsRadial = checkRadial(rhs,dxiPerp,a_block);
      h5.writePatch(dxiPerp,rhsRadial,"rhsViscRadial");
    }
  rhs *= (a_scale);
  a_Rhs += rhs;

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