template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::DivBPowell(
                        State& a_Rhs,
                        const State& a_USph,
                        T a_scale) const
{
  int a_block = this->block();
  State rhs(a_Rhs.box());
  rhs.setVal(0.0);
  PR_TIME("BoxOp_Euler::DivBPowell()");
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
  
    BoxData<T,1,MEM> adjDrDir = slice(adjDr,dir);  

    //Required mapping variables.
    Box bxface = WFace.box();
    BoxData<T,1,MEM> deta = CubedSphereShell::detA<T,MEM>(bxface,dxiPerp,dir); 
    Box bx = W.box();
    Box bxdir = bx.extrude(Point::Basis(dir));
    

    State fluxdir(bxface);
    forallInPlace([ ] PROTO_LAMBDA
      (Var<T,NUMCOMPS,MEM>& a_fluxout,
        const Var<T,NUMCOMPS,MEM>& a_Wsph,
        const Var<T,1,MEM>& a_deta,
        int a_dir)
      {  
        for (int comp = 0; comp < NUMCOMPS;comp++)
        {                  
          a_fluxout(comp) = a_deta(0)*a_Wsph(iBX + a_dir);
        }
      },fluxdir,WFace,deta,dir);
  
    for (int comp = 0; comp < NUMCOMPS;comp++)
    {
      BoxData<T,1,MEM> advfluxcomp = slice(fluxdir,comp);
      BoxData<T,1,MEM> advfluxprod = faceProduct(advfluxcomp,adjDrDir,dir);
      advfluxprod.copyTo(advfluxcomp);
    }

    // update RHS.
    BoxData<T,NUMCOMPS,MEM> rhsIncr;
    T dxInv = -1.0/(this->dx()[dir]);
    rhsIncr = Stencil<T>::FluxDivergence(dir)(fluxdir,dxInv);

    rhs += rhsIncr;
    outputbox &= rhsIncr.box();   
  }
  if (a_block == BLOCKTEST)
  {
    HDF5Handler h5;
    BoxData<T,NUMCOMPS,MEM> rhsRadial = checkRadial(rhs,dxiPerp,a_block);
    h5.writePatch(dxiPerp,rhsRadial,"rhsDivBRadial");
  }
  rhs *= (a_scale);

  // We have calculated -div(B)*dt in rhs. Need to multiply other Powell terms to get the final RHS.
  Array<Array<uint,DIM>,6> permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
  Array<Array<int,DIM>,6> sign = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}}; 
  Point high = bx.high();
  Point low = bx.low();
  high[0] = low[0];
  Box bx0(low,high);
  BoxData<T,DIM,MEM,DIM> A_matrix(bx);
  T offseta = half;
  T offsetb = half;
  forallInPlace_p(f_Amatrix,bx0,A_matrix,permute[a_block],sign[a_block],
                  dxiPerp,offseta,offsetb);
  spreadSlice(A_matrix);

  BoxData<T,DIM,MEM> Vsph = slice<T,NUMCOMPS,DIM,MEM>(W,iVX);
  BoxData<T,DIM,MEM> Vcart = Operator::_cellMatrixProductAB(A_matrix,Vsph,A_matrix,Vsph);

  BoxData<T,DIM,MEM> Bsph = slice<T,NUMCOMPS,DIM,MEM>(W,iBX);
  BoxData<T,DIM,MEM> Bcart = Operator::_cellMatrixProductAB(A_matrix,Bsph,A_matrix,Bsph);

  forallInPlace([ ] PROTO_LAMBDA
      (Var<T,NUMCOMPS,MEM>& a_rhs,
        const Var<T,DIM,MEM>& a_Vcart,
        const Var<T,DIM,MEM>& a_Bcart)
      {  
        T VdotB = a_Vcart(0)*a_Bcart(0) + a_Vcart(1)*a_Bcart(1) + a_Vcart(2)*a_Bcart(2);                
        a_rhs(iRHO) *= 0;
        a_rhs(iMOMX) *= a_Bcart(0)/4.0/M_PI;
        a_rhs(iMOMY) *= a_Bcart(1)/4.0/M_PI;
        a_rhs(iMOMZ) *= a_Bcart(2)/4.0/M_PI;
        a_rhs(iE) *= VdotB/4.0/M_PI;
        a_rhs(iBX) *= a_Vcart(0);
        a_rhs(iBY) *= a_Vcart(1);
        a_rhs(iBZ) *= a_Vcart(2);
        
      },rhs,Vcart,Bcart);

  a_Rhs += rhs;
}