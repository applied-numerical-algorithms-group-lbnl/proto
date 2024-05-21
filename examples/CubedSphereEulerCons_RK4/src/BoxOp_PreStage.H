#include "MHDReader.H"
#include "MHD_Constants.H"
extern MHDReader BC_global;
template <typename T, MemType MEM>
PROTO_KERNEL_START
void f_radialBCs_F(
    Point a_pt,
    Var<T, NUMCOMPS, MEM> &a_USph,
    Var<T> &a_radius,
    T a_dxradius,
    T a_gamma,
    int a_nradius)
{
  // Compute spherical BCs.
  T p0 = 1.0;
  T rho0 = 1.0;
  T eps = 0.1;
  T amplitude;
  T arg = (1.0 * a_pt[0] + .5 - 1.0 * a_nradius / 2) * a_dxradius / 2.;
  if (abs(arg) < .25)
  {
    amplitude = eps * pow(cos(2 * M_PI * arg), 6);
  }
  else
  {
    amplitude = 0.;
  }
  T rho = rho0 + amplitude * rho0;
  T p = p0 * pow(rho / rho0, a_gamma);
  T ur = amplitude * sqrt(a_gamma * p0 / rho0) / rho0;
  a_USph(iRHO) = rho;
  a_USph(iMOMX) = 0.0;//ur * rho;
  a_USph(iMOMY) = 0.0;
  a_USph(iMOMZ) = 0.0;
  a_USph(iBX) = 0.0;
  a_USph(iBY) = 0.0;
  a_USph(iBZ) = 0.0;
  T umag = sqrt(a_USph(iMOMX) * a_USph(iMOMX) + a_USph(iMOMY) * a_USph(iMOMY) + a_USph(iMOMZ) * a_USph(iMOMZ)) / a_USph(iRHO);
  T Bmag = sqrt(a_USph(iBX) * a_USph(iBX) + a_USph(iBY) * a_USph(iBY) + a_USph(iBZ) * a_USph(iBZ));
  a_USph(iE) = p / (a_gamma - 1.0) + 0.5 * a_USph(iRHO) * umag * umag + Bmag * Bmag/8.0/M_PI;
}
PROTO_KERNEL_END(f_radialBCs_F, f_radialBCs)


template <typename T, MemType MEM>
PROTO_KERNEL_START
void f_radialoutflowIBCs_F(
    Point a_pt,
    Var<T, NUMCOMPS, MEM> &a_USph,
    Var<T> &a_radius,
    Var<T, DIM> &a_X_cart,
    T a_gamma,
    int a_nradius)
{
  T r = sqrt(a_X_cart(0) * a_X_cart(0) + a_X_cart(1) * a_X_cart(1) + a_X_cart(2) * a_X_cart(2));
  T theta = acos(a_X_cart(2) / r);
  // Compute spherical BCs.
  double r_BC = CUBED_SPHERE_SHELL_R0; //in cm
  T p0 = 1.0e-7; // p near 21.5 c_SR is about 1e-7 dyne/cm2
  T rho0 = 700*c_MP; // rho at 21.5 c_SR is about 700/cm3
  T u0 = 500.0*1e5; // v at 21.5 c_SR is about 500 km/s
  T rho = rho0*pow(a_radius(0)/r_BC, -2.0);
  T p = p0 * pow(rho / rho0, a_gamma);
  T ur = u0;
  // T ur = u0 - u0/10 * sin(theta);
  a_USph(iRHO) = rho;
  a_USph(iMOMX) = ur * rho;
  a_USph(iMOMY) = 0.0;
  a_USph(iMOMZ) = 0.0;
  a_USph(iBX) = 0.0;
  a_USph(iBY) = 0.0;
  a_USph(iBZ) = 0.0;
  T umag = sqrt(a_USph(iMOMX) * a_USph(iMOMX) + a_USph(iMOMY) * a_USph(iMOMY) + a_USph(iMOMZ) * a_USph(iMOMZ)) / a_USph(iRHO);
  T Bmag = sqrt(a_USph(iBX) * a_USph(iBX) + a_USph(iBY) * a_USph(iBY) + a_USph(iBZ) * a_USph(iBZ));
  a_USph(iE) = p / (a_gamma - 1.0) + 0.5 * a_USph(iRHO) * umag * umag + Bmag * Bmag/8.0/M_PI;
}
PROTO_KERNEL_END(f_radialoutflowIBCs_F, f_radialoutflowIBCs)

// Radial boundary conditions here.
template<
  typename T,
  template<MemType> typename MAP,
  MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::
PreStagePatch(BoxData<T,NUMCOMPS,MEM>& a_USphStage,
              BoxData<T,NUMCOMPS,MEM>& a_BC,
              BoxData<T,1,MEM>& a_dVolr,
              Box a_blockBox,
              T a_time,
              T a_dt,
              int a_kstage)
{
  int a_block = this->block();
  Box bxGhosted = a_USphStage.box();
  BoxData<double> radius(a_dVolr.box());
  BoxData<double, DIM, HOST> Dr(a_dVolr.box());
  BoxData<double, DIM, HOST> adjDr(a_dVolr.box());
  radialMetrics(radius, Dr, adjDr, a_dVolr, Dr.box());
  Array<double, DIM> dx;
  dx = this->dx();
  double dxradius = dx[0];
  int thickness = std::round(1.0 / dxradius);
  // Set radial boundary condtions in radial ghost cells.
  int init_condition_type = ParseInputs::get_init_condition_type();
  // Outer BCs
  if (a_blockBox.high()[0] < bxGhosted.high()[0])
  {
    Point low = bxGhosted.low();
    low[0] = a_blockBox.high()[0] + 1;
    Box bdryBoxHigh(low, bxGhosted.high());
    if (init_condition_type == 0 || init_condition_type == 1){
      forallInPlace_p(f_radialBCs, bdryBoxHigh, a_USphStage, radius, dxradius, m_gamma, thickness);
    } 
    if (init_condition_type == 2 || init_condition_type == 3){
      // Copy last layer to ghost cells and scale by radius. // Extrapolation can make -ve P in ghost cells.
      Point source_low = Point(a_blockBox.high()[0],bxGhosted.low()[1],bxGhosted.low()[2]);
      Point source_high = Point(a_blockBox.high()[0],bxGhosted.high()[1],bxGhosted.high()[2]);
      Box sourceBox(source_low,source_high);
      State U_ghost(bdryBoxHigh);
      Scalar radius_ghost(bdryBoxHigh);
      Stencil<double> m_exterp_f_2nd;
      for (int i = 1; i <= NGHOST; i++ ) {
        a_USphStage.copyTo(U_ghost,sourceBox,Point::Basis(0)*(i));// Using shifting option of copyTo
        radius.copyTo(radius_ghost,sourceBox,Point::Basis(0)*(i));
      }
      U_ghost.copyTo(a_USphStage,bdryBoxHigh);
      //scale by radius.
      forallInPlace_p([ ] PROTO_LAMBDA
      (Point a_pt,
       Var<T, NUMCOMPS, MEM> &a_Usph,
       Var<T> &a_radius,
       Var<T> &a_radius_ghost)
      {  
        double rad = a_radius(0);
        double rad_last = a_radius_ghost(0);
        a_Usph(0) *= pow(rad_last/rad,2.0);
        a_Usph(1) *= pow(rad_last/rad,2.0);
        a_Usph(2) *= pow(rad_last/rad,2.0);
        a_Usph(3) *= pow(rad_last/rad,2.0);
        a_Usph(4) *= pow(rad_last/rad,2.0);
        a_Usph(5) *= pow(rad_last/rad,2.0);
        a_Usph(6) *= pow(rad_last/rad,1.0); 
        a_Usph(7) *= pow(rad_last/rad,1.0);
      },bdryBoxHigh,a_USphStage,radius,radius_ghost);
    }
  }
  // Inner BCs
  if (a_blockBox.low()[0] > bxGhosted.low()[0])
  {
    Point high = bxGhosted.high();
    high[0] = a_blockBox.low()[0] - 1;
    Box bdryBoxLow(bxGhosted.low(), high);
    if (init_condition_type == 0 || init_condition_type == 1){
      forallInPlace_p(f_radialBCs, bdryBoxLow, a_USphStage, radius, dxradius, m_gamma, thickness);
    }
    if (init_condition_type == 2){
      double half = 0.5;
      BoxData<double, DIM, HOST> XCart = forall_p<double,DIM,HOST>
      (f_cubedSphereMap3,radius.box(),radius,dx,half,half,a_block);  
      forallInPlace_p(f_radialoutflowIBCs, bdryBoxLow, a_USphStage, radius, XCart, m_gamma, thickness);
    }
    if (init_condition_type == 3){
      Point source_lo = Point(a_BC.box().low()[0],a_BC.box().low()[1],a_BC.box().low()[2]);
      Point source_hi = Point(a_BC.box().low()[0],a_BC.box().high()[1],a_BC.box().high()[2]);
      Box sourceBox(source_lo,source_hi);
      for (int i = 0; i <= NGHOST; i++){
        a_BC.copyTo(a_USphStage,sourceBox,Point::Basis(0)*(-i));
      }
      // Convert primitive spherical to conserved spherical variables. in correct units. scale by radius.
      forallInPlace_p([ ] PROTO_LAMBDA
      (Point a_pt,
       Var<T, NUMCOMPS, MEM> &a_USph,
       Var<T> &a_radius,
       T a_gamma)
      {  
        double r_BC = CUBED_SPHERE_SHELL_R0; //in cm
        double rad = a_radius(0);
        double rho = a_USph(0)*1.67262192e-24*pow(r_BC/rad,2.0);
        double u = a_USph(1);
        double v = a_USph(2);
        double w = a_USph(3);
        double p = a_USph(4)*1.0e-12*pow(r_BC/rad,2.0*a_gamma); // From picodyne to dyne 
        double Bx = a_USph(5)*1.0e-6*pow(r_BC/rad,2.0); // MicroGauss to Gauss
        double By = a_USph(6)*1.0e-6*pow(r_BC/rad,1.0);
        double Bz = a_USph(7)*1.0e-6*pow(r_BC/rad,1.0);
        double e = p/(a_gamma-1.0) + rho*(u*u+v*v+w*w)/2.0 + (Bx*Bx+By*By+Bz*Bz)/8.0/M_PI;
        a_USph(0) = rho; //rho
        a_USph(1) = rho*u; //Momentum-x
        a_USph(2) = rho*v; //Momentum-y
        a_USph(3) = rho*w; //Momentum-z
        a_USph(4) = e; //Energy
        a_USph(5) = Bx; //Bx
        a_USph(6) = By; //By
        a_USph(7) = Bz; //Bz
      },bdryBoxLow,a_USphStage,radius,m_gamma);
    }
  }
}


/// Exchange, block boundary conditions, transform from JU to USph. Set state variables in boundary registers.
template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::
preStageLevel(MBLevelBoxData<T,NUMCOMPS,MEM>& a_UStage,
              MBLevelBoxData<T,NUMCOMPS,MEM>& a_rhsTotal,
              MBLevelBoxData<T,1,MEM>& a_dVolrLev,
              MBInterpOp& a_iop,
              MBBoundaryRegister<T,NUMREG,MEM>& a_blockreg,
              T a_time,
              T a_dt,
              T a_weight,
              int a_kstage)
{
  // Transform JU -> USph.
  CubedSphereShell::consToSphInterpEuler(a_UStage,a_iop,a_dVolrLev, 4);
  // Compute viscous operators here, and increment a_rhsTotal.
  int linear_visc = ParseInputs::get_linear_viscosity();
  int divBcleaning = ParseInputs::get_divergence_cleaning();
  int init_condition_type = ParseInputs::get_init_condition_type();
  auto map = CubedSphereShell::Map(a_UStage.layout(), ghost());
  auto layout = a_UStage.layout();
  auto eulerOp = CubedSphereShell::Operator<BoxOp_EulerCubedSphere, double, HOST>(map);
  double gamma = ParseInputs::get_gamma();

  HDF5Handler h5;

  int rCoord = CUBED_SPHERE_SHELL_RADIAL_COORD;
  int thetaCoord = (rCoord + 1) % 3;
  int phiCoord = (rCoord + 2) % 3;
  MBLevelBoxData<double, 8, HOST> dstData(layout, Point::Basis(rCoord) + NGHOST*Point::Basis(thetaCoord) + NGHOST*Point::Basis(phiCoord));
  if (init_condition_type == 3) BC_global.BoxData_to_BC(dstData, map, a_time);

  // // Move from ghosted to non-ghosted data to easily see in VisIt
  // MBLevelBoxData<double, 8, HOST> dstData2(a_UStage.layout(), Point::Basis(rCoord) + NGHOST*Point::Basis(thetaCoord) + NGHOST*Point::Basis(phiCoord));
  // dstData2.setVal(0);
  // for (auto dit : dstData.layout())
  // {
  //   Point source_lo = Point(dstData[dit].box().low()[0],dstData[dit].box().low()[1],dstData[dit].box().low()[2]);
  //   Point source_hi = Point(dstData[dit].box().low()[0],dstData[dit].box().high()[1],dstData[dit].box().high()[2]);
  //   Box sourceBox(source_lo,source_hi);
  //   dstData[dit].copyTo(dstData2[dit],sourceBox,Point::Basis(0)*(1));
  // }
  // h5.writeMBLevel(map, dstData2, "DST_DATA2");


  for (auto dit : a_UStage.layout())
  {
    unsigned int block = layout.block(dit);
    Box blockBox = layout.getBlock(block).domain().box();
    auto &rhs_i = a_rhsTotal[dit];
    rhs_i.setVal(0.0);
    auto &USph_i = a_UStage[dit];
    eulerOp[dit].PreStagePatch(USph_i,dstData[dit],a_dVolrLev[dit],blockBox,a_time,a_dt,a_kstage);

    if (linear_visc == 1 && a_kstage == 0)
      {
        auto mbdata = a_blockreg.bounds(dit);
        Array<BoxData<double, NUMCOMPS>, DIM> fluxes;
        fluxes[0].define(rhs_i.box().extrude(0));
        fluxes[0].setToZero();
        fluxes[1].define(rhs_i.box().extrude(1));
        fluxes[1].setToZero();
        fluxes[2].define(rhs_i.box().extrude(2));
        fluxes[2].setToZero();
        eulerOp[dit].LinearVisc(rhs_i, fluxes, USph_i, a_dt);
        T one = 1.0;
        // eulerOp[dit].PostStagePatch(mbdata,fluxes,a_time,a_dt,one);
        //a_blockreg.increment(fluxes, dit, a_dt);
      }
    if (divBcleaning == 1 && a_kstage == 0) eulerOp[dit].DivBPowell(rhs_i, USph_i, a_dt);
  }
}
