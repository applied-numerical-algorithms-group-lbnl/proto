// #include <gtest/gtest.h>
#include "Proto.H"
#include "InputParser.H"
// #include "TestFunctions.H"
// #include "MBLevelMap_Shear.H"
// #include "MBLevelMap_XPointRigid.H"
// #include "Proto_CubedSphereShell.H"
#include "BoxOp_EulerCubedSphere.H"

void GetCmdLineArgumenti(int argc, const char **argv, const char *name, int *rtn)
{
  size_t len = strlen(name);
  for (int i = 1; i < argc; i += 2)
  {
    if (strcmp(argv[i] + 1, name) == 0)
    {
      *rtn = atoi(argv[i + 1]);
      std::cout << name << "="
                << " " << *rtn << std::endl;
      break;
    }
  }
}

PROTO_KERNEL_START
template <typename T, MemType MEM>
void f_radialInit_F(
    Point a_pt,
    Var<T, NUMCOMPS, MEM> &a_W,
    Var<T> &a_radius,
    T a_dxradius,
    T a_gamma,
    int a_nradius)
{
  // Compute spherical initial data.
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
  a_W(iRHO) = rho;
  a_W(iVX) = 0.;//ur;
  a_W(iVY) = 0.0;
  a_W(iVZ) = 0.0;
  a_W(iP) = p;
  a_W(iBX) = 0.0;
  a_W(iBY) = 0.0;
  a_W(iBZ) = 0.0;
}

PROTO_KERNEL_END(f_radialInit_F, f_radialInit)

PROTO_KERNEL_START
template <typename T, MemType MEM>
void f_nonradialInit_F(
    Point a_pt,
    Var<T, NUMCOMPS, MEM> &a_W,
    Var<T, DIM> &a_X_cart,
    T a_gamma,
    int a_nradius)
{
  // Compute spherical initial data.
  T p0 = 1.0;
  T rho0 = 1.0;
  T eps = 0.1;
  T amplitude;
  T arg = sqrt((a_X_cart(2)+0.5) * (a_X_cart(2)+0.5) + a_X_cart(0) * a_X_cart(0) + a_X_cart(1) * a_X_cart(1));
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
  a_W(iRHO) = rho;
  a_W(iVX) = 0.;//ur;
  a_W(iVY) = 0.0;
  a_W(iVZ) = 0.0;
  a_W(iP) = p;
  a_W(iBX) = 0.0;
  a_W(iBY) = 0.0;
  a_W(iBZ) = 0.0;
}
PROTO_KERNEL_END(f_nonradialInit_F, f_nonradialInit)

PROTO_KERNEL_START
template <typename T, MemType MEM>
void f_radialBCs_F(
    Point a_pt,
    Var<T, NUMCOMPS, MEM> &a_USph,
    Var<T> &a_radius,
    T a_dxradius,
    T a_gamma,
    int a_nradius)
{
  // Compute sphericl BCs.
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
  a_USph(iMOMX) = 0;//ur * rho;
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

void Write_W(MBLevelBoxData<double, NUMCOMPS, HOST> &a_JU,
             auto &eulerOp,
             MBInterpOp& a_iop,
             double thickness,
             int iter)
{
  Array<double, DIM> dx;
  double dxradius = 1.0 / thickness;
  double gamma = 5.0 / 3.0;
  HDF5Handler h5;
  typedef BoxOp_EulerCubedSphere<double, MBMap_CubedSphereShell, HOST> OP;
  auto layout = a_JU.layout();
  auto map = CubedSphereShell::Map(layout, OP::ghost());

  MBLevelBoxData<double, NUMCOMPS, HOST> USph(layout, OP::ghost());
  MBLevelBoxData<double, NUMCOMPS, HOST> JUTemp(layout, OP::ghost());
  MBLevelBoxData<double, NUMCOMPS, HOST> Wout(layout, OP::ghost());
  MBLevelBoxData<double, 1, HOST> dVolrLev(layout, OP::ghost() + Point::Basis(0, 2));
  MBLevelBoxData<double, 1, HOST> radiusLev(layout, OP::ghost() + 2);
  for (auto dit : a_JU.layout())
  {
    dx = eulerOp[dit].dx();
    BoxData<double, DIM, HOST> Dr(dVolrLev[dit].box());
    BoxData<double, DIM, HOST> adjDr(dVolrLev[dit].box());
    eulerOp[dit].radialMetrics(radiusLev[dit], Dr, adjDr, dVolrLev[dit], Dr.box());//, thickness);
    a_JU[dit].copyTo(JUTemp[dit]);
  }

  CubedSphereShell::consToSphInterpEuler(JUTemp, a_iop,dVolrLev, dx, 4);

  for (auto dit : a_JU.layout())
  {
    Box bxGhosted = JUTemp[dit].box();
    unsigned int block = layout.block(dit);
    // Set radial boundary condtions in radial ghost cells.
    Box blockBox = layout.getBlock(block).domain().box();
    if (blockBox.high()[0] < bxGhosted.high()[0])
    {
      Point low = bxGhosted.low();
      low[0] = blockBox.high()[0] + 1;
      Box bdryBoxHigh(low, bxGhosted.high());
      forallInPlace_p(f_radialBCs, bdryBoxHigh, JUTemp[dit], radiusLev[dit], dxradius, gamma, thickness);
    }
    if (blockBox.low()[0] > bxGhosted.low()[0])
    {
      Point high = bxGhosted.high();
      high[0] = blockBox.low()[0] + 1;
      Box bdryBoxLow(bxGhosted.low(), high);
      forallInPlace_p(f_radialBCs, bdryBoxLow, JUTemp[dit], radiusLev[dit], dxradius, gamma, thickness);
    }
  }

  for (auto dit : a_JU.layout())
  {
    State W;
    State WBar;
    eulerOp[dit].sphToPrim(W, WBar, JUTemp[dit], layout.block(dit));
    W.copyTo(Wout[dit]);
  }

  h5.writeMBLevel({"density","Vr","Vt","Vp","P","Br","Bt","Bp"}, map, Wout, "W_" + to_string(iter));
}

int main(int argc, char *argv[])
{
#ifdef PR_MPI
  MPI_Init(&argc, &argv);
#endif
  typedef BoxOp_EulerCubedSphere<double, MBMap_CubedSphereShell, HOST> OP;
  HDF5Handler h5;
  int domainSize = 32;
  int thickness = 32;
  InputArgs args;
  args.add("nsph", domainSize);
  args.add("nrad", thickness);
  args.parse(argc, argv);
  args.print();
  Array<double, DIM> offset = {0., 0., 0.};
  Array<double, DIM> exp = {1., 1., 1.};
  PR_TIMER_SETFILE(to_string(domainSize) + "_DIM" + to_string(DIM) //+ "_NProc" + to_string(numProc())
                   + "_CubeSphereTest.time.table");
  PR_TIMERS("MMBEuler");

  bool cullRadialGhost = true;
  bool use2DFootprint = true;
  int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD;
  // Array<Point, DIM+1> ghost;
  // ghost.fill(Point::Ones(5));
  Array<Array<uint, DIM>, 6> permute = {{2, 1, 0}, {2, 1, 0}, {1, 0, 2}, {0, 1, 2}, {1, 0, 2}, {0, 1, 2}};
  Array<Array<int, DIM>, 6> sign = {{-1, 1, 1}, {1, 1, -1}, {-1, 1, 1}, {1, 1, 1}, {1, -1, 1}, {-1, -1, 1}};
  auto domain =
      CubedSphereShell::Domain(domainSize, thickness, radialDir);
  Point boxSizeVect = Point::Ones(domainSize);
  boxSizeVect[radialDir] = thickness;
  MBDisjointBoxLayout layout(domain, boxSizeVect);

  int count = 0;
  for (auto dit : layout)
  {
    count++;
  }
  std::cout << "proc_id: " << procID() << ";      num boxes: " << count << std::endl;

  // initialize data and map
  auto map = CubedSphereShell::Map(layout, OP::ghost());
  MBLevelBoxData<double, NUMCOMPS, HOST> JU(layout, OP::ghost());
  MBLevelBoxData<double, NUMCOMPS, HOST> JU_Temp(layout, OP::ghost());
  MBLevelBoxData<double, NUMCOMPS, HOST> USph(layout, OP::ghost());
  MBLevelBoxData<double, NUMCOMPS, HOST> rhs(layout, Point::Zeros());
  MBLevelBoxData<double, 1, HOST> dVolrLev(layout, OP::ghost() + Point::Basis(0, 2));
  Array<double, DIM> dx;
  auto eulerOp = CubedSphereShell::Operator<BoxOp_EulerCubedSphere, double, HOST>(map);
  USph.setVal(0.);
  double dxradius = 1.0 / thickness;
  auto C2C = Stencil<double>::CornersToCells(4);
  // Set input solution.
  for (auto dit : layout)
  {
    dx = eulerOp[dit].dx();
    BoxData<double> radius(layout[dit].grow(6));
    double gamma = 5.0 / 3.0;
    BoxData<double, DIM, HOST> Dr(dVolrLev[dit].box());
    BoxData<double, DIM, HOST> adjDr(dVolrLev[dit].box());
    eulerOp[dit].radialMetrics(radius, Dr, adjDr, dVolrLev[dit], Dr.box());//, thickness);
    Box b_i = C2C.domain(layout[dit]).grow(6);
    BoxData<double, DIM> x_i(b_i.grow(Point::Ones()));
    auto block = layout.block(dit);
    auto &JU_i = JU[dit];
    BoxData<double, NUMCOMPS, HOST> JUTemp;
    BoxData<double, NUMCOMPS> WBar_i(JU_i.box());
    BoxData<double, NUMCOMPS, HOST> WPoint_i(JU_i.box());
    double half = 0.5;
    BoxData<double, DIM, HOST> XCart = forall_p<double,DIM,HOST>
    (f_cubedSphereMap3,radius.box(),radius,dx,half,half,block);  
    forallInPlace_p(f_nonradialInit, WPoint_i, XCart, gamma, thickness);
    // forallInPlace_p(f_radialInit, WPoint_i, radius, dxradius, gamma, thickness);
    eulerOp[dit].primToCons(JUTemp, WPoint_i, dVolrLev[dit], gamma, dx[2], block);
    JU_i.setVal(0.);
    JUTemp.copyTo(JU_i, layout[dit]);
  }

  MBInterpOp iop = CubedSphereShell::InterpOp<HOST>(JU.layout(),OP::ghost(),4);
  Write_W(JU, eulerOp, iop, thickness, 0);
  double dt = 0.01;
  double time = 0.0;
  for (int iter = 1; iter <= 5; iter++)
  {
    if (procID() == 0) cout << "iter = " << iter << endl;
    for (auto dit : layout)
    {
      JU[dit].copyTo(JU_Temp[dit]);
    }
    CubedSphereShell::consToSphInterpEuler(JU_Temp,iop, dVolrLev, dx, 4);
    int ghostTest = 6;
    for (auto dit : USph.layout())
    {
      PR_TIMERS("RHS Calculation");
      auto &rhs_i = rhs[dit];
      auto &USph_i = JU_Temp[dit];
      Box bx_i = layout[dit];
      Box bxGhosted = USph_i.box();
      BoxData<double> radius(layout[dit].grow(ghostTest));
      double gamma = 5.0 / 3.0;
      BoxData<double, DIM, HOST> Dr(dVolrLev[dit].box());
      BoxData<double, DIM, HOST> adjDr(dVolrLev[dit].box());
      unsigned int block = layout.block(dit);
      eulerOp[dit].radialMetrics(radius, Dr, adjDr, dVolrLev[dit], Dr.box());//, thickness);

      // Set radial boundary condtions in radial ghost cells.
      Box blockBox = layout.getBlock(block).domain().box();
      if (blockBox.high()[0] < bxGhosted.high()[0])
      {
        Point low = bxGhosted.low();
        low[0] = blockBox.high()[0] + 1;
        Box bdryBoxHigh(low, bxGhosted.high());
        forallInPlace_p(f_radialBCs, bdryBoxHigh, USph_i, radius, dxradius, gamma, thickness);
      }
      if (blockBox.low()[0] > bxGhosted.low()[0])
      {
        Point high = bxGhosted.high();
        high[0] = blockBox.low()[0] + 1;
        Box bdryBoxLow(bxGhosted.low(), high);
        forallInPlace_p(f_radialBCs, bdryBoxLow, USph_i, radius, dxradius, gamma, thickness);
      }
      int block_i = layout.block(dit);
      Array<BoxData<double, NUMCOMPS>, DIM> fluxes;
      fluxes[0].define(rhs_i.box().extrude(0));
      fluxes[1].define(rhs_i.box().extrude(1));
      fluxes[2].define(rhs_i.box().extrude(2));
      eulerOp[dit](rhs_i, fluxes, USph_i, block_i, 1.0);
      rhs_i *= dt;
      JU[dit] -= rhs_i;
      // CubedSphereShell::consToSphNGEuler(rhs_i,dVolrLev[dit],layout[dit],
      //                                        layout.getBlock(block).domain().box(),1.0/domainSize,
      //                                        layout.block(dit),4);

    }
    // h5.writeMBLevel({ }, map, rhs, "MBEulerCubedSphereRHSPass" + to_string(iter));
    time += dt;
    Write_W(JU, eulerOp, iop, thickness, iter);
  }
  PR_TIMER_REPORT();
#ifdef PR_MPI
  MPI_Finalize();
#endif
}
