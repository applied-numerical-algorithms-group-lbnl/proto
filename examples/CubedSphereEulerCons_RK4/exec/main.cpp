// #include <gtest/gtest.h>
#include "Proto.H"
#include "InputParser.H"
// #include "Lambdas.H"
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
#if 0
template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_initPatchData_(
                      Point                           a_pt,
                      Var<T,NUMCOMPS,MEM>&            a_U,
                      Var<T>&                         a_radius,
                      Var<T,DIM,MEM,DIM>&             a_A,
                      T                               a_dx,
                      T                               a_dxi0,
                      T                               a_gamma,
                      int                             a_nradius)
{
  // Compute Cartesian conserved quantitities from spherical data.
  T p0 = a_gamma;
  T rho0 = 1.0;
  T eps = .01;
  T amplitude;
  T arg = (a_pt[0] - a_nradius/2)/(a_nradius*1.0);
  if (abs(arg) < .5)
    {
      amplitude = eps*pow(cos(M_PI*arg),6);
    }else
    {
      amplitude = 0.;
    }
  T rho = rho0*(1.0 + amplitude);
  T p = p0*pow(rho/rho0,a_gamma);
  Array<T,DIM> vSpherical = {1.,0.,0.};
  vSpherical[0] *= amplitude;
  Array<T,DIM> vcart;
  T ke;
  for (int dim1 = 0; dim1 < DIM; dim1++)
    {
      vcart[dim1] = 0.;
      // Can set velocities here.
      for (int dim2 = 0; dim2 < DIM; dim2++)
        {
          vcart[dim1] += vSpherical[dim2]*a_A(dim1,dim2);
        }
      a_U(dim1 + 1) = vcart[dim1]*rho;
      ke += .5*vcart[dim1]*vcart[dim1];
    }
  a_U(0) = rho;
  a_U(NUMCOMPS-1) = p/(a_gamma-1) + rho*ke;
}
PROTO_KERNEL_END(f_initPatchData_, f_initPatchData)
#endif

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
  a_W(0) = rho;
  a_W(1) = ur;
  a_W(2) = 0.0;
  a_W(3) = 0.0;
  a_W(NUMCOMPS - 1) = p;
}
PROTO_KERNEL_END(f_radialInit_F, f_radialInit)
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
  ;
  T p = p0 * pow(rho / rho0, a_gamma);
  T ur = amplitude * sqrt(a_gamma * p0 / rho0) / rho0;
  a_USph(0) = rho;
  a_USph(1) = ur * rho;
  a_USph(2) = 0.0;
  a_USph(3) = 0.0;
  a_USph(NUMCOMPS - 1) = p / (a_gamma - 1.0);
}
PROTO_KERNEL_END(f_radialBCs_F, f_radialBCs)
int main(int argc, char *argv[])
{
#ifdef PR_MPI
  MPI_Init(&argc, &argv);
#endif
  typedef BoxOp_EulerCubedSphere<double, MBMap_CubedSphereShell, HOST> OP;
  HDF5Handler h5;
  int domainSize = 32;
  int thickness = 16;
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

  // initialize data and map
  auto map = CubedSphereShell::Map(layout, OP::ghost());
  MBLevelBoxData<double, NUMCOMPS, HOST> JU(layout, OP::ghost());
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
    int r_dir = CUBED_SPHERE_SHELL_RADIAL_COORD;
    double r0 = CUBED_SPHERE_SHELL_R0;
    double r1 = CUBED_SPHERE_SHELL_R1;
    double dr = (r1 - r0) / thickness;
    double dxi0 = dx[0];
    double gamma = 5.0 / 3.0;
    BoxData<double, DIM, HOST> Dr(dVolrLev[dit].box());
    BoxData<double, DIM, HOST> adjDr(dVolrLev[dit].box());
    eulerOp[dit].radialMetrics(radius, Dr, adjDr, dVolrLev[dit], Dr.box(), thickness);
    Box b_i = C2C.domain(layout[dit]).grow(6);
    BoxData<double, DIM> x_i(b_i.grow(Point::Ones()));
    auto block = layout.block(dit);
    auto &JU_i = JU[dit];
    BoxData<double, NUMCOMPS, HOST> JUTemp;
    BoxData<double, NUMCOMPS> WBar_i(JU_i.box());
    BoxData<double, NUMCOMPS, HOST> WPoint_i(JU_i.box());
    forallInPlace_p(f_radialInit, WPoint_i, radius, dxradius, gamma, thickness);
    eulerOp[dit].primToCons(JUTemp, WPoint_i, dVolrLev[dit], gamma, dx[2], block);
    JU_i.setVal(0.);
    JUTemp.copyTo(JU_i, layout[dit]);

    h5.writePatch(dx, JU_i, "JU:Block" + to_string(block));
  }
  cout << "Setup done" << endl;
  for (int ipass = 1; ipass < 2; ipass++)
  {
    if (ipass == 1)
    {
      cout << "performing evaluation for interpolated block boundary data." << endl;
      CubedSphereShell::consToSphInterpEuler(JU, dVolrLev, dx, 4);
    }
    int ghostTest = 6;
    for (auto dit : USph.layout())
    {
      PR_TIMERS("RHS Calculation");
      auto &rhs_i = rhs[dit];
      auto &USph_i = JU[dit];
      Box bx_i = layout[dit];
      Box bxGhosted = USph_i.box();
      BoxData<double> radius(layout[dit].grow(ghostTest));
      int r_dir = CUBED_SPHERE_SHELL_RADIAL_COORD;
      double r0 = CUBED_SPHERE_SHELL_R0;
      double r1 = CUBED_SPHERE_SHELL_R1;
      double dr = (r1 - r0) / thickness;
      double gamma = 5.0 / 3.0;
      BoxData<double, DIM, HOST> Dr(dVolrLev[dit].box());
      BoxData<double, DIM, HOST> adjDr(dVolrLev[dit].box());
      unsigned int block = layout.block(dit);
      eulerOp[dit].radialMetrics(radius, Dr, adjDr, dVolrLev[dit], Dr.box(), thickness);

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
      double dx = 1.0 / domainSize;
      BoxData<double, NUMCOMPS, HOST> Wfoo(layout[dit].grow(ghostTest));
      BoxData<double, NUMCOMPS, HOST> Utemp(USph_i.box());
      forallInPlace_p(f_radialInit, Wfoo, radius, dxradius, gamma, thickness);
      cout << "initial box" << Wfoo.box() << endl;
      eulerOp[dit].primToSph(Utemp, Wfoo, dVolrLev[dit], block_i);
      if (ipass == 0)
      {
        h5.writePatch(dx, Utemp,
                      "USphAnalyticPass" + to_string(ipass) + "Block" + to_string(block_i));
        eulerOp[dit](rhs_i, fluxes, Utemp, block_i, 1.0);
      }
      else
      {
        h5.writePatch(dx, USph_i,
                      "USphGhostPass" + to_string(ipass) + "Block" + to_string(block_i));
        eulerOp[dit](rhs_i, fluxes, USph_i, block_i, 1.0);
      }
      CubedSphereShell::consToSphNGEuler(rhs_i, dVolrLev[dit], layout[dit],
                                         layout.getBlock(block).domain().box(), dx,
                                         layout.block(dit), 4);
      Utemp -= USph_i;
      double maxpforce = rhs_i.absMax(1, 0, 0);
      BoxData<double, NUMCOMPS> rhs_coarse = Stencil<double>::AvgDown(2)(rhs_i);
      double maxpforceC = rhs_coarse.absMax(1, 0, 0);
      cout << std::fixed;
      cout << "pass = " << ipass << " , block = " << block_i << ": "
           << "absmax of rhs = "
           << maxpforce << " , "
           << "absmax of coarsened rhs = "
           << maxpforceC << endl;
      BoxData<double, NUMCOMPS, HOST> rhsSph(rhs_i.box());
      rhs_i.copyTo(rhsSph);
      CubedSphereShell::consToSphNGEuler(rhsSph, dVolrLev[dit], layout[dit],
                                         layout.getBlock(block).domain().box(), dx,
                                         layout.block(dit), 4);

      h5.writePatch(dx, rhsSph, "rhsPatchPass" + to_string(ipass) + "Block" + to_string(block_i));
    }
    h5.writeMBLevel({}, map, rhs, "MBEulerCubedSphereRHSPass" + to_string(ipass));
  }
  PR_TIMER_REPORT();
#ifdef PR_MPI
  MPI_Finalize();
#endif
}
