// #include <gtest/gtest.h>
#include "Proto.H"
#include "InputParser.H"
// #include "Lambdas.H"
// #include "MBLevelMap_Shear.H"
// #include "MBLevelMap_XPointRigid.H"
// #include "Proto_CubedSphereShell.H"
#include "BoxOp_EulerCubedSphere.H"

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
  a_W(iVX) = ur;
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
  int half_nr = a_nradius / 2;
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
PROTO_KERNEL_END(f_radialBCs_F, f_radialBCs)


int main(int argc, char *argv[])
{
#ifdef PR_MPI
  MPI_Init(&argc, &argv);
#endif
  
  int domainSize = 32;
  int thickness = 64;
  int iter_max = 2;
  double dt = 0.001;
  MBLevelBoxData<double, NUMCOMPS, HOST> U_conv_test[3]; 
  PR_TIMER_SETFILE(to_string(domainSize) + "_DIM" + to_string(DIM)
                 + "_CubeSphereTest.time.table");
  for (int lev=0; lev<3; lev++)
	{
    typedef BoxOp_EulerCubedSphere<double, MBMap_CubedSphereShell, HOST> OP;
    HDF5Handler h5;
    bool cullRadialGhost = true;
    bool use2DFootprint = true;
    int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD;
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
    double time = 0.0;
    for (int iter = 1; iter <= iter_max; iter++)
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
        auto &JU_i = JU[dit];
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
      }
      time += dt;
    }
    U_conv_test[lev].define(layout, {Point::Zeros(),Point::Zeros(),Point::Zeros(),Point::Zeros()});
    for (auto dit : layout)
    {
      // rhs[dit].copyTo(U_conv_test[lev][dit]);
      JU[dit].copyTo(U_conv_test[lev][dit]);
    }
    h5.writeMBLevel({}, map, U_conv_test[lev], "U_conv_test_" + to_string(lev));
    domainSize *= 2;
    thickness *= 2;
    if (procID() == 0) cout << "domainSize = " << domainSize/2 << " is done" << endl;
  }

  //Here, we perform the error calculations on a single patch.
	if(procID() == 0)
	{
		for (int varr = 0; varr < NUMCOMPS; varr++) {
					double ErrMax[2];
			for(int ilev=0; ilev<2; ilev++)
			{
				auto dit_lev=U_conv_test[ilev].begin();
				auto dit_levp1=U_conv_test[ilev+1].begin();

				BoxData<double,1> err=slice(U_conv_test[ilev][*dit_lev],varr);
				err-=Stencil<double>::AvgDown(2)(slice(U_conv_test[ilev+1][*dit_levp1],varr));
				ErrMax[ilev]=err.absMax();
				std::string filename="Comp_"+std::to_string(varr)+"_err_"+std::to_string(ilev);
        HDF5Handler h5;
        h5.writePatch({"err"}, 1, err, filename.c_str());					
				std::cout << "Lev: " << ilev << " , " << ErrMax[ilev] << std::endl;
			}
			double rate = log(abs(ErrMax[0]/ErrMax[1]))/log(2.0);
			std::cout << "order of accuracy for var " << varr << " = " << rate << std::endl;
		}
	}
PR_TIMER_REPORT();
#ifdef PR_MPI
  MPI_Finalize();
#endif
}
