#include "Proto.H"
#include "Inputs_Parsing.H"
#include "BoxOp_EulerCubedSphere.H"
template <typename T, MemType MEM>
BoxData<T,1,MEM> myRadialCoord(T a_dxi0, Box a_bx)
{
  int r_dir = CUBED_SPHERE_SHELL_RADIAL_COORD;
  double r0 = CUBED_SPHERE_SHELL_R0;
  double r1 = CUBED_SPHERE_SHELL_R1;
  T dxi0 = a_dxi0;
  double dr = (r1-r0)*dxi0;
  BoxData<T,1,MEM> radius =
    forall_p<T,1,MEM>([]PROTO_LAMBDA(Point a_pt,
                                     Var<double>& a_rad,
                                     double a_r0,
                                     double a_r1,
                                     double a_dr,
                                     double a_dxi0,
                                     int a_rdir)
                  {
                    // a_rad(0) = a_r0 + a_dr*a_pt[a_rdir];
                    // T rlow = a_rad(0);
                    // T rhigh = a_rad(0)+a_dr;
                    T C_rad = 1.0;
                    T R_t = (a_r1 - a_r0) / (exp(C_rad) - 1.0);
                    T etalow = a_pt[a_rdir]*a_dxi0;
                    T rlow = a_r0 + R_t * (exp(C_rad * etalow) - 1.0);                            // a_rad(0) = rlow;
                    a_rad(0) = rlow;                
                  },a_bx,r0,r1,dr,dxi0,r_dir);
  return radius;
}
void Write_W(MBLevelBoxData<double, NUMCOMPS, HOST> &a_JU,
             auto &eulerOp,
             MBInterpOp& a_iop,
             double thickness,
             int iter,
             double time = 0.0,
             double dt = 0.0)
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

  CubedSphereShell::consToSphInterpEuler(JUTemp, a_iop,dVolrLev, 4);

  for (auto dit : a_JU.layout())
  {
    int kstage = 0;
    eulerOp[dit].PreStagePatch(JUTemp[dit],dVolrLev[dit],time,dt,kstage);
  }

  for (auto dit : a_JU.layout())
  {
    State W;
    State WBar;
    eulerOp[dit].sphToPrim(W, WBar, JUTemp[dit], layout.block(dit));
    W.copyTo(Wout[dit]);
  }

  h5.writeMBLevel({"density","Vr","Vt","Vp","P","Br","Bt","Bp"}, map, Wout, "W_visc_" + to_string(iter));
}
template <typename T, MemType MEM>
void myXCart(BoxData<T,DIM,MEM>& a_xcart,
             // const BoxData<T,1,MEM>& a_radius,
             // Array<T,DIM> a_dxi,
             T a_offseta,
             T a_offsetb,  
             int a_block)
{
  Box bx = a_xcart.box();
  Box bx0 = bx.edge(Point::Basis(0,-1),1);
  auto permval = CubedSphereShell::perm(a_block);
  auto sign = CubedSphereShell::permsign(a_block);
  // DEBUG
  T one = 1.0;
  a_xcart.setToZero();
  forallInPlace_p([ ] PROTO_LAMBDA(Point a_pt,
                                   Var<T,DIM,MEM>& a_xcart,
                                   Array<uint,DIM> a_perm,
                                   Array<int,DIM> a_sign,
                                   T a_offseta,
                                   T a_offsetb,  
                                   T a_h)
                  {
                    T alpha = -.25*M_PI +
                      .5*M_PI*(a_pt[1]*a_h + a_offseta*a_h);
                    T beta = -.25*M_PI +
                        .5*M_PI*(a_pt[2]*a_h + a_offsetb*a_h);
                    T tana = tan(alpha);
                    T tanb = tan(beta);
                    T norm = 1.0/sqrt(1.0 + tana*tana + tanb*tanb); 
                    // a_xcart(0) = 1.0*a_pt[a_perm[0]];
                    // a_xcart(1) = 1.0*a_pt[a_perm[1]];                   
                    // a_xcart(2) = 1.0*a_pt[a_perm[2]];
                    a_xcart(0) = 1.0*a_pt[0];
                    a_xcart(1) = 1.0*a_pt[1];                   
                    a_xcart(2) = 1.0*a_pt[2];
                    // Array<T,DIM> xcart;
                    // if (a_pt[0] == 16)
  //                     {
  //                       //a_xcart(0) = 1.0*a_pt[1];
  //                       //a_xcart(1) = 1.0*a_pt[2];
  //                       xcart[0] = norm;
  //                       xcart[1] = tana*norm;
  //                       xcart[2] = tanb*norm;
  //                       a_xcart(0) = a_sign[0]*xcart[a_perm[0]];
  //                       a_xcart(1) = a_sign[1]*xcart[a_perm[1]];
  //                       a_xcart(2) = a_sign[2]*xcart[a_perm[2]];
  //                     }
                  },bx,a_xcart,permval,sign,a_offseta,a_offsetb,one);
}
int main(int argc, char *argv[])
{
#ifdef PR_MPI
  MPI_Init(&argc, &argv);
#endif
  ParseInputs::getInstance().parsenow(argc, argv);
  HDF5Handler h5;
  int domainSize = ParseInputs::get_domainSize();
  int thickness = ParseInputs::get_thickness();
  int max_iter = ParseInputs::get_max_iter();
  double dt = 0.01*ParseInputs::get_CFL();
  double dt_next = 0.0;
  int write_cadence = ParseInputs::get_write_cadence();
  int convTestType = ParseInputs::get_convTestType();
  int init_condition_type = ParseInputs::get_init_condition_type();
  Array<double, DIM> offset = {0., 0., 0.};
  Array<double, DIM> exp = {1., 1., 1.};
  MBLevelBoxData<double, NUMCOMPS, HOST> U_conv_test[3]; 
  PR_TIMER_SETFILE(to_string(domainSize) + "_DIM" + to_string(DIM) //+ "_NProc" + to_string(numProc())
                   + "_CubeSphereTest.time.table");
  PR_TIMERS("MMBEuler");
  int levmax = 1;
  if (convTestType == 0) levmax = 1;
  int lev = 0;
  typedef BoxOp_EulerCubedSphere<double, MBMap_CubedSphereShell, HOST> OP;
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
  MBLevelBoxData<double, NUMCOMPS, HOST> rhs_Temp(layout, Point::Zeros());
  MBLevelBoxData<double, NUMCOMPS, HOST> rhs(layout, Point::Zeros());
  MBLevelBoxData<double, 1, HOST> dVolrLev(layout, OP::ghost() + Point::Basis(0, 2));
  Array<double, DIM> dx;
  MBLevelBoxData<double,DIM, HOST> XCoords(layout, Point::Zeros());
  auto eulerOp = CubedSphereShell::Operator<BoxOp_EulerCubedSphere, double, HOST>(map);
  USph.setVal(0.);
  double dxradius = 1.0 / thickness;
  auto C2C = Stencil<double>::CornersToCells(4);
  // Set input solution.
  for (auto dit : layout)
    {
      dx = eulerOp[dit].dx();
      // BoxData<double> radius(layout[dit].grow(6));
      BoxData<double> radius(dVolrLev[dit].box());
      double gamma = 5.0 / 3.0;
      BoxData<double, DIM, HOST> Dr(dVolrLev[dit].box());
      BoxData<double, DIM, HOST> adjDr(dVolrLev[dit].box());
      eulerOp[dit].radialMetrics(radius, Dr, adjDr, dVolrLev[dit], Dr.box());
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
      //if (init_condition_type == 0) forallInPlace_p(f_radialInit, WPoint_i, radius, dxradius, gamma, thickness);
      //if (init_condition_type == 1) forallInPlace_p(f_nonradialInit, WPoint_i, XCart, gamma, thickness);
      //if (init_condition_type == 2) forallInPlace_p(f_radialoutflowInit, WPoint_i, radius, XCart, gamma, thickness);
      //eulerOp[dit].primToCons(JUTemp, WPoint_i, dVolrLev[dit], gamma, dx[2], block);
      //JU_i.setVal(0.);
      //JUTemp.copyTo(JU_i, layout[dit]);
    }

  MBInterpOp iop = CubedSphereShell::InterpOp<HOST>(JU.layout(),OP::ghost(),4);
  MBBoundaryRegister<double,DIM+1,HOST> blockreg;
  blockreg.define(layout,1,Point::Zeros());
  blockreg.clear();
  // Fill registers.
  {
    MBLevelBoxData<double, DIM+1, HOST> rhsTest(layout, Point::Zeros());
    HDF5Handler h5;
    auto domain = layout.domain();
    for (auto dit : layout)
      {
        Array<BoxData<double,DIM+1,HOST>,DIM> fluxes;
        Box bx = layout[dit];
        Box bxCorners = bx.extrude(Point::Ones());
        auto block = layout.block(dit);
        rhsTest[dit].setToZero();
        BoxData<double,1,HOST> radius =
          myRadialCoord<double,HOST>(dx[0],bxCorners);
        myXCart(XCoords[dit],0.,0.,block);
        for (int d = 0; d < DIM; d++)
          {
            fluxes[d].define(bxCorners);
            fluxes[d].setVal(2.5);              
          }
        blockreg.increment(fluxes,dit,1.0);     
      }
    h5.writeMBLevel({"X","Y","Z","block"}, map, XCoords, "XCoord");
    //blockreg.reflux(rhsTest);
    blockreg.exchange();
    for (auto dit : layout)
      {
        auto block = layout.block(dit);
            
        auto bx = layout[dit];
        auto& refluxOut = rhsTest[dit];         
        auto mbdata = blockreg.bounds(dit);
        BoxData<double,DIM+1,HOST> fluxLocalOut
          (bx.extrude(Point::Ones() - Point::Basis(0)));
        BoxData<double,DIM+1,HOST> fluxAdjOut
          (bx.extrude(Point::Ones() - Point::Basis(0)));
        // BoxData<double,DIM+1,HOST> refluxLocal(bx);
        // DEBUG
        fluxLocalOut.setVal(0.);
        fluxAdjOut.setVal(0.);
        for (auto it = mbdata.begin(); it != mbdata.end(); it++)
          {
            auto mb = *it;
                
            pair<int,Side::LoHiSide> direction = mb.direction(layout);
            int dir = direction.first;
            auto side = direction.second;
            auto R = mb.adjToLocal;
            auto Rinv = R.inverse();
            Point normal = Point::Basis(dir,side);
            Point adjNormal = Rinv(normal);
            adjNormal *= -1;
            Point adjBasis = adjNormal*adjNormal;
            int dp = adjNormal.dot(adjBasis);
            Box bxfaces = bx.extrude(Point::Basis(dir));
            Box bxedge = bxfaces.face(dir,side);
            BoxData<double,DIM+1,HOST> adjCopy(bxedge);
            if (dp*(2*side-1) == 1)
              {
                mb.adjData->shift(Point::Basis(dir));
                mb.adjData -> copyTo(adjCopy);
                mb.adjData->shift(Point::Basis(dir,-1));
              }
            else
              {
                mb.adjData -> copyTo(adjCopy);
              }
            if (adjCopy.absMax(0) == 100)
              { 
                cout << "block = " << block  << ", " << "direction = " << dir << ", " <<
                  "side = " << 2*side - 1 << ", " << ", dp = "<< dp <<
                  ", adj Box = " <<  mb.adjData -> box() << ", bxedge = " << bxedge << endl;
              }

            mb.localData -> copyTo(fluxLocalOut,bxedge);
            adjCopy.copyTo(fluxAdjOut,bxedge);      
            h5.writePatch({"X","Y","Z","block"}, fluxLocalOut, "fluxlocal"+to_string(block));
            h5.writePatch({"X","Y","Z","block"}, fluxAdjOut, "fluxadj"+to_string(block));             
          }
      }
  }  
  PR_TIMER_REPORT();
#ifdef PR_MPI
  MPI_Finalize();
#endif
}
