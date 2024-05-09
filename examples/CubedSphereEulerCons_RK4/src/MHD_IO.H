#include "Proto.H"
#include "Inputs_Parsing.H"
#include "BoxOp_EulerCubedSphere.H"

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
  double gamma = ParseInputs::get_gamma();
  string BC_file = ParseInputs::get_BC_file();
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
    eulerOp[dit].radialMetrics(radiusLev[dit], Dr, adjDr, dVolrLev[dit], Dr.box());
    a_JU[dit].copyTo(JUTemp[dit]);
  }

  CubedSphereShell::consToSphInterpEuler(JUTemp, a_iop,dVolrLev, 4);

  for (auto dit : a_JU.layout())
  {
    State W;
    // State WBar;
    // eulerOp[dit].sphToPrim(W, WBar, JUTemp[dit], layout.block(dit));
    // Second order for output purposes. This then doesn't need BCs.
    eulerOp[dit].sphToPrim_2O(W, JUTemp[dit], layout.block(dit)); 
    W.copyTo(Wout[dit]);
  }

  h5.writeMBLevel({"density","Vr","Vt","Vp","P","Br","Bt","Bp"}, map, Wout, "W_" + to_string(iter));
  // h5.writeMBLevel({"density","Momr","Momt","Momp","E","Br","Bt","Bp"}, map, JUTemp, "U_sph_" + to_string(iter));
}