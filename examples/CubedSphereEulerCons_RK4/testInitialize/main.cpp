#include "Proto.H"
#include "InputParser.H"
#include "BoxOp_EulerCubedSphere.H"
#include "MHD_IO.H"
#include "MHD_CME.H"
#include <chrono> // Used by timer

MHDReader BC_global;

int main(int argc, char *argv[])
{
  #ifdef PR_MPI
    MPI_Init(&argc, &argv);
  #endif
  ParseInputs::getInstance().parsenow(argc, argv);
  HDF5Handler h5;

  int domainSize = ParseInputs::get_domainSize();
  int thickness = ParseInputs::get_thickness();
  int boxSize_nonrad = ParseInputs::get_boxSize_nonrad();
  int boxSize_rad = ParseInputs::get_boxSize_rad();
  int max_iter = ParseInputs::get_max_iter();
  int temporal_order = ParseInputs::get_temporal_order();
  double gamma = ParseInputs::get_gamma();
  double dt = 0.0;
  double dt_next = 0.0;
  double time = 0.0;
  int restart_step = 0;
  int write_cadence = ParseInputs::get_write_cadence();
  int checkpoint_cadence = ParseInputs::get_checkpoint_cadence();
  int convTestType = ParseInputs::get_convTestType();
  int init_condition_type = ParseInputs::get_init_condition_type();
  int radial_refinement = ParseInputs::get_radial_refinement();
  int MBInterp_define = ParseInputs::get_MBInterp_define();
  string BC_file = ParseInputs::get_BC_file();
  string restart_file = ParseInputs::get_restart_file();
  if (init_condition_type == 3) BC_global.file_to_BoxData_vec(BC_file);
  
  double probe_cadence = ParseInputs::get_Probe_cadence();
  int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD;
  Array<double, DIM> offset = {0., 0., 0.};
  Array<double, DIM> exp = {1., 1., 1.};
  MBLevelBoxData<double, NUMCOMPS, HOST> U_conv_test[3]; 
  PR_TIMER_SETFILE(to_string(thickness) + "_DIM" + to_string(DIM) +"_"+to_string(boxSize_nonrad) + "_" + to_string(boxSize_rad)//+ "_NProc" + to_string(numProc())
                   + "_CubeSphereTest.time.table");
  PR_TIMERS("MMBEuler");
  int levmax = 3;
  auto domain =
    CubedSphereShell::Domain(domainSize, thickness, radialDir);
  // #if 0 //short test
  if (convTestType == 0) levmax = 1;
  if (convTestType == 4) levmax = 1;
  for (int lev=0; lev<levmax; lev++)
    {
      typedef BoxOp_EulerCubedSphere<double, MBMap_CubedSphereShell, HOST> OP;
      bool cullRadialGhost = true;
      bool use2DFootprint = true;
      
      Array<Array<uint, DIM>, 6> permute = {{1,2,0}, {1,2,0}, {1, 0, 2}, {0, 1, 2}, {1, 0, 2}, {0, 1, 2}};
      Array<Array<int, DIM>, 6> sign = {{1, -1, -1}, {1, 1, 1}, {1, -1, 1}, {1, 1, 1}, {-1, 1, 1}, {-1, -1, 1}}; 
      Point boxSizeVect = Point::Ones(boxSize_nonrad);
      boxSizeVect[radialDir] = boxSize_rad;
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
      MBLevelBoxData<double, NUMCOMPS, HOST> USph(layout, OP::ghost());
      MBLevelBoxData<double, NUMCOMPS, HOST> rhs(layout, Point::Zeros());
      MBLevelBoxData<double, 1, HOST> dVolrLev(layout, OP::ghost() + Point::Basis(0, 2));
      Array<double, DIM> dx;
      auto eulerOp = CubedSphereShell::Operator<BoxOp_EulerCubedSphere, double, HOST>(map);
      USph.setVal(0.);
      double dxradius = 1.0 / thickness;
      auto C2C = Stencil<double>::CornersToCells(4);

      int rCoord = CUBED_SPHERE_SHELL_RADIAL_COORD;
      int thetaCoord = (rCoord + 1) % 3;
      int phiCoord = (rCoord + 2) % 3;
      //Point ghst = Point::Ones(NGHOST);
      //MBLevelBoxData<double, 8, HOST> dstData(layout, ghst);
      MBLevelBoxData<double, 8, HOST> dstData(layout, Point::Basis(rCoord) + NGHOST*Point::Basis(thetaCoord) + NGHOST*Point::Basis(phiCoord));
      if (init_condition_type == 3) BC_global.BoxData_to_BC(dstData, map, time);

      // Point ghostForInterp = OP::ghost() + Point::Ones();
      //MBInterpOp iop = CubedSphereShell::InterpOpOld<HOST>(JU.layout(),
      //                                                  ghostForInterp,4);
      MBInterpOp iop;
      iop = CubedSphereShell::InterpOp<HOST>(JU.layout(),OP::ghost(),4);

      // Set input solution.
      for (auto dit : layout)
        {
          dx = eulerOp[dit].dx();
          BoxData<double> radius(dVolrLev[dit].box());
          BoxData<double, DIM, HOST> Dr(dVolrLev[dit].box());
          BoxData<double, DIM, HOST> adjDr(dVolrLev[dit].box());
          eulerOp[dit].radialMetrics(radius, Dr, adjDr, dVolrLev[dit], Dr.box());
        }
      
        {
        Reduction<double,Operation::Max,HOST> dtinv;
        dtinv.reset();
        MBLevelBoxData<double, NUMCOMPS+2, HOST> WSemiCME(layout, JU.ghost());
        MBLevelBoxData<double, NUMCOMPS, HOST> WSphCME(layout, JU.ghost());
        MBLevelBoxData<double, NUMCOMPS+2, HOST> WSemi(layout, JU.ghost());
        MBLevelBoxData<double, NUMCOMPS, HOST> WSph(layout, JU.ghost());      
        MBLevelBoxData<double, NUMCOMPS, HOST> Wout(layout,JU.ghost());
        MBLevelBoxData<double, NUMCOMPS+2, HOST> USemi(layout,JU.ghost());
        for (auto dit : layout)
          {
            dx = eulerOp[dit].dx();
            BoxData<double> radius(dVolrLev[dit].box());
            BoxData<double, DIM, HOST> Dr(dVolrLev[dit].box());
            BoxData<double, DIM, HOST> adjDr(dVolrLev[dit].box());
            eulerOp[dit].radialMetrics(radius, Dr, adjDr, dVolrLev[dit], Dr.box());
            auto block = layout.block(dit);
            auto &JU_i = JU[dit];
            double half = 0.5;
            BoxData<double, DIM, HOST> XCart = forall_p<double,DIM,HOST>
              (f_cubedSphereMap3,radius.box(),radius,dx,half,half,block);
            
            eulerOp[dit].initialize(WSph[dit], dstData[dit], radius,
                                    XCart, gamma, thickness, dx[2], block);
            // Increment WPoint_i with CME perturbation
            CubedSphereShell::
              WSphToWSemiPointwise(WSemi[dit],WSph[dit],dx[1],block);
            BoxData<double,NUMCOMPS,HOST> W_CME(WSemi[dit].box());
            W_CME.setVal(0.);
            forallInPlace_p(define_CME,W_CME,XCart);
            CubedSphereShell::
              WCartToWSemiPointwise(WSemiCME[dit],W_CME,dx[1],block);
            WSemi[dit] += WSemiCME[dit];
            CubedSphereShell::
              WSemiToUSemiPointwise<double,NUMCOMPS,HOST>(USemi[dit],WSemi[dit],dx[1],gamma,block);
          }
        CubedSphereShell::
          USemiSphPointwiseToJUAverage(JU,USemi,iop,dVolrLev,dx);
        Write_W(JU,eulerOp,iop,-1,0,0);
        h5.writeMBLevel({}, map,WSph,"WSphRef" );
        h5.writeMBLevel({}, map,WSemi,"WSemi" );
        h5.writeMBLevel({}, map,WSemiCME,"WSemiCME" );
        h5.writeMBLevel({}, map,USemi,"USemi" );
        }
 // initialization test only
    }
  PR_TIMER_REPORT();
#ifdef PR_MPI
  MPI_Finalize();
#endif
}

