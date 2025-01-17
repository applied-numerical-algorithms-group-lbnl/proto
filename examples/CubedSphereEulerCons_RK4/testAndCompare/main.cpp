#include "Proto.H"
#include "Inputs_Parsing.H"
#include "BoxOp_EulerCubedSphere.H"
#include "MHD_IO.H"
#include <chrono> // Used by timer
#define COMPARECOMPS NUMCOMPS+2
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
  string UOldInterp_file = ParseInputs::get_data_file_prefix();
  string UPreInterp_file = ParseInputs::get_checkpoint_file_prefix();
  if (init_condition_type == 3) BC_global.file_to_BoxData_vec(BC_file);
  
  double probe_cadence = ParseInputs::get_Probe_cadence();
  int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD; 
  PR_TIMER_SETFILE(to_string(thickness) + "_DIM" + to_string(DIM) +"_"+to_string(boxSize_nonrad) + "_" + to_string(boxSize_rad)//+ "_NProc" + to_string(numProc())
                   + "_testInterp.time.table");
  PR_TIMERS("MMBTestInterp");
  
  bool cullRadialGhost = true;
  bool use2DFootprint = true;
  auto domain =
    CubedSphereShell::Domain(domainSize, thickness, radialDir);
  typedef BoxOp_EulerCubedSphere<double, MBMap_CubedSphereShell, HOST> OP;
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
  
  Array<double, DIM> dx;
  auto eulerOp = CubedSphereShell::Operator<BoxOp_EulerCubedSphere, double, HOST>(map);
  double dxradius = 1.0 / thickness;
  //Point nghost = Point::Zeros();
  Point nghost = OP::ghost();
  MBLevelBoxData<double,COMPARECOMPS, HOST> USemi(layout,nghost);
  MBLevelBoxData<double, COMPARECOMPS, HOST> USemiOldInterp(layout,nghost);
  MBLevelBoxData<double, COMPARECOMPS, HOST> dU(layout,nghost);
  h5.readMBLevel(USemi, UPreInterp_file);
  h5.readMBLevel(USemiOldInterp, UOldInterp_file);
  MBInterpOp iop;
  iop = CubedSphereShell::InterpOp<HOST>(layout, OP::ghost(),4);
  USemi.exchange();  
  iop.apply(USemi,USemi);
  
  h5.writeMBLevel({}, map, USemiOldInterp, "USemiOldInterp");
  h5.writeMBLevel({}, map, USemi, "USemiNewInterp");
  for (auto dit : layout)
    {
      unsigned int block = layout.block(dit);
      h5.writePatch(1.0,USemiOldInterp[dit],"USemiOldInterp_"+to_string(block));
      h5.writePatch(1.0,USemi[dit],"USemiNewInterp_"+to_string(block));
      dU[dit].setVal(0.);
      dU[dit] += USemiOldInterp[dit];
      dU[dit] *= -1.;
      dU[dit] +=USemi[dit];
      forallInPlace([] PROTO_LAMBDA(Var<double,COMPARECOMPS,HOST>& du)
                    {
                      for (int comp = 0; comp < COMPARECOMPS; comp++)
                        {
                          du(comp) = abs(du(comp));
                        }
                    },dU[dit]);
      h5.writePatch(1.0,dU[dit],"dU_"+to_string(block));
    }
  h5.writeMBLevel({}, map, dU, "dU");
  PR_TIMER_REPORT();
#ifdef PR_MPI
  MPI_Finalize();
#endif
}
