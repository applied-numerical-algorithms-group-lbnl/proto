#include "Proto.H"
#include "Inputs_Parsing.H"
#include "BoxOp_EulerCubedSphere.H"
#include "MHD_IO.H"
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
  //auto eulerOp = CubedSphereShell::Operator<BoxOp_EulerCubedSphere, double, HOST>(map);
  double dxradius = 1.0 / thickness;
  
  MBInterpOp iop;
  iop = CubedSphereShell::InterpOp<HOST>(layout, OP::ghost(),4);

{
  Point targetPoint(90, 49, -7);
  MBIndex targetIndex = layout.find(Point::Zeros(), 3);
  MBDataPoint targetDataPoint(targetIndex, targetPoint, layout);
  auto pointOp = iop(targetDataPoint);
  if (procID() == 0)
  {
    std::cout << "target Point: " << pointOp.target() << std::endl;
    std::cout << "sources: " << std::endl;
    for (auto src : pointOp.sources())
    {
      std::cout << "\t" << src << std::endl;
    }
  }
}
{
  Point targetPoint(90, 50, -7);
  MBIndex targetIndex = layout.find(Point::Zeros(), 3);
  MBDataPoint targetDataPoint(targetIndex, targetPoint, layout);
  auto pointOp = iop(targetDataPoint);
  if (procID() == 0)
  {
    std::cout << "target Point: " << pointOp.target() << std::endl;
    std::cout << "sources: " << std::endl;
    for (auto src : pointOp.sources())
    {
      std::cout << "\t" << src << std::endl;
    }
  }
}
  MBLevelBoxData<double, 1, HOST> UComp(layout, OP::ghost());
  h5.readMBLevel(UComp, "UComp");
  h5.writeMBLevel({}, map, UComp, "UCompPre");
  int velycomp = 3;
  
  string strCounter = 
    "_" + to_string(restart_step+1) + "_" + to_string(0);
  for (auto dit : layout)
    {
      unsigned int block = layout.block(dit);
      h5.writePatch(1.0,UComp[dit],
                    "UCompPreBlock_"+to_string(block)+strCounter);
    }
  h5.writeMBLevelBounds({"data"}, UComp, "U0");
  UComp.exchange();
  h5.writeMBLevelBounds({"data"}, UComp, "U1");
  iop.apply(UComp,UComp);
  h5.writeMBLevelBounds({"data"}, UComp, "U2");

  h5.writeMBLevel({}, map, UComp, "UCompPost");
  for (auto dit : layout)
    {
      unsigned int block =layout.block(dit);
      h5.writePatch(1.0,UComp[dit],
                    "UCompPostBlock_"+to_string(block)+strCounter);
    }
  // #endif // end short test
  //#endif // end debug comment.
  PR_TIMER_REPORT();
#ifdef PR_MPI
  MPI_Finalize();
#endif
}
