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
  double dt = 0.01*ParseInputs::get_CFL();
  double dt_next = 0.0;
  double time = 0.0;
  int write_cadence = ParseInputs::get_write_cadence();
  int convTestType = ParseInputs::get_convTestType();
  int init_condition_type = ParseInputs::get_init_condition_type();
  string BC_file = ParseInputs::get_BC_file();
  BC_global.file_to_BoxData_vec(BC_file);
  
  double probe_cadence = ParseInputs::get_Probe_cadence();
  Array<double, DIM> offset = {0., 0., 0.};
  Array<double, DIM> exp = {1., 1., 1.};
  MBLevelBoxData<double, NUMCOMPS, HOST> U_conv_test[3]; 
  PR_TIMER_SETFILE(to_string(domainSize) + "_DIM" + to_string(DIM) //+ "_NProc" + to_string(numProc())
                   + "_CubeSphereTest.time.table");
  PR_TIMERS("MMBEuler");
  int levmax = 3;
  if (convTestType == 0) levmax = 1;
  for (int lev=0; lev<levmax; lev++)
	{
    typedef BoxOp_EulerCubedSphere<double, MBMap_CubedSphereShell, HOST> OP;
    bool cullRadialGhost = true;
    bool use2DFootprint = true;
    int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD;
    Array<Array<uint, DIM>, 6> permute = {{2, 1, 0}, {2, 1, 0}, {1, 0, 2}, {0, 1, 2}, {1, 0, 2}, {0, 1, 2}};
    Array<Array<int, DIM>, 6> sign = {{-1, 1, 1}, {1, 1, -1}, {-1, 1, 1}, {1, 1, 1}, {1, -1, 1}, {-1, -1, 1}};
    auto domain =
        CubedSphereShell::Domain(domainSize, thickness, radialDir);
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
    MBLevelBoxData<double, NUMCOMPS, HOST> JU_Temp(layout, OP::ghost());
    MBLevelBoxData<double, NUMCOMPS, HOST> USph(layout, OP::ghost());
    MBLevelBoxData<double, NUMCOMPS, HOST> rhs_Temp(layout, Point::Zeros());
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

    // Set input solution.
    for (auto dit : layout)
    {
      dx = eulerOp[dit].dx();
      BoxData<double> radius(dVolrLev[dit].box());
      BoxData<double, DIM, HOST> Dr(dVolrLev[dit].box());
      BoxData<double, DIM, HOST> adjDr(dVolrLev[dit].box());
      eulerOp[dit].radialMetrics(radius, Dr, adjDr, dVolrLev[dit], Dr.box());
      auto block = layout.block(dit);
      auto &JU_i = JU[dit];
      BoxData<double, NUMCOMPS, HOST> JUTemp;
      BoxData<double, NUMCOMPS, HOST> WPoint_i(JU_i.box());
      double half = 0.5;
      BoxData<double, DIM, HOST> XCart = forall_p<double,DIM,HOST>
      (f_cubedSphereMap3,radius.box(),radius,dx,half,half,block);  
      eulerOp[dit].initialize(WPoint_i, dstData[dit], radius, XCart, gamma, thickness);

      eulerOp[dit].primToCons(JUTemp, WPoint_i, dVolrLev[dit], gamma, dx[2], block);
      JU_i.setVal(0.);
      JUTemp.copyTo(JU_i, layout[dit]);
    }
    MBInterpOp iop = CubedSphereShell::InterpOp<HOST>(JU.layout(),OP::ghost(),4);
    MBLevelRK4<BoxOp_EulerCubedSphere, MBMap_CubedSphereShell, double> rk4(map, iop);
    Write_W(JU, eulerOp, iop, 0, time, dt_next);
    bool give_space_in_probe_file = true;
    double probe_cadence_temp = 0;
    for (int iter = 1; iter <= max_iter; iter++)
    {
      auto start = chrono::steady_clock::now();
      if (ParseInputs::get_convTestType() == 0){
        #ifdef PR_MPI
          MPI_Reduce(&dt_next, &dt, 1, MPI_DOUBLE, MPI_MIN, 0,MPI_COMM_WORLD);
          MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        #else
          dt = dt_next;
        #endif
        dt *= ParseInputs::get_CFL();
      }

      rk4.advance(JU, dt_next, dVolrLev, dt, time, temporal_order);
      time += dt;
      if (iter % write_cadence == 0) Write_W(JU, eulerOp, iop, iter, time, dt);
      auto end = chrono::steady_clock::now();
      
      int probe_cadence_new = floor(time/probe_cadence);
      if (probe_cadence_new > probe_cadence_temp || iter == 1){
        Probe(JU, map, eulerOp, iop, iter, time, dx, give_space_in_probe_file);
        give_space_in_probe_file = false;
        probe_cadence_temp = probe_cadence_new;
        if(procID() == 0) cout << "Probed" << endl;
      }
      

      if (procID() == 0) cout << "iter = " << iter << " dt = " << dt << " time = " << time  << " Time taken: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;
    }

    if (convTestType > 0){
      U_conv_test[lev].define(layout, {Point::Zeros(),Point::Zeros(),Point::Zeros(),Point::Zeros()});
      for (auto dit : layout)
      {
        JU[dit].copyTo(U_conv_test[lev][dit]);
      }
      h5.writeMBLevel({}, map, U_conv_test[lev], "U_conv_test_" + to_string(lev));
      
      domainSize *= 2;
      thickness *= 2;
      if (convTestType == 2){
        dt /= 2;
        max_iter *= 2;
      }
      if (procID() == 0 ) cout << "domainSize = " << domainSize/2 << " is done" << endl;
    }
  }

  if (convTestType > 0){
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
  }
  PR_TIMER_REPORT();
#ifdef PR_MPI
  MPI_Finalize();
#endif
}
