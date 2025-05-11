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
      
      Array<Array<uint, DIM>, 6> permute = {{2, 1, 0}, {2, 1, 0}, {1, 0, 2}, {0, 1, 2}, {1, 0, 2}, {0, 1, 2}};
      Array<Array<int, DIM>, 6> sign = {{-1, 1, 1}, {1, 1, -1}, {-1, 1, 1}, {1, 1, 1}, {1, -1, 1}, {-1, -1, 1}}; 
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
      if (restart_file.empty())
        {
        Reduction<double,Operation::Max,HOST> dtinv;
        dtinv.reset();
        MBLevelBoxData<double, NUMCOMPS, HOST> Wout(layout, Point::Zeros());
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
            eulerOp[dit].initialize(WPoint_i, dstData[dit], radius, XCart, gamma, thickness, dx[2], block);
            eulerOp[dit].dtInv(dtinv,WPoint_i);
            eulerOp[dit].primToCons(JUTemp, WPoint_i, dVolrLev[dit], gamma, dx[2], block);
            WPoint_i.copyTo(Wout[dit]);
            // if (procID()==5) h5.writePatch(1,Wout[dit],"W_CME0");
            JU_i.setVal(0.);
            JUTemp.copyTo(JU_i, layout[dit]);
          }
        } else {
          h5.readMBLevel(JU, restart_file);
          time = h5.time();
		      dt = h5.dt();   
          restart_step = h5.iter();
          if (procID() == 0) cout << "Restarting from step " << restart_step << " at time " << time << " with dt " << dt << endl;
        }

      OP::Insert_CME(JU,dVolrLev,iop,layout,time,dt,gamma);

      MBLevelRK4<BoxOp_EulerCubedSphere, MBMap_CubedSphereShell, double> rk4(map, iop);
      OP::P_floor(JU,dVolrLev,iop,layout,time,dt,gamma);
      Write_W(JU, eulerOp, iop, restart_step, time, dt);            
      {
        HDF5Handler h5;
        MBLevelBoxData<double, NUMCOMPS, HOST> JUTemp(JU.layout(), OP::ghost());
        JU.copyTo(JUTemp);
        CubedSphereShell::consToSphInterpEuler(JUTemp,iop,dVolrLev,4);
        
        for (auto dit : JUTemp.layout())
          {
            unsigned int block = layout.block(dit);
            Box blockBox = layout.getBlock(block).domain().box();           
            auto &USph_i = JUTemp[dit];
            eulerOp[dit].PreStagePatch(USph_i,JU[dit],dVolrLev[dit],blockBox,0.,0.,0);
          }
      }
    
  
    bool give_space_in_probe_file = true;
    double probe_cadence_temp = 0;
#if 0 // Begin debug comment.
    auto initialCons = CubedSphereShell::conservationSum(JU);
    for (int comp = 0; comp< 8; comp++)
      {
        if (procID()==0)
          cout << "initial conserved quantity, component " << comp << " = " << initialCons[comp] << endl;
      }
    double mass = initialCons[0];
    double energy =  initialCons[4];
    double momentum = sqrt(energy*mass);
    int discreteVol = domainSize*domainSize*thickness*6;
    if (procID() == 0)
      {
        cout << "mass = " << mass/discreteVol << ", energy = " << energy/discreteVol << ", momentum scale = " << momentum/discreteVol << endl;
      }
#endif // End debug comment.
    if (lev == 0)
      {
        double dtcfl1 = OP::dtCFL(JU,iop,dVolrLev);
        dt = 0.2*dtcfl1*ParseInputs::get_CFL();
      }
    if (convTestType > 2) max_iter = 1;
    
    for (int iter = restart_step + 1; iter <= max_iter; iter++)
    {
      auto start = chrono::steady_clock::now();
      if ((iter % ParseInputs::get_P_floor_cadence()) == 0) OP::P_floor(JU,dVolrLev,iop,layout,time,dt,gamma);
      if (convTestType == 0)
        {
          double dtcfl1 = OP::dtCFL(JU,iop,dVolrLev);
          dt = dtcfl1*ParseInputs::get_CFL();
        }
      if (convTestType < 3)
        {
          rk4.advance(JU, dVolrLev, dt, time, temporal_order);
          time += dt;
        }
#if 0
      else
        {
          PROTO_ASSERT(max_iter == 1,"trying to take more then one time step in applyOp test.");
          Array<double,8> opConsSums = OP::applyOp(JU,iop,dVolrLev);
          HDF5Handler h5;
          h5.writeMBLevel({}, map, JU, "applyOpTest" + to_string(lev));
          if (procID() == 0)
            {
              for (int comp = 0; comp < 8 ; comp++)
                {
                  //opConsSums[comp] *= dt;
                  if (comp == 0)
                    {
                      cout << "component:" << comp <<
                        ", scaled conserved mass = " <<
                        opConsSums[comp]/ mass  <<
                        endl;
                    }
                  else if (comp == 4)
                    {
                      cout << "component:" << comp <<
                        ", scaled conserved energy = " <<
                        opConsSums[comp] / energy <<
                        endl;
                    }
                  else
                    {
                      double norm = sqrt(energy/mass)*mass;
                      cout << "component:" << comp <<
                        ", scaled conserved quantity = " << opConsSums[comp] / norm <<  endl;
                    }
                }
            }
        }
#endif
      if (iter % write_cadence == 0)
        {
          Write_W(JU, eulerOp, iop, iter, time, dt);
          // Check conservation.
#if 0
          if (convTestType < 3)
            {
              Array<double,8> consRadial = rk4.getCons<8>();
              Array<double,8> consSum =
                CubedSphereShell::conservationSum(JU);        
              { 
                for (int comp = 0; comp < 8; comp++)
                  { 
                    consSum[comp] -= initialCons[comp];           
                  }
              }
              if (procID() == 0)
                {
                  double one = 1.0;
                  if (convTestType > 2) one = 0.;
                  for (int comp = 0; comp < 8; comp++)
                    {
                      if (comp == 0)
                        {
                          cout << "component:" << comp <<
                            ", change in scaled conserved mass = " <<
                            (consSum[comp] - consRadial[comp])/mass  <<
                            endl;
                        }
                      else if (comp == 4)
                        {
                          cout << "component:" << comp <<
                            ", change in scaled conserved energy = " <<
                            (consSum[comp] - consRadial[comp]) /energy <<
                            endl;
                        }
                      else
                        {
                          double norm = sqrt(energy/mass)*mass;
                          cout << "component:" << comp <<
                            ", change in scaled conserved quantity = " <<
                            (consSum[comp] - consRadial[comp]) / norm <<
                            endl;
                        }
                    }
                }
            }
#endif
        }
      // Checkpointing.
      if ((iter % checkpoint_cadence == 0) || (iter == max_iter))
        {
          Write_Checkpoint(JU, iter, restart_step, time, dt);
        }
     
      int probe_cadence_new = floor(time/probe_cadence);
      if (probe_cadence_new > probe_cadence_temp || iter == 1){
        Probe(JU, map, eulerOp, iop, iter, time, dx, give_space_in_probe_file);
        give_space_in_probe_file = false;
        probe_cadence_temp = probe_cadence_new;
      }
      
      auto end = chrono::steady_clock::now();
      if (procID() == 0) cout << "iter = " << iter << " dt = " << dt << " time = " << time  << " Time taken: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;
    }

      if ((convTestType > 0) && (convTestType != 4)) {
      U_conv_test[lev].define(layout, {Point::Zeros(),Point::Zeros(),Point::Zeros(),Point::Zeros()});
      for (auto dit : layout)
      {
        JU[dit].copyTo(U_conv_test[lev][dit]);
      }
      h5.writeMBLevel({}, map, U_conv_test[lev], "U_conv_test_" + to_string(lev));
      // int maxRadSize = 64;
      Point refRatio = 2*Point::Ones();
      if (radial_refinement)
        {
          refRatio = 2*Point::Ones() - Point::Basis(1) - Point::Basis(2);
          // maxRadSize = 1024;
          domainSize *= refRatio[1];
          boxSize_nonrad *= refRatio[1];
        } else {
          domainSize *= 2;
          boxSize_nonrad *= 2;
        }
      
      thickness *= 2;
      boxSize_rad *= 2;
      time = 0.;
      domain = domain.refine(refRatio); 
      // if (boxSize_nonrad > maxRadSize)
      //   {
      //     boxSize_nonrad = maxRadSize;
      //     if (procID() == 0) cout << "radial Box Size = " << boxSize_nonrad << endl;
      //   }
      if (convTestType == 2){
        dt /= 2;
        max_iter *= 2;
      }
      if (procID() == 0 ) cout << "thickness = " << thickness/2 << " is done" << endl;
    }
  }

    if((convTestType > 0) && (convTestType != 4))
      {
        //Here, we perform the error calculations on a single patch.
        // if(procID() == 0)
        Array<Array<double,8>,2> errmaxglobal;
        for(int lev=0; lev<2; lev++)
          {
            Array<Reduction<double,Operation::Max,HOST>,8> errmax;
            for (int comp = 0; comp < 8;comp++)
              {
                errmax[comp].reset();
              }
            auto layout = U_conv_test[lev].layout();
            MBLevelBoxData<double,8,HOST> err(layout,Point::Zeros());
            Point refRatio;
            if (radial_refinement)
              {
                refRatio = 2*Point::Ones() - Point::Basis(1) - Point::Basis(2);
              }
            else
              {
                refRatio = 2*Point::Ones();
              }          
            U_conv_test[lev + 1].coarsenTo(err,refRatio);
            
            for (auto dit : layout)
              {
                unsigned int block = layout.block(dit);
                Box blockBox = layout.getBlock(block).domain().box();
                Box radialBlockFaceHi = blockBox.face(0,Side::Hi,3);
                Box radialPatchFaceHi = radialBlockFaceHi&err[dit].box();
                Box radialBlockFaceLo = blockBox.face(0,Side::Lo,3);
                Box radialPatchFaceLo = radialBlockFaceLo&err[dit].box();
                err[dit] -= U_conv_test[lev][dit];
                for (int comp = 0; comp < 8; comp++)
                  {
                    auto errslice = slice(err[dit],comp);
                    if (!radialPatchFaceHi.empty())
                      {
                        errslice.setVal(0.,radialPatchFaceHi);
                      }
                    if (!radialPatchFaceLo.empty())
                      {
                        errslice.setVal(0.,radialPatchFaceLo);
                      }
                    double errPatch = errslice.absMax();
                    errmax[comp].reduce(&errPatch,1);
                  }
              }
            HDF5Handler h5;
            auto map = CubedSphereShell::Map(err);
            h5.writeMBLevel({"err"}, map, err,"Err"+to_string(lev));
            for (int comp = 0; comp < 8; comp++)
              {
                errmaxglobal[lev][comp] = errmax[comp].fetch();
                if (procID() == 0)
                  {
                    std::cout << "Lev = " << lev << ", component = "
                              << comp <<", error = " << errmaxglobal[lev][comp] << std::endl;
                  }
              }
          }
        for (int comp = 0; comp < 8 ; comp++)
          {
            double rate = log(abs(errmaxglobal[0][comp]/errmaxglobal[1][comp]))/log(2.0);
            if (procID() == 0) std::cout << "order of accuracy for var " << comp << " = " << rate << std::endl;
          }    
      }
// #endif // end short test
//#endif // end debug comment.
  PR_TIMER_REPORT();
#ifdef PR_MPI
  MPI_Finalize();
#endif
}
