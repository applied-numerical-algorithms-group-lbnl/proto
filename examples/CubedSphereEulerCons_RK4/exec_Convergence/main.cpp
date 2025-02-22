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
            JU_i.setVal(0.);
            JUTemp.copyTo(JU_i, layout[dit]);
          }
        } else {
          h5.readMBLevel(JU, restart_file);
          time = h5.time();
		      dt = h5.dt();   
          restart_step = h5.iter();
          if (procID() == 0) cout << "Restarting from step " << restart_step << " at time " << time << " with dt " << dt << endl;
          if (ParseInputs::get_CME_type() == 1){
            if (procID() == 0) cout << "Inserting CME" << endl;
            MBLevelBoxData<double, NUMCOMPS+2, HOST> Wout(layout, Point::Zeros());
            
            for (auto dit : layout)
            {
              dx = eulerOp[dit].dx();
              BoxData<double> radius(dVolrLev[dit].box());
              BoxData<double, DIM, HOST> Dr(dVolrLev[dit].box());
              BoxData<double, DIM, HOST> adjDr(dVolrLev[dit].box());
              eulerOp[dit].radialMetrics(radius, Dr, adjDr, dVolrLev[dit], Dr.box());
              auto block = layout.block(dit);
              auto &JU_i = JU[dit];

              double dxiPerp = dx[2];
              BoxData<double,NUMCOMPS+2,HOST> USemiSph(layout[dit]);
              auto& Wout_i = Wout[dit];
              USemiSph.setToZero();
              Box blockBox = layout.getBlock(block).domain().box();
              auto& dVolr_i = dVolrLev[dit];
              Box bx_i = layout[dit];
              CubedSphereShell::
              consToSemiSphNGEuler(USemiSph,JU_i,dVolr_i,bx_i,
                                  blockBox,dxiPerp,block,4);
              forallInPlace([ ] PROTO_LAMBDA
                (Var<double,NUMCOMPS+2,HOST>& a_Wout,
                  const Var<double,NUMCOMPS+2,HOST>& a_USemiSph,
                  double a_gamma,
                  const Var<double,1,HOST>& a_radius)
                {
                  a_Wout(0) = a_USemiSph(0);
                  a_Wout(1) = a_USemiSph(1) / a_Wout(0);
                  a_Wout(2) = a_USemiSph(2) / a_Wout(0);
                  a_Wout(3) = a_USemiSph(3) / a_Wout(0);
                  a_Wout(5) = a_USemiSph(5);
                  a_Wout(6) = a_USemiSph(6);
                  a_Wout(7) = a_USemiSph(7);
                  a_Wout(8) = a_USemiSph(8) / a_Wout(0);
                  a_Wout(9) = a_USemiSph(9);
                  double ke = a_Wout(0)*(a_Wout(1)*a_Wout(1)
                                      + a_Wout(2)*a_Wout(2)
                                      + a_Wout(3)*a_Wout(3) + a_Wout(8)*a_Wout(8))*.5;
                  double me = (a_Wout(5)*a_Wout(5)
                          + a_Wout(6)*a_Wout(6)
                          + a_Wout(7)*a_Wout(7) + a_Wout(9)*a_Wout(9))/(8.0*M_PI);
                  a_Wout(4) = (a_USemiSph(4) - ke - me)*(a_gamma - 1.0);
                }, Wout[dit], USemiSph, gamma, radius);


              BoxData<double, NUMCOMPS, HOST> JU_CME;
              BoxData<double, NUMCOMPS, HOST> W_CME(JU_i.box());
              W_CME.setVal(0.);
              double half = 0.5;
              BoxData<double, DIM, HOST> XCart = forall_p<double,DIM,HOST>
                (f_cubedSphereMap3,radius.box(),radius,dx,half,half,block);  
              define_CME_calc(W_CME, XCart); // Set CME values in cartesian coordinates

              // Calculate A_matrix
              Box bx = W_CME.box();
              
              Array<Array<uint,DIM>,6> permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
              Array<Array<int,DIM>,6> sign = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}}; 
              Point high = bx.high();
              Point low = bx.low();
              high[0] = low[0];
              Box bx0(low,high);
              BoxData<double ,DIM,HOST,DIM> A_matrix(bx);
              double offseta = half;
              double offsetb = half;
              forallInPlace_p(f_Amatrix,bx0,A_matrix,permute[block],sign[block],
                              dxiPerp,offseta,offsetb);
              spreadSlice(A_matrix);

              BoxData<double,DIM,HOST,DIM> A_matrix_inv(bx);
              A_matrix_inv.setToZero();
              forallInPlace(f_matinv3by3,A_matrix_inv,A_matrix);

              BoxData<double,DIM,HOST> V_CME_sph(bx);
              BoxData<double,DIM,HOST> B_CME_sph(bx);

              BoxData<double,DIM,HOST> V_CME_cart = slice<double,NUMCOMPS,DIM,HOST>(W_CME,iVX);
              BoxData<double,DIM,HOST> B_CME_cart = slice<double,NUMCOMPS,DIM,HOST>(W_CME,iBX);
              double one = 1.0;  
              forallInPlace(f_matVecProd,V_CME_sph,A_matrix_inv,V_CME_cart,one);
              forallInPlace(f_matVecProd,B_CME_sph,A_matrix_inv,B_CME_cart,one);

              forallInPlace([ ] PROTO_LAMBDA
                          (Var<double, NUMCOMPS, HOST> &a_a_W,
                          Var<double, DIM, HOST> &a_V_sph,
                          Var<double, DIM, HOST> &a_B_sph)
              {  
                a_a_W(iVX) = a_V_sph(0);
                a_a_W(iVY) = a_V_sph(1);
                a_a_W(iVZ) = a_V_sph(2);
                a_a_W(iBX) = a_B_sph(0);
                a_a_W(iBY) = a_B_sph(1);
                a_a_W(iBZ) = a_B_sph(2);
              },W_CME, V_CME_sph, B_CME_sph);

              eulerOp[dit].primToCons(JU_CME, W_CME, dVolrLev[dit], gamma, dx[2], block);

              forallInPlace([ ] PROTO_LAMBDA
                          (Var<double, NUMCOMPS, HOST> &a_JU,
                          Var<double, NUMCOMPS, HOST> &a_JU_CME,
                          Var<double, NUMCOMPS, HOST> &a_W_CME,
                          Var<double, NUMCOMPS+2, HOST> &a_Wout,
                          double a_gamma)
              { 
                if (a_W_CME(iRHO) != 0){
                  double p = a_Wout(iP);
                  double J = a_JU(iRHO)/a_Wout(iRHO);
                  a_JU(iRHO) += a_JU_CME(iRHO); 
                  a_JU(iMOMX) = a_JU(iRHO)*a_JU_CME(iMOMX)/a_JU_CME(iRHO);
                  a_JU(iMOMY) = a_JU(iRHO)*a_JU_CME(iMOMY)/a_JU_CME(iRHO);
                  a_JU(iMOMZ) = a_JU(iRHO)*a_JU_CME(iMOMZ)/a_JU_CME(iRHO);
                  a_JU(iBX) = a_JU_CME(iBX);
                  a_JU(iBY) = a_JU_CME(iBY);
                  a_JU(iBZ) = a_JU_CME(iBZ);
                  // e = p_from_JU/(gamma-1) + rho*(v^2)/2 + B^2/8pi
                  a_JU(iE) = J*p/(a_gamma-1) + 0.5*(a_JU(iMOMX)*a_JU(iMOMX) + a_JU(iMOMY)*a_JU(iMOMY) + a_JU(iMOMZ)*a_JU(iMOMZ))/a_JU(iRHO) + (a_JU(iBX)*a_JU(iBX) + a_JU(iBY)*a_JU(iBY) + a_JU(iBZ)*a_JU(iBZ))/J/8/c_PI;
                }
              },JU_i, JU_CME, W_CME, Wout[dit],gamma);
            }
          }
        }
      MBLevelRK4<BoxOp_EulerCubedSphere, MBMap_CubedSphereShell, double> rk4(map, iop);
    
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
        // auto map = CubedSphereShell::Map(JUTemp);
        // h5.writeMBLevel({}, map, JUTemp, "USphere");
      }
    
  //#if 0 // Begin debug comment.
    bool give_space_in_probe_file = true;
    double probe_cadence_temp = 0;
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
    if (lev == 0)
      {
        double dtcfl1 = OP::dtCFL(JU,iop,dVolrLev);
        // if (procID() == 0) cout << "dt_{CFL=1} = " << dtcfl1;
        dt = dtcfl1*ParseInputs::get_CFL();
        // if (procID() == 0) cout << " ,input CFl dt = " << dt << endl;
      }
    if (convTestType > 2) max_iter = 1;
    
    for (int iter = restart_step + 1; iter <= max_iter; iter++)
    {
      auto start = chrono::steady_clock::now();
      if (convTestType == 0)
        {
          double dtcfl1 = OP::dtCFL(JU,iop,dVolrLev);
          // if (procID() == 0) cout << "dt_{CFL=1} = " << dtcfl1;
          dt = dtcfl1*ParseInputs::get_CFL();
          // if (procID() == 0) cout << " ,input CFl dt = " << dt << endl;
        }
      if (convTestType < 3)
        {
          rk4.advance(JU, dVolrLev, dt, time, temporal_order);
          time += dt;
        }
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
      if (iter % write_cadence == 0)
        {
          Write_W(JU, eulerOp, iop, iter, time, dt);
          // Check conservation.
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
        }
      // Checkpointing.
      if ((iter % checkpoint_cadence == 0) || (iter == max_iter))
        {
          std::string check_file_name = ParseInputs::get_checkpoint_file_prefix() + "_" + std::to_string(iter);
          h5.setTime(time);
		      h5.setTimestep(dt);
          h5.setIter(iter);
          h5.writeMBLevel(JU, check_file_name);
          if (procID() == 0) cout << "Checkpointed at iter " << iter << endl;
          int iter_to_delete = iter - (checkpoint_cadence*ParseInputs::get_max_checkpoint_files());
          if (iter_to_delete > restart_step)
          {
            std::string filename_to_delete=ParseInputs::get_checkpoint_file_prefix()+"_"+std::to_string(iter_to_delete)+".hdf5";
            const char* str = filename_to_delete.c_str();
            if (procID() == 0) std::remove(str);
          }
        }
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

      if ((convTestType > 0) && (convTestType != 4)) {
      U_conv_test[lev].define(layout, {Point::Zeros(),Point::Zeros(),Point::Zeros(),Point::Zeros()});
      for (auto dit : layout)
      {
        JU[dit].copyTo(U_conv_test[lev][dit]);
      }
      h5.writeMBLevel({}, map, U_conv_test[lev], "U_conv_test_" + to_string(lev));
      int maxRadSize = 64;
      Point refRatio = 2*Point::Ones();
      if (radial_refinement)
        {
          refRatio = 2*Point::Ones() - Point::Basis(1) - Point::Basis(2);
          maxRadSize = 1024;
        }
      domainSize *= refRatio[1];
      boxSize_nonrad *= refRatio[1];
      thickness *= 2;
      boxSize_rad *= 2;
      time = 0.;
      domain = domain.refine(refRatio); 
      if (boxSize_nonrad > maxRadSize)
        {
          boxSize_nonrad = maxRadSize;
          if (procID() == 0) cout << "radial Box Size = " << boxSize_nonrad << endl;
        }
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
