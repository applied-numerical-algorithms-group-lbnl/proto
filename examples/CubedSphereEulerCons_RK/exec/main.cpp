//#include <gtest/gtest.h>
#include "Proto.H"
#include "InputParser.H"
#include "BoxOp_EulerCubedSphere.H"
#include "MHD_Input_Parsing.H"

Parsefrominputs inputs;

template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_radialInit_F(
                      Point                           a_pt,
                      Var<T,NUMCOMPS,MEM>&            a_W,
                      Var<T>&                         a_radius,
                      T                               a_dxradius,
                      T                               a_gamma,
                      int                             a_nradius)
{
  // Compute spherical initial data.
  T p0 = 1.0;
  T rho0 = 1.0;
  T eps = 0.1;
  T amplitude;
  T arg = (1.0*a_pt[0] + .5 - 1.0*a_nradius/2)*a_dxradius/2.;
  if (abs(arg) < .125)
    {
      amplitude = eps*pow(cos(2*M_PI*arg),16);
    }else
    {
      amplitude = 0.;
    }
  T rho = rho0 + amplitude*rho0;
  T p = p0*pow(rho/rho0,a_gamma);
  T ur =  amplitude*sqrt(a_gamma*p0/rho0)/rho0;
  a_W(iRHO) = rho;
  a_W(iVX) = 0.0;//ur;
  a_W(iVY) = 0.0;
  a_W(iVZ) = 0.0;
  a_W(iP) = p;
  a_W(iBX) = 0.0;
  a_W(iBY) = 0.0;
  a_W(iBZ) = 0.0;
}
PROTO_KERNEL_END(f_radialInit_F, f_radialInit)

// A function to initialize the data

// initialize data
void radialInit(MBLevelBoxData<double,NUMCOMPS,HOST>& a_JU,
                auto& eulerOp, 
                double domainSize,
                double thickness)
{ 
  typedef BoxOp_EulerCubedSphere<double, MBMap_CubedSphereShell, HOST> OP;
  auto layout = a_JU.layout();
  MBLevelBoxData<double, NUMCOMPS, HOST> WPoint(layout, OP::ghost()+Point::Ones());
  WPoint.setVal(0.);
  double dxradius = 1.0/thickness;
  double dx = 1.0/domainSize;
  for (auto dit : a_JU.layout())
  {      
    BoxData<double> radius(layout[dit].grow(6));
    BoxData<double,DIM,HOST> Dr(layout[dit].grow(6));
    BoxData<double,DIM,HOST> adjDr(layout[dit].grow(6));
    BoxData<double,1,HOST> dVolr(layout[dit].grow(6));
    eulerOp[dit].radialMetrics(radius,Dr,adjDr,dVolr,Dr.box(),thickness);
    auto block = layout.block(dit);
    auto& JU_i = a_JU[dit];
    auto& WPoint_i = WPoint[dit];
    BoxData<double,NUMCOMPS,HOST> JUTemp;
    forallInPlace_p(f_radialInit,WPoint_i,radius,dxradius,inputs.gamma,thickness);
    eulerOp[dit].primToCons(JUTemp,WPoint_i,dVolr,inputs.gamma,dx,block);
    JUTemp.copyTo(JU_i);    
  }
}

void Write_W(MBLevelBoxData<double,NUMCOMPS,HOST>& a_JU,
                auto& eulerOp, 
                double domainSize,
                double thickness,
                int iter)
{
  typedef BoxOp_EulerCubedSphere<double, MBMap_CubedSphereShell, HOST> OP;
  HDF5Handler h5;
  int ghostTest = 6;
  auto layout = a_JU.layout();
  MBLevelBoxData<double, NUMCOMPS, HOST> WNew_out(layout, Point::Zeros());
  double dx = 1.0/domainSize;
  auto map = CubedSphereShell::Map(layout,OP::ghost());
  for (auto dit :a_JU.layout())
  {
    BoxData<double> radius(layout[dit].grow(ghostTest));
    BoxData<double,DIM,HOST> Dr(layout[dit].grow(ghostTest));
    BoxData<double,DIM,HOST> adjDr(layout[dit].grow(ghostTest));
    BoxData<double,1,HOST> dVolr(layout[dit].grow(ghostTest));
    eulerOp[dit].radialMetrics(radius,Dr,adjDr,dVolr,Dr.box(),thickness);
    int block_i = layout.block(dit);
    auto& WNew_out_i = WNew_out[dit];
    BoxData<double,NUMCOMPS,HOST> WNewTemp;
    eulerOp[dit].consToPrim(WNewTemp,WNewTemp,a_JU[dit],dVolr,inputs.gamma,dx,block_i);
    WNewTemp.copyTo(WNew_out_i);
  }
  h5.writeMBLevel({ }, map, WNew_out, "W_"+to_string(iter));
}


int main(int argc, char* argv[])
{   
#ifdef PR_MPI
    MPI_Init (&argc, &argv);
#endif
  typedef BoxOp_EulerCubedSphere<double, MBMap_CubedSphereShell, HOST> OP;
  inputs.parsenow(argc, argv);
  #ifdef PR_MPI
	  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  HDF5Handler h5;
  int domainSize = inputs.domainSizey;
  int thickness = inputs.domainSizex;
  int boxSize_nonrad = inputs.BoxSize;
  int boxSize_rad = inputs.BoxSize;

  PR_TIMER_SETFILE(to_string(domainSize) 
                   + "_DIM" + to_string(DIM) //+ "_NProc" + to_string(numProc())
                   + "_CubeSphereTest.time.table");
  PR_TIMERS("MMBEuler");
  
  int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD; 
  auto domain = CubedSphereShell::Domain(domainSize, thickness, radialDir);
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
  auto map = CubedSphereShell::Map(layout,OP::ghost());
  MBLevelBoxData<double, NUMCOMPS, HOST> U(layout, OP::ghost());
  MBLevelBoxData<double, NUMCOMPS, HOST> rhs(layout, Point::Zeros());
  auto eulerOp = CubedSphereShell::Operator<BoxOp_EulerCubedSphere, double, HOST>(map);
  U.setVal(0.);  
  radialInit(U, eulerOp, domainSize, thickness);

  // FIXME: Commenting this out until the interpOp can be built without SVD failing
  //auto interpOp = CubedSphereShell::InterpOp<HOST>(layout, OP::ghost(), 4);
  //MBLevelRK4<BoxOp_EulerCubedSphere, MBMap_CubedSphereShell, double> rk4(map, interpOp);

  
  double dt = inputs.CFL;
  double time = 0.0;
  for (int iter = 0; iter <= inputs.maxStep; iter++)
  {
    if (procID() == 0) cout << "iter = " << iter << endl;
    // Write out the data
    if (iter % inputs.outputInterval == 0) Write_W(U, eulerOp, domainSize, thickness, iter);
    
    // time step
    U.exchange();
    for (auto dit :U.layout())
    {
      PR_TIMERS("RHS Calculation");
      eulerOp[dit](rhs[dit],U[dit],dt);
      U[dit] -= rhs[dit]; 
    }
    time += dt;
  }

  PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
}
