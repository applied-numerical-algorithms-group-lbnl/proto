//#include <gtest/gtest.h>
#include "Proto.H"
#include "InputParser.H"
#include "BoxOp_EulerCubedSphere.H"

 
void GetCmdLineArgumenti(int argc, const char** argv, const char* name, int* rtn)
{
    size_t len = strlen(name);
    for(int i=1; i<argc; i+=2)
    {
        if(strcmp(argv[i]+1,name) ==0)
        {
          *rtn = atoi(argv[i+1]);
            // std::cout<<name<<"="<<" "<<*rtn<<std::endl;
            break;
        }
    }
}
#if 0
template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_initPatchData_(
                      Point                           a_pt,
                      Var<T,NUMCOMPS,MEM>&            a_U,
                      Var<T>&                         a_radius,
                      Var<T,DIM,MEM,DIM>&             a_A,
                      T                               a_dx,
                      T                               a_dxi0,
                      T                               a_gamma,
                      int                             a_nradius)
{
  // Compute Cartesian conserved quantitities from spherical data.
  T p0 = a_gamma;
  T rho0 = 1.0;
  T eps = .01;
  T amplitude;
  T arg = (a_pt[0] - a_nradius/2)/(a_nradius*1.0);
  if (abs(arg) < .5)
    {
      amplitude = eps*pow(cos(M_PI*arg),6);
    }else
    {
      amplitude = 0.;
    }
  T rho = rho0*(1.0 + amplitude);
  T p = p0*pow(rho/rho0,a_gamma);
  Array<T,DIM> vSpherical = {1.,0.,0.};
  vSpherical[0] *= amplitude;
  Array<T,DIM> vcart;
  T ke;
  for (int dim1 = 0; dim1 < DIM; dim1++)
    {
      vcart[dim1] = 0.;
      // Can set velocities here.
      for (int dim2 = 0; dim2 < DIM; dim2++)
        {
          vcart[dim1] += vSpherical[dim2]*a_A(dim1,dim2);
        }
      a_U(dim1 + 1) = vcart[dim1]*rho;
      ke += .5*vcart[dim1]*vcart[dim1];
    }
  a_U(0) = rho;
  a_U(NUMCOMPS-1) = p/(a_gamma-1) + rho*ke;
}
PROTO_KERNEL_END(f_initPatchData_, f_initPatchData)
#endif

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
  T rho = rho0 + amplitude*rho0;//rho0*pow(1.0  - eps*arg,2);
  T p = p0*pow(rho/rho0,a_gamma);
  T ur =  amplitude*sqrt(a_gamma*p0/rho0)/rho0;
  //pow(1.0 - eps*arg,4)*sqrt(a_gamma*p/rho0)/rho0;
  a_W(0) = rho;
  a_W(1) = 0.0;//ur;
  a_W(2) = 0.0;
  a_W(3) = 0.0;
  a_W(NUMCOMPS-1) = p;
}
PROTO_KERNEL_END(f_radialInit_F, f_radialInit)

int main(int argc, char* argv[])
{   
#ifdef PR_MPI
    MPI_Init (&argc, &argv);
#endif
  HDF5Handler h5;
  int domainSize = 32;
  int thickness = 16;
  int boxSize_nonrad = 16;
  int boxSize_rad = 16;
  InputArgs args;
  args.add("nsph", domainSize);
  args.add("nrad",  thickness);
  args.parse(argc, argv);
  // args.print();
  Array<double,DIM> offset = {0.,0.,0.};
  Array<double,DIM> exp = {1.,1.,1.};
  PR_TIMER_SETFILE(to_string(domainSize) 
                   + "_DIM" + to_string(DIM) //+ "_NProc" + to_string(numProc())
                   + "_CubeSphereTest.time.table");
  PR_TIMERS("MMBEuler");
  
  bool cullRadialGhost = true;
  bool use2DFootprint = true;
  int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD;
  Array<Point, DIM+1> ghost;
  ghost.fill(Point::Ones(NGHOST));
   Array<Array<uint,DIM>,6> permute = {{2,1,0},{2,1,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
  Array<Array<int,DIM>,6> sign = {{-1,1,1},{1,1,-1},{-1,1,1},{1,1,1},{1,-1,1},{-1,-1,1}};     
  auto domain =
  CubedSphereShell::Domain(domainSize, thickness, radialDir);
  Point boxSizeVect = Point::Ones(boxSize_nonrad);
  boxSizeVect[radialDir] = boxSize_rad;
  MBDisjointBoxLayout layout(domain, boxSizeVect);

  // initialize data and map
  auto map = CubedSphereShell::Map(layout,ghost);
  MBLevelBoxData<double, NUMCOMPS, HOST> U(layout, ghost);
  MBLevelBoxData<double, NUMCOMPS, HOST> rhs(layout, Point::Zeros());
 
  // MBLevelRK4<BoxOp_EulerCubedSphere, MBMap_CubedSphereShell, double> rk4(map);

  auto eulerOp = CubedSphereShell::Operator<BoxOp_EulerCubedSphere, double, HOST>(map);
  U.setVal(0.);  
  double dx = 1.0/domainSize;
  double dxradius = 1.0/thickness;
  auto C2C = Stencil<double>::CornersToCells(4);
  MBLevelBoxData<double, NUMCOMPS, HOST> WPoint(layout, ghost+Point::Ones());
  WPoint.setVal(0.);


  MBLevelBoxData<double, NUMCOMPS, HOST> WNew_out(layout, Point::Zeros());


  for (auto dit : U.layout())
    {      
      BoxData<double> radius(layout[dit].grow(6));
      double gamma = 5.0/3.0;
      BoxData<double,DIM,HOST> Dr(layout[dit].grow(6));
      BoxData<double,DIM,HOST> adjDr(layout[dit].grow(6));
      BoxData<double,1,HOST> dVolr(layout[dit].grow(6));
      eulerOp[dit].radialMetrics(radius,Dr,adjDr,dVolr,Dr.box(),thickness);
      auto block = layout.block(dit);
      auto& JU_i = U[dit];
      auto& WPoint_i = WPoint[dit];
      BoxData<double,NUMCOMPS,HOST> JUTemp;
      forallInPlace_p(f_radialInit,WPoint_i,radius,dxradius,gamma,thickness);
      eulerOp[dit].primToCons(JUTemp,WPoint_i,dVolr,gamma,dx,block);
      JUTemp.copyTo(JU_i);    
    }


  int ghostTest = 6;
  double dt = 1e-3;
  double time = 0.0;
  for (int iter = 0; iter <= 2; iter++)
  {
    if (procID() == 0) cout << "iter = " << iter << endl;
    U.exchange();
    if (iter % 1 == 0)
    {
      for (auto dit :U.layout())
      {
        BoxData<double> radius(layout[dit].grow(ghostTest));
        double gamma = 5.0/3.0;
        BoxData<double,DIM,HOST> Dr(layout[dit].grow(ghostTest));
        BoxData<double,DIM,HOST> adjDr(layout[dit].grow(ghostTest));
        BoxData<double,1,HOST> dVolr(layout[dit].grow(ghostTest));
        eulerOp[dit].radialMetrics(radius,Dr,adjDr,dVolr,Dr.box(),thickness);
        int block_i = layout.block(dit);
        auto& WNew_out_i = WNew_out[dit];
        BoxData<double,NUMCOMPS,HOST> WNewTemp;
        eulerOp[dit].consToPrim(WNewTemp,WNewTemp,U[dit],dVolr,gamma,dx,block_i);
        WNewTemp.copyTo(WNew_out_i);
      }
      h5.writeMBLevel({ }, map, WNew_out, "W_"+to_string(iter));
    }
    

    for (auto dit :U.layout())
    {
      PR_TIMERS("RHS Calculation");
      auto& rhs_i = rhs[dit];
      auto& U_i = U[dit];
      double gamma = 5.0/3.0;
      Array<BoxData<double,NUMCOMPS>, DIM> fluxes;
      fluxes[0].define(rhs_i.box().extrude(0));
      fluxes[1].define(rhs_i.box().extrude(1));
      fluxes[2].define(rhs_i.box().extrude(2));
      eulerOp[dit](rhs_i,fluxes,U_i,1.0);
      rhs_i*=dt;
      U[dit]-=rhs_i; 
    }
    time += dt;
  }

  PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
}
