#include "MHD_Constants.H"
template <typename T, MemType MEM>
PROTO_KERNEL_START
void f_radialInit_F(
    Point a_pt,
    Var<T, NUMCOMPS, MEM> &a_W,
    Var<T> &a_radius,
    T a_gamma,
    int a_nradius)
{
  // Compute spherical initial data.
  T p0 = 1.0;
  T rho0 = 1.0;
  T eps = 0.1;
  T amplitude;
  T arg = (1.0 * a_pt[0] + .5 - 1.0 * a_nradius / 2) / a_nradius / 2.;
  if (abs(arg) < .25)
  {
    amplitude = eps * pow(cos(2 * M_PI * arg), 6);
  }
  else
  {
    amplitude = 0.;
  }
  T rho = rho0 + amplitude * rho0;
  T p = p0 * pow(rho / rho0, a_gamma);
  T ur = amplitude * sqrt(a_gamma * p0 / rho0) / rho0;
  a_W(iRHO) = rho;
  a_W(iVX) = 0.0; //ur;
  a_W(iVY) = 0.0;
  a_W(iVZ) = 0.0;
  a_W(iP) = p;
  a_W(iBX) = 0.0;
  a_W(iBY) = 0.0;
  a_W(iBZ) = 0.0;
}

PROTO_KERNEL_END(f_radialInit_F, f_radialInit)

template <typename T, MemType MEM>
PROTO_KERNEL_START
void f_nonradialInit_F(
    Point a_pt,
    Var<T, NUMCOMPS, MEM> &a_W,
    Var<T, DIM> &a_X_cart,
    T a_gamma,
    int a_nradius)
{
  // Compute spherical initial data.
  T p0 = 1.0;
  T rho0 = 1.0;
  T eps = 0.1;
  T amplitude;
  T arg = sqrt((a_X_cart(2)+0.7) * (a_X_cart(2)+0.7) + a_X_cart(0) * a_X_cart(0) + a_X_cart(1) * a_X_cart(1));
  if (abs(arg) < .25)
  {
    amplitude = eps * pow(cos(2 * M_PI * arg), 6);
  }
  else
  {
    amplitude = 0.;
  }
  T rho = rho0 + amplitude * rho0;
  T p = p0 * pow(rho / rho0, a_gamma);
  T ur = amplitude * sqrt(a_gamma * p0 / rho0) / rho0;
  a_W(iRHO) = rho;
  a_W(iVX) = 0.;//ur;
  a_W(iVY) = 0.0;
  a_W(iVZ) = 0.0;
  a_W(iP) = p;
  a_W(iBX) = 0.0;
  a_W(iBY) = 0.0;
  a_W(iBZ) = 0.0;
}
PROTO_KERNEL_END(f_nonradialInit_F, f_nonradialInit)

template <typename T, MemType MEM>
PROTO_KERNEL_START
void f_radialoutflowInit_F(
    Point a_pt,
    Var<T, NUMCOMPS, MEM> &a_W,
    Var<T> &a_radius,
    Var<T, DIM> &a_X_cart,
    T a_gamma,
    int a_nradius)
{
  // Compute spherical initial data.
  double r_BC = CUBED_SPHERE_SHELL_R0; //in cm
  T p0 = 1.0e-7; // p near 21.5 c_SR is about 1e-7 dyne/cm2
  T rho0 = 700*c_MP; // rho at 21.5 c_SR is about 700/cm3
  T u0 = 500.0*1e5; // v at 21.5 c_SR is about 500 km/s
  T rho = rho0*pow(a_radius(0)/r_BC, -2.0);
  T p = p0 * pow(rho / rho0, a_gamma);
  T ur = u0;
  a_W(iRHO) = rho;
  a_W(iVX) = ur;
  a_W(iVY) = 0.0;
  a_W(iVZ) = 0.0;
  a_W(iP) = p;
  a_W(iBX) = 0.0;
  a_W(iBY) = 0.0;
  a_W(iBZ) = 0.0;
}

PROTO_KERNEL_END(f_radialoutflowInit_F, f_radialoutflowInit)

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void BoxOp_EulerCubedSphere<T,MAP,MEM>::
initialize(State& a_W,
          const State& a_BC,
          const Scalar& a_radius,
          const Vector& a_Xcart,
          T a_gamma,
          int a_nradius) const
{
  int init_condition_type = ParseInputs::get_init_condition_type();
  if (init_condition_type == 0) forallInPlace_p(f_radialInit, a_W, a_radius, a_gamma, a_nradius);
  if (init_condition_type == 1) forallInPlace_p(f_nonradialInit, a_W, a_Xcart, a_gamma, a_nradius);
  if (init_condition_type == 2 || init_condition_type == 3) forallInPlace_p(f_radialoutflowInit, a_W, a_radius, a_Xcart, a_gamma, a_nradius);


  // The follwoing breaks if decompostion is there in radial direction. That's because BCs are available at -1.
  // if (init_condition_type == 3) {
  //   Point source_lo = Point(a_BC.box().low()[0],a_BC.box().low()[1],a_BC.box().low()[2]);
  //   Point source_hi = Point(a_BC.box().low()[0],a_BC.box().high()[1],a_BC.box().high()[2]);
  //   Box sourceBox(source_lo,source_hi);
  //   int low_index = a_W.box().low()[0];
  //   int high_index = a_W.box().high()[0];
  //   for (int i = low_index; i <= high_index; i++){
  //     a_BC.copyTo(a_W,sourceBox,Point::Basis(0)*(i+1));
  //   }
  //   // In primitive spherical. To correct units. scale by radius.
  //   forallInPlace([ ] PROTO_LAMBDA
  //   (Var<double, NUMCOMPS, HOST> &a_WSph,
  //   Var<double> &a_radius,
  //   double a_gamma)
  //   {  
  //     double r_BC = CUBED_SPHERE_SHELL_R0; //in cm
  //     double rad = a_radius(0);
  //     double rho = a_WSph(0)*c_MP*pow(r_BC/rad,2.0);
  //     double u = a_WSph(1);
  //     double v = a_WSph(2);
  //     double w = a_WSph(3);
  //     double p = a_WSph(4)*1.0e-12*pow(r_BC/rad,2.0*a_gamma); // From picodyne to dyne 
  //     double Bx = a_WSph(5)*1.0e-6*pow(r_BC/rad,2.0); // MicroGauss to Gauss
  //     double By = a_WSph(6)*1.0e-6*pow(r_BC/rad,1.0);
  //     double Bz = a_WSph(7)*1.0e-6*pow(r_BC/rad,1.0);
  //     a_WSph(0) = rho; //rho
  //     a_WSph(1) = u; //V-x
  //     a_WSph(2) = v; //V-y
  //     a_WSph(3) = w; //V-z
  //     a_WSph(4) = p; //Pressure
  //     a_WSph(5) = Bx; //Bx
  //     a_WSph(6) = By; //By
  //     a_WSph(7) = Bz; //Bz
  //   },a_W,a_radius,a_gamma);
  // }
}