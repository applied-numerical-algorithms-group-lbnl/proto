#include "MHD_Constants.H"
template <typename T, MemType MEM>
PROTO_KERNEL_START
void f_radialInit_F(
    Point a_pt,
    Var<T, NUMCOMPS, MEM> &a_W,
    Var<T> &a_radius,
    T a_gamma,
    T a_strength,
    int a_nradius)
{
  // Compute spherical initial data.
  T p0 = 1.0e-7; // p near 21.5 c_SR is about 1e-7 dyne/cm2
  T rho0 = 700*c_MP; // rho at 21.5 c_SR is about 700/cm3
  //T u0 = 500.0*1e5; // v at 21.5 c_SR is about 500 km/s
  T amplitude;
  T arg = (1.0 * a_pt[0] + .5 - 1.0 * a_nradius / 2) / a_nradius / 2.;
  if (arg < 0)
  {
    amplitude = a_strength * pow(cos(2 * M_PI * arg), 6);
  }
  else
  {
    amplitude = 0.;
  }
  T p = p0 + amplitude * p0;
  T W = sqrt(a_gamma*p0*rho0*(1.0 + ((a_gamma + 1.0)/(2.0*a_gamma))*(p-p0)/p0));
  T rho = 1.0/(1./rho0 - (p - p0)/(W*W));
  T ur = (p-p0)/W;
  a_W(iRHO) = rho;
  a_W(iVX) = ur;
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
  T p0 = 1.0e-7; // p near 21.5 c_SR is about 1e-7 dyne/cm2
  T rho0 = 700*c_MP; // rho at 21.5 c_SR is about 700/cm3
  T u0 = 500.0*1e5; // v at 21.5 c_SR is about 500 km/s
  T eps = 0.1;
  T amplitude;
  T arg = sqrt((a_X_cart(2)+0.7) * (a_X_cart(2)+0.7) + a_X_cart(0) * a_X_cart(0) + a_X_cart(1) * a_X_cart(1));
  // T arg = (1.0 * a_pt[0] + .5 - 1.0 * a_nradius / 2) / a_nradius / 2.;
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
          int a_nradius,
          T a_h,
          uint a_block) const
{
  int init_condition_type = ParseInputs::get_init_condition_type();
  T strength = ParseInputs::get_strength();
  if (init_condition_type == 0) forallInPlace_p(f_radialInit, a_W, a_radius, a_gamma, strength, a_nradius);
  if (init_condition_type == 1) forallInPlace_p(f_nonradialInit, a_W, a_Xcart, a_gamma, a_nradius);
  if (init_condition_type == 2) forallInPlace_p(f_radialoutflowInit, a_W, a_radius, a_Xcart, a_gamma, a_nradius);
  // if (init_condition_type == 3) forallInPlace_p(f_radialoutflowInit, a_W, a_radius, a_Xcart, a_gamma, a_nradius);


  // The following breaks if decompostion is there in radial direction. That's because BCs are available at -1.
  if (init_condition_type == 3) {
    Point source_lo = Point(a_BC.box().low()[0],a_BC.box().low()[1],a_BC.box().low()[2]);
    Point source_hi = Point(a_BC.box().low()[0],a_BC.box().high()[1],a_BC.box().high()[2]);
    Box sourceBox(source_lo,source_hi);
    int low_index = a_W.box().low()[0];
    int high_index = a_W.box().high()[0];
    for (int i = low_index; i <= high_index; i++){
      a_BC.copyTo(a_W,sourceBox,Point::Basis(0)*(i+1));
    }
    // In primitive spherical. To correct units. scale by radius.
    forallInPlace([ ] PROTO_LAMBDA
    (Var<double, NUMCOMPS, HOST> &a_WSph,
    Var<double> &a_radius,
    double a_gamma)
    {  
      double r_BC = CUBED_SPHERE_SHELL_R0; //in cm
      double rad = a_radius(0);
      double rho = a_WSph(0)*c_MP*pow(r_BC/rad,2.0);
      double u = a_WSph(1);
      double v = a_WSph(2);
      double w = a_WSph(3);
      double p = a_WSph(4)*1.0e-12*pow(r_BC/rad,2.0*a_gamma); // From picodyne to dyne 
      double Bx = a_WSph(5)*1.0e-6*pow(r_BC/rad,2.0); // MicroGauss to Gauss
      double By = a_WSph(6)*1.0e-6*pow(r_BC/rad,1.0);
      double Bz = a_WSph(7)*1.0e-6*pow(r_BC/rad,1.0);
      a_WSph(0) = rho; //rho
      a_WSph(1) = u; //V-x
      a_WSph(2) = v; //V-y
      a_WSph(3) = w; //V-z
      a_WSph(4) = p; //Pressure
      a_WSph(5) = Bx; //Bx
      a_WSph(6) = By; //By
      a_WSph(7) = Bz; //Bz
    },a_W,a_radius,a_gamma);
  }

  if (init_condition_type == 4){
    forallInPlace_p(f_radialInit, a_W, a_radius, a_gamma, strength, a_nradius);

    // Calculate A_matrix
    Box bx = a_W.box();
    double half = 0.5;
    double dxiPerp = a_h;
    Array<Array<uint,DIM>,6> permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
    Array<Array<int,DIM>,6> sign = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}}; 
    Point high = bx.high();
    Point low = bx.low();
    high[0] = low[0];
    Box bx0(low,high);
    BoxData<T,DIM,MEM,DIM> A_matrix(bx);
    T offseta = half;
    T offsetb = half;
    forallInPlace_p(f_Amatrix,bx0,A_matrix,permute[a_block],sign[a_block],
                    dxiPerp,offseta,offsetb);
    spreadSlice(A_matrix);

    BoxData<T,DIM,MEM,DIM> A_matrix_inv(bx);
    A_matrix_inv.setToZero();
    forallInPlace(f_matinv3by3,A_matrix_inv,A_matrix);

    // A constant U_cart
    BoxData<double, DIM> U_cart(bx);
    BoxData<double, DIM> W_sph(bx);
    forallInPlace_p([ ] PROTO_LAMBDA
    (Point a_pt,
     Var<double, DIM, HOST> &a_U_cart,
     double a_a_gamma,
     int a_a_nradius)
    {  
      double eps = 0.1;
      double amplitude;
      double arg = (1.0 * a_pt[0] + .5 - 1.0 * a_a_nradius / 2) / a_a_nradius / 2.;
      if (abs(arg) < .25)
      {
        amplitude = eps * pow(cos(2 * M_PI * arg), 6);
      }
      else
      {
        amplitude = 0.;
      }
      double p0 = 1.0e-7; // p near 21.5 c_SR is about 1e-7 dyne/cm2
      double rho0 = 700*c_MP; // rho at 21.5 c_SR is about 700/cm3
      a_U_cart(0) = 0.1*amplitude *sqrt(a_a_gamma * p0 / rho0);
      a_U_cart(1) = 0.0;
      a_U_cart(2) = 0.0;
    },U_cart, a_gamma, a_nradius);

    T one = 1.0;  
    forallInPlace(f_matVecProd,W_sph,A_matrix_inv,U_cart,one);

    forallInPlace([ ] PROTO_LAMBDA
    (Var<double, NUMCOMPS, HOST> &a_a_W,
     Var<double, DIM, HOST> &a_W_sph)
    {  
      a_a_W(iVX) += a_W_sph(0);
      a_a_W(iVY) += a_W_sph(1);
      a_a_W(iVZ) += a_W_sph(2);
    },a_W, W_sph);
  }

  
}
