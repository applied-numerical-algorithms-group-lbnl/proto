template<typename T,MemType MEM>
PROTO_KERNEL_START
void del2_W_c_calcF(Var<T,NUMCOMPS,MEM>& del2_W_c,
                    const Var<T,NUMCOMPS,MEM>& a_W_ave,
                    const Var<T,NUMCOMPS,MEM>& a_W_ave_behind,
                    const Var<T,NUMCOMPS,MEM>& a_W_ave_ahead)
{
  for (int i=0; i< NUMCOMPS; i++) {
    del2_W_c(i) = a_W_ave_behind(i) - 2.0*a_W_ave(i) + a_W_ave_ahead(i);
  }
}
PROTO_KERNEL_END(del2_W_c_calcF, del2_W_c_calc)

template<typename T,MemType MEM>
PROTO_KERNEL_START
void del3_W_calcF(Var<T,NUMCOMPS,MEM>& a_del3_W,
                  const Var<T,NUMCOMPS,MEM>& a_del2_W_c,
                  const Var<T,NUMCOMPS,MEM>& a_del2_W_c_behind)
{
  for (int i=0; i< NUMCOMPS; i++) {
    a_del3_W(i) = a_del2_W_c(i) - a_del2_W_c_behind(i);
  }
}
PROTO_KERNEL_END(del3_W_calcF, del3_W_calc)

const bool printlim = false;

template<typename T,MemType MEM>
PROTO_KERNEL_START
void limiter_calcF(
  const Point& a_pt,
  Var<T,NUMCOMPS,MEM>& a_W_ave_L_ahead_limited,
  Var<T,NUMCOMPS,MEM>& a_W_ave_R_limited,
  const Var<T,NUMCOMPS,MEM>& a_W_ave_L_ahead,
  const Var<T,NUMCOMPS,MEM>& a_W_ave_R,
  const Var<T,NUMCOMPS,MEM>& a_W_ave,
  const Var<T,NUMCOMPS,MEM>& a_W_ave_ahead,
  const Var<T,NUMCOMPS,MEM>& a_W_ave_ahead2,
  const Var<T,NUMCOMPS,MEM>& a_W_ave_behind,
  const Var<T,NUMCOMPS,MEM>& a_W_ave_behind2,
  const Var<T,NUMCOMPS,MEM>& a_del2_W_c,
  const Var<T,NUMCOMPS,MEM>& a_del2_W_c_ahead,
  const Var<T,NUMCOMPS,MEM>& a_del2_W_c_behind,
  const Var<T,NUMCOMPS,MEM>& a_del3_W,
  const Var<T,NUMCOMPS,MEM>& a_del3_W_L,
  const Var<T,NUMCOMPS,MEM>& a_del3_W_R,
  const Var<T,NUMCOMPS,MEM>& a_del3_W_behind,
  const Var<T,NUMCOMPS,MEM>& a_del3_W_ahead,
  const Var<T,NUMCOMPS,MEM>& a_del3_W_ahead2,
  int a_dir,
  const double a_dx,
  const double a_dy,
  const double a_dz)
{
  double rhs = 1.0 - 1.0e-12;
  double rhs_test, lhs_test, rhs_test2, lhs_test2,rhs_test_a, lhs_test_a, rhs_rho, lhs_rho;
  double a_del_W_f_m, a_del_W_f_p, a_del2_W_f, a_del3_W_C, a_del2_W_lim, a_rho_i, a_del3_W_min, a_del3_W_max;
  for (int i=0; i< NUMCOMPS; i++)
  {
    a_del_W_f_m = a_W_ave(i) - a_W_ave_R(i);
    a_del_W_f_p = a_W_ave_L_ahead(i) - a_W_ave(i);
    a_del2_W_f = 6.0*(a_W_ave_R(i) - 2.0*a_W_ave(i) + a_W_ave_L_ahead(i));
    a_del3_W_C = (a_del2_W_c_ahead(i) - a_del2_W_c_behind(i))/2.0;
    if ((a_del2_W_c_behind(i) >= 0.0) && (a_del2_W_c(i) >= 0.0) && (a_del2_W_c_ahead(i) >= 0.0) && (a_del2_W_f >= 0.0)) {
      a_del2_W_lim = std::min({std::abs(a_del2_W_f), 1.25*std::abs(a_del2_W_c_behind(i)),
                                1.25*std::abs(a_del2_W_c(i)), 1.25*std::abs(a_del2_W_c_ahead(i))});
    } else if ((a_del2_W_c_behind(i) < 0.0) && (a_del2_W_c(i) < 0.0) && (a_del2_W_c_ahead(i) < 0.0) && (a_del2_W_f < 0.0)) {
      a_del2_W_lim = -1.0*std::min({std::abs(a_del2_W_f), 1.25*std::abs(a_del2_W_c_behind(i)),
                                    1.25*std::abs(a_del2_W_c(i)), 1.25*std::abs(a_del2_W_c_ahead(i))});
    } else {
      a_del2_W_lim = 0.0;
    }
    rhs_rho = 1.0e-12*std::max({std::abs(a_W_ave_behind2(i)), std::abs(a_W_ave_behind(i)), std::abs(a_W_ave(i)), std::abs(a_W_ave_ahead(i)), std::abs(a_W_ave_ahead2(i))});
    lhs_rho = std::abs(a_del2_W_f);
    if (lhs_rho <= rhs_rho) {
      a_rho_i = 0.0;
    } else {
      a_rho_i = a_del2_W_lim/a_del2_W_f;
    }
    a_del3_W_min = std::min({a_del3_W_behind(i), a_del3_W(i), a_del3_W_ahead(i), a_del3_W_ahead2(i)});
    a_del3_W_max = std::max({a_del3_W_behind(i), a_del3_W(i), a_del3_W_ahead(i), a_del3_W_ahead2(i)});

    if (printlim)
    {
      if (i==1 && a_dir == 0){
        if (a_pt[0]  >= 16 && a_pt[0]  <= 25 && a_pt [1] == 45 ) {
          cout << setw(30) << setprecision(18) << a_pt[0]
                << setw(30) << setprecision(18) << a_pt[1]
              << setw(30) << setprecision(18) << a_W_ave(i)
                << setw(30) << setprecision(18) << a_W_ave_L_ahead(i)
                << setw(30) << setprecision(18) << a_W_ave_R(i)
                << setw(30) << setprecision(18) << a_del_W_f_p
                << setw(30) << setprecision(18) << a_del_W_f_m
                << setw(30) << setprecision(18) << a_del2_W_c(i)
                << setw(30) << setprecision(18) << a_del2_W_f
                << setw(30) << setprecision(18) << a_del3_W(i)
                << setw(30) << setprecision(18) << a_del2_W_lim
                << setw(30) << setprecision(18) << a_rho_i
                << setw(30) << setprecision(18) << a_del3_W_min
                << setw(30) << setprecision(18) << a_del3_W_max
                << setw(30) << setprecision(18) << a_del3_W_L(i)
                << setw(30) << setprecision(18) << a_del3_W_R(i)
                << setw(30) << setprecision(18) << a_del3_W_C
              << endl;
        }
      }
    }

    if ((a_del_W_f_m * a_del_W_f_p <= 0.0) || ((a_W_ave(i) - a_W_ave_behind2(i))*(a_W_ave_ahead2(i) - a_W_ave(i)) <= 0.0)) {
      if (a_rho_i < rhs) {
        rhs_test = a_del3_W_max - a_del3_W_min;
        lhs_test = 0.1*std::max({std::abs(a_del3_W_min), std::abs(a_del3_W_max)});
        lhs_test_a = std::min({std::abs(a_del3_W_L(i)),std::abs(a_del3_W_R(i))})*1.1;
        rhs_test_a = std::abs(a_del3_W_C);
        if (lhs_test <= rhs_test){
        //if (!(lhs_test_a > rhs_test_a && a_del3_W_L(i)*a_del3_W_R(i) > 0.)) {
          rhs_test2 = 2.0*std::abs(a_del_W_f_m);
          lhs_test2 = std::abs(a_del_W_f_p);
          if (a_del_W_f_m * a_del_W_f_p < 0.0) {
            a_W_ave_L_ahead_limited(i) = a_W_ave(i) + a_rho_i*a_del_W_f_p;
          } else if (lhs_test2 >= rhs_test2) {
            a_W_ave_L_ahead_limited(i) = a_W_ave(i) + 2.0*(1.0-a_rho_i)*a_del_W_f_m + a_rho_i*a_del_W_f_p;
          }
          rhs_test2 = 2.0*std::abs(a_del_W_f_p);
          lhs_test2 = std::abs(a_del_W_f_m);
          if (a_del_W_f_m * a_del_W_f_p < 0.0) {
            a_W_ave_R_limited(i) = a_W_ave(i) - a_rho_i*a_del_W_f_m;
          } else if (lhs_test2 >= rhs_test2) {
            a_W_ave_R_limited(i) = a_W_ave(i) - 2.0*(1.0-a_rho_i)*a_del_W_f_p - a_rho_i*a_del_W_f_m;
          }

        }
      }
    } else {
      rhs_test2 = 2.0*std::abs(a_del_W_f_m);
      lhs_test2 = std::abs(a_del_W_f_p);
      if (lhs_test2 >= rhs_test2) {
        a_W_ave_L_ahead_limited(i) = a_W_ave(i) + 2.0*a_del_W_f_m;
      }
      rhs_test2 = 2.0*std::abs(a_del_W_f_p);
      lhs_test2 = std::abs(a_del_W_f_m);
      if (lhs_test2 >= rhs_test2) {
        a_W_ave_R_limited(i) = a_W_ave(i) - 2.0*a_del_W_f_p;
      }
    }
  }
}
PROTO_KERNEL_END(limiter_calcF, limiter_calc)


template<typename T,MemType MEM>
PROTO_KERNEL_START
void eta_tilde_d_calcF(Var<double,1>& a_eta_tilde_d,
                        const Var<T,NUMCOMPS,MEM>& a_W_ave,
                        const Var<T,NUMCOMPS,MEM>& a_W_ave_ahead,
                        const Var<T,NUMCOMPS,MEM>& a_W_ave_ahead2,
                        const Var<T,NUMCOMPS,MEM>& a_W_ave_behind,
                        const Var<T,NUMCOMPS,MEM>& a_W_ave_behind2,
                        int a_d)
{
  double z1=0.85, z0=0.75, delta = 0.33, si;
  double p, p_ahead, p_ahead2, p_behind, p_behind2;

  p = a_W_ave(iP) + (a_W_ave(iBX)*a_W_ave(iBX)+a_W_ave(iBY)*a_W_ave(iBY)+a_W_ave(iBZ)*a_W_ave(iBZ))/8.0/M_PI;
  p_ahead = a_W_ave_ahead(iP) + (a_W_ave_ahead(iBX)*a_W_ave_ahead(iBX)+a_W_ave_ahead(iBY)*a_W_ave_ahead(iBY)+a_W_ave_ahead(iBZ)*a_W_ave_ahead(iBZ))/8.0/M_PI;
  p_ahead2 = a_W_ave_ahead2(iP) + (a_W_ave_ahead2(iBX)*a_W_ave_ahead2(iBX)+a_W_ave_ahead2(iBY)*a_W_ave_ahead2(iBY)+a_W_ave_ahead2(iBZ)*a_W_ave_ahead2(iBZ))/8.0/M_PI;
  p_behind = a_W_ave_behind(iP) + (a_W_ave_behind(iBX)*a_W_ave_behind(iBX)+a_W_ave_behind(iBY)*a_W_ave_behind(iBY)+a_W_ave_behind(iBZ)*a_W_ave_behind(iBZ))/8.0/M_PI;
  p_behind2 = a_W_ave_behind2(iP) + (a_W_ave_behind2(iBX)*a_W_ave_behind2(iBX)+a_W_ave_behind2(iBY)*a_W_ave_behind2(iBY)+a_W_ave_behind2(iBZ)*a_W_ave_behind2(iBZ))/8.0/M_PI;

  double arg = std::abs(p_ahead-p_behind)/std::abs(p_ahead2-p_behind2);
  if (arg > z1) {si = 0.0;}
  if (arg < z0) {si = 1.0;}
  if ((arg<=z1) && (arg>=z0)) {si = 1.0 - (arg-z0)/(z1-z0);}

  double lhs1 = a_W_ave_behind(1+a_d)-a_W_ave_ahead(1+a_d);
  double lhs2 = std::abs(p_ahead-p_behind)/std::min({p_ahead,p_behind});
  if ((lhs1 > 0.0) && (lhs2 > delta)) {
    a_eta_tilde_d(0) = si;
  } else {
    a_eta_tilde_d(0) = 0.0;
  }
}
PROTO_KERNEL_END(eta_tilde_d_calcF, eta_tilde_d_calc)

template<typename T,MemType MEM>
PROTO_KERNEL_START
void eta_d_calcF(Var<T,1,MEM>& a_eta_d,
                  const Var<T,1,MEM>& a_eta_tilde_d,
                  const Var<T,1,MEM>& a_eta_tilde_d_ahead,
                  const Var<T,1,MEM>& a_eta_tilde_d_behind)
{
  a_eta_d(0) = std::min({a_eta_tilde_d(0),a_eta_tilde_d_ahead(0),a_eta_tilde_d_behind(0)});
}
PROTO_KERNEL_END(eta_d_calcF, eta_d_calc)

template<typename T,MemType MEM>
PROTO_KERNEL_START
void eta_calcF(Var<T,1,MEM>& a_eta,
                const Var<T,1,MEM>& a_eta_d,
                const Var<T,1,MEM>& a_eta_d_old)
{
  a_eta(0) = std::min({a_eta_d(0),a_eta_d_old(0)});
}
PROTO_KERNEL_END(eta_calcF, eta_calc)

template<typename T,MemType MEM>
PROTO_KERNEL_START
void Flat_calcF(Var<T,NUMCOMPS,MEM>& a_flattened,
                const Var<T,NUMCOMPS,MEM>& a_not_flattened,
                const Var<T,NUMCOMPS,MEM>& a_W_ave,
                const Var<double,1>& a_eta,
                int a_dir)
{
  for (int i=0; i< NUMCOMPS; i++) {
    a_flattened(i) = (1.0-a_eta(0))*a_not_flattened(i) + a_eta(0)*a_W_ave(i);       //This works
  }
}
PROTO_KERNEL_END(Flat_calcF, Flat_calc)

template<typename T,MemType MEM>
PROTO_KERNEL_START
void compare_calcF(const Point& a_pt,
        Var<T,NUMCOMPS,MEM>& a_W_diff,
                const Var<T,NUMCOMPS,MEM>& a_W_ave_R_limited,
                const Var<T,NUMCOMPS,MEM>& a_W_ave)
{
  for (int i=0; i< NUMCOMPS; i++) {
    if (a_W_ave_R_limited(i) != a_W_ave(i)) cout << a_pt[0] << a_pt[1] << endl;
  }
}
PROTO_KERNEL_END(compare_calcF, compare_calc)


template<
  typename T,
  template<MemType> typename MAP,
  MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::Limiters_4O(BoxData<double,NUMCOMPS>& a_W_ave_L_lim_flat,
                  BoxData<double,NUMCOMPS>& a_W_ave_R_lim_flat,
                  BoxData<double,NUMCOMPS>& a_W_ave_L,
                  BoxData<double,NUMCOMPS>& a_W_ave_R,
                  BoxData<double,NUMCOMPS>& a_W_ave,
                  BoxData<double,NUMCOMPS>& a_W_bar,
                  const int a_d,
                  const double a_dx,
                  const double a_dy,
                  const double a_dz) const
{
  static Stencil<double> m_copy;
  m_copy = 1.0*Shift(Point::Zeros());
  //Limiter Starts here
  State W_ave_L_ahead = alias(a_W_ave_L,Point::Basis(a_d)*(-1));
  //State W_ave_L_ahead_limited = alias(W_ave_L_ahead);  // This works, but won't be able to compare diffrences as if W_ave_L_ahead_limited changes, so does W_ave_L_ahead
  State W_ave_L_ahead_limited = m_copy(W_ave_L_ahead);
  //State W_ave_R_limited = alias(a_W_ave_R);
  State W_ave_R_limited = m_copy(a_W_ave_R);
  int limiter_apply = 1;
  int slope_flattening_apply = 1;
  if (limiter_apply == 1) {
    State W_ave_ahead = alias(a_W_ave,Point::Basis(a_d)*(-1));
    State W_ave_ahead2 = alias(a_W_ave,Point::Basis(a_d)*(-2));
    State W_ave_behind = alias(a_W_ave,Point::Basis(a_d)*(1));
    State W_ave_behind2 = alias(a_W_ave,Point::Basis(a_d)*(2));
    State del2_W_c  = forall<double,NUMCOMPS,MEM>( del2_W_c_calc, a_W_ave, W_ave_behind, W_ave_ahead);
    State del2_W_c_ahead = alias(del2_W_c,Point::Basis(a_d)*(-1));
    State del2_W_c_behind = alias(del2_W_c,Point::Basis(a_d)*(1));
    State del3_W = forall<double,NUMCOMPS,MEM>( del3_W_calc, del2_W_c, del2_W_c_behind);
    State del3_W_L = forall<double,NUMCOMPS,MEM>( del3_W_calc, del2_W_c, del2_W_c_behind);
    State del3_W_R = forall<double,NUMCOMPS,MEM>( del3_W_calc, del2_W_c_ahead, del2_W_c);
    State del3_W_behind = alias(del3_W,Point::Basis(a_d)*(1));
    State del3_W_ahead = alias(del3_W,Point::Basis(a_d)*(-1));
    State del3_W_ahead2 = alias(del3_W,Point::Basis(a_d)*(-2));
    forallInPlace_p( limiter_calc,W_ave_L_ahead_limited, W_ave_R_limited, W_ave_L_ahead, a_W_ave_R,
                      a_W_ave,W_ave_ahead, W_ave_ahead2, W_ave_behind, W_ave_behind2, del2_W_c, del2_W_c_ahead, del2_W_c_behind, del3_W,
                      del3_W_L, del3_W_R, del3_W_behind, del3_W_ahead, del3_W_ahead2, a_d, a_dx, a_dy, a_dz);
  }
  Box dbx0 = W_ave_L_ahead_limited.box();
  Box dbx1 = W_ave_L_ahead.box();
  State W_diff(dbx0);
  //forallInPlace_p(compare_calc, W_diff, W_ave_R_limited, a_W_ave_R);
  // Slope flattening starts here
  if (slope_flattening_apply == 1) {
    Scalar eta;
    Scalar eta_old;
    for (int d = 0; d < DIM; d++)
    {
      State W_bar_ahead = alias(a_W_bar,Point::Basis(d)*(-1));
      State W_bar_ahead2 = alias(a_W_bar,Point::Basis(d)*(-2));
      State W_bar_behind = alias(a_W_bar,Point::Basis(d)*(1));
      State W_bar_behind2 = alias(a_W_bar,Point::Basis(d)*(2));
      Scalar eta_tilde_d = forall<double,1,MEM>( eta_tilde_d_calc, a_W_bar, W_bar_ahead, W_bar_ahead2,
                                            W_bar_behind, W_bar_behind2, d);
      Scalar eta_tilde_d_ahead = alias(eta_tilde_d,Point::Basis(d)*(-1));
      Scalar eta_tilde_d_behind = alias(eta_tilde_d,Point::Basis(d)*(1));
      Scalar eta_d = forall<double,1,MEM>( eta_d_calc, eta_tilde_d, eta_tilde_d_ahead, eta_tilde_d_behind);
      if (d>0) {
        eta = forall<double,1,MEM>(eta_calc, eta_d, eta_old);
      } else {
        eta = alias(eta_d);
      }
      eta_old = alias(eta);
    }
    State W_ave_L_ahead_lim_flat = forall<double,NUMCOMPS,MEM>(Flat_calc, W_ave_L_ahead_limited, a_W_ave, eta, a_d);
    a_W_ave_R_lim_flat = forall<double,NUMCOMPS,MEM>(Flat_calc, W_ave_R_limited, a_W_ave, eta, a_d);
    a_W_ave_L_lim_flat = alias(W_ave_L_ahead_lim_flat,Point::Basis(a_d)*(1));
    //Slope flattening ends here
  } else {
    State W_ave_L_ahead_lim_flat = alias(W_ave_L_ahead_limited);
    a_W_ave_R_lim_flat = alias(W_ave_R_limited);
    a_W_ave_L_lim_flat = alias(W_ave_L_ahead_lim_flat,Point::Basis(a_d)*(1));
  }
}

template<typename T,MemType MEM>
PROTO_KERNEL_START
void minmod_calcF(const Point& a_pt,
                  Var<T,NUMCOMPS,MEM>& W_low_ahead_limited,
                  Var<T,NUMCOMPS,MEM>& W_high_limited,
                  const Var<T,NUMCOMPS,MEM>& a_W,
                  const Var<T,NUMCOMPS,MEM>& a_dWL,
                  const Var<T,NUMCOMPS,MEM>& a_dWR,
                  Var<T,DIM,MEM>& a_dx_sph,
                  Var<T,DIM,MEM>& a_dx_sph_behind,
                  Var<T,DIM,MEM>& a_dx_sph_ahead,
                  int a_d)
{
  double c1     = 1.0/(a_dx_sph(a_d) + a_dx_sph_ahead(a_d));
  double c2     = 1.0/(a_dx_sph(a_d) + a_dx_sph_behind(a_d));

  for (int i = 0; i< NUMCOMPS; i++){ 
    double WP     = a_W(i);
    double dWR    = a_dWR(i);
    double dWL    = a_dWL(i);

    // SLOPE  = dx*minmod( dWR*c1, dWL*c2 )

    // double signA = /abs(dWR*c1);
    double minmod;
    double A = dWR*c1;
    double B = dWL*c2;
    if (A*B < 0){
      minmod = 0;
    } else if (abs(A) < abs(B)){
      minmod = A;
    } else {
      minmod = B;
    }

    // double minmod = signA*max(0.0,min(abs(dWR*c1),signA*dWL*c2));
    double SLOPE = a_dx_sph(a_d)*minmod;

    W_low_ahead_limited(i) = WP + SLOPE;
    W_high_limited(i) = WP - SLOPE;

    // Reduce order if density or p become negative at faces
    // if (i==0 || i==4){
    if (i==4){
      if (W_low_ahead_limited(i) < 0) W_low_ahead_limited(i) = 0.;//WP;
      if (W_high_limited(i) < 0) W_high_limited(i) = 0.;//WP;
    }
  }
}
PROTO_KERNEL_END(minmod_calcF, minmod_calc)




template<
  typename T,
  template<MemType> typename MAP,
  MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::Limiters_minmod
(BoxData<double,NUMCOMPS>& a_W_low_lim,
 BoxData<double,NUMCOMPS>& a_W_high_lim,
 BoxData<double,NUMCOMPS>& a_W,
 BoxData<double,DIM>& a_dx_sph,
 const int a_d) const
{
  static Stencil<double> m_difference;
  m_difference = -1.0*Shift(Point::Zeros()) + 1.0*Shift(Point::Basis(a_d)); 
  static Stencil<double> m_copy;
  m_copy = 1.0*Shift(Point::Zeros());

  State W_low_ahead = alias(a_W,Point::Basis(a_d)*(-1));
  State W_low_ahead_limited = m_copy(W_low_ahead);
  State W_high_limited = m_copy(a_W);
  int limiter_apply = 1;
  if (limiter_apply == 1) {
    State dWR = m_difference(a_W);
    State dWL = alias(dWR,Point::Basis(a_d)*(1));
    BoxData<double,DIM> a_dx_sph_behind = alias(a_dx_sph,Point::Basis(a_d)*(1));
    BoxData<double,DIM> a_dx_sph_ahead = alias(a_dx_sph,Point::Basis(a_d)*(-1));
    forallInPlace_p(minmod_calc,W_low_ahead_limited, W_high_limited, a_W, dWL,
                      dWR, a_dx_sph, a_dx_sph_behind, a_dx_sph_ahead, a_d);
  }
  a_W_high_lim = alias(W_high_limited);
  a_W_low_lim = alias(W_low_ahead_limited,Point::Basis(a_d)*(1));		
}
