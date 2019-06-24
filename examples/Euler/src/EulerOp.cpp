#include "Proto.H"
#include "EulerOp.H"
#include "Proto_Timer.H"


//double EulerOp::s_gamma = 1.4;
//double EulerOp::s_dx = 1.0;

typedef BoxData<double,1,1,1> Scalar;
typedef BoxData<double,NUMCOMPS,1,1> Vector;

namespace EulerOp {

  PROTO_KERNEL_START
  void 
  consToPrimF(State&         a_W, 
              const State&   a_U,
              double         a_gamma)
  {
    double rho = a_U(0);
    double v2 = 0.0;
    double gamma = a_gamma;
    a_W(0) = rho;
    
    for (int i = 1; i <= DIM; i++)
      {
        double v;
        v = a_U(i) / rho;
        
        a_W(i) = v;
        v2 += v*v;
      }
    
    a_W(NUMCOMPS-1) = (a_U(NUMCOMPS-1) - .5 * rho * v2) * (gamma - 1.0);
  }
  PROTO_KERNEL_END(consToPrimF, consToPrim)

  PROTO_KERNEL_START
  void upwindStateF(State& a_out,
                    const State& a_low,
                    const State& a_high,
                    int   a_dir,
                    double a_gamma)
  {
    const double& rhol = a_low(0);
    const double& rhor = a_high(0);
    const double& ul = a_low(a_dir+1);
    const double& ur = a_high(a_dir+1);
    const double& pl = a_low(NUMCOMPS-1);
    const double& pr = a_high(NUMCOMPS-1);
    double gamma = a_gamma;
    double rhobar = (rhol + rhor)*.5;
    double pbar = (pl + pr)*.5;
    double ubar = (ul + ur)*.5;
    double cbar = sqrt(gamma*pbar/rhobar);
    double pstar = (pl + pr)*.5 + rhobar*cbar*(ul - ur)*.5;
    double ustar = (ul + ur)*.5 + (pl - pr)/(2*rhobar*cbar);
    int sign;
    if (ustar > 0) 
      {
        sign = -1;
        for (int icomp = 0;icomp < NUMCOMPS;icomp++)
          {
            a_out(icomp) = a_low(icomp);
          }
      }
    else
      {
        sign = 1;
        for (int icomp = 0;icomp < NUMCOMPS;icomp++)
          {
            a_out(icomp) = a_high(icomp);
          }
      }
    if (cbar + sign*ubar > 0)
      {
        a_out(0) += (pstar - a_out(NUMCOMPS-1))/(cbar*cbar);
        a_out(a_dir+1) = ustar;
        a_out(NUMCOMPS-1) = pstar;
      }
  }
  PROTO_KERNEL_END(upwindStateF, upwindState)

  PROTO_KERNEL_START
  void getFluxF(State& a_F, const State& a_W, 
                int    a_d,
                double a_gamma)
  {
    double F0 = a_W(a_d+1)*a_W(0);
    double W2 = 0.0;
    double gamma = a_gamma;

    a_F(0) = F0;

    for (int d = 1; d <= DIM; d++)
      {
        double Wd = a_W(d);

        a_F(d) = Wd*F0;
        W2 += Wd*Wd;
      }

    a_F(a_d+1) += a_W(NUMCOMPS-1);
    a_F(NUMCOMPS-1) = gamma/(gamma - 1.0) * a_W(a_d+1) * a_W(NUMCOMPS-1) + 0.5 * F0 * W2;
  }
  PROTO_KERNEL_END(getFluxF, getFlux)

  PROTO_KERNEL_START
  void waveSpeedBoundF(Var<double,1>& a_speed,
                       const State& a_W,
                       double       a_gamma)
  {
    double gamma = a_gamma;
    a_speed(0) = DIM*sqrt(gamma*a_W(NUMCOMPS-1)/a_W(0));
    for (int dir = 1 ; dir <= DIM; dir++)
      {
        a_speed(0) += a_W(dir);
      }
  }
  PROTO_KERNEL_END(waveSpeedBoundF, waveSpeedBound)



  double step(BoxData<double,NUMCOMPS>& a_Rhs,
                             const BoxData<double,NUMCOMPS>& a_U,
                             const Box& a_rangeBox)
  {
    static Stencil<double> m_laplacian;
    static Stencil<double> m_deconvolve;
    static Stencil<double> m_laplacian_f[DIM];
    static Stencil<double> m_deconvolve_f[DIM];
    static Stencil<double> m_interp_H[DIM];
    static Stencil<double> m_interp_L[DIM];
    static Stencil<double> m_divergence[DIM];
    static bool initialized = false;
    if(!initialized)
    {
      m_laplacian = Stencil<double>::Laplacian();
      m_deconvolve = (-1.0/24.0)*m_laplacian + (1.0)*Shift(Point::Zeros());
      for (int dir = 0; dir < DIM; dir++)
      {
        m_laplacian_f[dir] = Stencil<double>::LaplacianFace(dir);
        m_deconvolve_f[dir] = (-1.0/24.0)*m_laplacian_f[dir] + 1.0*Shift(Point::Zeros());
        m_interp_H[dir] = Stencil<double>::CellToEdgeH(dir);
        m_interp_L[dir] = Stencil<double>::CellToEdgeL(dir);
        m_divergence[dir] = Stencil<double>::FluxDivergence(dir);
      }
      initialized =  true;
    }


    using namespace std;
    PR_TIME("EulerOp::operator");
    a_Rhs.setVal(0.0);

    double gamma = s_gamma;
    double retval;
    //PR_TIME("EulerOp::operator::W_bar");
    Vector W_bar = forall<double,NUMCOMPS>(consToPrim,a_U, gamma);
    //PR_TIME("EulerOp::operator::U");
    Vector U = m_deconvolve(a_U);
    //PR_TIME("EulerOp::operator::W");
    Vector W = forall<double,NUMCOMPS>(consToPrim,U, gamma);
    Scalar umax = forall<double>(waveSpeedBound,a_rangeBox,W, gamma);
    retval = umax.absMax();
    //PR_TIME("EulerOp::operator::W_ave");
    Vector W_ave = m_laplacian(W_bar,1.0/24.0);
    W_ave += W;
    for (int d = 0; d < DIM; d++)
      {
        //PR_TIME("EulerOp::operator::W_ave_f::interpolation");
        Vector W_ave_low = m_interp_L[d](W_ave);
        Vector W_ave_high = m_interp_H[d](W_ave);
        //PR_TIME("EulerOp::operator::W_ave_f::upwinding");
        Vector W_ave_f = forall<double,NUMCOMPS>(
                                                 upwindState,W_ave_low, W_ave_high,d,  gamma);
#if DIM>1
        //PR_TIME("EulerOp::operator::F_bar_f");
        Vector F_bar_f = forall<double,NUMCOMPS>(getFlux, W_ave_f, d,  gamma);
#endif
        //PR_TIME("EulerOp::operator::W_f");
#if DIM>1
        Vector W_f = m_deconvolve_f[d](W_ave_f);
#else
        Vector W_f = W_ave_f;
#endif
        //PR_TIME("EulerOp::operator::F_ave");
        Vector F_ave_f = forall<double,NUMCOMPS>(getFlux, W_f, d, gamma);
#if DIM>1
        // F_bar_f *= (1./24.);
        F_ave_f += m_laplacian_f[d](F_bar_f,1.0/24.0);
#endif
        //PR_TIME("EulerOp::operator::minusDivF");
        a_Rhs += m_divergence[d](F_ave_f);
      }
    //PR_TIME("EulerOp::operator::RHS*=-1.0/dx");
    a_Rhs *= -1./s_dx;
    return retval;
  }
}
