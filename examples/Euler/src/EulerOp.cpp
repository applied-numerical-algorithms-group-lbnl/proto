#include "Proto.H"
#include "EulerOp.H"
#include "Proto_Timer.H"

//#define DATAFLOW_ON 0
#if DATAFLOW_ON > 0
#include "Proto_DataFlow.H"
#endif

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

  Stencil<double> m_laplacian;
  Stencil<double> m_deconvolve;
  Stencil<double> m_laplacian_f[DIM];
  Stencil<double> m_deconvolve_f[DIM];
  Stencil<double> m_interp_H[DIM];
  Stencil<double> m_interp_L[DIM];
  Stencil<double> m_divergence[DIM];

  bool s_init() // init function in place of a constructor
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
    return true;
  }

  double step(BoxData<double,NUMCOMPS>& a_Rhs,
              const BoxData<double,NUMCOMPS>& a_U,
              const Box& a_rangeBox)
  {
    static bool initFlag=s_init();
    using namespace std;
    PR_TIME("EulerOp::operator");
    a_Rhs.setVal(0.0);

    double gamma = s_gamma;
    double retval;

#if DATAFLOW_ON > 0
    DataFlowFactory& fac = DataFlowFactory::get();
    fac.init("euler_step", "retval", "d", "i"); //typeid(retval).name());
    fac.constants({"N", "C", "G", "D"}, {NUMCELLS, NUMCOMPS, NGHOST, DIM});
    fac.newSpace<double,NUMCOMPS>("rhs", a_Rhs);
#endif
    //PR_TIME("EulerOp::operator::W_bar");
    Vector W_bar = forall<double,NUMCOMPS>(consToPrim,a_U, gamma);
#if DATAFLOW_ON > 0
    fac.newComp<double,NUMCOMPS>("consToPrim", {"U"}, "W_bar", W_bar, consToPrim, a_U, gamma);
#endif
    //PR_TIME("EulerOp::operator::U");
    Vector U = m_deconvolve(a_U);
#if DATAFLOW_ON > 0
    fac.newComp<double,NUMCOMPS>("deconvolve", "u", m_deconvolve, a_U, 1.0, U);
#endif
    //PR_TIME("EulerOp::operator::W");
    Vector W = forall<double,NUMCOMPS>(consToPrim,U, gamma);
#if DATAFLOW_ON > 0
    fac.newComp<double,NUMCOMPS>("consToPrim", {"u"}, "W", W, consToPrim, U, gamma);
#endif
    Scalar umax = forall<double>(waveSpeedBound,a_rangeBox,W, gamma);
#if DATAFLOW_ON > 0
    fac.newComp<double>("waveSpeedBound", {"W"}, "umax", umax, waveSpeedBound, a_rangeBox, W, gamma);
#endif
    retval = umax.absMax();
#if DATAFLOW_ON > 0
    fac.newComp<double>("absMax", "retval", "absmax(", retval, umax);
#endif

    //PR_TIME("EulerOp::operator::W_ave");
    Vector W_ave = m_laplacian(W_bar,1.0/24.0);
#if DATAFLOW_ON > 0
    fac.newComp<double,NUMCOMPS>("laplacian", "W_ave", m_laplacian, W_bar, 1.0/24.0, W_ave);
#endif
    W_ave += W;
#if DATAFLOW_ON > 0
    fac.newComp<double,NUMCOMPS>("increment", "W_ave2", "+=", W_ave, W);
#endif

    for (int d = 0; d < DIM; d++)
      {
        //PR_TIME("EulerOp::operator::W_ave_f::interpolation");
        Vector W_ave_low = m_interp_L[d](W_ave);
#if DATAFLOW_ON > 0
        string dim = "d" + to_string(d+1);
        fac.newComp<double,NUMCOMPS>("interpL_" + dim, "W_aveL_" + dim, m_interp_L[d], W_ave, 1.0, W_ave_low);
#endif
        Vector W_ave_high = m_interp_H[d](W_ave);
#if DATAFLOW_ON > 0
        fac.newComp<double,NUMCOMPS>("interpH_" + dim, "W_aveH_" + dim, m_interp_H[d], W_ave, 1.0, W_ave_high);
#endif
        //PR_TIME("EulerOp::operator::W_ave_f::upwinding");
        Vector W_ave_f = forall<double,NUMCOMPS>(upwindState,W_ave_low, W_ave_high, d, gamma);
#if DATAFLOW_ON > 0
        fac.newComp<double,NUMCOMPS>("upwindState", {"W_aveL_" + dim, "W_aveH_" + dim}, "W_ave_f_" + dim,
                                     W_ave_f, upwindState, W_ave_low, W_ave_high, d, gamma);
#endif

#if DIM>1
        //PR_TIME("EulerOp::operator::F_bar_f");
        Vector F_bar_f = forall<double,NUMCOMPS>(getFlux, W_ave_f, d, gamma);
#if DATAFLOW_ON > 0
        fac.newComp<double,NUMCOMPS>("getFlux", {"W_ave_f_" + dim}, "F_bar_f_" + dim, F_bar_f, getFlux, W_ave_f, d, gamma);
#endif
#endif
        //PR_TIME("EulerOp::operator::W_f");
#if DIM>1
        Vector W_f = m_deconvolve_f[d](W_ave_f);
#if DATAFLOW_ON > 0
        fac.newComp<double,NUMCOMPS>("deconvolve_f_" + dim, "W_f_" + dim, m_deconvolve_f[d], W_ave_f, 1.0, W_f);
#endif
#else
        Vector W_f = W_ave_f;
#endif
        //PR_TIME("EulerOp::operator::F_ave");
        Vector F_ave_f = forall<double,NUMCOMPS>(getFlux, W_f, d, gamma);
#if DATAFLOW_ON > 0
        fac.newComp<double,NUMCOMPS>("getFlux", {"W_f_" + dim}, "F_ave_f_" + dim, F_ave_f, getFlux, W_f, d, gamma);
#endif

#if DIM>1
        F_bar_f *= (1./24.);
#if DATAFLOW_ON > 0
        fac.newComp<double,NUMCOMPS>("smul_" + dim, "F_bar2_f_" + dim, "*=", F_bar_f, 1./24.);
#endif
#if DATAFLOW_ON > 0
        Vector F_lap_f = m_laplacian_f[d](F_bar_f,1.0/24.0);
        fac.newComp<double,NUMCOMPS>("lap_f_" + dim, "F_lap_f_" + dim, m_laplacian_f[d], F_bar_f, 1.0/24.0, F_lap_f);
        F_ave_f += F_lap_f;
        fac.newComp<double,NUMCOMPS>("inc_f_" + dim, "F_ave2_f_" + dim, "+=", F_ave_f, F_lap_f);
#else
        F_ave_f += m_laplacian_f[d](F_bar_f,1.0/24.0);
#endif
#endif
        //PR_TIME("EulerOp::operator::minusDivF");
#if DATAFLOW_ON > 0
        Vector F_div_f = m_divergence[d](F_ave_f);
        fac.newComp<double,NUMCOMPS>("div_f_" + dim, "F_div_f_" + dim, m_divergence[d], F_ave_f, 1.0, F_div_f);
        a_Rhs += F_div_f;
        fac.newComp<double,NUMCOMPS>("inc_rhs_" + dim, "rhs_" + dim, "+=", a_Rhs, F_div_f);
#else
        a_Rhs += m_divergence[d](F_ave_f);
#endif
      }
    //PR_TIME("EulerOp::operator::RHS*=-1.0/dx");
    a_Rhs *= -1./s_dx;
#if DATAFLOW_ON > 0
    fac.newComp<double,NUMCOMPS>("muldx", "rhs", "*=", a_Rhs, -1./s_dx);

    // Fuse Commands
    pdfg::fuse({"consToPrim1", "deconvolve", "consToPrim2", "waveSpeedBound1", "absMax"});
    pdfg::fuse({"laplacian", "increment", "interpL_d1", "interpH_d1", "interpL_d2", "interpH_d2"});
    pdfg::fuse({"upwindState1", "getFlux1", "smul_d1"});
    pdfg::fuse({"deconvolve_f_d1", "getFlux2"});
    pdfg::fuse({"lap_f_d1", "inc_f_d1", "div_f_d1", "inc_rhs_d1"});
    pdfg::fuse({"upwindState2", "getFlux3", "smul_d2"});
    pdfg::fuse({"deconvolve_f_d2", "getFlux4"});
#if DIM<3
    pdfg::fuse({"lap_f_d2", "inc_f_d2", "div_f_d2", "inc_rhs_d2", "muldx"});
#else
    pdfg::fuse({"lap_f_d2", "inc_f_d2", "div_f_d2", "inc_rhs_d2"});
    pdfg::fuse({"interpL_d3", "interpH_d3"});
    pdfg::fuse({"upwindState3", "getFlux5", "smul_d3"});
    pdfg::fuse({"deconvolve_f_d3", "getFlux6"});
    pdfg::fuse({"lap_f_d3", "inc_f_d3", "div_f_d3", "inc_rhs_d3", "muldx"});
#endif

    pdfg::perfmodel();
    fac.print("out/euler_step.json");
    fac.codegen("out/euler_step.h");
#endif
    return retval;
  }
}
