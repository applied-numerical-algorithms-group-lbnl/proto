#include "Proto.H"
#include "EulerOp.H"
//#include "CommonTemplates.H"
#include "Proto_Timer.H"

typedef BoxData<double> Scalar;
typedef BoxData<double,NUMCOMPS> Vector;

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

        double outval = a_out(0) + (pstar - a_out(NUMCOMPS-1))/(cbar*cbar);
        if (cbar + sign*ubar > 0)
        {
            a_out(0) = outval;
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

    PROTO_KERNEL_START 
    unsigned int InitializeStateF(State& a_U,
            const V& a_x,
            const double a_gamma)
    {
        double rho0 = a_gamma;
        double p0 = 1.;
        double umag = 0.;
        double rho = rho0;
        rho += .01*rho0*sin(2*2*PI*a_x(0));
        double p = p0*pow(rho/rho0,a_gamma);
        a_U(0) = rho;
        double c0 = sqrt(a_gamma*p0/rho0);
        double c = sqrt(a_gamma*p/rho);
        umag = 2*(c-c0)/(a_gamma-1.);
        a_U(1) = rho*umag;
        //NOTE: This assumes that NUMCOMPS=DIM+2
        for(int dir=2; dir<=DIM; dir++)
            a_U(dir)=0.0;
        double ke = 0.5*umag*umag;
        a_U(NUMCOMPS-1) = p/(a_gamma-1.0) + rho*ke;
        return 0;
    }
    PROTO_KERNEL_END(InitializeStateF, InitializeState)

        //=================================================================================================
    PROTO_KERNEL_START
    void iotaFuncF(Point           & a_p,
            V               & a_X,
            double            a_h)
    {
        for (int ii = 0; ii < DIM; ii++)
        {
            a_X(ii) = a_p[ii]*a_h + 0.5*a_h;
        }
    }
    PROTO_KERNEL_END(iotaFuncF,iotaFunc)

        void initializeState(BoxData<double,NUMCOMPS>& a_state,
                const double a_dx,
                const double a_gamma)
        {
            Box dbx0=a_state.box();
            Box dbx = dbx0.grow(NGHOST);
            Box dbx1 = dbx.grow(1);
            BoxData<double,NUMCOMPS> UBig(dbx1);
            BoxData<double,DIM> x(dbx1);
            forallInPlace_p(iotaFunc, dbx1, x, a_dx);
            forallInPlace(InitializeState,dbx1,UBig,x,a_gamma);
            Stencil<double> Lap2nd = Stencil<double>::Laplacian();
            a_state |= Lap2nd(UBig,dbx,1.0/24.0); 
            a_state += UBig;
        }

    void step(BoxData<double,NUMCOMPS>& a_Rhs,
            const BoxData<double,NUMCOMPS>& a_U,
            const Box& a_rangeBox,
            const double a_dx,
            const double a_gamma,
            Reduction<double, Abs>& a_Rxn,
            bool a_computeMaxWaveSpeed,
            bool a_callBCs)
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
                m_interp_H[dir] = Stencil<double>::CellToFaceH(dir);
                m_interp_L[dir] = Stencil<double>::CellToFaceL(dir);
                m_divergence[dir] = Stencil<double>::FluxDivergence(dir);
            }
            initialized =  true;
        }


        using namespace std;
        PR_TIME("EulerOp::operator");
        a_Rhs.setVal(0.0);

        unsigned long long int sqrtnum = 10;  //this one is just a guess
        unsigned long long int ctoprmnum  = 4*DIM + 5;
        unsigned long long int upwindnum  = sqrtnum + 25;
        unsigned long long int getfluxnum = 9 + DIM*3;
        unsigned long long int wavespdnum = sqrtnum +3 + DIM;

        auto idOp = (1.0)*Shift(Point::Zeros());
        double gamma = a_gamma;
        //PR_TIME("EulerOp::operator::W_bar");
        if(a_callBCs)
        {
            BoxData<double, NUMCOMPS>& castU = const_cast<BoxData<double, NUMCOMPS> &>(a_U);
            int nghost = a_rangeBox.low()[0] - castU.box().low()[0];
            for(int idir = 0; idir < DIM; idir++)
            {
                protocommon::enforceSGBoundaryConditions(castU, nghost, idir);
            }
        }
        Vector W_bar = forallOp<double,NUMCOMPS>(ctoprmnum, "consToPrim", consToPrim,a_U, gamma); // 
        //PR_TIME("EulerOp::operator::U");
        Vector U = m_deconvolve(a_U);
        //PR_TIME("EulerOp::operator::W");
        Vector W    = forallOp<double,NUMCOMPS>(ctoprmnum, "consToPrim",consToPrim,U, gamma);
        if(a_computeMaxWaveSpeed)
        {
            Scalar umax = forallOp<double>(wavespdnum, "wavespeed", waveSpeedBound,a_rangeBox,W, gamma);

            umax.absMax(a_Rxn);
        }

        //PR_TIME("EulerOp::operator::W_ave");
        Vector W_ave = m_laplacian(W_bar,1.0/24.0);
        {
            PR_TIME("EulerOp::operator::W_ave +=");
            W_ave += idOp(W);
        }
        for (int d = 0; d < DIM; d++)
        {
            PR_TIME("EulerOp::operator::dimLoop");
            //PR_TIME("EulerOp::operator::W_ave_f::interpolation");
            Vector W_ave_low = m_interp_L[d](W_ave);
            Vector W_ave_high = m_interp_H[d](W_ave);
            //PR_TIME("EulerOp::operator::W_ave_f::upwinding");
            Vector W_ave_f = forallOp<double,NUMCOMPS>(upwindnum, "upwind",
                    upwindState,W_ave_low, W_ave_high,d,  gamma);
#if DIM>1
            //PR_TIME("EulerOp::operator::F_bar_f");
            Vector F_bar_f = forallOp<double,NUMCOMPS>(getfluxnum, "getflux", getFlux, W_ave_f, d,  gamma);
#endif
            //PR_TIME("EulerOp::operator::W_f");
#if DIM>1
            Vector W_f = m_deconvolve_f[d](W_ave_f);
#else
            Vector W_f = W_ave_f;
#endif
            //PR_TIME("EulerOp::operator::F_ave");
            Vector F_ave_f = forallOp<double,NUMCOMPS>(getfluxnum, "getflux", getFlux, W_f, d, gamma);
#if DIM>1
            // F_bar_f *= (1./24.);
            F_ave_f += m_laplacian_f[d](F_bar_f,1.0/24.0);
#endif
            //PR_TIME("EulerOp::operator::minusDivF");
            a_Rhs += m_divergence[d](F_ave_f);
        }

        //PR_TIME("EulerOp::operator::RHS*=-1.0/dx");
        a_Rhs *= -1./a_dx;
    }
}
