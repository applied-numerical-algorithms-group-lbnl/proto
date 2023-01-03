#pragma once
#ifndef _BOX_OP_EULER_
#define _BOX_OP_EULER_

#include "Proto.H"

#define NUMCOMPS DIM+2

using namespace Proto;

typedef BoxData<double> Scalar;
typedef BoxData<double, NUMCOMPS> Vector;
//State: [rho, G0, G1, ..., E]
// Gi = rho*vi
// E = p/(gamma-1) + 0.5*rho*|v|^2

template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_consToPrim_(
        Var<T, NUMCOMPS, MEM>&          a_W, 
        const Var<T, NUMCOMPS, MEM>&  a_U,
        double                          a_gamma)
{
    double rho = a_U(0);
    double v2 = 0.0;
    a_W(0) = rho;

    for (int i = 1; i <= DIM; i++)
    {
        double v;
        v = a_U(i) / rho;

        a_W(i) = v;
        v2 += v*v;
    }

    a_W(NUMCOMPS-1) = (a_U(NUMCOMPS-1) - .5 * rho * v2) * (a_gamma - 1.0);
}
PROTO_KERNEL_END(f_consToPrim_, f_consToPrim)

template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_upwindState_(
        Var<T, NUMCOMPS, MEM>&       a_out,
        const Var<T, NUMCOMPS, MEM>& a_low,
        const Var<T, NUMCOMPS, MEM>& a_high,
        int                          a_dir,
        double                       a_gamma)
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
PROTO_KERNEL_END(f_upwindState_, f_upwindState)
    
template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_getFlux_(
        Var<T, NUMCOMPS, MEM>&       a_F,
        const Var<T, NUMCOMPS, MEM>& a_W, 
        int                          a_dir,
        double                       a_gamma)
{
    double F0 = a_W(a_dir+1)*a_W(0);
    double W2 = 0.0;
    double gamma = a_gamma;

    a_F(0) = F0;

    for (int d = 1; d <= DIM; d++)
    {
        double Wd = a_W(d);

        a_F(d) = Wd*F0;
        W2 += Wd*Wd;
    }

    a_F(a_dir+1) += a_W(NUMCOMPS-1);
    a_F(NUMCOMPS-1) = gamma/(gamma - 1.0) * a_W(a_dir+1) * a_W(NUMCOMPS-1) + 0.5 * F0 * W2;
}
PROTO_KERNEL_END(f_getFlux_, f_getFlux)

template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_waveSpeedBound_(Var<double,1>& a_speed,
        const Var<T, NUMCOMPS, MEM>& a_W,
        double       a_gamma)
{
    a_speed(0) = DIM*sqrt(a_gamma*a_W(NUMCOMPS-1)/a_W(0));
    for (int dir = 1 ; dir <= DIM; dir++)
    {
        a_speed(0) += a_W(dir);
    }
}
PROTO_KERNEL_END(f_waveSpeedBound_, f_waveSpeedBound)

template<typename T, MemType MEM = MEMTYPE_DEFAULT>
class BoxOp_Euler : public BoxOp<T, NUMCOMPS, 1, MEM>
{
    public:

    T gamma = 1.4;
    mutable T umax;

    // How many ghost cells does the operator need from the state variables
    inline static Point ghost() { return Point::Ones(4);}
    
    // How many ghost cells does the operator need from the auxiliary variables
    inline static Point auxGhost() { return Point::Zeros();}
    
    // What is the intended order of accuracy of the operator
    inline static constexpr int order() { return 4; }
    
    inline BoxOp_Euler() : BoxOp<T, NUMCOMPS, 1, MEM>() {}
    inline BoxOp_Euler(Array<double, DIM> a_dx) : BoxOp<T, NUMCOMPS, 1, MEM>(a_dx) {}
    inline BoxOp_Euler(double a_dx) : BoxOp<T, NUMCOMPS, 1, MEM>(a_dx) {}
    
    // Initialization
    inline void init()
    {
        for (int dir = 0; dir < DIM; dir++)
        {
            m_interp_H[dir] = Stencil<double>::CellToFaceH(dir);
            m_interp_L[dir] = Stencil<double>::CellToFaceL(dir);
            m_divergence[dir] = Stencil<double>::FluxDivergence(dir);
            m_laplacian_f[dir] = Stencil<double>::LaplacianFace(dir);
        }
    }

    // Helper Function
    inline void computeFlux(
            BoxData<T, NUMCOMPS>& a_flux,
            const BoxData<T, NUMCOMPS>& a_W_ave,
            int a_dir) const
    {
        PR_TIME("BoxOp_Euler::computeFlux");
        unsigned long long int sqrtnum = 10;  //this one is just a guess
        unsigned long long int getfluxnum = 9 + DIM*3;
        unsigned long long int upwindnum  = sqrtnum + 25;
        Vector W_ave_L = m_interp_L[a_dir](a_W_ave); 
        Vector W_ave_H = m_interp_H[a_dir](a_W_ave); 
        Vector W_ave_f = forallOp<double,NUMCOMPS>(
                upwindnum, "upwind", f_upwindState, W_ave_L, W_ave_H, a_dir, gamma);
#if DIM>1
        Vector F_bar_f = forallOp<double,NUMCOMPS>(
                getfluxnum, "getflux", f_getFlux, W_ave_f, a_dir,  gamma);
        Vector W_f = Operator::deconvolveFace(W_ave_f, a_dir);
#else
        Vector W_f = W_ave_f;
#endif
        a_flux = forallOp<double,NUMCOMPS>(
                getfluxnum, "getflux", f_getFlux, W_f, a_dir, gamma);
#if DIM>1
        a_flux += m_laplacian_f[a_dir](F_bar_f, 1.0/24.0);
#endif
    }
   
    // Flux Definition
    inline void flux(
            BoxData<T, NUMCOMPS>& a_flux,
            const BoxData<T, NUMCOMPS>& a_U,
            int a_dir) const
    {
        PR_TIME("BoxOp_Euler::flux");
        // FLOP COUNTS
        unsigned long long int ctoprmnum  = 4*DIM + 5;
        
        Vector W_bar = forallOp<double, NUMCOMPS>(
            ctoprmnum, "consToPrim", f_consToPrim, a_U, gamma);
        Vector U = Operator::deconvolve(a_U);
        Vector W = forallOp<double, NUMCOMPS>(
            ctoprmnum, "consToPrim", f_consToPrim, U, gamma);
        Vector W_ave = Operator::_convolve(W, W_bar);
        computeFlux(a_flux, W_ave, a_dir);
    }
    
    // Apply Operator
    inline void operator()(
            BoxData<T, NUMCOMPS>&                   a_Rhs,
            Array<BoxData<T, NUMCOMPS>, DIM>&  a_fluxes,
            const BoxData<T, NUMCOMPS>&             a_U,
            T                                       a_scale = 1.0) const
    {
        T dx = this->dx()[0];
        PR_TIME("BoxOp_Euler::operator()");
        // FLOP COUNTS
        unsigned long long int sqrtnum = 10;  //this one is just a guess
        unsigned long long int ctoprmnum  = 4*DIM + 5;
        unsigned long long int wavespdnum = sqrtnum +3 + DIM;
    
        // COMPUTE W_AVE
        a_Rhs.setVal(0.0);    
        Vector W_bar = forallOp<double, NUMCOMPS>(
            ctoprmnum, "consToPrim", f_consToPrim, a_U, gamma);
        Vector U = Operator::deconvolve(a_U);
        Vector W = forallOp<double, NUMCOMPS>(
            ctoprmnum, "consToPrim", f_consToPrim, U, gamma);
        Vector W_ave = Operator::_convolve(W, W_bar);
        
        // COMPUTE MAX WAVE SPEED
        Box rangeBox = a_U.box().grow(-ghost());
        Scalar uabs = forallOp<double>(
            wavespdnum, "wavespeed", f_waveSpeedBound, rangeBox, W, gamma);
        umax = uabs.absMax();

        // COMPUTE DIV FLUXES
        for (int dir = 0; dir < DIM; dir++)
        {
            computeFlux(a_fluxes[dir], W_ave, dir);
            a_Rhs += m_divergence[dir](a_fluxes[dir]);
        }
        a_Rhs *= (-1.0*a_scale / dx); //Assuming isotropic grid spacing
    }


    private:
    Stencil<T> m_interp_H[DIM];
    Stencil<T> m_interp_L[DIM];
    Stencil<T> m_divergence[DIM];
    Stencil<T> m_laplacian_f[DIM];
};

#endif //end include guard
