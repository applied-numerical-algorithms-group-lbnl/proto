#ifndef _EULEROP_H_
#define _EULEROP_H_

#include "Proto.H"

constexpr unsigned int NUMCOMPS=DIM+2;
#define NGHOST 4

using namespace Proto;

typedef Var<double,NUMCOMPS> State;

class EulerOp
{
public:
    static double s_dx;
    static double s_gamma;

    /// Default constructor.
    EulerOp();

    ///  Destructor.
    ~EulerOp(){};

    /// Compute the Right-Hand side, conforming to the RK4 interface.
    double operator()(BoxData<double,NUMCOMPS>& a_Rhs,
                      const BoxData<double,NUMCOMPS>& a_U,
                      const Bx& a_dbx0);
    /// Pointwise Functions. 

private:
    Stencil<double> m_laplacian;
    Stencil<double> m_deconvolve;
    Stencil<double> m_laplacian_f[DIM];
    Stencil<double> m_deconvolve_f[DIM];
    Stencil<double> m_interp_H[DIM];
    Stencil<double> m_interp_L[DIM];
    Stencil<double> m_divergence[DIM];
};
#endif //end include guard
