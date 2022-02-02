#pragma once
#include "BaseOp.H"

/// Second Order Laplacian Operator
/**
    Implements Red-Black Gauss Seidel relaxation if GSRB is true.
*/
class LaplaceOp : public BaseOp<1, 2> {
    public:
        typedef Proto::BoxData<Real> patch; 
        typedef Proto::Var<Real> var; 

        inline LaplaceOp(){}
        /// Init
        /**
            Initializes the relaxation parameter as: \f$ \lambda = -\frac{h^2}{4*DIM}\f$. The values is doubled if GSRB = true
            Boundary interpolation uses a 3rd order accurate quadratic interpolant.
        */
        inline void init();
        
        /// Flux
        /**
            Computes the flux as \f$ F_i \equiv \frac{1}{h}(\phi_i - \phi_{i-1})\f$
        */
        inline void flux(patch& a_flux, const patch& a_phi, int a_dir,
                Real a_dx, Real a_scale = 1.0) const;
#ifdef GSRB
        /// Relax
        /**
            Overrides the default Jacobi iteration with Red - Black Gauss Seidel iteration.
        */
        inline void relax(LevelData<FArrayBox>& a_phi,
            const LevelData<FArrayBox>& a_R,
            int a_n) const;
#endif
    private:
#ifdef GSRB
        Proto::Stencil<Real> m_L;                   ///< Stencil of the full Laplace operator. Equivalent to the 2*DIM+1 point operator.
        Proto::Stencil<Real> m_GSRB[ipow<2>(DIM)];  ///< Set of Stencils used in GSRB algorithm
        Proto::Stencil<Real> m_GSRB_I[ipow<2>(DIM)];///< Set of Stencils used in GSRB algorithm
#endif
};

void LaplaceOp::init()
{
    this->m_lambda = -(this->m_dx*this->m_dx)/(4.0*DIM);
    
    std::vector<Proto::Point> shifts;
    auto K = Proto::Box::Cube(3).shift(Proto::Point::Ones(-1));
    for (auto iter = K.begin(); iter != K.end(); ++iter)
    {
        shifts.push_back(*iter);
    }
    this->buildInterpStencils(shifts, 2);
    
#ifdef GSRB
    (this->m_lambda) *= 2.0;
    m_L = Proto::Stencil<Real>::Laplacian();
    auto kernel = Proto::Box::Cube(2);
    for (auto iter = kernel.begin(); iter != kernel.end(); ++iter)
    {
        int index = kernel.index(*iter);
        m_GSRB[index] = m_L;
        m_GSRB[index] *= (-1.0/(this->m_dx*this->m_dx))*Proto::Shift(*iter);
        m_GSRB[index].srcRatio() = Proto::Point::Ones(2);
        m_GSRB[index].destRatio() = Proto::Point::Ones(2);
        m_GSRB[index].destShift() = *iter;
        
        m_GSRB_I[index] = 1.0*Proto::Shift(*iter);
        m_GSRB_I[index].srcRatio() = Proto::Point::Ones(2);
        m_GSRB_I[index].destRatio() = Proto::Point::Ones(2);
        m_GSRB_I[index].destShift() = *iter;
    }
#endif
}

void LaplaceOp::flux(patch& a_flux,
        const patch& a_phi,
        int a_dir,
        Real a_dx,
        Real a_scale) const
{
    a_flux = this->m_flux[a_dir](a_phi, a_scale/a_dx);
}
#ifdef GSRB
template<class FArrayBox>
void LaplaceOp<FArrayBox>::relax(LevelData<FArrayBox>& a_phi,
        const LevelData<FArrayBox>& a_R,
        int a_n) const
{
    auto layout = a_phi.disjointBoxLayout();
    auto iter = a_phi.dataIterator();
    auto K = Proto::Box::Cube(2);
    for (int ii = 0; ii < a_n; ii++)
    {
        for (int icolor = 0; icolor < 2; icolor++)
        {
            a_phi.exchange();
            for (iter.reset(); iter.ok(); ++iter)
            {
                patch phi = a_phi[iter];
                patch rhs = a_R[iter];
                patch temp(layout[iter]);
                for (auto kiter = K.begin(); kiter != K.end(); ++kiter)
                {
                    int color = 0;
                    for (int ii = 0; ii < DIM; ii++)
                    {
                        color += (*kiter)[ii];
                    }
                    color = color % 2;
                    if (color == icolor)
                    {
                        int index = K.index(*kiter);
                        temp |= m_GSRB[index](phi);
                        temp += m_GSRB_I[index](rhs);
                        phi += m_GSRB_I[index](temp, this->m_lambda);
                    }
                }
            }
        }
    }// end iterations
}
#endif
