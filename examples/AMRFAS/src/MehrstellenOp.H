#pragma once
#include "BaseOp.H"

/// Fourth Order Mehrstellen Operator
/**
    Implements the 4th order accurate Mehrstellen Laplacian. This is equivalent to the 9 point Laplacian in 2D or the 19 point Laplacian in 3D.
    This operator is not implemented for any other value of DIM.
    
    MehrstellenOp assumes that the user has properly added the Mehrstellen correction of \f$ \frac{1}{12}\nabla_2\rho \f$ to the right-hand side.
    The operator MehrstellenCorrectionOp is provided to assist in computing this term conservatively. 
*/
class MehrstellenOp : public BaseOp<1, 2> {
    public:
        typedef Proto::BoxData<Real> patch; 
        typedef Proto::Var<Real> var; 

        inline MehrstellenOp()
        {
            PROTO_ASSERT(DIM >= 2 && DIM <= 3, "MehrstellenOp is not defined for DIM = %i. Use DIM = 2 or 3.", DIM);
        }
        /// Init
        /**
            Initializes the relaxation parameter as \f$ \lambda = -\frac{3h^2}{20}\f$ in 2D or \f$ \lambda = - \frac{15h^2}{128}\f$ in 3D. The values is doubled if GSRB = true
            Boundary interpolation uses a 5th order accurate quartic interpolant.
        */
        inline void init();
        /// Flux
        /**
            Computes the flux \f$ F_i = \frac{1}{h}(f_{i} - f_{i-1}) \f$ where \f$ f_i \equiv \phi_i + \nabla_2^\perp\phi \f$
        */
        inline void flux(patch& a_flux, const patch& a_phi, int a_dir,
                Real a_dx, Real a_scale = 1.0) const;
#ifdef GSRB
        /// Relax
        /**
            Overrides the default Jacobi iteration with 4-Color Gauss Seidel iteration.
        */
        inline void relax(LevelData<FArrayBox>& a_phi,
            const LevelData<FArrayBox>& a_R,
            int a_n) const;
#endif
    private:
        Proto::Stencil<Real> m_Ld2[DIM];            ///< Stencils for computing the 2nd order perpendicular Laplacian
#ifdef GSRB
        Proto::Stencil<Real> m_L;                   ///< Stencil for computing the 9 or 19 point Mehrstellen operator
        Proto::Stencil<Real> m_GSRB[ipow<2>(DIM)];  ///< Stencils for 4-Color Gauss Seidel
        Proto::Stencil<Real> m_GSRB_I[ipow<2>(DIM)];///< Stencils for 4-Color Gauss Seidel
#endif
};

void MehrstellenOp::init()
{
    PROTO_ASSERT(DIM == 2 || DIM == 3,
    "MehrstellenOp not defined for DIM = %i",DIM);
#if DIM==2
    this->m_lambda = -0.5*(this->m_dx*this->m_dx)*(6.0/20.0);  //9-Point
#elif DIM==3
    this->m_lambda = -0.5*(this->m_dx*this->m_dx)*(15.0/64.0); //19-Point
#endif
    for (int ii = 0; ii < DIM; ii++)
    {
        m_Ld2[ii] = Proto::Stencil<Real>::LaplacianFace(ii,2);
    }
    this->buildInterpStencils(3,
            Proto::Box::Cube(5).shift(Proto::Point::Ones(-2)),
            4);
#ifdef GSRB
    (this->m_lambda) *= 2.0;
#if DIM==2
    m_L = Proto::Stencil<Real>::Laplacian_9();
#elif DIM==3
    m_L = Proto::Stencil<Real>::Laplacian_19();
#endif
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

void MehrstellenOp::flux(patch& a_flux,
        const patch& a_phi,
        int a_dir,
        Real a_dx,
        Real a_scale) const
{
    PR_TIME("MehrstellenOp::flux");
    //F(phi) = 1/h * ( flux_d(phi) + h^2/12*L2_d(flux_d(phi));
    auto S =  1.0*Proto::Shift::Zeros() + (1.0/12.0) * m_Ld2[a_dir];
    
    patch f0 = S(a_phi);
    a_flux = this->m_flux[a_dir](f0,a_scale/a_dx);
}

#ifdef GSRB
void MehrstellenOp::relax(LevelData<FArrayBox>& a_phi,
        const LevelData<FArrayBox>& a_R,
        int a_n) const
{
    PR_TIME("MehrstellenOp::relax");
    auto layout = a_phi.disjointBoxLayout();
    auto iter = a_phi.dataIterator();
    auto K = Proto::Box::Cube(2);
    for (int ii = 0; ii < a_n; ii++)
    {
        for (int icolor = 0; icolor < 4; icolor++)
        {
            a_phi.exchange();
            for (iter.reset(); iter.ok(); ++iter)
            {
                patch phi = a_phi[iter];
                patch rhs = a_R[iter];
                patch temp(layout[iter]);
                
                temp |= m_GSRB[icolor](phi);
                temp += m_GSRB_I[icolor](rhs);
                phi += m_GSRB_I[icolor](temp, this->m_lambda);
                #if DIM==3
                temp |= m_GSRB[7 - icolor](phi);
                temp += m_GSRB_I[7 - icolor](rhs);
                phi += m_GSRB_I[7 - icolor](temp, this->m_lambda);
                #endif
            }
        }
    }// end iterations
}
#endif
