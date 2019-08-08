#pragma once
#include "BaseOp.H"

class MehrstellenCorrectionOp : public BaseOp<1,1> {
    public: 
        typedef Proto::BoxData<Real> patch; 
        typedef Proto::Var<Real> var;
        inline void init();
        inline void flux(patch& a_flux, const patch& a_phi, int a_dir, 
            Real a_dx, Real a_scale = 1.0) const;
};

void MehrstellenCorrectionOp::init()
{
    this->m_lambda = -1;
    this->m_interpBC.push_back(Proto::InterpStencil<Real>::Quadratic(this->m_refRatio));
}

void MehrstellenCorrectionOp::flux(patch& a_flux,
        const patch& a_phi,
        int a_dir,
        Real a_dx,
        Real a_scale) const
{
    //  F = h^2/12 * 1/h(phi_high - phi_low) = h/12 * (phi_high - phi_low)
    a_flux = this->m_flux[a_dir](a_phi, a_dx/12.0);
}
