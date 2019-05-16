#pragma once
#include "BaseOp.H"

template<class DATA>
class MehrstellenOp : public BaseOp<1, 2, DATA> {
    public:
        typedef Proto::BoxData<Real> patch; 
        typedef Proto::Var<Real> var; 

        inline MehrstellenOp()
        {
            //std::cout << "Created a MehrstellenOp" << std::endl;
        }
        inline void init();
        inline patch flux(const patch& a_phi, int a_dir,
                Real a_dx, Real a_scale = 1.0) const;
        inline void operator()(patch& a_rhs, patch& a_phi,
                Real a_dx, Real a_scale = 1.0) const;    
        inline void operator()(LevelData<DATA>& a_rhs,
                LevelData<DATA>& a_phi,
                Real a_dx, Real a_scale = 1.0) const;
        inline void precondition(
            std::vector<std::shared_ptr<LevelData<DATA>>> a_rho);
#ifdef GSRB
        inline void relax(LevelData<DATA>& a_phi,
            const LevelData<DATA>& a_R,
            int a_n) const;
#endif
    private:
        Proto::Stencil<Real> m_Ld2[DIM];
#ifdef GSRB
        Proto::Stencil<Real> m_L;
        Proto::Stencil<Real> m_GSRB[ipow<2>(DIM)];
        Proto::Stencil<Real> m_GSRB_I[ipow<2>(DIM)];
#endif
};

template<class DATA>
void MehrstellenOp<DATA>::init()
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

template<class DATA>
Proto::BoxData<Real> MehrstellenOp<DATA>::flux(const patch& a_phi,
        int a_dir,
        Real a_dx,
        Real a_scale) const
{
    //F(phi) = 1/h * ( flux_d(phi) + h^2/12*L2_d(flux_d(phi));
    auto S =  1.0*Proto::Shift::Zeros() + (1.0/12.0) * m_Ld2[a_dir];
    
    patch f0 = S(a_phi);
    patch f1 = this->m_flux[a_dir](f0,a_scale/a_dx);
    return f1;
}

#ifdef GSRB
template<class DATA>
void MehrstellenOp<DATA>::relax(LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_R,
        int a_n) const
{
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
