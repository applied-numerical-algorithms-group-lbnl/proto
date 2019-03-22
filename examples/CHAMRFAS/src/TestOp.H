#include "BaseOp.H"

#define GSRB TRUE

template<class DATA>
class TestOp : public BaseOp<1, DATA> {
    public:
        typedef Proto::BoxData<Real> patch; 
        typedef Proto::Var<Real> var; 

        inline TestOp(){}
        inline void init();
        inline Proto::Point ghost() const {return Proto::Point::Ones();}
        inline patch flux(const patch& a_phi, int a_dir,
                Real a_dx, Real a_scale = 1.0) const;
        inline void operator()(patch& a_rhs, const patch& a_phi,
                Real a_dx, Real a_scale = 1.0) const;    
        inline void operator()(LevelData<DATA>& a_rhs,
                const LevelData<DATA>& a_phi,
                Real a_dx, Real a_scale = 1.0) const;
#ifdef GSRB
        inline void relax(LevelData<DATA>& a_phi,
            const LevelData<DATA>& a_R,
            int a_n) const;
#endif
    private:
        Proto::Stencil<Real> m_flux[DIM];
        Proto::Stencil<Real> m_div[DIM]; 
        Proto::Stencil<Real> m_L;
        Proto::Stencil<Real> m_GSRB[ipow<2>(DIM)];
        Proto::Stencil<Real> m_GSRB_I[ipow<2>(DIM)];
};

template<class DATA>
void TestOp<DATA>::init()
{
    this->m_lambda = -(this->m_dx*this->m_dx)/(4.0*DIM);
#ifdef GSRB
    (this->m_lambda) *= 2.0;
#endif
#if DIM==2
    //this->m_lambda = -0.5*(this->m_dx*this->m_dx)*(6.0/20.0);  //9-Point
#elif DIM==3
    //this->m_lambda = -0.5*(this->m_dx*this->m_dx)*(15.0/64.0); //27-Point
#endif
    for (int ii = 0; ii < DIM; ii++)
    {
        m_flux[ii] = 1.0*Proto::Shift::Zeros()
                   - 1.0*Proto::Shift::Basis(ii,-1);
        m_div[ii] = 1.0*Proto::Shift::Basis(ii,1)
                  - 1.0*Proto::Shift::Zeros();
    }
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
}

template<class DATA>
Proto::BoxData<Real> TestOp<DATA>::flux(const patch& a_phi,
        int a_dir,
        Real a_dx,
        Real a_scale) const
{
    Proto::BoxData<Real> flux = m_flux[a_dir](a_phi, a_scale/a_dx);
    return flux;
}

template<class DATA>
void TestOp<DATA>::operator()(patch& a_rhs, const patch& a_phi,
        Real a_dx, Real a_scale) const
{
    for (int ii = 0; ii < DIM; ii++)
    {
        patch flux = m_flux[ii](a_phi, 1.0/a_dx);
        if (ii == 0)
        {
            a_rhs |= m_div[ii](flux, a_scale/a_dx);
        } else {
            a_rhs += m_div[ii](flux, a_scale/a_dx);
        }
    }
}

template<class DATA>
void TestOp<DATA>::operator()(LevelData<DATA>& a_rhs,
        const LevelData<DATA>& a_phi,
        Real a_dx, Real a_scale) const
{
    auto iter = a_phi.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        patch phi = a_phi[iter];
        patch rhs = a_rhs[iter];
        (*this)(rhs, phi, a_dx, a_scale);
    }
}
#ifdef GSRB
template<class DATA>
void TestOp<DATA>::relax(LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_R,
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
