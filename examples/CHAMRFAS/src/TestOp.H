#include "BaseOp.H"

template<class DATA>
class TestOp : public BaseOp<1, DATA> {
    public:
        typedef Proto::BoxData<Real> patch; 

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
    private:
        Proto::Stencil<Real> m_flux[DIM];
        Proto::Stencil<Real> m_div[DIM]; 
};

template<class DATA>
void TestOp<DATA>::init()
{
    this->m_lambda = -(this->m_dx*this->m_dx)/(4.0*DIM);
#if DIM == 2
    //this->m_lambda = -0.5*(this->m_dx*this->m_dx)*(6.0/20.0);  //9-Point
#elif DIM == 3
    //this->m_lambda = -0.5*(this->m_dx*this->m_dx)*(15.0/64.0); //27-Point
#endif
    for (int ii = 0; ii < DIM; ii++)
    {
        m_flux[ii] = 1.0*Proto::Shift::Zeros()
                   - 1.0*Proto::Shift::Basis(ii,-1);
        m_div[ii] = 1.0*Proto::Shift::Basis(ii,1)
                  - 1.0*Proto::Shift::Zeros();
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
