#include "EulerLevelDataRK4.H"

EulerLevelDataState::EulerLevelDataState()
{}

EulerLevelDataState::~EulerLevelDataState()
{}

EulerLevelDataState::EulerLevelDataState(DisjointBoxLayout& dbl):
    m_dbl(dbl),
    m_U(m_dbl,Point::Zero())
{}

void EulerLevelDataState::increment(const EulerLevelDataDX& a_DX)
{
    for(DataIterator dit=a_DX.m_DU.begin(); *dit!=dit.end(); ++dit) {
        m_U[*dit]+=a_DX.m_DU[*dit];
    }
}

EulerLevelDataDX::EulerLevelDataDX()
{}

EulerLevelDataDX::~EulerLevelDataDX()
{}

void EulerLevelDataDX::init(EulerLevelDataState& a_State)
{
    m_dbl=a_State.m_dbl;
    m_DU.define(m_dbl,Point::Zero());
    m_DU.setToZero();
}

void EulerLevelDataDX::increment(double& a_weight, const EulerLevelDataDX& a_incr)
{
    for(DataIterator dit=m_DU.begin(); *dit!=dit.end(); ++dit) {
        const BoxData<double,NCOMPS>& incr=a_incr.m_DU[*dit];
        BoxData<double,NCOMPS> temp(incr.box());
        incr.copyTo(temp);
        temp*=a_weight;
        m_DU[*dit]+=temp;
    }
}

void EulerLevelDataDX::operator*=(const double& a_weight)
{
    for(DataIterator dit=m_DU.begin(); *dit!=dit.end(); ++dit) {
        m_DU[*dit]*=a_weight;
    }
}

EulerLevelDataRK4Op::EulerLevelDataRK4Op()
{}

EulerLevelDataRK4Op::~EulerLevelDataRK4Op()
{}

void EulerLevelDataRK4Op::operator()(EulerLevelDataDX& a_DX,
                                     double a_time,
                                     double a_dt,
                                     EulerLevelDataState& a_State)
{
}



