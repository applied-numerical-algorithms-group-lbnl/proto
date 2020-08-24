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
        const BoxData<double,NUMCOMPS>& incr=a_incr.m_DU[*dit];
        BoxData<double,NUMCOMPS> temp(incr.box());
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
    (a_DX.m_DU).setToZero();
    LevelBoxData<double,NUMCOMPS> new_state(a_State.m_dbl,Point::Ones(NGHOST));
    (a_State.m_U).copyTo(new_state);
    for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
        new_state[*dit]+=(a_DX.m_DU)[*dit];
    }
    new_state.exchange();
    for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
        Reduction<double> rxn; //Dummy: not used
        //Set the last two arguments to false so as not to call routines that would don't work in parallel yet
        EulerOp::step(a_DX.m_DU[*dit],new_state[*dit],a_State.m_U[*dit].box(), rxn, false, false);
    }
}



