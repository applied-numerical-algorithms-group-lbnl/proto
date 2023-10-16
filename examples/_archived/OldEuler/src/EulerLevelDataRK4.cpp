#include "EulerLevelDataRK4.H"

EulerLevelDataState::EulerLevelDataState()
{}

EulerLevelDataState::~EulerLevelDataState()
{}

EulerLevelDataState::EulerLevelDataState(const ProblemDomain& a_probDom,
        const Point& a_boxSize,
        const double a_dx,
        const double a_gamma):
    m_dx(a_dx),
    m_gamma(a_gamma),
    m_dbl(a_probDom, a_boxSize)//,
    //m_U(m_dbl,Point::Zero())
{
    m_U.define(m_dbl,Point::Zero());
}

void EulerLevelDataState::increment(const EulerLevelDataDX& a_DX)
{
    for (auto dit : a_DX.m_DU)
    {
        m_U[dit]+=a_DX.m_DU[dit];
    }
    /*
       for(DataIterator dit=a_DX.m_DU.begin(); *dit!=dit.end(); ++dit) {
       m_U[*dit]+=a_DX.m_DU[*dit];
       }
     */
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
    //for(DataIterator dit=m_DU.begin(); *dit!=dit.end(); ++dit) {
    for(auto dit : m_DU) {
        const BoxData<double,NUMCOMPS>& incr=a_incr.m_DU[dit];
        BoxData<double,NUMCOMPS> temp(incr.box());
        incr.copyTo(temp);
        temp*=a_weight;
        m_DU[dit]+=temp;
    }
}

void EulerLevelDataDX::operator*=(const double& a_weight)
{
    //for(DataIterator dit=m_DU.begin(); *dit!=dit.end(); ++dit) {
    for (auto dit : m_DU){
        m_DU[dit]*=a_weight;
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
    PR_TIME("EulerRK4::operator");
    auto idOp = (1.0)*Shift(Point::Zeros());
    LevelBoxData<double,NUMCOMPS> new_state(a_State.m_dbl,Point::Ones(NGHOST));
    (a_State.m_U).copyTo(new_state);
    
    for (auto dit : new_state) {
        new_state[dit]+=idOp((a_DX.m_DU)[dit]);
    }
    new_state.exchange();
    for (auto dit : new_state) {
        Reduction<double, Abs> rxn; //Dummy: not used
        //Set the last two arguments to false so as not to call routines that would don't work in parallel yet
        EulerOp::step(a_DX.m_DU[dit],new_state[dit],a_State.m_U[dit].box(), a_State.m_dx, a_State.m_gamma, rxn, false, false);
    }
    a_DX*=a_dt;
}



