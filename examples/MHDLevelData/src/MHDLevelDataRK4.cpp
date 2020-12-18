#include "MHDLevelDataRK4.H"

MHDLevelDataState::MHDLevelDataState()
{}

MHDLevelDataState::~MHDLevelDataState()
{}

MHDLevelDataState::MHDLevelDataState(const ProblemDomain& a_probDom,
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

void MHDLevelDataState::increment(const MHDLevelDataDX& a_DX)
{
    for(DataIterator dit=a_DX.m_DU.begin(); *dit!=dit.end(); ++dit) {
        m_U[*dit]+=a_DX.m_DU[*dit];
    }
}

MHDLevelDataDX::MHDLevelDataDX()
{}

MHDLevelDataDX::~MHDLevelDataDX()
{}

void MHDLevelDataDX::init(MHDLevelDataState& a_State)
{
    m_dbl=a_State.m_dbl;
    m_DU.define(m_dbl,Point::Zero());
    m_DU.setToZero();
}

void MHDLevelDataDX::increment(double& a_weight, const MHDLevelDataDX& a_incr)
{
    for(DataIterator dit=m_DU.begin(); *dit!=dit.end(); ++dit) {
        const BoxData<double,NUMCOMPS>& incr=a_incr.m_DU[*dit];
        BoxData<double,NUMCOMPS> temp(incr.box());
        incr.copyTo(temp);
        temp*=a_weight;
        m_DU[*dit]+=temp;
    }
}

void MHDLevelDataDX::operator*=(const double& a_weight)
{
    for(DataIterator dit=m_DU.begin(); *dit!=dit.end(); ++dit) {
        m_DU[*dit]*=a_weight;
    }
}

MHDLevelDataRK4Op::MHDLevelDataRK4Op()
{}

MHDLevelDataRK4Op::~MHDLevelDataRK4Op()
{}

void MHDLevelDataRK4Op::operator()(MHDLevelDataDX& a_DX,
                                     double a_time,
                                     double a_dt,
                                     MHDLevelDataState& a_State)
{
    LevelBoxData<double,NUMCOMPS> new_state(a_State.m_dbl,Point::Ones(NGHOST));
    (a_State.m_U).copyTo(new_state);
    for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
        new_state[*dit]+=(a_DX.m_DU)[*dit];
    }
    new_state.exchange();
    for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
        Reduction<double> rxn; //Dummy: not used
        //Set the last two arguments to false so as not to call routines that would don't work in parallel yet
        MHDOp::step(a_DX.m_DU[*dit],new_state[*dit],a_State.m_U[*dit].box(), a_State.m_dx, a_State.m_gamma, rxn, false, false);
    }
    a_DX*=a_dt;
}



MHDLevelDataEulerOp::MHDLevelDataEulerOp()
{}

MHDLevelDataEulerOp::~MHDLevelDataEulerOp()
{}

void MHDLevelDataEulerOp::operator()(MHDLevelDataDX& a_DX,
                                     double a_time,
                                     double a_dt,
                                     MHDLevelDataState& a_State)
{
    LevelBoxData<double,NUMCOMPS> new_state(a_State.m_dbl,Point::Ones(NGHOST));
    (a_State.m_U).copyTo(new_state);
    for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
        new_state[*dit]+=(a_DX.m_DU)[*dit];
    }
    new_state.exchange();
    for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
        Reduction<double> rxn; //Dummy: not used
        //Set the last two arguments to false so as not to call routines that would don't work in parallel yet
        MHDOp::step2(a_DX.m_DU[*dit],new_state[*dit],a_State.m_U[*dit].box(), a_State.m_dx, a_State.m_gamma, rxn, false, false);
    }
    a_DX*=a_dt;
}


MHDLevelDataViscosityOp::MHDLevelDataViscosityOp()
{}

MHDLevelDataViscosityOp::~MHDLevelDataViscosityOp()
{}

void MHDLevelDataViscosityOp::operator()(MHDLevelDataDX& a_DX,
                                     double a_time,
                                     double a_dt,
                                     MHDLevelDataState& a_State)
{
    LevelBoxData<double,NUMCOMPS> new_state(a_State.m_dbl,Point::Ones(NGHOST));
    (a_State.m_U).copyTo(new_state);
    for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
        new_state[*dit]+=(a_DX.m_DU)[*dit];
    }
    new_state.exchange();
    for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
        Reduction<double> rxn; //Dummy: not used
        //Set the last two arguments to false so as not to call routines that would don't work in parallel yet
        MHDOp::step3(a_DX.m_DU[*dit],new_state[*dit],a_State.m_U[*dit].box(), a_State.m_dx, a_State.m_gamma, rxn, false, false);
    }
    a_DX*=a_dt;
}