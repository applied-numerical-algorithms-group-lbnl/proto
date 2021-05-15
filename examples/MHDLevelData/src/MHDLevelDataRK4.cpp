#include "MHDLevelDataRK4.H"
#include "MHD_Set_Boundary_Values.H"
#include "MHD_Mapping.H"

extern int LowBoundType;
extern int HighBoundType;

MHDLevelDataState::MHDLevelDataState()
{
}

MHDLevelDataState::~MHDLevelDataState()
{
}

MHDLevelDataState::MHDLevelDataState(const ProblemDomain& a_probDom,
                                     const Point& a_boxSize,
                                     const double a_dx,
                                     const double a_gamma) :
	m_dx(a_dx),
	m_gamma(a_gamma),
	m_dbl(a_probDom, a_boxSize),//,
	m_probDom(a_probDom)
	//m_U(m_dbl,Point::Zero())
{
	m_U.define(m_dbl,Point::Zero());
	m_U_old.define(m_dbl,Point::Zero());
}

void MHDLevelDataState::increment(const MHDLevelDataDX& a_DX)
{
	m_U_old.setToZero();
	for(DataIterator dit=a_DX.m_DU.begin(); *dit!=dit.end(); ++dit) {
		m_U_old[*dit]+=m_U[*dit]; // This is needed in artificial viscosity calculation
		m_U[*dit]+=a_DX.m_DU[*dit];
	}
}

MHDLevelDataDX::MHDLevelDataDX()
{
}

MHDLevelDataDX::~MHDLevelDataDX()
{
}

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
{
}

MHDLevelDataRK4Op::~MHDLevelDataRK4Op()
{
}

void MHDLevelDataRK4Op::operator()(MHDLevelDataDX& a_DX,
                                   double a_time,
                                   double a_dt,
                                   MHDLevelDataState& a_State)
{
	LevelBoxData<double,NUMCOMPS> new_state(a_State.m_dbl,Point::Ones(NGHOST));
	LevelBoxData<double,1> Jacobian_ave(a_State.m_dbl,Point::Ones(NGHOST));
	(a_State.m_U).copyTo(new_state);
	for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
		new_state[*dit]+=(a_DX.m_DU)[*dit];
		MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave[*dit],a_State.m_dx,a_State.m_U[*dit].box());
	}
	new_state.exchange();
	Jacobian_ave.exchange();
	for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
		if (LowBoundType != 0 || HighBoundType != 0) {
			MHD_Set_Boundary_Values::Set_Jacobian_Values(Jacobian_ave[*dit],a_State.m_U[*dit].box(),a_State.m_probDom,a_State.m_dx, a_State.m_gamma, LowBoundType,HighBoundType);
			MHD_Set_Boundary_Values::Set_Boundary_Values(new_state[*dit],a_State.m_U[*dit].box(),a_State.m_probDom,a_State.m_dx, a_State.m_gamma,Jacobian_ave[*dit], LowBoundType,HighBoundType);
		}
		Reduction<double> rxn; //Dummy: not used
		//Set the last two arguments to false so as not to call routines that would don't work in parallel yet
		MHDOp::step(a_DX.m_DU[*dit],new_state[*dit],a_State.m_U[*dit].box(), a_State.m_dx, a_State.m_gamma, rxn,Jacobian_ave[*dit], false, false);
	}
	a_DX*=a_dt;
}



MHDLevelDatadivBOp::MHDLevelDatadivBOp()
{
}

MHDLevelDatadivBOp::~MHDLevelDatadivBOp()
{
}

void MHDLevelDatadivBOp::operator()(MHDLevelDataDX& a_DX,
                                    double a_time,
                                    double a_dt,
                                    MHDLevelDataState& a_State)
{
	LevelBoxData<double,NUMCOMPS> new_state(a_State.m_dbl,Point::Ones(NGHOST));
	LevelBoxData<double,1> Jacobian_ave(a_State.m_dbl,Point::Ones(NGHOST));
	(a_State.m_U).copyTo(new_state);

	for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
		new_state[*dit]+=(a_DX.m_DU)[*dit];
		MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave[*dit],a_State.m_dx,a_State.m_U[*dit].box());
	}
	new_state.exchange();
	Jacobian_ave.exchange();
	for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
		if (LowBoundType != 0 || HighBoundType != 0) {
			MHD_Set_Boundary_Values::Set_Jacobian_Values(Jacobian_ave[*dit],a_State.m_U[*dit].box(),a_State.m_probDom,a_State.m_dx, a_State.m_gamma, LowBoundType,HighBoundType);
			MHD_Set_Boundary_Values::Set_Boundary_Values(new_state[*dit],a_State.m_U[*dit].box(),a_State.m_probDom,a_State.m_dx, a_State.m_gamma,Jacobian_ave[*dit], LowBoundType,HighBoundType);
		}

		Reduction<double> rxn; //Dummy: not used
		//Set the last two arguments to false so as not to call routines that would don't work in parallel yet
		MHD_divB_Cleaning::step(a_DX.m_DU[*dit],new_state[*dit],a_State.m_U[*dit].box(), a_State.m_dx, a_State.m_gamma, rxn, Jacobian_ave[*dit], false, false);
	}
	a_DX*=a_dt;
}


MHDLevelDataViscosityOp::MHDLevelDataViscosityOp()
{
}

MHDLevelDataViscosityOp::~MHDLevelDataViscosityOp()
{
}

void MHDLevelDataViscosityOp::operator()(MHDLevelDataDX& a_DX,
                                         double a_time,
                                         double a_dt,
                                         MHDLevelDataState& a_State)
{
	LevelBoxData<double,NUMCOMPS> new_state(a_State.m_dbl,Point::Ones(NGHOST));
	LevelBoxData<double,1> Jacobian_ave(a_State.m_dbl,Point::Ones(NGHOST));
	//(a_State.m_U).copyTo(new_state);
	(a_State.m_U_old).copyTo(new_state);
	for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
		new_state[*dit]+=(a_DX.m_DU)[*dit];
		MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave[*dit],a_State.m_dx,a_State.m_U[*dit].box());
	}
	new_state.exchange();
	Jacobian_ave.exchange();
	for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
		if (LowBoundType != 0 || HighBoundType != 0) {
			MHD_Set_Boundary_Values::Set_Jacobian_Values(Jacobian_ave[*dit],a_State.m_U[*dit].box(),a_State.m_probDom,a_State.m_dx, a_State.m_gamma, LowBoundType,HighBoundType);
			MHD_Set_Boundary_Values::Set_Boundary_Values(new_state[*dit],a_State.m_U[*dit].box(),a_State.m_probDom,a_State.m_dx, a_State.m_gamma,Jacobian_ave[*dit], LowBoundType,HighBoundType);
		}
		Reduction<double> rxn; //Dummy: not used
		//Set the last two arguments to false so as not to call routines that would don't work in parallel yet
		MHD_Artificial_Viscosity::step(a_DX.m_DU[*dit],new_state[*dit],a_State.m_U[*dit].box(), a_State.m_dx, a_State.m_gamma, rxn, Jacobian_ave[*dit], false, false);
	}
	a_DX*=a_dt;
}
