#include "MHDLevelDataRK4.H"
#include "MHD_Set_Boundary_Values.H"
#include "MHD_Mapping.H"
#include "MHD_Output_Writer.H"

extern int LowBoundType;
extern int HighBoundType;
extern int grid_type_global;
//extern int domainSizex;
//extern int domainSizey;
//extern int domainSizez;
extern int BoxSize;
extern bool pole_correction;

MHDLevelDataState::MHDLevelDataState()
{
}

MHDLevelDataState::~MHDLevelDataState()
{
}

MHDLevelDataState::MHDLevelDataState(const ProblemDomain& a_probDom,
                                     const Point& a_boxSize,
                                     const double a_dx,
                                     const double a_dy,
                                     const double a_dz,
                                     const double a_gamma) :
    m_dx(a_dx),
    m_dy(a_dy),
    m_dz(a_dz),
	m_gamma(a_gamma),
	m_dbl(a_probDom, a_boxSize),
	m_probDom(a_probDom)
	//m_U(m_dbl,Point::Zero())
{
	m_U.define(m_dbl,Point::Zero());
	m_U_old.define(m_dbl,Point::Zero());
	m_Jacobian_ave.define(m_dbl,Point::Ones(NGHOST));
	m_N_ave_f.define(m_dbl,Point::Ones(NGHOST));
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
                                   MHDLevelDataState& a_State
								   )
{



	LevelBoxData<double,NUMCOMPS> new_state(a_State.m_dbl,Point::Ones(NGHOST));
	//LevelBoxData<double,1> Jacobian_ave(a_State.m_dbl,Point::Ones(NGHOST));
    auto idOp = (1.0)*Shift(Point::Zeros());
	(a_State.m_U).copyTo(new_state);
	for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
		//new_state[*dit]+=(a_DX.m_DU)[*dit];
        new_state[*dit]+=idOp((a_DX.m_DU)[*dit]);  //Phil found doing this improves performance
		//MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave[*dit],a_State.m_dx,a_State.m_dy,a_State.m_dz,a_State.m_U[*dit].box());
	}

	//(a_State.m_Jacobian_ave).exchange();
    new_state.exchange();
#if DIM == 3
	if (grid_type_global == 2 && pole_correction){
        LevelBoxData<double,NUMCOMPS> U_pole_long;
        LevelBoxData<double,NUMCOMPS> U_pole2_temp(a_State.m_dbl,Point::Ones(NGHOST));
        LevelBoxData<double,NUMCOMPS> U_pole_long_ahead;
        LevelBoxData<double,NUMCOMPS> rotateddata_pole(a_State.m_dbl,Point::Ones(NGHOST));
        int domainSizex = a_State.m_probDom.box().high()[0] + 1;
        int domainSizey = a_State.m_probDom.box().high()[1] + 1;
        U_pole_long.define(DisjointBoxLayout(a_State.m_probDom,Point(domainSizex, domainSizey, BoxSize)), {{0,domainSizey/2,0}});
        U_pole_long_ahead.define(DisjointBoxLayout(a_State.m_probDom,Point(domainSizex, domainSizey, BoxSize)), {{0,0,0}});


        static Stencil<double> m_right_shift;
        m_right_shift = 1.0*Shift(Point::Basis(2)*(-NGHOST));
		for(DataIterator dit2=new_state.begin(); *dit2!=dit2.end(); ++dit2){
			 U_pole2_temp[*dit2] = m_right_shift(new_state[*dit2]);
		}
        (U_pole2_temp).copyTo(U_pole_long);
		U_pole_long.exchange();
		static Stencil<double> m_ahead_shift;
		m_ahead_shift = 1.0*Shift(Point::Basis(1)*(domainSizey/2));
		for(DataIterator dit2=U_pole_long.begin(); *dit2!=dit2.end(); ++dit2){
			 U_pole_long_ahead[*dit2] = m_ahead_shift(U_pole_long[*dit2]);
		}
		(U_pole_long_ahead).copyTo(rotateddata_pole);
        rotateddata_pole.exchange();


        for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
            MHD_Set_Boundary_Values::Set_Zaxis_Values(new_state[*dit],a_State.m_probDom, rotateddata_pole[*dit]);
        }
	}
#endif

    for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
        //MHD_Set_Boundary_Values::Set_Zaxis_Values(new_state[*dit],a_State.m_U[*dit].box(),a_State.m_probDom, rotateddata[*dit]);
		if (LowBoundType != 0 || HighBoundType != 0) {
			MHD_Set_Boundary_Values::Set_Jacobian_Values((a_State.m_Jacobian_ave)[*dit],a_State.m_U[*dit].box(),a_State.m_probDom,a_State.m_dx,a_State.m_dy,a_State.m_dz, a_State.m_gamma, LowBoundType,HighBoundType);
			MHD_Set_Boundary_Values::Set_Boundary_Values(new_state[*dit],a_State.m_U[*dit].box(),a_State.m_probDom,a_State.m_dx,a_State.m_dy,a_State.m_dz, a_State.m_gamma,(a_State.m_Jacobian_ave)[*dit], LowBoundType,HighBoundType);
		}
    }



    for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
		Reduction<double> rxn; //Dummy: not used
		//Set the last two arguments to false so as not to call routines that would don't work in parallel yet
		MHDOp::step(a_DX.m_DU[*dit],new_state[*dit],a_State.m_U[*dit].box(), a_State.m_dx, a_State.m_dy, a_State.m_dz, a_State.m_gamma, rxn,(a_State.m_Jacobian_ave)[*dit],(a_State.m_N_ave_f)[*dit], false, false);
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
    auto idOp = (1.0)*Shift(Point::Zeros());
	//LevelBoxData<double,1> Jacobian_ave(a_State.m_dbl,Point::Ones(NGHOST));
	(a_State.m_U).copyTo(new_state);

	for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
        //new_state[*dit]+=(a_DX.m_DU)[*dit];
        new_state[*dit]+=idOp((a_DX.m_DU)[*dit]);  //Phil found doing this improves performance
		//MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave[*dit],a_State.m_dx,a_State.m_dy,a_State.m_dz,a_State.m_U[*dit].box());
	}
    new_state.exchange();
	//(a_State.m_Jacobian_ave).exchange();
#if DIM == 3
	if (grid_type_global == 2 && pole_correction){
        LevelBoxData<double,NUMCOMPS> U_pole_long;
        LevelBoxData<double,NUMCOMPS> U_pole2_temp(a_State.m_dbl,Point::Ones(NGHOST));
        LevelBoxData<double,NUMCOMPS> U_pole_long_ahead;
        LevelBoxData<double,NUMCOMPS> rotateddata_pole(a_State.m_dbl,Point::Ones(NGHOST));
        int domainSizex = a_State.m_probDom.box().high()[0] + 1;
        int domainSizey = a_State.m_probDom.box().high()[1] + 1;
        U_pole_long.define(DisjointBoxLayout(a_State.m_probDom,Point(domainSizex, domainSizey, BoxSize)), {{0,domainSizey/2,0}});
        U_pole_long_ahead.define(DisjointBoxLayout(a_State.m_probDom,Point(domainSizex, domainSizey, BoxSize)), {{0,0,0}});


        static Stencil<double> m_right_shift;
        m_right_shift = 1.0*Shift(Point::Basis(2)*(-NGHOST));
		for(DataIterator dit2=new_state.begin(); *dit2!=dit2.end(); ++dit2){
			 U_pole2_temp[*dit2] = m_right_shift(new_state[*dit2]);
		}
        (U_pole2_temp).copyTo(U_pole_long);
		U_pole_long.exchange();
		static Stencil<double> m_ahead_shift;
		m_ahead_shift = 1.0*Shift(Point::Basis(1)*(domainSizey/2));
		for(DataIterator dit2=U_pole_long.begin(); *dit2!=dit2.end(); ++dit2){
			 U_pole_long_ahead[*dit2] = m_ahead_shift(U_pole_long[*dit2]);
		}
		(U_pole_long_ahead).copyTo(rotateddata_pole);
        rotateddata_pole.exchange();


        for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
            MHD_Set_Boundary_Values::Set_Zaxis_Values(new_state[*dit],a_State.m_probDom, rotateddata_pole[*dit]);
        }
	}
#endif
	for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
		if (LowBoundType != 0 || HighBoundType != 0) {
			MHD_Set_Boundary_Values::Set_Jacobian_Values((a_State.m_Jacobian_ave)[*dit],a_State.m_U[*dit].box(),a_State.m_probDom,a_State.m_dx,a_State.m_dy,a_State.m_dz, a_State.m_gamma, LowBoundType,HighBoundType);
			MHD_Set_Boundary_Values::Set_Boundary_Values(new_state[*dit],a_State.m_U[*dit].box(),a_State.m_probDom,a_State.m_dx,a_State.m_dy,a_State.m_dz, a_State.m_gamma,(a_State.m_Jacobian_ave)[*dit], LowBoundType,HighBoundType);
		}
    }



    for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
		Reduction<double> rxn; //Dummy: not used
		//Set the last two arguments to false so as not to call routines that would don't work in parallel yet
		MHD_divB_Cleaning::step(a_DX.m_DU[*dit],new_state[*dit],a_State.m_U[*dit].box(), a_State.m_dx, a_State.m_dy, a_State.m_dz, a_State.m_gamma, rxn, (a_State.m_Jacobian_ave)[*dit],(a_State.m_N_ave_f)[*dit], false, false);
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
	//LevelBoxData<double,1> Jacobian_ave(a_State.m_dbl,Point::Ones(NGHOST));
	//(a_State.m_U).copyTo(new_state);
    auto idOp = (1.0)*Shift(Point::Zeros());
	(a_State.m_U_old).copyTo(new_state);
	for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
        //new_state[*dit]+=(a_DX.m_DU)[*dit];
        new_state[*dit]+=idOp((a_DX.m_DU)[*dit]);  //Phil found doing this improves performance
		//MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave[*dit],a_State.m_dx,a_State.m_dy,a_State.m_dz,a_State.m_U[*dit].box());
	}
	new_state.exchange();
	//(a_State.m_Jacobian_ave).exchange();
#if DIM == 3
	if (grid_type_global == 2 && pole_correction){
        LevelBoxData<double,NUMCOMPS> U_pole_long;
        LevelBoxData<double,NUMCOMPS> U_pole2_temp(a_State.m_dbl,Point::Ones(NGHOST));
        LevelBoxData<double,NUMCOMPS> U_pole_long_ahead;
        LevelBoxData<double,NUMCOMPS> rotateddata_pole(a_State.m_dbl,Point::Ones(NGHOST));
        int domainSizex = a_State.m_probDom.box().high()[0] + 1;
        int domainSizey = a_State.m_probDom.box().high()[1] + 1;
        U_pole_long.define(DisjointBoxLayout(a_State.m_probDom,Point(domainSizex, domainSizey, BoxSize)), {{0,domainSizey/2,0}});
        U_pole_long_ahead.define(DisjointBoxLayout(a_State.m_probDom,Point(domainSizex, domainSizey, BoxSize)), {{0,0,0}});


        static Stencil<double> m_right_shift;
        m_right_shift = 1.0*Shift(Point::Basis(2)*(-NGHOST));
		for(DataIterator dit2=new_state.begin(); *dit2!=dit2.end(); ++dit2){
			 U_pole2_temp[*dit2] = m_right_shift(new_state[*dit2]);
		}
        (U_pole2_temp).copyTo(U_pole_long);
		U_pole_long.exchange();
		static Stencil<double> m_ahead_shift;
		m_ahead_shift = 1.0*Shift(Point::Basis(1)*(domainSizey/2));
		for(DataIterator dit2=U_pole_long.begin(); *dit2!=dit2.end(); ++dit2){
			 U_pole_long_ahead[*dit2] = m_ahead_shift(U_pole_long[*dit2]);
		}
		(U_pole_long_ahead).copyTo(rotateddata_pole);
        rotateddata_pole.exchange();


        for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
            MHD_Set_Boundary_Values::Set_Zaxis_Values(new_state[*dit],a_State.m_probDom, rotateddata_pole[*dit]);
        }
	}
#endif
	for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
		if (LowBoundType != 0 || HighBoundType != 0) {
			MHD_Set_Boundary_Values::Set_Jacobian_Values((a_State.m_Jacobian_ave)[*dit],a_State.m_U[*dit].box(),a_State.m_probDom,a_State.m_dx,a_State.m_dy,a_State.m_dz, a_State.m_gamma, LowBoundType,HighBoundType);
			MHD_Set_Boundary_Values::Set_Boundary_Values(new_state[*dit],a_State.m_U[*dit].box(),a_State.m_probDom,a_State.m_dx,a_State.m_dy,a_State.m_dz, a_State.m_gamma,(a_State.m_Jacobian_ave)[*dit], LowBoundType,HighBoundType);
		}
    }


    for(DataIterator dit=new_state.begin(); *dit!=dit.end(); ++dit) {
		Reduction<double> rxn; //Dummy: not used
		//Set the last two arguments to false so as not to call routines that would don't work in parallel yet
		MHD_Artificial_Viscosity::step(a_DX.m_DU[*dit],new_state[*dit],a_State.m_U[*dit].box(), a_State.m_dx, a_State.m_dy, a_State.m_dz, a_State.m_gamma, rxn, (a_State.m_Jacobian_ave)[*dit],(a_State.m_N_ave_f)[*dit], false, false);
	}
	a_DX*=a_dt;
}
