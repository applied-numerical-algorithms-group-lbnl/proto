#include "Proto.H"
#include "MHD_Set_Boundary_Values.H"
#include "CommonTemplates.H"
#include "Proto_Timer.H"
#include "Proto_WriteBoxData.H"
#include "Proto_LevelBoxData.H"
#include "Proto_ProblemDomain.H"
#include "MHDOp.H"
#include "MHD_Mapping.H"
#include "MHD_Initialize.H"
#include "MHD_Input_Parsing.H"
typedef BoxData<double,1,HOST> Scalar;
typedef BoxData<double,NUMCOMPS,HOST> Vector;

extern Parsefrominputs inputs;

namespace MHD_Set_Boundary_Values {

	PROTO_KERNEL_START
	void dot_pro_calcFF(State& a_dot_pro,
	                    const Var<double,1>& a_d_perp_N_s,
	                    const State& a_d_perp_F)
	{
		for (int i=0; i< NUMCOMPS; i++) {
			a_dot_pro(i) = (a_d_perp_N_s(0)*a_d_perp_F(i));
		}
	}
	PROTO_KERNEL_END(dot_pro_calcFF, dot_pro_calcF)


	PROTO_KERNEL_START
	void iotaFuncF(Point           & a_p,
	               V               & a_X,
	               const double a_dx,
	               const double a_dy,
	               const double a_dz)
	{
		for (int ii = 0; ii < DIM; ii++)
		{
			double dxd;
			if (ii == 0) dxd = a_dx;
			if (ii == 1) dxd = a_dy;
			if (ii == 2) dxd = a_dz;
			a_X(ii) = a_p[ii]*dxd + 0.5*dxd;
		}
	}
	PROTO_KERNEL_END(iotaFuncF,iotaFunc)




	void Set_Boundary_Values(BoxData<double,NUMCOMPS>& a_JU,
	                         const Box& a_rangeBox,
	                         const ProblemDomain& a_probDomain,
	                         const double a_dx,
	                         const double a_dy,
	                         const double a_dz,
	                         const double a_gamma,
	                         const BoxData<double,1>& Jacobian_ave,
							 BoxData<double,DIM*DIM>& a_detAA_avg,
							 BoxData<double,DIM*DIM>& a_detAA_inv_avg,
	                  		 BoxData<double,1>& a_r2rdot_avg,
	                 		 BoxData<double,1>& a_detA_avg,
	                         const int a_lowBCtype,
	                         const int a_highBCtype)
	{


		Box dbx0 = a_JU.box();
		//Filling ghost cells for low side of dir == 0 here. This will be the inner boundary in r direction once we map to polar, spherical, or cubed sphere grids.
		if (dbx0.low()[0] < a_probDomain.box().low()[0] && a_lowBCtype == 1) {
			Point ghost_low = dbx0.low();
#if DIM == 2
			Point ghost_high = Point(a_probDomain.box().low()[0]-1, dbx0.high()[1]);
#endif
#if DIM == 3
			Point ghost_high = Point(a_probDomain.box().low()[0]-1, dbx0.high()[1], dbx0.high()[2]);
#endif
			Box BoundBox(ghost_low,ghost_high);
			Vector a_U(BoundBox), a_JU_ghost(BoundBox);
			BoxData<double,DIM> eta(BoundBox);
			forallInPlace_p(iotaFunc, BoundBox, eta, a_dx, a_dy, a_dz);
			BoxData<double,DIM> x(BoundBox);
			MHD_Mapping::eta_to_x_calc(x,eta);
			MHD_Initialize::InitializeStatecalc(a_U,x,eta,a_gamma,BoundBox);
			
			if (inputs.grid_type_global == 2 && inputs.initialize_in_spherical_coords == 1){
				MHD_Mapping::U_Sph_ave_to_JU_calc_func(a_JU_ghost, a_U, a_detAA_avg, a_r2rdot_avg, a_detA_avg);
			} else {
				a_JU_ghost = forall<double,NUMCOMPS>(dot_pro_calcF, Jacobian_ave, a_U);
			}
			// a_JU_ghost = forall<double,NUMCOMPS>(dot_pro_calcF, Jacobian_ave, a_U);
			a_JU_ghost.copyTo(a_JU,BoundBox);
		}


		//Filling ghost cells for high side of dir == 0 here. This will be the outer boundary in r direction once we map to polar, spherical, or cubed sphere grids.
		if (dbx0.high()[0] > a_probDomain.box().high()[0] && a_highBCtype == 1) {
#if DIM == 2
			Point ghost_low = Point(a_probDomain.box().high()[0]+1, dbx0.low()[1]);
#endif
#if DIM == 3
			Point ghost_low = Point(a_probDomain.box().high()[0]+1, dbx0.low()[1], dbx0.low()[2]);
#endif
			Point ghost_high = dbx0.high();
			// cout << ghost_low[0] << " " << ghost_low[1] << " " << ghost_low[2] << " " << ghost_high[0] << " " << ghost_high[1] << " " << ghost_high[2] << endl;
			Box BoundBox(ghost_low,ghost_high);
			Vector a_U(BoundBox), a_JU_ghost(BoundBox);;
			BoxData<double,DIM> eta(BoundBox);
			forallInPlace_p(iotaFunc, BoundBox, eta, a_dx, a_dy, a_dz);
			BoxData<double,DIM> x(BoundBox);
			MHD_Mapping::eta_to_x_calc(x,eta);
			MHD_Initialize::InitializeStatecalc(a_U,x,eta,a_gamma,BoundBox);
			if (inputs.grid_type_global == 2 && inputs.initialize_in_spherical_coords == 1){
				MHD_Mapping::U_Sph_ave_to_JU_calc_func(a_JU_ghost, a_U, a_detAA_avg, a_r2rdot_avg, a_detA_avg);
			} else {
				a_JU_ghost = forall<double,NUMCOMPS>(dot_pro_calcF, Jacobian_ave, a_U);
			}
			// a_JU_ghost = forall<double,NUMCOMPS>(dot_pro_calcF, Jacobian_ave, a_U);
			a_JU_ghost.copyTo(a_JU,BoundBox);
		}



		//Filling ghost cells for low side of dir == 0 here, to make an open boundary. Ghost cells are filled by innermost cells of problem domain. This will be the inner boundary in r direction once we map to polar, spherical, or cubed sphere grids.
		if (dbx0.low()[0] < a_probDomain.box().low()[0] && a_lowBCtype == 2) {
#if DIM == 2
			Point source_low = Point(a_probDomain.box().low()[0],a_JU.box().low()[1]);
			Point source_high = Point(a_probDomain.box().low()[0],a_JU.box().high()[1]);
#endif
#if DIM == 3
			Point source_low = Point(a_probDomain.box().low()[0],a_JU.box().low()[1],a_JU.box().low()[2]);
			Point source_high = Point(a_probDomain.box().low()[0],a_JU.box().high()[1],a_JU.box().high()[2]);
#endif
			Box sourceBox(source_low,source_high);
			Box sourceBox1 = sourceBox.grow(1);

			Vector a_U(sourceBox1);

			// if (inputs.grid_type_global == 2 && inputs.initialize_in_spherical_coords == 1){
			// 	MHD_Mapping::JU_to_U_Sph_ave_calc_func(a_U, a_JU, a_detAA_inv_avg, a_r2rdot_avg, a_detA_avg, false);
			// } else {
			// 	MHD_Mapping::JU_to_U_2ndOrdercalc(a_U, a_JU, Jacobian_ave, sourceBox1);
			// }
			MHD_Mapping::JU_to_U_2ndOrdercalc(a_U, a_JU, Jacobian_ave, sourceBox1);
			
			Point ghost_low = dbx0.low();
#if DIM == 2
			Point ghost_high = Point(a_probDomain.box().low()[0]-1, dbx0.high()[1]);
#endif
#if DIM == 3
			Point ghost_high = Point(a_probDomain.box().low()[0]-1, dbx0.high()[1], dbx0.high()[2]);
#endif
			Box BoundBox(ghost_low,ghost_high);
			// Box BoundBox1 = BoundBox.grow(1);
			Box BoundBox1 = BoundBox.grow(0);
			Vector a_U_ghost(BoundBox), a_JU_ghost(BoundBox);
			Scalar Jacobian_ave2(BoundBox1);
			if (inputs.grid_type_global == 2){
				MHD_Mapping::Jacobian_ave_sph_calc_func(Jacobian_ave2,a_dx, a_dy, a_dz);
			} else {
				MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave2,a_dx, a_dy, a_dz,BoundBox1);
			}
			for (int i = 1; i <= NGHOST; i++ ) {
				a_U.copyTo(a_U_ghost,sourceBox,Point::Basis(0)*(-i));// Using shifting option of copyTo
			}
			// if (inputs.grid_type_global == 2 && inputs.initialize_in_spherical_coords == 1){
			// 	MHD_Mapping::U_Sph_ave_to_JU_calc_func(a_JU_ghost, a_U_ghost, a_detAA_avg, a_r2rdot_avg, a_detA_avg);
			// } else {
			// 	a_JU_ghost = forall<double,NUMCOMPS>(dot_pro_calcF, Jacobian_ave2, a_U_ghost);
			// }
			a_JU_ghost = forall<double,NUMCOMPS>(dot_pro_calcF, Jacobian_ave2, a_U_ghost);
			a_JU_ghost.copyTo(a_JU,BoundBox);
		}





		//Filling ghost cells for high side of dir == 0 here, to make an open boundary. Ghost cells are filled by outermost cells of problem domain. This will be the outer boundary in r direction once we map to polar, spherical, or cubed sphere grids.
		if (dbx0.high()[0] > a_probDomain.box().high()[0] && a_highBCtype == 2) {
#if DIM == 2
			Point source_low = Point(a_probDomain.box().high()[0],a_JU.box().low()[1]);
			Point source_high = Point(a_probDomain.box().high()[0],a_JU.box().high()[1]);
#endif
#if DIM == 3
			Point source_low = Point(a_probDomain.box().high()[0],a_JU.box().low()[1],a_JU.box().low()[2]);
			Point source_high = Point(a_probDomain.box().high()[0],a_JU.box().high()[1],a_JU.box().high()[2]);
#endif
			Box sourceBox(source_low,source_high);
			Box sourceBox1 = sourceBox.grow(1);

			Vector a_U(sourceBox1);
			
			if (inputs.grid_type_global == 2 && inputs.initialize_in_spherical_coords == 1){
				MHD_Mapping::JU_to_U_Sph_ave_calc_func(a_U, a_JU, a_detAA_inv_avg, a_r2rdot_avg, a_detA_avg, false);
			} else {
				MHD_Mapping::JU_to_U_2ndOrdercalc(a_U, a_JU, Jacobian_ave, sourceBox1);
			}
			// MHD_Mapping::JU_to_U_2ndOrdercalc(a_U, a_JU, Jacobian_ave, sourceBox1);

#if DIM == 2
			Point ghost_low = Point(a_probDomain.box().high()[0]+1, dbx0.low()[1]);
#endif
#if DIM == 3
			Point ghost_low = Point(a_probDomain.box().high()[0]+1, dbx0.low()[1], dbx0.low()[2]);
#endif
			Point ghost_high = dbx0.high();
			Box BoundBox(ghost_low,ghost_high);
			// Box BoundBox1 = BoundBox.grow(1);
			Box BoundBox1 = BoundBox.grow(0);
			Vector a_U_ghost(BoundBox), a_JU_ghost(BoundBox);
			Scalar Jacobian_ave2(BoundBox1);
			if (inputs.grid_type_global == 2){
				MHD_Mapping::Jacobian_ave_sph_calc_func(Jacobian_ave2,a_dx, a_dy, a_dz);
			} else {
				MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave2,a_dx, a_dy, a_dz,BoundBox1);
			}

			for (int i = 1; i <= NGHOST; i++ ) {
				a_U.copyTo(a_U_ghost,sourceBox,Point::Basis(0)*(i));// Using shifting option of copyTo
			}
			if (inputs.grid_type_global == 2 && inputs.initialize_in_spherical_coords == 1){
				MHD_Mapping::U_Sph_ave_to_JU_calc_func(a_JU_ghost, a_U_ghost, a_detAA_avg, a_r2rdot_avg, a_detA_avg);
			} else {
				a_JU_ghost = forall<double,NUMCOMPS>(dot_pro_calcF, Jacobian_ave2, a_U_ghost);
			}
			// a_JU_ghost = forall<double,NUMCOMPS>(dot_pro_calcF, Jacobian_ave2, a_U_ghost);
			a_JU_ghost.copyTo(a_JU,BoundBox);
		}


	}


	void Set_Jacobian_Values(BoxData<double,1>& a_Jacobian_ave,
	                         const Box& a_rangeBox,
	                         const ProblemDomain& a_probDomain,
	                         const double a_dx,
	                         const double a_dy,
	                         const double a_dz,
	                         const double a_gamma,
	                         const int a_lowBCtype,
	                         const int a_highBCtype)
	{


		Box dbx0 = a_Jacobian_ave.box();
		//Filling Jacobian values in ghost cells for low side of dir == 0 here. This will be the inner boundary in r direction once we map to polar, spherical, or cubed sphere grids.
		if (dbx0.low()[0] < a_probDomain.box().low()[0] && a_lowBCtype != 0) {
			Point ghost_low = dbx0.low();
#if DIM == 2
			Point ghost_high = Point(a_probDomain.box().low()[0]-1, dbx0.high()[1]);
#endif
#if DIM == 3
			Point ghost_high = Point(a_probDomain.box().low()[0]-1, dbx0.high()[1], dbx0.high()[2]);
#endif
			Box BoundBox(ghost_low,ghost_high);
			Box BoundBox1 = BoundBox.grow(1);
			Scalar Jacobian_ave2(BoundBox1);
			if (inputs.grid_type_global == 2){
				MHD_Mapping::Jacobian_ave_sph_calc_func(Jacobian_ave2,a_dx, a_dy, a_dz);
			} else {
				MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave2,a_dx, a_dy, a_dz,BoundBox1);
			}
			Jacobian_ave2.copyTo(a_Jacobian_ave,BoundBox);
		}

		//Filling Jacobian values in ghost cells for high side of dir == 0 here. This will be the outer boundary in r direction once we map to polar, spherical, or cubed sphere grids.
		if (dbx0.high()[0] > a_probDomain.box().high()[0] && a_highBCtype != 0) {
#if DIM == 2
			Point ghost_low = Point(a_probDomain.box().high()[0]+1, dbx0.low()[1]);
#endif
#if DIM == 3
			Point ghost_low = Point(a_probDomain.box().high()[0]+1, dbx0.low()[1], dbx0.low()[2]);
#endif
			Point ghost_high = dbx0.high();
			Box BoundBox(ghost_low,ghost_high);
			Box BoundBox1 = BoundBox.grow(1);
			Scalar Jacobian_ave2(BoundBox1);
			if (inputs.grid_type_global == 2){
				MHD_Mapping::Jacobian_ave_sph_calc_func(Jacobian_ave2,a_dx, a_dy, a_dz);
			} else {
				MHD_Mapping::Jacobian_Ave_calc(Jacobian_ave2,a_dx, a_dy, a_dz,BoundBox1);
			}
			Jacobian_ave2.copyTo(a_Jacobian_ave,BoundBox);
		}
	}


	void Set_Zaxis_Values(BoxData<double,NUMCOMPS>& a_JU,
						  const ProblemDomain& a_probDomain,
	                      const BoxData<double,NUMCOMPS>& a_JU_pole)
	{

		Box dbx0 = a_JU.box();
		//Filling Jacobian values in ghost cells for theta = 0 degrees here.
		if (dbx0.low()[1] < a_probDomain.box().low()[1]) {
			Point ghost_low = dbx0.low();
			Point ghost_high = Point(dbx0.high()[0], a_probDomain.box().low()[1]-1, dbx0.high()[2] );
			
			Box BoundBox(ghost_low,ghost_high);
			Vector a_JU_ghost(BoundBox);
			for (int i = 1; i <= NGHOST; i++ ) {
				Point source_low = Point(a_JU.box().low()[0],a_probDomain.box().low()[1]+NGHOST+i-1, a_JU.box().low()[2]);
				Point source_high = Point(a_JU.box().high()[0], a_probDomain.box().low()[1]+NGHOST+i-1, a_JU.box().high()[2]);
				Box sourceBox(source_low,source_high);
				a_JU_pole.copyTo(a_JU_ghost,sourceBox,Point::Basis(1)*(-(2*i-1+NGHOST)));// Using shifting option of copyTo
			}
			a_JU_ghost.copyTo(a_JU,BoundBox);
		}

		//Filling Jacobian values in ghost cells for theta = 180 degrees here.
		if (dbx0.high()[1] > a_probDomain.box().high()[1]) {
			Point ghost_low = Point(dbx0.low()[0], a_probDomain.box().high()[1]+1, dbx0.low()[2]);
			Point ghost_high = dbx0.high();
			Box BoundBox(ghost_low,ghost_high);
			Vector a_JU_ghost(BoundBox);
			for (int i = 1; i <= NGHOST; i++ ) {
				Point source_low = Point(a_JU.box().low()[0], a_probDomain.box().high()[1]-i+1+NGHOST, a_JU.box().low()[2]);
				Point source_high = Point(a_JU.box().high()[0], a_probDomain.box().high()[1]-i+1+NGHOST, a_JU.box().high()[2]);
				Box sourceBox(source_low,source_high);
				a_JU_pole.copyTo(a_JU_ghost,sourceBox,Point::Basis(1)*((2*i-(NGHOST+1))));// Using shifting option of copyTo
			}
			a_JU_ghost.copyTo(a_JU,BoundBox);
		}


	}
}
