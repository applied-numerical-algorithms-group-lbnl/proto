#pragma once

#include "Proto.H"

#include "MHDOp.H"

using namespace Proto;

namespace MHD_Artificial_Viscosity {

// Used to implement artificial viscosity
	void step(BoxData<double,NUMCOMPS>& a_Rhs,
	          const BoxData<double,NUMCOMPS>& a_U,
	          const Box& a_rangeBox,
	          const double a_dx,
	          const double a_dy,
	          const double a_dz,
	          const double a_gamma,
	          Reduction<double>& a_Rxn,
	          BoxData<double,1>& a_Jacobian_ave,
			  BoxData<double,DIM*DIM>& a_N_ave_f,
	          bool a_computeMaxWaveSpeed,
	          bool a_callBCs);

	void step_spherical(BoxData<double,NUMCOMPS>& a_Rhs,
	                    const BoxData<double,NUMCOMPS>& a_JU,
	                    const Box& a_rangeBox,
	                    const double a_dx,
	                    const double a_dy,
	                    const double a_dz,
	                    const double a_gamma,
	                    Reduction<double>& a_Rxn,
	                    BoxData<double,1>& a_Jacobian_ave,
	                    BoxData<double,DIM*DIM>& a_N_ave_f,
						BoxData<double,DIM*DIM>& a_A_1_avg,
						BoxData<double,DIM*DIM>& a_A_2_avg,
						BoxData<double,DIM*DIM>& a_A_3_avg,
	                    BoxData<double,DIM*DIM>& a_detAA_avg,
	                    BoxData<double,DIM*DIM>& a_detAA_inv_avg,
	                    BoxData<double,1>& a_r2rdot_avg,
	                    BoxData<double,1>& a_detA_avg,
	                    BoxData<double,1>& a_r2detA_1_avg,
	                    BoxData<double,DIM*DIM>& a_r2detAA_1_avg,
	                    BoxData<double,DIM>& a_r2detAn_1_avg,
	                    BoxData<double,1>& a_rrdotdetA_2_avg,
	                    BoxData<double,DIM*DIM>& a_rrdotdetAA_2_avg,
	                    BoxData<double,DIM>& a_rrdotd3ncn_2_avg,
	                    BoxData<double,1>& a_rrdotdetA_3_avg,
	                    BoxData<double,DIM*DIM>& a_rrdotdetAA_3_avg,
	                    BoxData<double,DIM>& a_rrdotncd2n_3_avg,
	                    bool a_computeMaxWaveSpeed,
	                    bool a_callBCs);
}
