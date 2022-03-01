#pragma once

#include "Proto.H"

#include "MHDOp.H"

using namespace Proto;

namespace MHD_divB_Cleaning {
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
}
