#pragma once

#include "Proto.H"
#include "MHDOp.H"
//TODO: These variables need to go into a common file or be defined in a more accessible way
//constexpr unsigned int NUMCOMPS=2*DIM+2;
typedef Var<double,NUMCOMPS> State;


namespace MHD_Limiters {
	void MHD_Limiters(BoxData<double,NUMCOMPS>& a_W_ave_low_lim_flat,
	                  BoxData<double,NUMCOMPS>& a_W_ave_high_lim_flat,
	                  BoxData<double,NUMCOMPS>& a_W_ave_low,
	                  BoxData<double,NUMCOMPS>& a_W_ave_high,
	                  BoxData<double,NUMCOMPS>& a_W_ave,
	                  BoxData<double,NUMCOMPS>& a_W_bar,
	                  const int a_d,
	                  const double a_dx,
	                  const double a_dy,
	                  const double a_dz);
}
