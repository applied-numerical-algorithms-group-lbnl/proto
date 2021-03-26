#pragma once

#include "Proto.H"
#include "MHDOp.H"
//TODO: These variables need to go into a common file or be defined in a more accessible way
//constexpr unsigned int NUMCOMPS=2*DIM+2;
typedef Var<double,NUMCOMPS> State;


namespace MHD_Limiters{
	void MHD_Limiters(BoxData<double,NUMCOMPS>& a_W_ave_low_lim_flat,
				  BoxData<double,NUMCOMPS>& a_W_ave_high_lim_flat,
				  const BoxData<double,NUMCOMPS>& a_W_ave_low,
				  const BoxData<double,NUMCOMPS>& a_W_ave_high,
				  const BoxData<double,NUMCOMPS>& a_W_ave,
				  const int a_d,
                  const double a_dx);
}