#pragma once

#include "Proto.H"
#include "MHDOp.H"
//TODO: These variables need to go into a common file or be defined in a more accessible way
//constexpr unsigned int NUMCOMPS=2*DIM+2;
typedef Var<double,NUMCOMPS> State;


namespace MHD_Mapping {
	void Jacobian_Ave_calc(BoxData<double,1>& a_Jacobian_ave,
	                       const double a_dx,
	                       const Box& a_dbx0);

	void JU_to_U_calc(BoxData<double,NUMCOMPS>& a_U_demapped,
	                  const BoxData<double,NUMCOMPS>& a_a_U,
	                  BoxData<double,1>& a_Jacobian_ave,
	                  const Box& a_dbx0);

	void JU_to_U_2ndOrdercalc(BoxData<double,NUMCOMPS>& a_U_demapped,
	                          const BoxData<double,NUMCOMPS>& a_a_U,
	                          const BoxData<double,1>& a_Jacobian_ave,
	                          const Box& a_dbx0);

	void JU_to_W_calc(BoxData<double,NUMCOMPS>& a_W,
	                  const BoxData<double,NUMCOMPS>& a_JU,
	                  const double a_dx,
	                  const double a_gamma);

	void JU_to_W_bar_calc(BoxData<double,NUMCOMPS>& a_W_bar,
	                      const BoxData<double,NUMCOMPS>& a_JU,
	                      const double a_dx,
	                      const double a_gamma);


	void N_ave_f_calc_func(BoxData<double,1>& a_N_ave_f,
	                       const int a_s,
	                       const int a_d,
	                       const double a_dx);
	void eta_calc(BoxData<double,DIM>& a_eta,
	              const Box& a_bx,
	              const double a_dx);

	void eta_to_x_calc(BoxData<double,DIM>& a_x,
	                   const BoxData<double,DIM>& a_eta);

	void phys_coords_calc(BoxData<double,DIM>& a_x,
	                      const Box& dbx1,
	                      const double a_dx);

	void out_data_calc(BoxData<double,NUMCOMPS+DIM>& a_out_data,
	                   const BoxData<double,DIM>& a_phys_coords,
	                   const BoxData<double,NUMCOMPS>& a_W);
}
