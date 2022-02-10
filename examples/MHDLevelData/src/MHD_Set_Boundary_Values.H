#pragma once

#include "Proto.H"
#include "MHDOp.H"

namespace MHD_Set_Boundary_Values {

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
	                         const int a_highBCtype);

	void Set_Jacobian_Values(BoxData<double,1>& a_Jacobian_ave,
	                         const Box& a_rangeBox,
	                         const ProblemDomain& a_probDomain,
	                         const double a_dx,
	                         const double a_dy,
	                         const double a_dz,
	                         const double a_gamma,
	                         const int a_lowBCtype,
	                         const int a_highBCtype);

	void Set_Zaxis_Values(BoxData<double,NUMCOMPS>& a_JU,
 						  const ProblemDomain& a_probDomain,
 	                      const BoxData<double,NUMCOMPS>& a_JU_pole);

}
