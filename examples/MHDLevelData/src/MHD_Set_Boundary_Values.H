#pragma once

#include "Proto.H"
#include "MHDOp.H"

namespace MHD_Set_Boundary_Values {

	void Set_Boundary_Values(BoxData<double,NUMCOMPS>& a_JU,
	                         const Box& a_rangeBox,
	                         const ProblemDomain& a_probDomain,
	                         const double a_dx,
	                         const double a_gamma,
	                         const BoxData<double,1>& Jacobian_ave,
	                         const int a_lowBCtype,
	                         const int a_highBCtype);

	void Set_Jacobian_Values(BoxData<double,1>& a_Jacobian_ave,
	                         const Box& a_rangeBox,
	                         const ProblemDomain& a_probDomain,
	                         const double a_dx,
	                         const double a_gamma,
	                         const int a_lowBCtype,
	                         const int a_highBCtype);

}
