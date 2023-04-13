#pragma once

#include "Proto.H"
#include "MHDOp.H"

using namespace Proto;

namespace MHD_CFL
{
	void Min_dt_calc_func(double& a_dt,
	                    const BoxData<double,NUMCOMPS>& a_W_ave,
                        const double a_dx,
                        const double a_dy,
                        const double a_dz,
	                    const double a_gamma);						 					 		 
}
