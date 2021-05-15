#pragma once

#include "Proto.H"
#include "MHDOp.H"

using namespace Proto;

namespace MHD_Initialize
{
	void initializeState(BoxData<double,NUMCOMPS>& a_state,
	                     const double a_dx,
	                     const double a_gamma);
}
