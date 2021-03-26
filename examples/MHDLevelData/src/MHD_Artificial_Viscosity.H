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
            const double a_gamma,
            Reduction<double>& a_Rxn,
            bool a_computeMaxWaveSpeed,
            bool a_callBCs);
}