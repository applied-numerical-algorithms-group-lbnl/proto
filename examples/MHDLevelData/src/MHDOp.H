#pragma once

#include "Proto.H"

//TODO: These variables need to go into a common file or be defined in a more accessible way
constexpr unsigned int NUMCOMPS=2*DIM+2;
#define NGHOST 8
#define PI 3.141592653589793

using namespace Proto;

typedef Var<double,DIM> V;
typedef Var<double,NUMCOMPS> State;

namespace MHDOp
{

    void initializeState(BoxData<double,NUMCOMPS>& a_state,
                         const double a_dx,
                         const double a_gamma);

    // Compute the right-hand side for the MHD equations (i.e. div(flux)).
    // Note: a_Rhs is NOT multiplied by the time step, which may be required
    // by the RK4 implementation.
	void step(BoxData<double,NUMCOMPS>& a_Rhs,
            const BoxData<double,NUMCOMPS>& a_U,
            const Box& a_dbx0,
            const double a_dx,
            const double a_gamma,
            Reduction<double>& a_Rxn,
            bool a_computeMaxWaveSpeed=true,
            bool a_callBCs=true);
			
	void step2(BoxData<double,NUMCOMPS>& a_Rhs,
            const BoxData<double,NUMCOMPS>& a_U,
            const Box& a_dbx0,
            const double a_dx,
            const double a_gamma,
            Reduction<double>& a_Rxn,
            bool a_computeMaxWaveSpeed=true,
            bool a_callBCs=true);

    void step3(BoxData<double,NUMCOMPS>& a_Rhs,
            const BoxData<double,NUMCOMPS>& a_U,
            const Box& a_dbx0,
            const double a_dx,
            const double a_gamma,
            Reduction<double>& a_Rxn,
            bool a_computeMaxWaveSpeed=true,
            bool a_callBCs=true);			
  
    // If one is so motivated, one can expose internal pointwise functions.
    // It is not necessary for correctness, these are all internal functions
    // used by the "step" function above.  
	PROTO_KERNEL_START
	void consToPrimF(State&         a_W, 
             const State&   a_U,
             double         a_gamma);

	PROTO_KERNEL_START
	void upwindStateF(State& a_out,
                   const State& a_low,
                   const State& a_high,
                   int   a_dir,
                   double a_gamma);


}

