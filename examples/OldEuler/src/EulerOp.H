#pragma once

#include "Proto.H"
#include "CommonTemplates.H"

//TODO: These variables need to go into a common file or be defined in a more accessible way
constexpr unsigned int NUMCOMPS=DIM+2;
#define NGHOST 4
#define PI 3.141592653589793

using namespace Proto;

typedef Var<double,DIM> V;
typedef Var<double,NUMCOMPS> State;

namespace EulerOp
{
// Initialize state based on the following equations:
// rho = rho0+0.01*rho0*sin(4*PI*x(0))
// p = p0*pow(rho/rho0,gamma)
// c0 = sqrt(gamma*p0/rho0)
// c = sqrt(gamma*p/rho)
// umag = 2*(c-c0)/(gamma-1)
// state(0)=rho
// state(1)=rho*umag
// state(NUMCOMPS-1)=p/(gamma-1)+0.5*rho*umag*umag
// The state in all other directions is set to 0.
// In the implementation, rho0=gamma, p0=1.
    void initializeState(BoxData<double,NUMCOMPS>& a_state,
                         const double a_dx,
                         const double a_gamma);

  /// Compute the right-hand side for the Euler equations (i.e. div(flux)).
    // Note: a_Rhs is NOT multiplied by the time step, which may be required
// by the RK4 implementation.
  void step(BoxData<double,NUMCOMPS>& a_Rhs,
            const BoxData<double,NUMCOMPS>& a_U,
            const Box& a_dbx0,
            const double a_dx,
            const double a_gamma,
            Reduction<double, Abs>& a_Rxn,
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

