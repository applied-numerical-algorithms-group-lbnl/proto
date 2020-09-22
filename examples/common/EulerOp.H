#pragma once

#include "Proto.H"

constexpr unsigned int NUMCOMPS=DIM+2;
#define NGHOST 4

using namespace Proto;

typedef Var<double,NUMCOMPS> State;

namespace EulerOp
{
  /// Compute the Right-Hand side, conforming to the RK4 interface.
  void step(BoxData<double,NUMCOMPS>& a_Rhs,
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

