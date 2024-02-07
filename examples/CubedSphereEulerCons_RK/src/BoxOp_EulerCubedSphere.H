#pragma once
#ifndef _BOX_OP_EULER_
#define _BOX_OP_EULER_

#include "Proto.H"
//#include "Proto_CubedSphereShell.H"
//#include "Proto_CubedSphereKernels.H"

#define NUMCOMPS DIM+2
#define NGHOST 5

using namespace Proto;
typedef BoxData<double> Scalar;
typedef BoxData<double, NUMCOMPS> State;
typedef BoxData<double, DIM> Vector;
#define BLOCKTEST -1

//State: [rho, G0, G1, ..., E]
// Gi = rho*vi
// E = p/(gamma-1) + 0.5*rho*|v|^2
//template<typename T, MemType MEM>
PROTO_KERNEL_START
void
f_thresholdF(
             Var<short>& a_tags,
             Var<double, NUMCOMPS>& a_U)
{
  double thresh = 1.001;
  if (a_U(0) > thresh) {a_tags(0) = 1;}
  else {a_tags(0) = 0;};
}
PROTO_KERNEL_END(f_thresholdF, f_threshold);


template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM = MEMTYPE_DEFAULT>
class BoxOp_EulerCubedSphere : public MBBoxOp<T, MAP, NUMCOMPS, 0, MEM>
{
public:
  using MBBoxOp<T, MAP, NUMCOMPS,0,MEM>::MBBoxOp;

  T m_gamma = 5.0/3.0;
  mutable T umax;

  // How many ghost cells does the operator need from the state variables
  inline static Point ghost() { return Point::Ones(4);}
    
  // How many ghost cells does the operator need from the auxiliary variables
  inline static Point auxGhost() { return Point::Zeros();}
    
  // What is the intended order of accuracy of the operator
  inline static constexpr int order() { return 4; }

  // Initialization
  inline void init()
  {
  };  
  
  // Apply BCs by filling ghost cells in stage values. For Euler, this is done by calling
  // exchange. For the MHD code, it will be more complicated.
  // The interface is very provisional. We expect it to evolve as we d more real problems.
  // Apply Operator
  inline void operator()
      (
       State&                          a_Rhs,
       Array<State, DIM>&              a_fluxes,
       const State&                    a_JU,
       const BoxData<T,DIM,MEM>&       a_Dr,
       const BoxData<T,DIM,MEM>&       a_adjDr,
       const BoxData<T,1,MEM>&         a_dVolr,
       T                               a_h,
       int                             a_block,
       T                               a_scale = 1.0) const;

  // Same as the above but conforms to Proto's API
  inline void operator()
      (
       State&                          a_Rhs,
       Array<State, DIM>&              a_fluxes,
       const State&                    a_JU,
       T                               a_scale = 1.0) const;
  
  inline void primToFlux(
                         State& a_physFlux,
                         const State& a_WFace,
                         const State& a_WFaceBar,
                         const BoxData<T,1,MEM>& a_adjDdir,
                         T a_h,
                         T a_gamma,
                         uint a_block,
                         uint a_centering) const;
  
  inline void consToPrim(
                         State& a_W,
                         State& a_WBar,
                         const State& a_JU,
                         const BoxData<T,1,MEM>& a_dVolr,
                         T a_gamma,
                         T a_h,
                         uint a_block) const;
  
  inline void primToCons(
                         State& a_JU,
                         const State& a_WPoint,
                         const BoxData<T>& a_dVolr,
                         T a_gamma,
                         T a_h,
                         uint a_block) const;
  
  inline BoxData<T,1,MEM> kineticEnergyW(
                                         const State& a_WPoint,                    
                                         T a_h,
                                         uint a_block,
                                         uint a_centering) const;

  inline BoxData<T,NUMCOMPS,MEM> checkRadial(
                                             State& a_inputs,
                                             T a_h,
                                             unsigned int a_block) const;
  inline void radialMetrics(
                            BoxData<T, 1, MEM>& a_radius,
                            BoxData<T,DIM,MEM>& a_Diagr,
                            BoxData<T,DIM,MEM>& a_adjDr,
                            BoxData<T, 1, MEM>& a_dVolr,
                            Box a_bx,
                            int a_nradius) const;
private:
  
};
#include "BoxOp_EulerCubedSphereImplem.H"

#ifdef PR_AMR
/* TODO: this doesn't belong here and it doesn't compile
static inline void generateTags(
                                TagData& a_tags,
                                BoxData<T, NUMCOMPS>& a_state)
{
  forallInPlace(f_threshold, a_tags, a_state);
}
*/
#endif


#endif //end include guard
