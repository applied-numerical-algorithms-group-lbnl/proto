#pragma once
#ifndef _BOX_OP_EULER_
#define _BOX_OP_EULER_

#include "Proto.H"
#include "Inputs_Parsing.H"

#define NUMCOMPS 2*DIM+2
#define NUMREG NUMCOMPS
#define NGHOST 7
#if DIM == 2
#define iRHO 0
#define iVX 1
#define iMOMX 1
#define iVY 2
#define iMOMY 2
#define iE 3
#define iP 3
#define iBX 4
#define iBY 5
#if TURB == 1
	#define iRHOZ2 6
	#define iZ2 6
	#define iRHOZ2SIGMA 7
	#define iSIGMA 7
	#define iRHOLAMBDA 8
	#define iLAMBDA 8
#endif
#endif

#if DIM == 3
#define iRHO 0
#define iVX 1
#define iMOMX 1
#define iVY 2
#define iMOMY 2
#define iVZ 3
#define iMOMZ 3
#define iE 4
#define iP 4
#define iBX 5
#define iBY 6
#define iBZ 7
#if TURB == 1
	#define iRHOZ2 8
	#define iZ2 8
	#define iRHOZ2SIGMA 9
	#define iSIGMA 9
	#define iRHOLAMBDA 10
	#define iLAMBDA 10
#endif
#endif

using namespace Proto;
typedef BoxData<double> Scalar;
typedef BoxData<double, NUMCOMPS> State;
typedef BoxData<double, DIM> Vector;
#define BLOCKTEST -1

//State: [rho, G0, G1, ..., E]
// Gi = rho*vi
// E = p/(gamma-1) + 0.5*rho*|v|^2
//template<typename T, MemType MEM>
template<typename T,MemType MEM>
PROTO_KERNEL_START
void f_sphToPrim_F(
                   Var<T,NUMCOMPS,MEM>& a_wpoint,
                   Var<T,NUMCOMPS,MEM>& a_upoint)                           
{
  Var<T,NUMCOMPS,MEM> wpoint;
  a_wpoint(iRHO) = a_upoint(iRHO);
  a_wpoint(iVX) =
    a_upoint(iMOMX)/a_upoint(iRHO);
  a_wpoint(iVY) =
    a_upoint(iMOMY)/a_upoint(iRHO);
  a_wpoint(iVZ) =
    a_upoint(iMOMZ)/a_upoint(iRHO);
  a_wpoint(iE) = a_upoint(iP);  
  a_wpoint(iBX) = a_upoint(iBX);
  a_wpoint(iBY) = a_upoint(iBY);
  a_wpoint(iBZ) = a_upoint(iBZ);
  
}
PROTO_KERNEL_END(f_sphToPrim_F,f_sphToPrim)
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

  T m_gamma = ParseInputs::get_gamma();
  mutable T umax;

  // How many ghost cells does the operator need from the state variables
  inline static Point ghost() { return Point::Ones(NGHOST);}
    
  // How many ghost cells does the operator need from the auxiliary variables
  inline static Point auxGhost() { return Point::Zeros();}
    
  // What is the intended order of accuracy of the operator
  inline static constexpr int order() { return 4; }

  // Initialization
  inline void initialize( State& a_W,
                        const State& a_BC,
                        const Scalar& a_radius,
                        const Vector& a_Xcart,
                        T a_gamma,
                        int a_nradius) const;

  // New factorization of the operator starts here.
  // New operator() - calls RadialMetrics from inside the operator.
  
   
  
  inline void operator()
  (
   State&                          a_Rhs,
   Array<State, DIM>&              a_fluxes,
   T&                              a_dt_next,
   const State&                    a_USph,
   T                               a_scale = 1.0) const;
  
  inline void primToFlux(
                         State& a_physFlux,
                         const State& a_WFace,
                         const State& a_WFaceBar,
                         const BoxData<T,1,MEM>& a_adjDdir,
                         uint a_block,
                         uint a_centering) const;
  // Implemented.
  inline void sphToPrim(
                         State& a_W,
                         State& a_WBar,
                         const State& a_USph,
                         uint a_block) const;

  // Implemented.
  inline void sphToPrim_2O(
                         State& a_W,
                         const State& a_USph,
                         uint a_block) const;
  
  // Implemented.
  inline void primToSph(
                         State& a_USph,
                         const State& a_WPoint,
                         const BoxData<T,1,MEM>& a_dVolr,
                         uint a_block) const;
  // Implemented.
  inline BoxData<T,NUMCOMPS,MEM> consToSph(State& a_input,
                                           Box a_bx,
                                          unsigned int a_block) const;

  /// Compute conserved quantities only on valid cells, given point values
  /// of primitive variables. Typically used for initialization. 
  inline void primToCons(
                         State& a_JU,
                         const State& a_WPoint,
                         const BoxData<T>& a_dVolr,
                         T a_gamma,
                         T a_h,
                         uint a_block) const;
  
  // Implemented.
  inline BoxData<T,NUMCOMPS,MEM> sphToCart(State& a_input,
                                           Box a_bx,
                                         unsigned int a_block) const;

  //Implemented.
  inline void cartToSph(
                        State& a_USph,
                        const State& a_JU,
                        const BoxData<T,1,MEM>& a_dVolr,
                        uint a_block) const;

  inline BoxData<T,1,MEM> kineticEnergyW(
                                         const State& a_WPoint,
                                         uint a_block,
                                         uint a_centering) const;

  inline BoxData<T,1,MEM> magneticEnergyW(
                                         const State& a_WPoint,
                                         uint a_block,
                                         uint a_centering) const;

  inline BoxData<T,1,MEM> vdotbW(
                                const State& a_WPoint,
                                uint a_block,
                                uint a_centering) const;                                       

  inline void radialMetrics(
                            BoxData<T, 1, MEM>& a_radius,
                            BoxData<T,DIM,MEM>& a_Diagr,
                            BoxData<T,DIM,MEM>& a_adjDr,
                            BoxData<T, 1, MEM>& a_dVolr,
                            Box a_bx) const;

  inline void primToMappedFlux(
                              State& a_fluxdir,
                              const State& a_physFlux,
                              const State& a_WFace,
                              const State& a_WFaceBar,
                              const BoxData<T,1,MEM>& a_adjDrDir,
                              const BoxData<T,DIM,MEM,DIM>& a_deta_a,
                              const BoxData<T,1,MEM>& a_deta,
                              const BoxData<T,DIM,MEM>& a_ntrow,
                              uint a_block,
                              uint a_dir) const; 

  inline void primToMappedRusanovFlux(
                              State& a_fluxdir,
                              const State& a_WFace_L,
                              const State& a_WFaceBar_L,
                              const State& a_WFace_R,
                              const State& a_WFaceBar_R,
                              const BoxData<T,1,MEM>& a_adjDdir,
                              const BoxData<T,DIM,MEM,DIM>& a_deta_a,
                              const BoxData<T,1,MEM>& a_deta,
                              const BoxData<T,DIM,MEM,DIM>& a_A,
                              const BoxData<T,DIM,MEM>& a_ntrow,
                              uint a_block,
                              uint a_dir) const;

  // End new factorization of operator().
  // Linear Viscocity starts here.

  inline void LinearVisc(
                        State& a_Rhs,
                        Array<BoxData<double, NUMCOMPS>, DIM>& a_fluxes,
                        const State& a_USph,
                        T a_scale = 1.0) const;
  

  inline void primToMappedViscFlux(
                                  State& a_fluxdir,
                                  const State& a_physFlux,
                                  const State& a_WFace,
                                  const State& a_WFaceBar,
                                  const BoxData<T,1,MEM>& a_adjDdir,
                                  const BoxData<T,DIM,MEM,DIM>& a_deta_a,
                                  const BoxData<T,1,MEM>& a_deta,
                                  const BoxData<T,DIM,MEM,DIM>& a_A,
                                  const BoxData<T,DIM,MEM>& a_ntrow,
                                  uint a_block,
                                  uint a_dir) const;

  // End Linear Viscosity.


  // DivB cleaning starts here.

  inline void DivBPowell(
                        State& a_Rhs,
                        const State& a_USph,
                        T a_scale = 1.0) const;

  // End DivB cleaning.

  // Begin Limiter.
  inline void Limiters_4O(BoxData<double,NUMCOMPS>& a_W_ave_L_lim_flat,
                              BoxData<double,NUMCOMPS>& a_W_ave_R_lim_flat,
                              BoxData<double,NUMCOMPS>& a_W_ave_L,
                              BoxData<double,NUMCOMPS>& a_W_ave_R,
                              BoxData<double,NUMCOMPS>& a_W_ave,
                              BoxData<double,NUMCOMPS>& a_W_bar,
                              const int a_d,
                              const double a_dx,
                              const double a_dy,
                              const double a_dz) const;

  inline void Limiters_minmod(BoxData<double,NUMCOMPS>& a_W_low_lim,
                  BoxData<double,NUMCOMPS>& a_W_high_lim,
                  BoxData<double,NUMCOMPS>& a_W,
                  BoxData<double,DIM>& a_dx_sph,
                  const int a_d) const;                   

  // End Limiter.
  
  /// Local boundary conditions here. 
  inline void
  PreStagePatch(BoxData<T,NUMCOMPS,MEM>& a_UStage,
                BoxData<T,NUMCOMPS,MEM>& a_BC,
                BoxData<T,1,MEM>& a_dVolr,
                Box a_blockBox,
                T a_time,
                T a_dt,
                int a_kstage);
  
  /// Exchange, block boundary conditions, transform from JU to USph. 
  // Implemented.
  static inline void
  preStageLevel(MBLevelBoxData<T,NUMCOMPS,MEM>& a_UStage,
                MBLevelBoxData<T,NUMCOMPS,MEM>& a_rhsTotal,               
                MBLevelBoxData<T,1,MEM>& a_dVolrLev,
                MBInterpOp& a_iop,
                MBBoundaryRegister<T,NUMREG,MEM>& a_blockreg,
                T a_time,
                T a_dt,
                T a_weight,
                int a_kstage);

  /// Compute fluxes, and RHS.
  static inline void
  StageLevel(MBLevelBoxData<T,NUMCOMPS,MEM>& a_rhs, 
             T& a_dt_next,
             MBLevelBoxData<T,NUMCOMPS,MEM>& a_UStage,
             MBBoundaryRegister<T,NUMREG,MEM>& a_blockreg,
             T a_time,
             T a_dt,
             T a_weight,
             int a_kstage);

  /// Fill boundary registers. Need to add boundary register argument.
  // Implemented.
  inline void
  PostStagePatch(
                 std::vector<MBBoundaryData<T,NUMREG,MEM> >& a_mbdata,
                 Array<BoxData<T,NUMCOMPS,MEM>,DIM>& a_fluxes, 
                 T a_time,
                 T a_dt,
                 T a_weight);
                 

  /// Update solution, impose flux matching.
  static inline void
  postStepLevel(MBLevelBoxData<T,NUMCOMPS,MEM>& a_rhsTotal,
                MBBoundaryRegister<T,NUMREG,MEM,PR_CELL>& a_blockReg);

  inline BoxData<T,NUMCOMPS,MEM> checkRadial(
                                             State& a_inputs,
                                             T a_h,
                                             unsigned int a_block) const;

  inline void radialGhost(State& a_USph,
                          BoxData<T,1,MEM>& a_radius,
                          T a_dxi0,
                          Box a_ghostedBox,
                          Box a_blockBox);
  
  // static inline void
  // preStageLevel(MBLevelBoxData<T,NUMCOMPS,HOST>& UStage,
  //               const MBInterpOp& a_iop,
  //               const MBLevelBoxData<T,NUMCOMPS,MEM>& a_dVolrLev);
  
  static inline void
  postStageLevel(MBLevelBoxData<T,NUMCOMPS,MEM>& a_rhsStage,
                 MBBoundaryRegister<T,NUMREG,MEM,PR_CELL>& a_blockReg,
                 const MBLevelBoxData<T,NUMCOMPS,MEM>& a_dVolrLev);
private:
  
};
#include "BoxOp_Utilities.H"
#include "BoxOp_Viscocity.H"
#include "BoxOp_DivBcleaning.H"
#include "BoxOp_Limiters.H"
#include "BoxOp_PreStage.H"
#include "BoxOp_EulerCubedSphereImplem.H"
#include "BoxOp_Initialize.H"

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
