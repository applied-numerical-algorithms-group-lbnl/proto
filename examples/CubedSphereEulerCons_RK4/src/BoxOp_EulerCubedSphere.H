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
#define VEL00 1.E5
#define RHO00  1.67262192e-24/100.
#define P00 1.E12
#define B00 1.E6
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
  
  a_wpoint(iRHO) = max(a_upoint(iRHO),RHO00);
  a_wpoint(iVX) =
    a_upoint(iMOMX)/a_wpoint(iRHO);
  a_wpoint(iVY) =
    a_upoint(iMOMY)/a_wpoint(iRHO);
  a_wpoint(iVZ) =
    a_upoint(iMOMZ)/a_wpoint(iRHO);
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
template<typename T,MemType MEM>
PROTO_KERNEL_START
void f_speedmax_F
(
 Var<T,1,MEM>& a_speedmax,
 const Var<T,NUMCOMPS,MEM>& a_W,
 const Var<T,DIM,MEM,DIM>& a_amat,
 const Var<T,DIM,MEM>& a_adjA,
 const Var<T,1,MEM>& a_adjDrDir,
 const Var<T,1,MEM>& a_detA,
 const Var<T,1,MEM>& a_dVolr,
 const Var<T,DIM,MEM>& a_Dr,
 const T a_dxiPerp,
 const T a_dxiR,
 const int a_d,
 const T a_gamma)
{ 
  double dx_d, dt_d;
  if (a_d > 0){
    T adjA_mag = sqrt(a_adjA(0)*a_adjA(0) + a_adjA(1)*a_adjA(1) + a_adjA(2)*a_adjA(2));
    dx_d = a_dVolr(0) * a_detA(0) * a_dxiPerp / a_adjDrDir(0) / adjA_mag;
  } 
  if (a_d == 0){
    dx_d = a_Dr(0) * a_dxiR;
  }   
  vecType velSph;
  vecType velCart;
  vecType bSph;
  vecType bCart;

  velSph[0] = a_W(iVX);       
  velSph[1] = a_W(iVY);
  velSph[2] = a_W(iVZ);
  velCart[0] = 0.;       
  velCart[1] = 0.;
  velCart[2] = 0.;
      
  bSph[0] = a_W(iBX);       
  bSph[1] = a_W(iBY);
  bSph[2] = a_W(iBZ);
  bCart[0] = 0.;       
  bCart[1] = 0.;
  bCart[2] = 0.;
  T half = .5;
  T one = 1.0;
  T v2 = 0.;
  T v1 = 0.;
  T b2 = 0.;
  for (int dir0 = 0; dir0 < DIM; dir0++)
    {
      for (int dir1 = 0; dir1 < DIM; dir1++)
        {
          velCart[dir0] += a_amat(dir0,dir1)*velSph[dir1];
          bCart[dir0] += a_amat(dir0,dir1)*bSph[dir1];
        }
      v2 += velCart[dir0]*velCart[dir0];
      v1 += abs(velCart[dir0]);
      b2 += bCart[dir0]*bCart[dir0];
    }
  T rho = a_W(iRHO);
  T B_mag = sqrt(b2);
  T p = a_W(iP);
  if (p < 0.0) p = 0.0;
  T af = sqrt(a_gamma*p/rho + b2/4.0/M_PI/rho);  
  a_speedmax(0) += (af + sqrt(v2))/dx_d;
}
PROTO_KERNEL_END(f_speedmax_F, f_speedmax);
template<typename T,MemType MEM>
PROTO_KERNEL_START
void f_csq_F
(
 Var<T,1,MEM>& a_csq,
 const Var<T,NUMCOMPS,MEM>& a_W,
 const Var<T,DIM,MEM,DIM>& a_amat,
 const T a_gamma)
{ 
  vecType velSph;
  vecType velCart;
  vecType bSph;
  vecType bCart;

  velSph[0] = a_W(iVX);       
  velSph[1] = a_W(iVY);
  velSph[2] = a_W(iVZ);
  velCart[0] = 0.;       
  velCart[1] = 0.;
  velCart[2] = 0.;
      
  bSph[0] = a_W(iBX);       
  bSph[1] = a_W(iBY);
  bSph[2] = a_W(iBZ);
  bCart[0] = 0.;       
  bCart[1] = 0.;
  bCart[2] = 0.;
  T half = .5;
  T one = 1.0;
  T v2 = 0.;
  T v1 = 0.;
  T b2 = 0.;
  for (int dir0 = 0; dir0 < DIM; dir0++)
    {
      for (int dir1 = 0; dir1 < DIM; dir1++)
        {
          velCart[dir0] += a_amat(dir0,dir1)*velSph[dir1];
          bCart[dir0] += a_amat(dir0,dir1)*bSph[dir1];
        }
      v2 += velCart[dir0]*velCart[dir0];
      v1 += abs(velCart[dir0]);
      b2 += bCart[dir0]*bCart[dir0];
    }
  T rho = a_W(iRHO);
  T B_mag = sqrt(b2);
  T p = a_W(iP);
  if (p < 0.0) p = 0.0;
  a_csq(0) = a_gamma*p/rho + b2/4.0/M_PI/rho;  
}
PROTO_KERNEL_END(f_csq_F, f_csq);
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
                        int a_nradius,
                        T a_h,
                        uint a_block) const;

  // New factorization of the operator starts here.
  // New operator() - calls RadialMetrics from inside the operator.

  inline void operator()
  (
   State&                          a_Rhs,
   State&                          a_Rhs_by_J,
   Array<State, DIM>&              a_fluxes,
   const State&                    a_USph,
   T                               a_dt,
   T                               a_scale = 1.0,
   int                             a_counter = -1) const;
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

  inline void dtInv(Reduction<T,Operation::Max,MEM>& a_dtInv,const State& a_W) const;
   
   inline static T dtCFL(
                         MBLevelBoxData<T,NUMCOMPS,MEM>& a_JU,
                         MBInterpOp& a_iop,
                         MBLevelBoxData<T,1,MEM>& a_dVolr);
  
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
  // Viscosity starts here.

  inline void
  viscFlux(
           State& a_Rhs,
           Array<BoxData<double, NUMCOMPS>, DIM>& a_fluxes,
           const State& a_USph,
           T a_scale,
           T a_dt) const;
  inline void LinearVisc(
                        State& a_Rhs,
                        Array<BoxData<double, NUMCOMPS>, DIM>& a_fluxes,
                        const State& a_USph,
                        T a_scale = 1.0) const;

  inline void
  LinearViscCoef(
                 BoxData<T,1,MEM>& a_LNViscCoef,
                 const State& a_WFace,
                 const BoxData<T,1,MEM>& a_deta,
                 const BoxData<T,DIM,MEM,DIM>& a_A,
                 const BoxData<T,DIM,MEM>& a_ntrow,
                 uint a_block,
                 uint a_dir) const;
  inline void
  NonLinearViscCoef(
                    BoxData<T,1,MEM>& a_NLViscCoef,
                    const State& a_W,
                    const BoxData<T,DIM,MEM>& a_adjDr,
                    const BoxData<T,DIM,MEM,DIM>& a_A,
                    const BoxData<T,DIM,MEM>& a_ntrow,
                    const BoxData<T,1,MEM>& a_dVolr,
                    const BoxData<T,1,MEM>& a_deta,
                    Array<T,DIM> a_dxi,
                    uint a_block,
                    uint a_dir) const;

  inline void
  primToMappedViscFlux(
                       State& a_fluxdir,
                       const State& a_physFlux,
                       const State& a_WFace,
                       const BoxData<T,1,MEM>& a_adjDdir,
                       const BoxData<T,DIM,MEM,DIM>& a_deta_a,
                       const BoxData<T,1,MEM>& a_deta,
                       const BoxData<T,DIM,MEM,DIM>& a_A,
                       const BoxData<T,DIM,MEM>& a_ntrow,
                       uint a_block,
                       uint a_dir) const;

  // End Viscosity.

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

// Begin CME
  static inline void
  Insert_CME(MBLevelBoxData<T, NUMCOMPS, MEM> &a_JU,
    MBLevelBoxData<T, 1, MEM> &a_dVolrLev,
    MBInterpOp &a_iop,
    MBDisjointBoxLayout &a_layout,
    double a_time,
    double a_dt,
    double a_gamma);
// End CME

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
                Array<Reduction<T,Operation::Sum,MEM>,NUMCOMPS>& a_consSum,
                T a_time,
                T a_dt,
                T a_weight,
                int a_kstage);

  /// Compute fluxes, and RHS.
  static inline void
  StageLevel(MBLevelBoxData<T,NUMCOMPS,MEM>& a_rhs,
             MBLevelBoxData<T,NUMCOMPS,MEM>& a_rhs_by_J,
             MBLevelBoxData<T,NUMCOMPS,MEM>& a_UStage,
             MBBoundaryRegister<T,NUMREG,MEM>& a_blockreg,
             Array<Reduction<T,Operation::Sum,MEM>,NUMCOMPS>& a_consSum,
             T a_time,
             T a_dt,
             T a_weight,
             int a_kstage);

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
  
  /// Apply operator. Usede for convergence testing.
  static Array<T,NUMCOMPS>
  applyOp(MBLevelBoxData<T,NUMCOMPS,MEM>& a_JU,
          MBInterpOp& a_iop,
          MBLevelBoxData<T,1,MEM>& a_dVolrLev);
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
  int m_count = 0;
  
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
