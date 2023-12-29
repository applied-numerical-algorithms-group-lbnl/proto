#pragma once
#ifndef _BOX_OP_EULER_
#define _BOX_OP_EULER_

#include "Proto.H"
#include "Proto_CubedSphereShell.H"

#define NUMCOMPS DIM+2

using namespace Proto;
typedef BoxData<double> Scalar;
typedef BoxData<double, NUMCOMPS> State;
typedef BoxData<double, DIM> Vector;

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
template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_consToPrim_(
        Var<T, NUMCOMPS, MEM>&          a_W, 
        const Var<T, NUMCOMPS, MEM>&    a_U,
        double                          a_gamma)
{
    double rho = a_U(0);
    double v2 = 0.0;
    a_W(0) = rho;

    for (int i = 1; i <= DIM; i++)
    {
        double v;
        v = a_U(i) / rho;

        a_W(i) = v;
        v2 += v*v;
    }

    a_W(NUMCOMPS-1) = (a_U(NUMCOMPS-1) - .5 * rho * v2) * (a_gamma - 1.0);
    
}
PROTO_KERNEL_END(f_consToPrim_, f_consToPrim);

template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_geomA_(
                     Point               a_pt,
                     Var<T,DIM,MEM,DIM>& a_A,
                     Var<T>&             a_detA,
                     Array<uint,DIM>&    a_permute,
                     Array<int,DIM>&     a_sign,                     
                     T                   a_dx,
                     int                 a_dir)
{
  
  // Compute geometric quantities at face centers.
  Array<T,DIM> x;

  for (int dir = 1; dir < DIM; dir++)
    {
      x[dir] = a_pt[a_permute[dir]]*a_dx - .5*M_PI;
    }
  T alpha,beta;
  switch (a_dir)
    {
    case 0:
      alpha = x[1] + .5*a_dx;
      beta = x[2] + .5*a_dx;
      break;
    case 1:
      alpha = x[1] + .5*a_dx;
      beta = x[2];
      break;
    case 2:
      alpha = x[1];
      beta = x[2] + .5*a_dx;
      break;
    }
  T seca = 1./cos(alpha);
  T secb = 1./cos(beta);
  T tana = tan(alpha);
  T tanb = tan(beta);
  T norm = 1./sqrt(1. + tana*tana + tanb*tanb);
  T norm3 = norm*norm*norm;
  a_detA(0) = seca*seca*secb*secb*norm;
  
  a_A(a_permute[0],0) = norm;
  a_A(a_permute[1],0) = -seca*seca*tana*norm3;
  a_A(a_permute[2],0) = secb*secb*tanb*norm3;
      
  a_A(a_permute[0],1) = tana*norm;
  a_A(a_permute[1],1) = seca*seca*(tanb*tanb + 1)*norm3;
  a_A(a_permute[2],1) = -secb*secb*tana*tanb*norm3;
      
  a_A(a_permute[0],2) = tanb*norm;
  a_A(a_permute[1],2) = -seca*seca*tana*tanb*norm3;
  a_A(a_permute[2],2) = secb*secb*(tana*tana + 1)*norm3;

  for (int dir0 = 0; dir0 < DIM; dir0++)
    {
      for (int dir1 = 0; dir1 < DIM; dir1++)
        {
          a_A(a_permute[dir0],dir1) *= a_sign[dir0];
        }
    }
}
PROTO_KERNEL_END(f_geomA_, f_geomA)
template<typename T,MemType MEM>
void primToFlux(
                State& a_physFlux,
                const State& a_WFacePoint,
                T a_gamma,
                uint a_block,
                uint a_centering)
{
  // Compute fluxes at face centers.
  auto pointflux = forall<T,NUMCOMPS>
    ([ ] PROTO_LAMBDA
     (
      Var<T,NUMCOMPS>& a_flux,
      const Var<T,NUMCOMPS>& a_wface,
      const T a_gamma,
      uint a_block,
      uint a_centering)
     {
       T velnorm = a_wface(a_centering+1);
       a_flux(0) = a_wface(0)*velnorm;
       T veltan1 = a_wface((a_centering+1)%DIM + 1);
       T veltan2 = a_wface((a_centering+2)%DIM + 1);
       a_flux(a_centering+1) = a_flux(0)*velnorm;
       a_flux((a_centering+1)%DIM + 1) = a_flux(0)*veltan1;
       a_flux((a_centering+2)%DIM + 1) = a_flux(0)*veltan2;
       T ke = a_wface(0)*(velnorm*velnorm + veltan1*veltan1 + veltan2*veltan2)/2;
       a_flux(NUMCOMPS - 1) =
       velnorm*(ke + a_wface(NUMCOMPS-1)*a_gamma/(a_gamma - 1.0));
     },a_WFacePoint,a_gamma,a_block,a_centering);

  // Compute Average fluxes from pointwise values.
  a_physFlux = (1.0*Shift(Point::Zeros()) +
                (1.0/24.0)*Stencil<T>::LaplacianFace(a_centering))(pointflux);
}
template<typename T,MemType MEM>
void consToPrim(
                State& a_W,
                State& a_WBar,
                const State& a_JU,
                const BoxData<T,1,MEM>& a_dVolr,
                T a_gamma,
                T a_h,
                uint a_block)
{
  // Get radial conserved variables.
  // Momentum.

  // Geometric data.
  Box bx = a_JU.box();
  BoxData<T,1,MEM> J(bx.grow(1));
  scaledDetA(J,a_dVolr,bx,a_h,a_block,0,true);
  BoxData<T,1,MEM> Jinv(bx);
  forallInPlace([] PROTO_LAMBDA
                (Var<T,1,MEM>& a_jinv,
                 const Var<T,1,MEM>& a_j) {a_jinv(0) = 1.0/a_j(0);},J,Jinv);
  
  BoxData<T,DIM,MEM,DIM> dVolrDetA_A(bx);
  scaledDetA_A(dVolrDetA_A,a_dVolr,bx,a_h,a_block,0);
  BoxData<T,DIM,MEM,DIM> invcoef(bx);
  forallInPlace(f_matinv3by3,invcoef,dVolrDetA_A);
   // DEBUG
  a_JU.copyTo(a_W);
  return;
  // Set up storage for radial variables.
  BoxData<T,NUMCOMPS,MEM> UrBar(bx);
  a_WBar.define(bx);
  Box bx4 = bx.grow(-1);
  BoxData<T,NUMCOMPS,MEM> Ur(bx4);
  a_W.define(bx);

  //Momentum, radial coordinates, second and fourth order.
  BoxData<T,DIM,MEM> UrMomBar = slice<T,NUMCOMPS,DIM,MEM>(UrBar,1);
  BoxData<T,DIM,MEM> UrMom = slice<T,NUMCOMPS,DIM,MEM>(Ur,1);
  BoxData<T,DIM,MEM> JUMom = slice<T,NUMCOMPS,DIM,MEM>(a_JU,1);
  //JUMom.copyTo(UrMom);
  //cellMatrixQuotientMatAndInv<T,DIM,MEM>(UrMom,UrMomBar,dVolrDetA_A,invcoef,JUMom);
  
  //Density, energy; second and fourth order.
  BoxData<T,1,MEM> rho = slice(Ur,0);
  BoxData<T,1,MEM> rhoE = slice(Ur,NUMCOMPS-1);
  BoxData<T,1,MEM> rhoBar = slice(UrBar,0);
  BoxData<T,1,MEM> rhoEBar = slice(UrBar,NUMCOMPS-1);
  BoxData<T,1,MEM> JRho = slice(a_JU,1);
  BoxData<T,1,MEM> JRhoE = slice(a_JU,NUMCOMPS-1);
   
  //cellMatrixQuotientMatAndInv(rho,rhoBar,J,Jinv,JRho);
  //cellMatrixQuotientMatAndInv(rhoE,rhoEBar,J,Jinv,JRhoE);
 
  // transform fron conserved variables to primitive variables.
  //Operator::deconvolve(Ur,UrBar);
  //forallInPlace(f_consToPrim,a_W,Ur,a_gamma);
  //forallInPlace(f_consToPrim,a_WBar,UrBar,a_gamma);
  //a_W = Operator::_convolve(a_W,a_WBar);
}
template<typename T,MemType MEM>
void primToCons(
                State& a_JU,
                const State& a_WPoint,
                const BoxData<T>& a_dVolr,
                T a_gamma,
                T a_h,
                uint a_block)
{
  // Compute fluxes at face centers.
  Box bx = a_JU.box().grow(1);
  
  // Compute Cartesian UPoint.
  BoxData<T,NUMCOMPS,MEM> UPoint =
    forall_p<T,NUMCOMPS,MEM>([ ] PROTO_LAMBDA
                (Point a_pt,
                 Var<T,NUMCOMPS,MEM>& a_upoint,
                 const Var<T,NUMCOMPS>& a_wpoint,
                 T a_gamma,
                 T a_h,
                 int a_block)
                {
                  vecType velSph;
                  vecType velCart;
                  //cout << a_pt << " , "<< a_wpoint(1) << " , " << a_wpoint(2) << " , " << a_wpoint(3) << endl;
                  velSph[0] = a_wpoint(1);       
                  velSph[1] = a_wpoint(2);
                  velSph[2] = a_wpoint(3);
                  matType amat,amat2;
                  T half = .5;
                  T one = 1.0;
                  f_cofAmatrix_(a_pt,amat,a_h,half,half,a_block);
                  f_matVecProd_<T,DIM,MEM,DIM>(velCart,amat,velSph,one);
                  T ke = (velCart[0]*velCart[0] + velCart[1]*velCart[1] + velCart[2]*velCart[2])/2;              
                  a_upoint(0) = a_wpoint(0);
                  a_upoint(1) = a_wpoint(0)*velSph[0];
                  a_upoint(2) = a_wpoint(0)*velSph[1];
                  a_upoint(3) = a_wpoint(0)*velSph[2];
                  a_upoint(4) = a_wpoint(4)/(a_gamma- 1.0) + a_wpoint(0)*ke;
                },bx,a_WPoint,a_gamma,a_h,a_block);
  BoxData<T,NUMCOMPS,MEM> U = (1.0*Shift(Point::Zeros())
                               + (1.0/24.0)*Stencil<T>::Laplacian())(UPoint);
  // U (spherical) -> JU (cartesian).
  // Momentum.
  //a_JU.define(U.box().grow(-1));
  BoxData<T,DIM,MEM> JUMom = slice<T,NUMCOMPS,DIM,MEM>(a_JU,1);
  BoxData<T,DIM,MEM> UMom = slice<T,NUMCOMPS,DIM,MEM>(U,1);
  BoxData<T,DIM,MEM,DIM> dVolrDetA_A(bx);
  scaledDetA_A(dVolrDetA_A,a_dVolr,bx,a_h,a_block,0,true);
#if 0
  unsigned int col = 0;
  forallInPlace_p
    ([ ] PROTO_LAMBDA
     (Point a_pt,
      Var<T,DIM,MEM,DIM>& a_mat,
      Var<T,DIM,MEM>& a_vec,
      unsigned int a_col)
     {
       Array<T,DIM> position;
       for (int ind = 0; ind < DIM; ind++)
         {
           a_vec(ind) = a_mat(ind,a_col);
         }
     },dVolrDetA_A,JUMom,col);
#endif
  
  Operator::matrixProductATB<T,3,3,3,1,1,HOST>(JUMom,dVolrDetA_A,UMom);
  // Density, energy.
  BoxData<T> J(U.box());
  scaledDetA(J,a_dVolr,bx,a_h,a_block,0,true);
  //DEBUG
  //J.setVal(1.0*a_block); 
  BoxData<T,1,MEM> Jrho = slice(a_JU,0);
  BoxData<T,1,MEM> rho = slice(U,0);
  Operator::cellProduct(Jrho,J,rho);
  
  //J.copyTo(Jrho);
  BoxData<T,1,MEM> JrhoE = slice(a_JU,NUMCOMPS-1);
  BoxData<T,1,MEM> rhoE = slice(U,NUMCOMPS-1);
  Operator::cellProduct(JrhoE,J,rhoE);
  
}
template<typename T, MemType MEM = MEMTYPE_DEFAULT>
class BoxOp_EulerCubedSphere : public BoxOp<T, NUMCOMPS, 1, MEM>
{
    public:
    using BoxOp<T,NUMCOMPS,1,MEM>::BoxOp;

    T m_gamma = 5.0/3.0;
  Array<Array<uint,DIM>,6> m_permute = {{2,1,0},{2,1,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
  Array<Array<int,DIM>,6> m_sign = {{-1,1,1},{1,1,-1},{-1,1,1},{1,1,1},{1,-1,1},{-1,-1,1}};
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
   T                               a_dx,
   int                             a_block,
   T                               a_scale = 1.0) const
  {       
    PR_TIME("BoxOp_Euler::operator()");        
    // Primitive Variables.

    // Transform to primitive variables, spherical velocities.
    Box bx = a_JU.box();
    State W;
    State WBar;
    consToPrim(W,WBar,a_JU,a_dVolr,m_gamma,a_dx,a_block);
    a_JU.copyTo(a_Rhs);
    // Compute Fluxes, rhs.
    // a_Rhs.setVal(0.0);
    // Loop over directions.
    
    for (int dir = 0; dir < 1; dir++)
      {
#if 0
        PR_TIME("BoxOp_Euler::operator()::directional flux calculation");        
        State WFace = Stencil<T>::CellToFace(dir)(W);
        State WFaceBar = Stencil<T>::CellToFace(dir,Side::Lo,2)(WBar);
        // Fluxes on dir face.

        // Physical fluxes.
        // Deconvolve W on the face.        
        State WFacePoint=
        Operator::_deconvolveFace(WFace,WFaceBar,dir);
     
        // Average of advective fluxes in spherical variables at the
        // face center, given point values of W.
        BoxData<T,NUMCOMPS> physFlux;
        primToFlux<T,MEM>(physFlux,WFacePoint,m_gamma,a_block,dir);
        State fluxdir(a_Rhs.box().extrude(Point::Ones() - Point::Zeros()));

        // Advective fluxes for Cartesian variables.
        // Scalar fluxes.
        // Geometric quantities.
        BoxData<T,1,MEM> adjDdetA(fluxdir.box().grow(Point::Ones() - Point::Basis(dir)));
        BoxData<T,1,MEM> adjDrDir = slice(a_adjDr,dir);
        scaledDetA(adjDdetA,adjDrDir,bx,a_dx,a_block,dir);
        BoxData<T,1,MEM> fluxRho = slice(fluxdir,0);
        BoxData<T,1,MEM> fluxRhoE = slice(fluxdir,NUMCOMPS-1);
        Operator::faceProduct(fluxRho,slice(physFlux,0),adjDdetA,dir);
        Operator::faceProduct(fluxRhoE,slice(physFlux,NUMCOMPS-1),adjDdetA,dir);
        
        // Advective fluxes for momentum.
        BoxData<T,DIM,MEM,DIM> adjDdetA_A(fluxdir.box().grow(Point::Ones() - Point::Basis(dir)));         
        scaledDetA_A(adjDdetA_A,adjDrDir,bx,a_dx,a_block,dir);
        BoxData<T,DIM,MEM> fluxdirMom = slice<T,NUMCOMPS,DIM,MEM>(fluxdir,1);
        BoxData<T,DIM,MEM> physFluxMom = slice<T,NUMCOMPS,DIM,MEM>(physFlux,1);
        BoxData<T,DIM,MEM> fluxMom;
        fluxMom = Operator::_faceMatrixProductATB
          (adjDdetA_A,physFluxMom,adjDdetA_A,physFluxMom,dir);
        fluxMom.copyTo(fluxdirMom);
        // Momentum fluxes - pressure forces.
        BoxData<T,1,MEM> pslice = slice(WFace,NUMCOMPS-1);
        Box bx = W.box();
        Box bxdir = bx.extrude(Point::Basis(dir));
        BoxData<T,DIM,MEM> NTrow(bxdir);
        NTRowCS<T,MEM>(NTrow,adjDrDir,fluxdir.box(),a_dx,a_block,dir);
        BoxData<T,DIM,MEM> tmp =
          Operator::_faceMatrixProductAB(NTrow,pslice,NTrow,pslice,dir);
        fluxMom += tmp;

        // Save fluxes, update RHS.
        fluxdir.copyTo(a_fluxes[dir]);
        a_Rhs += Stencil<T>::FluxDivergence(dir)(a_fluxes[dir]);
#endif
      }
    a_Rhs *= (a_scale); 
  };
#ifdef PR_AMR
static inline void generateTags(
                                TagData& a_tags,
                                BoxData<T, NUMCOMPS>& a_state)
{
  forallInPlace(f_threshold, a_tags, a_state);
}
#endif
private:
 
};
template<typename T, MemType MEM>
inline void radialMetrics(
                          BoxData<T, 1, MEM>& a_radius,
                          BoxData<T,DIM,MEM>& a_Diagr,
                          BoxData<T,DIM,MEM>& a_adjDr,
                          BoxData<T, 1, MEM>& a_dVolr,
                          Box a_bx,
                          int a_nradius)
{
  int r_dir = CUBED_SPHERE_SHELL_RADIAL_COORD;
  double r0 = CUBED_SPHERE_SHELL_R0;
  double r1 = CUBED_SPHERE_SHELL_R1;
  double dr = (r1-r0)/a_nradius;
  forallInPlace_p([]PROTO_LAMBDA(Point a_pt,
                                 Var<double>& a_rad,
                                 Var<double,DIM>& a_diagr,
                                 Var<double,DIM>& a_adjdiagr,
                                 Var<double>& a_dVolr,
                                 double a_r0,
                                 double a_r1,
                                 double a_dr,
                                 int a_rdir)
                  {
                    a_rad(0) = a_r0 + a_dr*a_pt[a_rdir];
                    T rlow = a_rad(0);
                    T rhigh = a_rad(0)+a_dr;
                    a_rad(0) = rlow;
                    a_dVolr(0) = (rhigh*rhigh*rhigh - rlow*rlow*rlow)/3.0;
                    a_diagr(0) = (rhigh - rlow)/a_dr;
                    a_diagr(1) = (rhigh + rlow)/2.0;
                    a_diagr(2) = (rhigh + rlow)/2.0;
                    a_adjdiagr(0) = rlow*rlow;
                    a_adjdiagr(1) = .5*(rhigh*rhigh - rlow*rlow)/a_dr;
                    a_adjdiagr(2) = .5*(rhigh*rhigh - rlow*rlow)/a_dr;
                  },a_bx,a_radius,a_Diagr,a_adjDr,a_dVolr,r0,r1,dr,r_dir);
}
#endif //end include guard
