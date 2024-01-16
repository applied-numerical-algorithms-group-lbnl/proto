#pragma once
#ifndef _BOX_OP_EULER_
#define _BOX_OP_EULER_

#include "Proto.H"
#include "Proto_CubedSphereShell.H"
#include "Proto_CubedSphereKernels.H"

#define NUMCOMPS DIM+2

using namespace Proto;
typedef BoxData<double> Scalar;
typedef BoxData<double, NUMCOMPS> State;
typedef BoxData<double, DIM> Vector;
#define blocktest -1

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

// Cmpute kinetic energy / mass.
template<typename T,MemType MEM>
BoxData<T,1,MEM> kineticEnergyW(
                                const State& a_WPoint,                                
                                T a_h,
                                uint a_block)
{
  PR_TIMERS("kineticEnergyW");
  Array<Array<uint,DIM>,6>
     permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
  Array<Array<int,DIM>,6>
    sign = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};  
  // Compute fluxes at face centers.
  Box bx = a_WPoint.box();

  // Compute Cartesian UPoint.
  T half = .5;
  BoxData<T,DIM,MEM,DIM> Amat =
    forall_p<T,DIM,MEM,DIM>(f_Amatrix,bx,permute[a_block],sign[a_block],a_h,half,half);
  BoxData<T,1,MEM> ke =
    forall<T,1,MEM>([ ] PROTO_LAMBDA
                    (                
                     Var<T,1,MEM> a_ke,
                     const Var<T,DIM,MEM,DIM>& a_amat,
                     const Var<T,NUMCOMPS>& a_wpoint)
                    {
                      vecType velSph;
                      vecType velCart;
                  
                      velSph[0] = a_wpoint(1);       
                      velSph[1] = a_wpoint(2);
                      velSph[2] = a_wpoint(3);                  
                      T half = .5;
                      T one = 1.0;
                      a_ke(0) = 0.;
                      for (int dir0 = 0; dir0 < DIM; dir0++)
                        {
                          velCart[dir0] = 0.;
                          for (int dir1 = 0; dir1 < DIM; dir1++)
                            {
                              velCart[dir0] += a_amat(dir0,dir1)*velSph[dir1];
                            }
                          a_ke(0) += velCart[dir0]*velCart[dir0]/2;
                        }
                      a_ke(0) *= a_wpoint(0);
                    },Amat,a_WPoint);
  return ke;
}
template<typename T,MemType MEM>
void primToFlux(
                State& a_physFlux,
                const State& a_WFace,
                const State& a_WFaceBar,
                const BoxData<T,1,MEM>& a_adjDdir,
                T a_h,
                T a_gamma,
                uint a_block,
                uint a_centering)
{
  // Compute fluxes at face centers.
  PR_TIMERS("primToFlux");
  Stencil<T> deconvolveStencilRadial = (1.0*Shift(Point::Zeros())
                                    - (1.0/24.0)*Stencil<T>::Derivative(2,0,2));
  // Deconvolve W to obtain point values.
  BoxData<T,NUMCOMPS,MEM> WFacePoint = _deconvolveFace(a_WFace,a_WFaceBar,a_centering);
  BoxData<T,1,MEM> kePoint = kineticEnergyW<T,MEM>(WFacePoint,a_h,a_block);
  BoxData<T,1,MEM> adjDdirPoint;
  if (a_centering != 0)
    {
      adjDdirPoint = deconvolveStencilRadial(a_adjDdir);
    }
  else
    {
      adjDdirPoint.define(a_adjDdir.box());
      a_adjDdir.copyTo(adjDdirPoint);
    }
  State fluxPoint = forall<T,NUMCOMPS>
    ([ ] PROTO_LAMBDA
     (
      Var<T,NUMCOMPS>& a_flux,
      const Var<T,NUMCOMPS>& a_wface,      
      const Var<T,1>& a_ke,
      const Var<T>& a_adjddir,
      const T a_gamma,
      uint a_block,
      uint a_centering)
     {
       T velnorm = a_wface(a_centering+1);
       a_flux(0) = a_wface(0)*velnorm*a_adjddir(0);
       T veltan1 = a_wface((a_centering+1)%DIM + 1);
       T veltan2 = a_wface((a_centering+2)%DIM + 1);
       a_flux(a_centering+1) = a_flux(0)*velnorm;
       a_flux((a_centering+1)%DIM + 1) = a_flux(0)*veltan1;
       a_flux((a_centering+2)%DIM + 1) = a_flux(0)*veltan2;
       a_flux(NUMCOMPS - 1) =
       (a_ke(0) + a_wface(NUMCOMPS-1)*a_gamma/(a_gamma - 1.0))*a_adjddir(0)*velnorm;
     },WFacePoint,kePoint,adjDdirPoint,a_h,a_block,a_centering);

  // Compute Average fluxes from pointwise values.
  a_physFlux = (1.0*Shift(Point::Zeros()) +
                (1.0/24.0)*Stencil<T>::LaplacianFace(a_centering))(fluxPoint);
}
template<typename T,MemType MEM>
BoxData<T,NUMCOMPS,MEM> checkRadial(State& a_inputs,
                                    T a_h,
                                    unsigned int a_block)
{
  BoxData<T,NUMCOMPS,MEM> outputs(a_inputs.box());
  outputs.setToZero();
  Box bx = a_inputs.box();
  BoxData<T,1,MEM> deta = detA<T,MEM>(bx.grow(2),a_h,0); 
  BoxData<T,1,MEM> detaInv = forall<T,1,MEM>
    ([ ] PROTO_LAMBDA
     (Var<T,1,MEM>& a_output, const Var<T,1,MEM>& a_deta)
     {
       a_output(0) = 1.0/a_deta(0);
     }
     ,deta);
  BoxData<T,1,MEM> rho = slice(outputs,0);
  BoxData<T,1,MEM> rhoIn = slice(a_inputs,0);
  rho += rhoIn*detaInv;
  BoxData<T,1,MEM> rhoE = slice(outputs,NUMCOMPS-1);
  BoxData<T,1,MEM> rhoEIn = slice(a_inputs,NUMCOMPS-1);
  rhoE += rhoEIn*detaInv;
  
  BoxData<T,DIM,MEM,DIM> deta_amat = detA_A<T,MEM>(bx.grow(2),a_h,a_block,0); 
  BoxData<T,DIM,MEM,DIM> invcoef(bx);
  invcoef.setToZero();
  forallInPlace(f_matinv3by3,invcoef,deta_amat);
  BoxData<T,DIM,MEM> Mom = slice<T,NUMCOMPS,DIM,MEM>(outputs,1);
  BoxData<T,DIM,MEM> MomIn = slice<T,NUMCOMPS,DIM,MEM>(a_inputs,1);
  BoxData<T,DIM,MEM> tmpmom(outputs.box());
  T one = 1.0;  
  forallInPlace(f_matVecProd,tmpmom,invcoef,MomIn,one);
  tmpmom.copyTo(Mom);
  return outputs;
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
  PR_TIMERS("ConsToPrim");
  a_W.define(a_JU.box().grow(-1));
  a_WBar.define(a_JU.box());
  a_W.setToZero();
  a_WBar.setToZero();
  
  Stencil<T> convolveStencil = (1.0*Shift(Point::Zeros())
                                + (1.0/24.0)*Stencil<T>::Laplacian());
  Stencil<T> deconvolveStencil = (1.0*Shift(Point::Zeros())
                                  - (1.0/24.0)*Stencil<T>::Laplacian());
  Stencil<T> convolveStencilRadial = (1.0*Shift(Point::Zeros())
                                      + (1.0/24.0)*Stencil<T>::Derivative(2,0,2));
  Stencil<T> deconvolveStencilRadial = (1.0*Shift(Point::Zeros())
                                        - (1.0/24.0)*Stencil<T>::Derivative(2,0,2));
  Box bx = a_JU.box();
  BoxData<T,1,MEM> deta = detA<T,MEM>(bx.grow(2),a_h,0); 
  BoxData<T,1,MEM> detaInv = forall<T,1,MEM>
    ([ ] PROTO_LAMBDA
     (Var<T,1,MEM>& a_output, const Var<T,1,MEM>& a_deta)
     {
       a_output(0) = 1.0/a_deta(0);
     }
     ,deta);
  
  BoxData<T,DIM,MEM,DIM> deta_amat = detA_A<T,MEM>(bx.grow(2),a_h,a_block,0); 
  BoxData<T,DIM,MEM,DIM> invcoef(bx);
  invcoef.setToZero();
  forallInPlace(f_matinv3by3,invcoef,deta_amat);
 
  // Set up storage for radial variables.
  a_W.define(bx);
  Box bx4 = bx.grow(-1);
  BoxData<T,NUMCOMPS,MEM> U(a_JU.box());
  // JU -> V_r*(rho,rho*vel_cart,rho*E) = JU / det(A);
  {
    PR_TIMERS("ConsToRadialCons");
    BoxData<T,1,MEM> tmp(a_JU.box()); 
    BoxData<T,1,MEM> JRho = slice(a_JU,0);
    BoxData<T,1,MEM> rho = slice(U,0);
    tmp = JRho*detaInv;
    tmp.copyTo(rho);
    BoxData<T,1,MEM> JRhoE = slice(a_JU,NUMCOMPS-1);
    BoxData<T,1,MEM> rhoE = slice(U,NUMCOMPS-1);
    tmp = JRhoE*detaInv;
    tmp.copyTo(rhoE);
    // V_r*rho*v_cart -> V_r*rho*v_sphere = (detA A)^{-1} V_r*rho*v_cart;
    BoxData<T,DIM,MEM> JUMom = slice<T,NUMCOMPS,DIM,MEM>(a_JU,1);
    BoxData<T,DIM,MEM> UMom = slice<T,NUMCOMPS,DIM,MEM>(U,1);
    BoxData<T,DIM,MEM> tmpmom(JUMom.box());
    T one = 1.0;  
    forallInPlace(f_matVecProd,tmpmom,invcoef,JUMom,one);
    tmpmom.copyTo(UMom);
    // Deconvolve-pointvise-convolve to go from volr rho Ur to W.
    // divide by dVolr and deconvolve.
    BoxData<T> tmp2(U.box());
    for (int comp = 0; comp < NUMCOMPS;comp++)
      {      
        BoxData<T> Ucomp = slice(U,comp);
        tmp2 = cellQuotient(Ucomp,a_dVolr);
        tmp2.copyTo(Ucomp);
    }
  }
  {
    PR_TIMERS("Cons Spherical to W, Wbar");
    BoxData<T,NUMCOMPS,MEM> UPoint = deconvolveStencil(U);    
    // Pointwise UrPoint -> WPoint.
    BoxData<T,NUMCOMPS,MEM> WPoint =
      forall<T,NUMCOMPS,MEM>([ ] PROTO_LAMBDA
                             (Var<T,NUMCOMPS,MEM>& a_wpoint,
                              Var<T,NUMCOMPS,MEM>& a_upoint)                           
                             {
                               Var<T,NUMCOMPS,MEM> wpoint;
                               a_wpoint(0) = a_upoint(0);
                               a_wpoint(1) = a_upoint(1)/a_upoint(0);
                               a_wpoint(2) = a_upoint(2)/a_upoint(0);
                               a_wpoint(3) = a_upoint(3)/a_upoint(0);
                               a_wpoint(4) = a_upoint(4);                             
                             },UPoint);

  
    // Compute pressure.
    BoxData<T,1,MEM> kePoint = kineticEnergyW<T,MEM>(WPoint,a_h,a_block); 
    forallInPlace([ ] PROTO_LAMBDA
                  (Var<T,NUMCOMPS>& a_wpoint,
                   const Var<T,1,MEM>& a_ke,
                   T a_gamma)
                  {
                    a_wpoint(NUMCOMPS-1) = (a_wpoint(NUMCOMPS-1) - a_ke(0))*(a_gamma - 1.0);
                  },WPoint,kePoint,a_gamma);
    if (a_block == blocktest)
          {
            HDF5Handler h5;
            BoxData<T,NUMCOMPS,MEM> tmpPlot(WPoint.box().grow(-2));
            WPoint.copyTo(tmpPlot);
            h5.writePatch(a_h,tmpPlot,"ConsToPrim::WPoint");
          }
    a_W = convolveStencil(WPoint);
    WPoint.copyTo(a_WBar);
  }
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
  PR_TIMERS("primToCons");
  Array<Array<uint,DIM>,6>
     permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
  Array<Array<int,DIM>,6>
    sign = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};
  Stencil<T> convolveStencil = (1.0*Shift(Point::Zeros())
                                + (1.0/24.0)*Stencil<T>::Laplacian());
  Stencil<T> deconvolveStencil = (1.0*Shift(Point::Zeros())
                                  - (1.0/24.0)*Stencil<T>::Laplacian());
  Stencil<T> convolveStencilRadial = (1.0*Shift(Point::Zeros())
                                + (1.0/24.0)*Stencil<T>::Derivative(2,0,2));
  Stencil<T> deconvolveStencilRadial = (1.0*Shift(Point::Zeros())
                                        - (1.0/24.0)*Stencil<T>::Derivative(2,0,2));
  Box bx = a_WPoint.box().grow(2);
  
  // Primitive variables -> V_r*(rho,vel_sph,rho*E).                                     
  BoxData<T,1,MEM> ke = kineticEnergyW<T,MEM>(a_WPoint,a_h,a_block);
  BoxData<T,NUMCOMPS,MEM> UPoint =
    forall<T,NUMCOMPS,MEM>([ ] PROTO_LAMBDA
                           (Var<T,NUMCOMPS,MEM>& a_upoint,
                            const Var<T,NUMCOMPS>& a_wpoint,
                            const Var<T,1,MEM>& a_ke,
                            T a_gamma)
                           {
                             a_upoint(0) = a_wpoint(0);
                             a_upoint(1) = a_wpoint(1)*a_upoint(0);
                             a_upoint(2) = a_wpoint(2)*a_upoint(0);
                             a_upoint(3) = a_wpoint(3)*a_upoint(0);
                             a_upoint(4) = a_wpoint(4)/(a_gamma - 1.0) + a_ke(0);
                           },a_WPoint,ke,a_gamma);
  
  BoxData<T,NUMCOMPS,MEM> U = convolveStencil(UPoint);

  // U -> U*dVolr.
  
  for (int comp = 0; comp < NUMCOMPS;comp++)
    {
      BoxData<T> tmp;
      BoxData<T> Ucomp = slice(U,comp);
      tmp = a_dVolr*Ucomp;
      tmp.copyTo(Ucomp);
    }

  // spherical velocities -> Cartesian velocities.
  a_JU.define(U.box().grow(-1));
  BoxData<T,DIM,MEM> JUMom = slice<T,NUMCOMPS,DIM,MEM>(a_JU,1);
  BoxData<T,DIM,MEM> UMom = slice<T,NUMCOMPS,DIM,MEM>(U,1);
  BoxData<T,DIM,MEM,DIM> amat = detA_A<T,MEM>(bx,a_h,a_block,0);
  
  T one = 1.0;
  forallInPlace(f_matVecProd,JUMom,amat,UMom,one);
  JUMom.copyTo(UMom);

  // scale by det(A)

  BoxData<T,1,MEM> deta = detA<T,MEM>(bx,a_h,0);
  {
    BoxData<T,1,MEM> JUcomp = slice(a_JU,0);
    BoxData<T,1,MEM> Ucomp = slice(U,0);
    BoxData<T,1,MEM> tmp;
    tmp = Ucomp * deta;
    tmp.copyTo(JUcomp);
  }
  {
    BoxData<T,1,MEM> JUcomp = slice(a_JU,NUMCOMPS-1);
    BoxData<T,1,MEM> Ucomp = slice(U,NUMCOMPS-1);
    BoxData<T,1,MEM> tmp;
    tmp = Ucomp * deta;
    tmp.copyTo(JUcomp);
  }
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
   T                               a_h,
   int                             a_block,
   T                               a_scale = 1.0) const
  {
    PR_TIME("BoxOp_Euler::operator()");        
    // Primitive Variables.

    // Transform to primitive variables, spherical velocities.
    Box bx = a_JU.box();
    State W;
    State WBar;
    if (a_block == blocktest)
          {
            HDF5Handler h5;           
            h5.writePatch(a_h,a_JU,"JU");
          }
    consToPrim(W,WBar,a_JU,a_dVolr,m_gamma,a_h,a_block);
    if (a_block == blocktest)
          {
            HDF5Handler h5;           
            h5.writePatch(a_h,W,"W");
            h5.writePatch(a_h,WBar,"WBar");
          }
    // Compute Fluxes, rhs.
    a_Rhs.setVal(0.0);
    // Loop over directions.
    
    for (int dir = 0; dir < DIM; dir++)
      {
        PR_TIME("BoxOp_Euler::operator()::directional flux calculation");        
        State WFace = Stencil<T>::CellToFace(dir)(W);
        State WFaceBar = Stencil<T>::CellToFace(dir,Side::Lo,2)(WBar);
        if (a_block == blocktest)
          {
            HDF5Handler h5;
            BoxData<T,NUMCOMPS,MEM> tmpPlot(a_Rhs.box().extrude(Point::Basis(dir)));
            WFace.copyTo(tmpPlot);
            h5.writePatch(a_h,tmpPlot,"WFace"+std::to_string(dir));
            WFaceBar.copyTo(tmpPlot);
            h5.writePatch(a_h,tmpPlot,"WFaceBar"+std::to_string(dir));
          }
        // Fluxes on dir face.

        // Physical fluxes.
     
        // Average of advective fluxes in spherical variables at the
        // face center, given point values of W.
        BoxData<T,NUMCOMPS> physFlux;
        BoxData<T,1,MEM> adjDrDir = slice(a_adjDr,dir);  
        primToFlux<T,MEM>(physFlux,WFace,WFaceBar,adjDrDir,a_h,m_gamma,a_block,dir);
        Box bxface = physFlux.box(); 
        State fluxdir(bxface);
        if (a_block == blocktest)
          {
            HDF5Handler h5;
            BoxData<T,NUMCOMPS,MEM> tmpPlot(a_Rhs.box().extrude(Point::Basis(dir)));
            physFlux.copyTo(tmpPlot);
            h5.writePatch(a_h,tmpPlot,"physFlux"+std::to_string(dir));
          }
        // Advective fluxes for Cartesian variables.
        // Scalar fluxes.
        // Geometric quantities.
        
        // Advective fluxes for Cartesian momentum.
        T one = 1.0;
        BoxData<T,DIM,MEM,DIM> deta_a = detA_A<T,MEM>(bxface,a_h,a_block,dir);
        BoxData<T,DIM,MEM> fluxdirMom = slice<T,NUMCOMPS,DIM,MEM>(fluxdir,1);
        BoxData<T,DIM,MEM> physFluxMom = slice<T,NUMCOMPS,DIM,MEM>(physFlux,1);
        forallInPlace(f_matVecProd,fluxdirMom,deta_a,physFluxMom,one);
        //BoxData<T,DIM,MEM> tmp =
        //  Operator::_faceMatrixProductAB(deta_a,physFluxMom,deta_a,physFluxMom,dir);
        //tmp.copyTo(fluxdirMom);
        // Scale all advective fluxes by det(A).
        BoxData<T,1,MEM> deta = detA<T,MEM>(bxface,a_h,dir);
        
        forallInPlace([ ] PROTO_LAMBDA
                      (Var<T,NUMCOMPS,MEM>& a_advflux,
                       const Var<T,NUMCOMPS,MEM>& a_fluxin,
                       const Var<T,1,MEM>& a_deta,
                       const Var<T,1,MEM>& a_adjddir)
                      {                    
                            a_advflux(0) = a_deta(0)*a_fluxin(0);
                            a_advflux(NUMCOMPS-1) = a_deta(0)*a_fluxin(NUMCOMPS-1);
                            for (int comp = 0; comp < NUMCOMPS;comp++)
                              {
                                a_advflux(comp) *= a_adjddir(0);
                              }
                      },fluxdir,physFlux,deta,adjDrDir);
        if (a_block == blocktest)
          {
            HDF5Handler h5;
            BoxData<T,NUMCOMPS,MEM> tmpPlot(a_Rhs.box().extrude(Point::Basis(dir)));
            fluxdir.copyTo(tmpPlot);
            h5.writePatch(a_h,tmpPlot,"advectiveFlux"+std::to_string(dir));
            h5.writePatch(a_h,adjDrDir,"AdjDrDir"+std::to_string(dir));
          }
        // Momentum fluxes - pressure forces.
        {
          BoxData<T,1,MEM> pslice = slice(WFace,NUMCOMPS-1);
          Box bx = W.box();
          Box bxdir = bx.extrude(Point::Basis(dir));
          BoxData<T,DIM,MEM> ntrow = adjA<T,MEM>(bx,a_h,a_block,dir);
          BoxData<T,DIM,MEM> tmp2(ntrow.box());
          BoxData<T,1,MEM> tmp3 = faceProduct(adjDrDir,pslice,dir);
        for (int dirvec = 0; dirvec < DIM; dirvec++)
          {
            BoxData<T> ntslice = slice(ntrow,dirvec);          
            ntslice *= tmp3;
          }
       
        fluxdirMom += ntrow;
         if (a_block == blocktest)
          {
            HDF5Handler h5;
            BoxData<T,NUMCOMPS,MEM> tmpPlot(a_Rhs.box().extrude(Point::Basis(dir)));
            fluxdir.copyTo(tmpPlot);
            h5.writePatch(a_h,tmpPlot,"totalFlux"+std::to_string(dir));
            if (dir < DIM)
              {
                BoxData<T,NUMCOMPS,MEM> tmpPlot2
                  = checkRadial<T,MEM>(tmpPlot,a_h,a_block);
                h5.writePatch(a_h,tmpPlot2,"radialTotalFlux"+std::to_string(dir));
              }           
          }
        // Save fluxes, update RHS.
         fluxdir.copyTo(a_fluxes[dir]);
         BoxData<T,NUMCOMPS,MEM> rhsIncr(a_Rhs.box());
         rhsIncr.setToZero();
         T hdirInv;
         hdirInv = 1.0*a_Rhs.box().size(dir);         
         rhsIncr += Stencil<T>::FluxDivergence(dir)(fluxdir,hdirInv);
         a_Rhs += rhsIncr;       
        }
      }
    if (a_block == blocktest)
      {
        HDF5Handler h5;
        BoxData<T,NUMCOMPS,MEM> tmpPlot
          = checkRadial<T,MEM>(a_Rhs,a_h,a_block);
        h5.writePatch(a_h,tmpPlot,"radialRHS");
        h5.writePatch(a_h,a_Rhs,"RHS");
      }                 
    a_Rhs *= (a_scale);
  };
  private:
  
};
#ifdef PR_AMR
static inline void generateTags(
                                TagData& a_tags,
                                BoxData<T, NUMCOMPS>& a_state)
{
  forallInPlace(f_threshold, a_tags, a_state);
}
#endif

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
                    a_dVolr(0) = (rhigh*rhigh*rhigh - rlow*rlow*rlow)/3.0/a_dr;
                    a_diagr(0) = (rhigh - rlow)/a_dr;
                    a_diagr(1) = (rhigh + rlow)/2.0;
                    a_diagr(2) = (rhigh + rlow)/2.0;
                    a_adjdiagr(0) = rlow*rlow;
                    a_adjdiagr(1) = .5*(rhigh*rhigh - rlow*rlow)/a_dr;
                    a_adjdiagr(2) = .5*(rhigh*rhigh - rlow*rlow)/a_dr;
                  },a_bx,a_radius,a_Diagr,a_adjDr,a_dVolr,r0,r1,dr,r_dir);
}
#endif //end include guard
