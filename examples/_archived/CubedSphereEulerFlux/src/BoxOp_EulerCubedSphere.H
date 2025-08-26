#pragma once
#ifndef _BOX_OP_EULER_
#define _BOX_OP_EULER_

#include "Proto.H"

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
};
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
PROTO_KERNEL_END(f_consToPrim_, f_consToPrim)

template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_geomFace_(
                 Point               a_pt,                               
                 Array<uint,DIM>&    a_permute,
                 Array<int,DIM>&     a_sign,
                 Var<T,DIM,MEM,DIM>& a_A,
                 Var<T,DIM,MEM,DIM>& a_adjA,
                 Var<T>&             a_detA,
                 T                   a_dx,
                 int                 a_dir)
{
  
  // Compute geometric quantities at face centers.
  Array<T,DIM> x;

  for (int dir = 1; dir < DIM; dir++)
    {
      x[dir] = a_pt[dir]*a_dx - .5*M_PI;
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
      
  a_adjA(a_permute[0],0) = a_detA(0);
  a_adjA(a_permute[1],0) = tana*a_detA(0);
  a_adjA(a_permute[2],0) = tanb*a_detA(0);
      
  a_adjA(a_permute[0],1) = seca*seca*tanb*norm;
  a_adjA(a_permute[1],1) = 0.;
  a_adjA(a_permute[2],1) = -seca*seca*norm;
      
  a_adjA(a_permute[0],2) = tana*secb*secb*norm;
  a_adjA(a_permute[1],2) = -secb*secb*norm;
  a_adjA(a_permute[2],2) = 0.;         

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
          a_adjA(a_permute[dir0],dir1) *= a_sign[dir0];
        }
    }
}
PROTO_KERNEL_END(f_geomFace_, f_geomFace)
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
template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_CSFlux_(
                 Var<T,NUMCOMPS>&    a_flux,
                 Var<T>&             a_adjD,
                 Var<T,DIM,MEM,DIM>& a_A,
                 Var<T,DIM,MEM,DIM>& a_adjA,
                 Var<T>&             a_detA,
                 Var<T,NUMCOMPS>&    a_Wface,
                 T                   a_gamma,
                 int                 a_dir)
{
      
  // Mass, energy flux.
  T rho = a_Wface(0);
  T p = a_Wface(NUMCOMPS-1);
  T velnorm = a_Wface(a_dir+1);
  T veltan1 = a_Wface((a_dir + 1)%DIM + 1);
  T veltan2 = a_Wface((a_dir + 2)%DIM + 1);
  a_flux(0) = a_detA(0)*a_adjD(0)*rho;
  T engPlusWork = a_gamma*p/(a_gamma-1.0)
                + rho*(velnorm*velnorm + veltan1*veltan1 + veltan2*veltan2)/2;
  a_flux(NUMCOMPS-1) = a_detA(0)*a_adjD(0)*engPlusWork;
  
  // Momentum flux.
  Array<T,DIM> vcart;
  for (int dim1 = 0; dim1 < DIM; dim1++)
    {
      vcart[dim1] = 0.;
      for (int dim2 = 0; dim2 < DIM; dim2++)
        {
          vcart[dim1] += a_A(dim1,dim2)*a_Wface(dim2+1);
        }
      a_flux(dim1+1) = a_detA(0)*a_adjD(0)*vcart[dim1]*rho;
      //Pressure forces.
      a_flux(dim1+1) += a_adjD(0)*a_adjA(dim1,a_dir)*p;
    }
}

PROTO_KERNEL_END(f_CSFlux_, f_CSFlux)
template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_consToPrimSphere_(
                         Var<T,NUMCOMPS>&          a_W,
                         const Var<T,DIM,MEM,DIM>& a_A,
                         const Var<T,NUMCOMPS>&    a_U,
                         const Array<uint,DIM>&    a_permute,
                         const Array<int,DIM>&     a_sign,
                         T                         a_gamma)
{
  // Primitive variables, Cartesian Velocities.    
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
    // Transform Cartesian velocities to spherical velocities: v_cyl = A^T.v_cart
    Array<T,DIM> vCyl;
    for (int dim1 = 0; dim1 < DIM; dim1++)
      {
        vCyl[dim1] = 0.;
        for (int dim2 = 0; dim2 < DIM; dim2++)
          {           
            vCyl[dim1] += a_A(dim2,dim1)*a_W(dim2+1);
          }
      }
    for (int dim = 0 ;dim < DIM ; dim++)
      {
        a_W(dim+1) = vCyl[dim];
      }
}
PROTO_KERNEL_END(f_consToPrimSphere_, f_consToPrimSphere)   
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
    }
   
    // Apply BCs by filling ghost cells in stage values. For Euler, this is done by calling
    // exchange. For the MHD code, it will be more complicated.
    // The interface is very provisional. We expect it to evolve as we d more real problems.
    // Apply Operator
  inline void operator()(
                         State&                          a_Rhs,
                         Array<State, DIM>&              a_fluxes,
                         const State&                    a_U,
                         BoxData<T>&                     a_radius,
                         BoxData<T>&                     a_rdot,
                         int                             a_block,
                         T                               a_dx,
                         T                               a_scale = 1.0) const
  {       
    PR_TIME("BoxOp_Euler::operator()");        
    // Primitive Variables.

    State U = Operator::deconvolve(a_U);
        
    // Transform to primitive variables, spherical velocities.
    
    Box bx = a_U.box();
    Point baseHigh = bx.high();
    Point baseLow = bx.low();
    baseHigh[0] = baseLow[0];
    Box basebx(baseHigh,baseLow);
    BoxData<T,DIM,MEM,DIM> A(bx);
    BoxData<T> detA(bx);
    forallInPlace_p(f_geomA,basebx,A,detA,m_permute[a_block],m_sign[a_block],a_dx,0);
    Operator::spreadSlice(A);
    Operator::spreadSlice(detA);
    
    State W = forall<T,NUMCOMPS>
      (f_consToPrimSphere,A,U,m_permute[a_block],m_sign[a_block],m_gamma);
    State W_bar = forall<T,NUMCOMPS>
      (f_consToPrimSphere,A,a_U,m_permute[a_block],m_sign[a_block],m_gamma);  
        
    // Transform to spherical velocities (update in place).

    // MAX WAVE SPEED
    //Box rangeBox = a_U.box().grow(-ghost());
    //Scalar uabs = forall<double>(f_waveSpeedBound, rangeBox, W, gamma);
    //umax = uabs.absMax();

    // Compute Fluxes, rhs.
    a_Rhs.setVal(0.0);
    // Loop over directions.
    
    for (int dir = 0; dir < DIM; dir++)
      {
        PR_TIME("BoxOp_Euler::operator()::directional flux calculation");        
        State WFace = Stencil<T>::CellToFace(dir)(W);
          
        // Fluxes on dir face.

        // Geometric Data.
          
        Box bx = WFace.box();
        Point baseHigh = bx.high();
        Point baseLow = bx.low();
        baseHigh[0] = baseLow[0];
        Box basebx(baseHigh,baseLow);
        BoxData<T,DIM,MEM,DIM> A(bx);
        BoxData<T,DIM,MEM,DIM> adjA(bx);
        BoxData<T> detA(bx);
        forallInPlace_p(f_geomFace,basebx,
                        m_permute[a_block],m_sign[a_block],A,adjA,detA,a_dx,dir);
        Operator::spreadSlice(A);
        Operator::spreadSlice(adjA);
        Operator::spreadSlice(detA);
        BoxData<T> faceRadius = Stencil<T>::faceToCell(dir)(a_radius);
        BoxData<T> adjD;
        if (dir == 0)
          {
            adjD = faceRadius*faceRadius;
          }
        else
          {
            BoxData<T> faceRdot = Stencil<T>::faceToCell(dir)(a_rdot);
            adjD = faceRadius*faceRdot;
          }
        
        // Physical fluxes.
        State& fluxdir = a_fluxes[dir];
        State flux = forall<T,NUMCOMPS>(f_CSFlux,adjD,A,adjA,detA,WFace,m_gamma,dir);
        fluxdir = Operator::convolveFace(flux,dir);    

        a_Rhs += Stencil<T>::FluxDivergence(dir)(a_fluxes[dir]);
      }
    a_Rhs *= (a_scale); 
  }
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
inline void radiusRDot(
                       BoxData<T, 1, MEM>& a_radius,
                       BoxData<T, 1, MEM>& a_rDot,
                       int a_nradius)
{
  int r_dir = CUBED_SPHERE_SHELL_RADIAL_COORD;
  double r0 = CUBED_SPHERE_SHELL_R0;
  double r1 = CUBED_SPHERE_SHELL_R1;
  double dr = (r1-r0)/a_nradius;
  forallInPlace_p([]PROTO_LAMBDA(Point a_pt,
                                 Var<double>& a_rad,
                                 Var<double>& a_rdot,
                                 double a_r0,
                                 double a_r1,
                                 double a_dr,
                                 int a_rdir)
                  {
                    a_rad(0) = a_r0 + a_dr*a_pt[a_rdir];
                    a_rdot(0) = a_r1 - a_r0;
                  },a_radius,a_rDot,r0,r1,dr,r_dir);
};
#endif //end include guard
