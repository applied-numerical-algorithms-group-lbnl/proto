#ifndef _MULTIGRID_H_
#define _MULTIGRID_H_
#include "Proto.H"
using namespace Proto;
#include "MakeLaplaceStencils.H"
template<typename T>
PROTO_KERNEL_START
void f_bc_F(
            Point a_pt,
            Var<T>& a_phi,
            const Var<T>& a_phishift,
            T a_dx,
            T a_dirichlet,
            Point a_ptshift)
{
  T dist = 0.;
  for (int dir = 0; dir < DIM; dir++)
    {
      T xloc = (a_pt[dir] + a_ptshift[dir]
                - std::copysign(.5,a_ptshift[dir]*1.0))*a_dx;
      dist += xloc*xloc;
    }
  dist = sqrt(dist);
  a_phi(0) = a_dirichlet*dist  + a_phishift(0);
  //a_phi(0) = a_phishift(0);
}
PROTO_KERNEL_END(f_bc_F,f_bc)
template<typename T, MemType MEM>
class Multigrid
{
public:
  Multigrid(
            DisjointBoxLayout a_layout,
            int a_solverType,
            int a_numlevels);
  
  void mgRelax(LevelBoxData<T,1,MEM>& a_phi,
               const LevelBoxData<T,1,MEM>& a_rho,
               T a_h,
               int a_level);
  
  void pointRelax(LevelBoxData<T,1,MEM>& a_phi,
                  const LevelBoxData<T,1,MEM>& a_rho,
                  T a_h,
                  int a_p,
                  int a_level);
  
  void residual(LevelBoxData<T,1,MEM>& a_res,
                LevelBoxData<T,1,MEM>& a_phi, 
                const LevelBoxData<T,1,MEM>& a_rho,
                T a_h);

  T resnorm(LevelBoxData<T,1,MEM>& a_phi, 
            const LevelBoxData<T,1,MEM>& a_rho,
            T a_h);
  
  void average(BoxData<T>& a_rhoCoarse,
               const BoxData<T>& a_rho,
               Box a_bxCoarseValid);
  
  void interpolate(BoxData<T>& a_phi, 
                   const BoxData<T>& a_phiCoarse);
  
  void setGhost(BoxData<T,1,MEM>& a_phi,
                T a_h,
                ProblemDomain domain,
                Box a_bx,
                int a_order = 2,
                int a_ghost = 1
                );
  
  int diamStencil(){return m_ghost;};
  Stencil<T> smoother(){return m_smoother;};
  T diagCoefLaplacian(){return m_diagCoefLaplacian;};
  T diagCoefSmoother(){return m_diagCoefSmooth;};
  
private:
  int m_p;
  int m_pbottom;
  Stencil<T> m_laplacian;
  Stencil<T> m_smoother;
  T m_diagCoefLaplacian;
  T m_diagCoefSmooth;
  int m_ghost;
  int m_numLevels;
  int m_order;
  Array<Stencil<T>,DIM> m_sten42low,m_sten2low,m_sten42high,m_sten2high;
  std::vector<Stencil<T>> m_sten4;
  BoxData<int> m_sten41ij;
  T m_dirichlet2 = 2.0;
  T m_dirichlet42 = 64.0/5.0;
  T m_dirichlet41 = 16.0/5.0;
  //T m_dirichlet41 = 2.0;
  
  std::vector<DisjointBoxLayout> m_layouts;
};
#include "MultigridImplem.H"
#endif
