#ifndef _SGMultigrid_H_
#define _SGMultigrid_H_
#include <cmath>
#include <memory>
#include "Proto.H"

#if DIM==2
#define MG_NUM_COLORS 4
#else
#define MG_NUM_COLORS 8
#endif

using namespace Proto;
using std::array;
/// Multigrid: solve a FV discretization of Poisson's equation on a rectangle.  periodic boundary conditions.
class SGMultigridLevel
{
public:

  ///
  SGMultigridLevel(const double & a_alpha,
                   const double & a_beta,
                   const double & a_dx,
                   const Box    & a_domain);

  ///make coarse from fine 
  SGMultigridLevel(const SGMultigridLevel& a_finerLevel);

  
  ///  
  ~SGMultigridLevel()
  {
  }

  /// 
  void residual(BoxData<double, 1>       & a_res,
                const BoxData<double, 1> & a_phi,
                const BoxData<double, 1> & a_rhs);


  /// 
  void relax(BoxData<double, 1>      & a_phi,
             const BoxData<double, 1>& a_rhs);


  /// average down residual to next coarser level.
  void restrictResidual(BoxData<double, 1>      & a_resc,
                        const BoxData<double, 1>& a_res);

  /// Piecewise constant interpolation of coarse correction to increment fine solution.
  void prolongIncrement(BoxData<double, 1>      & a_phiFine,
                        const BoxData<double, 1>& a_deltaCoarse);


  /// 
  void applyOp(BoxData<double, 1>       & a_lph,
               const BoxData<double, 1> & a_phi);

  /// Multigrid v-cycle.
  void vCycle(BoxData<double, 1>      & a_phi,
              const BoxData<double, 1>& a_rhs);

  void
  enforceBoundaryConditions(BoxData<double, 1>& a_phi, int a_ghost = 1);

  std::shared_ptr<SGMultigridLevel> m_coarser;
private:
  void defineStencils();
  void defineCoarserObjects();
  void getMultiColors();

  /// weak construction introduces unnecessary complications
  SGMultigridLevel()
  {
  }

  bool   m_hasCoarser; //not strictly necessary--could check the pointer
  double m_alpha;     
  double m_beta;       
  double m_dx;        
  Box     m_domain;       
  double m_lambda; //relaxation coeff

  array<Point,MG_NUM_COLORS>  s_colors;
  Stencil<double> m_negoperator; //always need -Lap(phi)
  Stencil<double> m_restrict;

  array<Stencil<double>, MG_NUM_COLORS> m_relaxOpPhi;
  array<Stencil<double>, MG_NUM_COLORS> m_updateOpPhi;
  array<Stencil<double>, MG_NUM_COLORS> m_relaxOpRhs;
  array<Stencil<double>, MG_NUM_COLORS> m_prolong;

  BoxData<double, 1> m_resid;
  BoxData<double, 1> m_residC;
  BoxData<double, 1> m_deltaC; 
};


///class that outsiders actually call
class SGMultigrid
{
public:
  static int  s_numSmoothDown ;
  static int  s_numSmoothUp   ;

  /// 
  SGMultigrid(const double & a_alpha,
              const double & a_beta,
              const double & a_dx,
              const Box     & a_box);

  /// 
  void residual(BoxData<double, 1>       & a_res,
                const BoxData<double, 1> & a_phi,
                const BoxData<double, 1> & a_rhs);

  /// 
  void applyOp(BoxData<double, 1>       & a_lph,
               const BoxData<double, 1> & a_phi);

  /// Multigrid v-cycle.
  void vCycle(BoxData<double, 1>      & a_phi,
              const BoxData<double, 1>& a_rhs);

  /// average down residual to next coarser level.
  void restrictResidual(BoxData<double, 1>      & a_resc,
                        const BoxData<double, 1>& a_res)
  {
    m_finest->m_coarser->restrictResidual(a_resc, a_res);
  }

  /// Piecewise constant interpolation of coarse correction to increment fine solution.
  void prolongIncrement(BoxData<double, 1>      & a_phiFine,
                        const BoxData<double, 1>& a_deltaCoarse)
  {
    m_finest->m_coarser->prolongIncrement(a_phiFine, a_deltaCoarse);
  }

  /// 
  void relax(BoxData<double, 1>      & a_phi,
             const BoxData<double, 1>& a_rhs)
  {
    m_finest->relax(a_phi, a_rhs);
  }
  ///
  void
  enforceBoundaryConditions(BoxData<double, 1>& a_phi)
  {
    m_finest->enforceBoundaryConditions(a_phi);
  }
private:
  /// weak construction introduces unnecessary complications
  SGMultigrid()
  {
  }

  std::shared_ptr<SGMultigridLevel> m_finest;
};
#endif
