#ifndef _LEVELMULTIGRID_H_
#define _LEVELMULTIGRID_H_
#include <cmath>
#include <cstring>
#include <memory>
#include "Proto.H"
#include "Proto_LevelBoxData.H"
#include "Proto_DisjointBoxLayout.H"
#define MAXBOXSIZE 64
using namespace std;
using namespace Proto;

/// Multigrid: solve a FV discretization of Poisson's equation on a rectangle.
class LevelMultigrid
{
public:
  /// Not assuming strong construction, so we have constructor and define.
  LevelMultigrid(
            const Box& a_box,
            double a_dx,
            int a_level
            );
  void define(
         const Box& a_box,
         double a_dx,
         int a_level
         );
  /// Compute coarse residual.
  void coarseResidual(
                      LevelBoxData<double >& a_resc,
                      LevelBoxData<double >& a_phi,
                      LevelBoxData<double >&a_rhs
                );
  /// Compute max norm of the residual.
  double resnorm(
                LevelData<BoxData<double > >& a_phi,
                LevelData<BoxData<double > >& a_rhs
                );
  /// Point relaxation. Currently point Jacobi.
  void pointRelax(
                  LevelBoxData<double >& a_phi,
                  LevelBoxData<double >& a_rhs,
                  int a_numIter
                  );
  /// average down residual to next coarser level. Assumes exact averaging on control volumes.
  void avgDown(
               LevelBoxData<double >& a_resc,
               LevelBoxData<double >& a_res
               );
  /// Piecewise constant interpolation of coarse correction to increment fine solution.
  void fineInterp(
                  LevelBoxData<BoxData<double > >& a_phiFine,
                  LevelBoxData<BoxData<double > >& a_deltaCoarse
                  );
  /// Multigrid v-cycle.
  void vCycle(
              LevelBoxData<BoxData<double > >& a_phi,
              LevelBoxData<BoxData<double > >& a_rhs
              );
private:
  LevelBoxData<double > m_localCoarse;
  LevelBoxData<double > m_resc;
  LevelBoxData<double > m_delta;
  Box m_box;
  DisjointBoxLayout m_dbl;
  int m_domainSize;
  shared_ptr<LevelMultigrid > m_coarsePtr;
  double m_dx;
  double m_lambda;
  int m_level;
  int m_preRelax = 2*DIM;
  int m_postRelax = 2*DIM;
  int m_bottomRelax = 10;
  InterpStencil<double> m_fineInterp = InterpStencil<double>::PiecewiseConstant(Point::Ones()*2); 
};
#endif
