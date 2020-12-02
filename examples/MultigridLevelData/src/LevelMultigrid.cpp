#include "LevelMultigrid.H"
#include "Proto_Timer.H"
#include "Proto_WriteBoxData.H"
using namespace std;

#ifdef PROTO_BRICK
#include "brick.h"
#include "bricksetup.h"
#include "multiarray.h"

namespace {

#if DIM == 3
#define SCRIPT "laplacian3d.py"
#else
#define SCRIPT "laplacian2d.py"
#endif

  void optimized_laplacian(Proto_brick &src,
                           Proto_brick &dest,
                           int brick_idx,
                           bool initToZero,
                           vector<double> &coefs) {
    brick(SCRIPT, BVEC, (BDIM), (BFOLD), brick_idx);
  }
}
#endif

static auto laplacian = Stencil<double>::Laplacian();

LevelMultigrid::LevelMultigrid(
                     const Box& a_bx,
                     double a_dx,
                     int a_level
                     )
{
  this->define(a_bx,a_dx,a_level);
};
/// Define function. a_bx is the domain, a_dxthe mesh spacing, and 
/// a_level is the number of multigrid levels (a_level = 0 gives you point relaxation.).
void
LevelMultigrid::define(
                  const Box& a_bx,
                  double a_dx,
                  int a_level
                  )
{
  m_box = a_bx;
  array<bool,DIM> per;
  per.fill(true);
  int size = m_box.size(0);
  m_dbl = DisjointBoxLayout(a_bx,min(MAXBOXSIZE,size),per);
  m_level = a_level;
  m_dx = a_dx;
  m_lambda = m_dx*m_dx/(4*DIM);

  if (m_level > 0)
    {
      // Set up next coarser level.
      // m_resc, m_delta are created with their own DBL set up to have the largest
      // patch size possible, while m_localCoarse is defined on a box-by-box coarsening of m_dbl.

      DisjointBoxLayout blCoarseLocal = m_dbl;
      blCoarseLocal.coarsen(2);
      DisjointBoxLayout blCoarse(a_bx.coarsen(2),min(MAXBOXSIZE,size/2),per);
      m_resc.define(blCoarse,Point::Ones(BSIZE));
      m_delta.define(blCoarse,Point::Ones(BSIZE));
      m_localCoarse.define(blCoarseLocal,Point::Zeros());
      // Pointer to next-coarser multigrid.
      m_coarsePtr = 
        shared_ptr<LevelMultigrid >(new LevelMultigrid(a_bx.coarsen(2),2*m_dx,m_level-1));
    }

#ifdef PROTO_BRICK
  laplacian.host_optimized = (void*)&optimized_laplacian;
#endif
};
void
LevelMultigrid::coarseResidual(
                    LevelData<BoxData<double > >& a_resc,
                    LevelData<BoxData<double > >& a_phi,
                    LevelData<BoxData<double > >& a_rhs
                    )
{
  PR_TIMERS("residual");
  a_phi.exchange();
  double hsqinv = 1./(m_dx*m_dx);
  for (int i = 0;i < m_dbl.size();i++)
    {
      BoxData<double>& phi = a_phi[i];
      BoxData<double>& rhs = a_rhs[i];
      BoxData<double>& rescLocal = m_localCoarse[i];
      BoxData<double> res(m_dbl[i]);
      res.setVal(0.);
      res += rhs;
      res += Stencil<double>::Laplacian()(phi,-hsqinv);
      rescLocal |= Stencil<double>::AvgDown(2)(res);
    }
  m_localCoarse.copyTo(a_resc);
};
double
LevelMultigrid::resnorm(
                    LevelData<BoxData<double > >& a_phi,
                    LevelData<BoxData<double > >& a_rhs
                    )
{
  PR_TIMERS("resnorm");  
  a_phi.exchange();
  double hsqinv = 1./(m_dx*m_dx);
  double maxnorm = 0.;

  for (int i = 0;i < m_dbl.size();i++)
    {
      BoxData<double>& phi = a_phi[i];
      BoxData<double>& rhs = a_rhs[i];
      BoxData<double> res(m_dbl[i]);
      res.setVal(0.);
      res -= rhs;
      res += laplacian(phi,hsqinv);
      maxnorm = max(res.absMax(),maxnorm);
    }
  return maxnorm;
};
void
LevelMultigrid::pointRelax(
                      LevelData<BoxData<double > >& a_phi,
                      LevelData<BoxData<double > >& a_rhs,
                      int a_numIter
                      ) 
{
  PR_TIMERS("relax");
  
  double wgt = 1.0/(4*DIM);
  auto phi_update = [this](Var<double>& phi,Var<double>& temp, Var<double>& rhs) {
    phi(0) += temp(0) + (-m_lambda) * rhs(0);
  };
  for (int iter = 0; iter < a_numIter;iter++)
    {
    ///if (iter % BSIZE == 0)
      a_phi.exchange();
      auto diag = (-m_lambda)*Shift(Point::Zeros());
      for (int i = 0;i < m_dbl.size();i++)
        { 
          BoxData<double>& phi = a_phi[i];
          BoxData<double>& rhs = a_rhs[i];
          auto bx = m_dbl[i];
          BoxData<double> temp(bx);

          laplacian.apply(phi, temp, bx, true, 1./(4.*DIM));
          forallInPlace(phi_update,bx,phi,temp,rhs);
        }
    }
}
void
LevelMultigrid::fineInterp(
                   LevelData<BoxData<double > >& a_phi,
                   LevelData<BoxData<double > >& a_delta
                   )
{
  PR_TIMERS("fineInterp");
  a_delta.copyTo(m_localCoarse);
  
  for (int i = 0;i < m_dbl.size();i++)
    {
      BoxData<double>& phi = a_phi[i];
      BoxData<double>& delta = m_localCoarse[i];
      Box K(Point::Zeros(),Point::Ones());
      for (auto itker = K.begin();!itker.done();++itker)
        {
          phi += m_fineInterp(*itker)(delta,m_dbl[i].coarsen(2));
        }
    }
};
void 
LevelMultigrid::vCycle(
                  LevelData<BoxData<double > >& a_phi,
                  LevelData<BoxData<double > >& a_rhs
                  )
{
  PR_TIMERS("vcycle");
  a_rhs.exchange();
  if (m_level > 0) 
    {
      pointRelax(a_phi,a_rhs,m_preRelax);
      coarseResidual(m_resc,a_phi,a_rhs);
      m_delta.setToZero(); 
      m_coarsePtr->vCycle(m_delta,m_resc);
      fineInterp(a_phi,m_delta);
      pointRelax(a_phi,a_rhs,m_postRelax);
    }
  else
    {
      pointRelax(a_phi,a_rhs,m_bottomRelax);
    }
};
