#include "LevelMultigrid.H"
#include "Proto_Timer.H"
#include "Proto_WriteBoxData.H"
using namespace std;

LevelMultigrid::LevelMultigrid(
                     const Box& a_bx,
                     double a_dx,
                     int a_level
                     )
{
  this->define(a_bx,a_dx,a_level);
};
/// Define function. a_bx is the domain, a_dx the mesh spacing, and 
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
  m_dbl = DisjointBoxLayout(ProblemDomain(a_bx,per),min(MAXBOXSIZE,size)*Point::Ones());
  m_level = a_level;
  m_dx = a_dx;
  m_lambda = m_dx*m_dx/(4*DIM);

  if (m_level > 0)
    {
      // Set up next coarser level.
      // m_resc, m_delta are created with their own DBL set up to have the largest
      // patch size possible, while m_localCoarse is defined on a box-by-box coarsening of m_dbl.

      DisjointBoxLayout blCoarseLocal = m_dbl;
      blCoarseLocal.coarsen(2*Point::Ones());
      DisjointBoxLayout blCoarse(ProblemDomain(a_bx.coarsen(2),per),
                                 min(MAXBOXSIZE,size/2)*Point::Ones());
      m_resc.define(blCoarse,Point::Zeros());      
      m_delta.define(blCoarse,Point::Ones());
      m_localCoarse.define(blCoarseLocal,Point::Zeros());
      // Pointer to next-coarser multigrid.
      m_coarsePtr = 
        shared_ptr<LevelMultigrid >(new LevelMultigrid(a_bx.coarsen(2),2*m_dx,m_level-1));
    }
};
void
LevelMultigrid::coarseResidual(
                    LevelBoxData<double >& a_resc,
                    LevelBoxData<double >& a_phi,
                    LevelBoxData<double >& a_rhs
                    )
{
  PR_TIMERS("residual");
  a_phi.exchange();
  double hsqinv = 1./(m_dx*m_dx);
  for (auto dit=a_phi.begin();*dit != dit.end();++dit)
    {
      BoxData<double>& phi = a_phi[*dit];
      BoxData<double>& rhs = a_rhs[*dit];
      BoxData<double>& rescLocal = m_localCoarse[*dit];
      BoxData<double> res(m_dbl[*dit]);
      res.setVal(0.);
      res += rhs;
      res += Stencil<double>::Laplacian()(phi,-hsqinv);
      rescLocal |= Stencil<double>::AvgDown(2)(res);
    }
  m_localCoarse.copyTo(a_resc);
};
double
LevelMultigrid::resnorm(
                    LevelBoxData<double >& a_phi,
                    LevelBoxData<double >& a_rhs
                    )
{
  PR_TIMERS("resnorm");  
  a_phi.exchange();
  double hsqinv = 1./(m_dx*m_dx);
  double maxnorm = 0.;
  for (auto dit=a_phi.begin();*dit != dit.end();++dit)
    {
      BoxData<double>& phi = a_phi[*dit];
      BoxData<double>& rhs = a_rhs[*dit];
      BoxData<double> res(m_dbl[*dit]);
      res.setVal(0.);
      res -= rhs;
      res += Stencil<double>::Laplacian()(phi,hsqinv);
      maxnorm = max(res.absMax(),maxnorm);
    }
  return maxnorm;
};
void
LevelMultigrid::pointRelax(
                      LevelBoxData<double >& a_phi,
                      LevelBoxData<double >& a_rhs,
                      int a_numIter
                      ) 
{
  PR_TIMERS("relax");
  
  double wgt = 1.0/(4*DIM);

  for (int iter = 0; iter < a_numIter;iter++)
    {
      a_phi.exchange();
      auto diag = (-m_lambda)*Shift(Point::Zeros());
      for (auto dit=a_phi.begin();*dit != dit.end();++dit)
        { 
          BoxData<double>& phi = a_phi[*dit];
          BoxData<double>& rhs = a_rhs[*dit];
          BoxData<double> temp = Stencil<double>::Laplacian()(phi,1./(4.*DIM));
          temp += diag(rhs);
          phi+= temp;
        }
    }
}
void
LevelMultigrid::fineInterp(
                   LevelBoxData<double >& a_phi,
                   LevelBoxData<double >& a_delta
                   )
{
  PR_TIMERS("fineInterp");
  a_delta.copyTo(m_localCoarse);
  
  for (auto dit=a_phi.begin();*dit != dit.end();++dit)
    {
      BoxData<double>& phi = a_phi[*dit];
      BoxData<double>& delta = m_localCoarse[*dit];
      Box K(Point::Zeros(),Point::Ones());
      for (auto itker = K.begin();!itker.done();++itker)
        {
          phi += m_fineInterp(*itker)(delta,m_dbl[*dit].coarsen(2));
        }
    }
};
void 
LevelMultigrid::vCycle(
                  LevelBoxData<double >& a_phi,
                  LevelBoxData<double >& a_rhs
                  )
{
  PR_TIMERS("vcycle");  
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
