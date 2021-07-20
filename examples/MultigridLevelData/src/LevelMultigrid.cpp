#include "LevelMultigrid.H"
#include "Proto_Timer.H"
#include "Proto_WriteBoxData.H"
using namespace std;

LevelMultigrid::LevelMultigrid(
                     const DisjointBoxLayout& a_dbl,
                     double a_dx,
                     int a_level
                     )
{
  this->define(a_dbl,a_dx,a_level);
};
/// Define function. a_bx is the domain, a_dx the mesh spacing, and
/// a_level is the number of multigrid levels (a_level = 0 gives you point relaxation.).
void
LevelMultigrid::define(
                  const DisjointBoxLayout& a_dbl,
                  double a_dx,
                  int a_level
                  )
{
  m_dbl = a_dbl;

  //cout << "Multigrid dbl:\n" << endl;
  //cout << m_dbl << endl;
  m_level = a_level;
  m_dx = a_dx;
  m_lambda = m_dx*m_dx/(4*DIM);
  Point boxsize = m_dbl.boxSize();
  ProblemDomain pd = m_dbl.problemDomain();
  // cout << "in multigrid::define - ProblemDomain = " << m_dbl.problemDomain() << endl;
  if (m_level > 0)
    {
      // Set up next coarser level.
      // m_resc, m_delta are created with their own DBL set up to have the largest
      // patch size possible, while m_localCoarse is defined on a box-by-box coarsening of m_dbl.

      DisjointBoxLayout dblCoarseLocal = m_dbl.coarsen(2*Point::Ones());
      //cout << "dblCoarseLocal: \n"<< dblCoarseLocal << endl;
      DisjointBoxLayout dblCoarse;
      ProblemDomain pdCoarse = pd.coarsen(2*Point::Ones());
      //cout << "checking pd, pdCoarse, Boxsize:\n" <<
      // pd <<  "\n size = " << pd.size() << "\n" << pdCoarse << "\n size = "<< pdCoarse.size() << "\n boxsize = " << boxsize << endl;
      //cout << "boxes : " << pd.box() << " , " << pdCoarse.box() << endl;
      if (pdCoarse.size()%boxsize == Point::Zeros())
        {
          dblCoarse = DisjointBoxLayout(pdCoarse,boxsize);
        }
      else
        {
          dblCoarse =  DisjointBoxLayout(pdCoarse,pdCoarse.size());
        }
      m_resc.define(dblCoarse,Point::Zeros());
      m_delta.define(dblCoarse,Point::Ones());
      m_localCoarse.define(dblCoarseLocal,Point::Zeros());
      // Pointer to next-coarser multigrid.
      m_coarsePtr =
      shared_ptr<LevelMultigrid >(new LevelMultigrid(dblCoarse,2*m_dx,m_level-1));
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
  //cout << "local coarsened DBL: " << endl;
  //cout << m_localCoarse.getDBL() << endl;
  //cout << "fine DBL: " << endl;
  //cout << a_phi.getDBL() << endl;
  for (auto dit=a_phi.begin();*dit != dit.end();++dit)
    {
      BoxData<double>& phi = a_phi[*dit];
      BoxData<double>& rhs = a_rhs[*dit];
      BoxData<double>& rescLocal = m_localCoarse[*dit];
      BoxData<double> res(dit.box());
      res.setVal(0.);
      res += rhs;
      res += Stencil<double>::Laplacian()(phi,-hsqinv);
      rescLocal |= Stencil<double>::AvgDown(2)(res);
    }
  //cout << "coarseResidual: entering copyTo" << endl;
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
      BoxData<double> res(dit.box());

      res.setVal(0.);
      res -= rhs;
      res += Stencil<double>::Laplacian()(phi,hsqinv);
      //Reduction doesn't seem to work properly for me
      m_rxn.reset();
      res.absMax(m_rxn);
      maxnorm = max(m_rxn.fetch(),maxnorm);
      //maxnorm = max(res.absMax(),maxnorm);
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
      auto id = (1.0)*Shift(Point::Zeros());
      for (auto dit=a_phi.begin();*dit != dit.end();++dit)
        {
          BoxData<double>& phi = a_phi[*dit];
          BoxData<double>& rhs = a_rhs[*dit];
          BoxData<double> temp = Stencil<double>::Laplacian()(phi,wgt);
          temp += diag(rhs);
          phi+= id(temp);
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
          phi += m_fineInterp(*itker)(delta,dit.box().coarsen(2));
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
      //cout << "entering pointrelax" << endl;
      pointRelax(a_phi,a_rhs,m_preRelax);
      //cout << "entering coarseResidual" << endl;
      coarseResidual(m_resc,a_phi,a_rhs);
      m_delta.setToZero();
      //cout << "entering vCycle" << endl;
      m_coarsePtr->vCycle(m_delta,m_resc);
      //cout << "entering fineInterp" << endl;
      fineInterp(a_phi,m_delta);
      //cout << "entering pointRelax" << endl;
      pointRelax(a_phi,a_rhs,m_postRelax);
    }
  else
    {
      //cout << "entering pointRelax - bottom" << endl;
      pointRelax(a_phi,a_rhs,m_bottomRelax);
    }
};
