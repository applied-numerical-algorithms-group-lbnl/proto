#pragma once

#include "Proto.H"
#include "LevelData.H"
#include "AMRUtils.H"
#include "DisjointBoxLayout.H"
#include <cmath>

#define NUMCOMPS 1

using namespace Proto;

int fileNum = 0;
typedef BoxData<Real, NUMCOMPS> BD;
typedef Var<Real, NUMCOMPS> V;
typedef BaseFab<Real> FAB;

template<class DATA>
class TestOp
{
public:

     
    inline TestOp(){}
    inline TestOp(const DisjointBoxLayout& a_layout, Real a_dx){define(a_layout, a_dx);};
    inline void define(const DisjointBoxLayout& a_layout, Real a_dx);

    static constexpr int numcomps(){return NUMCOMPS;}

    inline Point ghost() const {return m_L.ghost();}  
    inline Real dx() const {return m_dx;};
    
    // Multigrid stuff
    Real residual(LevelData<DATA>& a_res, LevelData<DATA>& a_phi, const LevelData<DATA>& a_R) const;
    void relax(LevelData<DATA>& a_phi, LevelData<DATA>& a_res, const LevelData<DATA>& a_R, int a_n) const;
    void coarsen(LevelData<DATA>& a_phiC, LevelData<DATA>& a_phi);
    void coarseRhs(LevelData<DATA>& a_RC, LevelData<DATA>& a_phiC, const LevelData<DATA>& a_res);
    void fineCorrection(LevelData<DATA>& a_phi, LevelData<DATA>& a_phiC, const LevelData<DATA>& a_phiC0);
    
    // Initialization
    inline static void initialCondition(LevelData<DATA>& a_phi, Real dx);
    inline static void forcing(LevelData<DATA>& a_R, Real dx);
    inline static void solution(LevelData<DATA>& a_soln, Real dx);

private:
    Stencil<Real> m_L;              //actual operator
    Stencil<Real> m_avg;            //averaging stencil
    InterpStencil<Real> m_interp;   //interpolation stencil
    Real m_dx;                      
    Real m_lambda;

    LevelData<DATA> m_coarseTemp;   //temp for coarsening / refining in multigrid

};

template<class DATA>
void TestOp<DATA>::define(const DisjointBoxLayout& a_layout, Real a_dx)
{
    #if DIM == 2
    m_L = Stencil<Real>::Laplacian_9();
    m_lambda = -0.5*(a_dx*a_dx)*(6.0/20.0);
    #elif DIM == 3
    m_L = Stencil<Real>::Laplacian_27();
    m_lambda = -0.5*(a_dx*a_dx)*(15.0/64.0);
    #endif          //undefined otherwise for the time being
    m_dx = a_dx;
    m_avg = Stencil<Real>::AvgDown(MG_REFRATIO);
    m_interp = InterpStencil<Real>::PiecewiseConstant(Point::Ones(MG_REFRATIO));

    DisjointBoxLayout coarseLayout;
    coarsen_dbl(coarseLayout, a_layout, MG_REFRATIO);
    m_coarseTemp.define(coarseLayout, NUMCOMPS);
}

template<class DATA>
Real TestOp<DATA>::residual(LevelData<DATA>& a_res,
                            LevelData<DATA>& a_phi,
                            const LevelData<DATA>& a_R) const
{
  CH_TIME("TestOp::residual");
  a_phi.exchange();
  auto iter = a_phi.getBoxes().dataIterator();
  Real resnorm = 0.0;
  for (iter.begin(); iter.ok(); ++iter)
  {
    BD phi = a_phi[iter()];
    BD res = a_res[iter()];
    BD R = a_R[iter()];
    //Bx b = a_phi.getBoxes()[iter()];
    res |= m_L(phi,1.0/(m_dx*m_dx));
    // use forall to fuse subtraction and reduction ops
    forallInPlace(
    [](V& a_lphi, const V& a_R)//, Param<Real>& a_resnorm)
    {
      for (int comp = 0; comp < NUMCOMPS; comp++)
      {
        a_lphi(comp) = a_R(comp) - a_lphi(comp);
        //a_resnorm() = max(abs(a_lphi(comp)), a_resnorm());
      }
    }, res, R);//,resnorm);
  }
  //return resnorm;
  return 0;
}

template<class DATA>
void TestOp<DATA>::relax(LevelData<DATA>& a_phi,
                         LevelData<DATA>& a_res,
                         const LevelData<DATA>& a_R, int a_n) const
{
    //TODO: make sure LD inputs all have same DBL
  CH_TIME("TestOp::relax");
  auto iter = a_phi.getBoxes().dataIterator();
  for (int ii = 0; ii < a_n; ii++)
    {
      a_phi.exchange();
      for (iter.begin(); iter.ok(); ++iter)
      {
        BD phi = a_phi[iter()];
        BD res = a_res[iter()];
        BD R = a_R[iter()];
        res |= m_L(phi, 1.0/(m_dx*m_dx));
        double lambda = m_lambda;
        forallInPlace(
        [lambda](V& a_phi, V& a_lphi, const V& a_R)// Param<Real>& lambda)
        { 
          for (int comp = 0; comp < NUMCOMPS; comp++)
          {
            a_lphi(comp) = a_R(comp) - a_lphi(comp);
            a_phi(comp) += (a_lphi(comp)*lambda);
          }
        }, phi, res, R);//, m_lambda);
      }
    }
}

template<class DATA>
void TestOp<DATA>::coarsen(LevelData<DATA>& a_phiC,
                           LevelData<DATA>& a_phi)
{
    //TODO: make sure LD inputs all have same DBL
  CH_TIME("TestOp::coarsen");
  a_phi.exchange();
  auto iter = a_phi.getBoxes().dataIterator();
  for (iter.begin(); iter.ok(); ++iter)
    {
      BD phi = a_phi[iter()];
      BD temp = m_coarseTemp[iter()];
      //Bx b = m_coarseTemp.getBoxes()[iter()];
      temp |= m_avg(phi);
    }
  m_coarseTemp.copyTo(a_phiC);
}


template<class DATA>
void TestOp<DATA>::coarseRhs(LevelData<DATA>& a_RC,
                             LevelData<DATA>& a_phiC,
                             const LevelData<DATA>& a_res)
{
  CH_TIME("TestOp::coarseRhs");
  auto fiter = a_res.getBoxes().dataIterator();
  for (fiter.begin(); fiter.ok(); ++fiter)
  {
    BD res = a_res[fiter()];
    BD temp = m_coarseTemp[fiter()];
    //Bx b = m_coarseTemp.getBoxes()[fiter()];
    temp |= m_avg(res);
  }
  m_coarseTemp.copyTo(a_RC);  
  a_phiC.exchange(); 
  auto citer = a_RC.getBoxes().dataIterator();
  
  for (citer.begin(); citer.ok(); ++citer)
  {
    BD R = a_RC[citer()];
    BD phiC = a_phiC[citer()];
    //Bx b = a_RC.getBoxes()[citer()];
    double dx = m_dx*MG_REFRATIO;
    R += m_L(phiC,1.0/(dx*dx));
  }
}

template<class DATA>
void TestOp<DATA>::fineCorrection(LevelData<DATA>& a_phi,
                                  LevelData<DATA>& a_phiC,
                                  const LevelData<DATA>& a_phiC0)
{
  CH_TIME("TestOp::fineCorrection");
  auto citer = a_phiC.dataIterator();
  for (citer.begin(); citer.ok(); ++citer)
  {
    BD phiC = a_phiC[citer()];
    BD phiC0 = a_phiC0[citer()];
    phiC -= phiC0; //phiC is now a correction
  }
  a_phiC.copyTo(m_coarseTemp);
  
  auto fiter = a_phi.getBoxes().dataIterator();
  for (fiter.begin(); fiter.ok(); ++fiter)
  {
    BD phi = a_phi[fiter()];
    BD temp = m_coarseTemp[fiter()];
    //Bx b = m_coarseTemp.getBoxes()[fiter()];
    phi += m_interp(temp);
  }
}
template<class DATA>
void TestOp<DATA>::initialCondition(LevelData<DATA>& a_phi, Real dx)
{
    auto iter = a_phi.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        BD phi = a_phi[iter()];
        forallInPlace_p([dx](Point p, V& phi)
        {
            double x = p[0]*dx + 0.5*dx + .25*M_PI;
            phi(0) = -cos(x) + 0.1*cos(x);
            phi(0) = 0.;
        },phi);
    }
}
template<class DATA>
void TestOp<DATA>::forcing(LevelData<DATA>& a_R, Real dx)
{
    auto iter = a_R.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        BD R = a_R[iter()];
        forallInPlace_p([dx](Point p, V& R)
        {
            double x = p[0]*dx + 0.5*dx;
            R(0) = cos(x);
        },R);
    }
}
template<class DATA>
void TestOp<DATA>::solution(LevelData<DATA>& a_soln, Real dx)
{
    auto iter = a_soln.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        BD soln = a_soln[iter()];
        forallInPlace_p([dx](Point p, V& soln)
        {
            double x = p[0]*dx + 0.5*dx;
            soln(0) = -cos(x);
        },soln);
    }
}

