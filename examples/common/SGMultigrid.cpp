#include "SGMultigrid.H"
#include "Proto.H"
#include "Proto_Timer.H"
#include "CommonTemplates.H"


using namespace std;
using namespace Proto;


int  SGMultigrid::s_numSmoothDown = 2;
int  SGMultigrid::s_numSmoothUp   = 2;
bool SGMultigrid::s_usePointJacoby = true;

/****/
void
SGMultigridLevel::
applyOp(BoxData<double, 1>       & a_lph,
        const BoxData<double, 1> & a_phi)
                    
{
  PR_TIME("sgmg::applyop");
  BoxData<double,1>& phi = const_cast<BoxData<double, 1>& >(a_phi);
  enforceBoundaryConditions(phi);
  a_lph |= m_negoperator(phi, m_domain);
  a_lph *= -1.0;
}
/****/
void
SGMultigrid::
applyOp(BoxData<double, 1>       & a_lph,
        const BoxData<double, 1> & a_phi)
{
  return m_finest->applyOp(a_lph, a_phi);
}
/***/
SGMultigrid::
SGMultigrid(const double & a_alpha,
            const double & a_beta,
            const double & a_dx,
            const Bx     & a_domain)
{
  m_finest = std::shared_ptr<SGMultigridLevel>(new SGMultigridLevel(a_alpha, a_beta, a_dx, a_domain));
}
/***/
void
SGMultigrid::
residual(BoxData<double, 1>       & a_res,
         const BoxData<double, 1> & a_phi,
         const BoxData<double, 1> & a_rhs)
{
  PR_TIME("sgmg::resid");
  return m_finest->residual(a_res, a_phi, a_rhs);
}
/***/
void
SGMultigrid::
vCycle(BoxData<double, 1>       & a_phi,
       const BoxData<double, 1> & a_rhs)
{
  PR_TIME("sgmg::vcycle");
  return m_finest->vCycle(a_phi, a_rhs);
}
/***/
void
SGMultigridLevel::
getMultiColors()
{

#if DIM==2
  s_colors[0] = Point::Zeros();//(0,0)
  s_colors[1] = Point::Ones();//(1,1)
  s_colors[2] = Point::Zeros() + Point::Basis(1);//(0,1)
  s_colors[3] = Point::Zeros() + Point::Basis(0);//(1,0)
#elif DIM==3
  s_colors[0] = Point::Zeros();//(0,0,0)
  s_colors[1] = Point::Zeros() + Point::Basis(0) + Point::Basis(1);//(1,1,0)
  s_colors[2] = Point::Zeros() + Point::Basis(1) + Point::Basis(2);//(0,1,1)
  s_colors[3] = Point::Zeros() + Point::Basis(0) + Point::Basis(2);//(1,0,1)
  s_colors[4] = Point::Zeros() + Point::Basis(1);//(0,1,0)
  s_colors[5] = Point::Zeros() + Point::Basis(0);//(1,0,0)
  s_colors[6] = Point::Zeros() + Point::Basis(2);//(0,0,1)
  s_colors[7] = Point::Ones();//(1,1,1)
#else
  compiler_error_this_code_is_only_written_for_dim_2_or_3();
#endif
}
/***/
SGMultigridLevel::
SGMultigridLevel(const double & a_alpha,
                 const double & a_beta,
                 const double & a_dx,
                 const Bx     & a_domain)
{
  m_alpha     = a_alpha;
  m_beta      = a_beta;
  m_dx        = a_dx;
  m_domain    = a_domain;
  m_resid.define(m_domain);

  defineStencils();

  defineCoarserObjects();
}
/***/
void
SGMultigridLevel::
defineCoarserObjects()
{
  PR_TIME("sgmglevel::defineCoarser");
  if(m_domain.coarsenable(4))
  {
    Bx coardom = m_domain.coarsen(2);
    Bx growdom = coardom.grow(1);
    m_residC.define(coardom);    
    m_deltaC.define(growdom);
    m_coarser = std::shared_ptr<SGMultigridLevel>(new SGMultigridLevel(*this));
    m_hasCoarser = true;
  }
  else
  {
    m_hasCoarser= false;
  }
}
/***/
SGMultigridLevel::
SGMultigridLevel(const SGMultigridLevel& a_finerLevel)
{
  PR_TIME("sgmglevel::constructor");
  m_alpha     = a_finerLevel.m_alpha;
  m_beta      = a_finerLevel.m_beta;
  m_dx        = 2*a_finerLevel.m_dx;
  m_domain    = a_finerLevel.m_domain.coarsen(2);
  m_resid.define(m_domain);
  defineStencils();
  defineCoarserObjects();

}
/***/
void
SGMultigridLevel::
defineStencils()
{
  PR_TIME("sgmglevel::definestencils");
  getMultiColors();

  //always need -lapl(phi) so store that.
  m_negoperator = (-m_alpha)*Shift(Point::Zeros()) + (-m_beta/(m_dx*m_dx))*(Stencil<double>::Laplacian(2));

  double safety = 1.0;
  double diag = m_alpha + (m_beta*(-2.*DIM)/(m_dx*m_dx));
  m_lambda = safety/diag;

      
  m_restrict = Stencil<double>();
  double numpts = double(MG_NUM_COLORS);
  for(int icolor = 0; icolor < MG_NUM_COLORS; icolor++)
  {
    m_relaxOpPhi[icolor] = (m_lambda)*m_negoperator;
    m_relaxOpPhi[icolor] *= (1.0)*Shift(s_colors[icolor]);

    m_relaxOpRhs[icolor] = (m_lambda)*Shift(s_colors[icolor]);

    m_relaxOpPhi[icolor].destRatio() = Point::Ones(2);
    m_relaxOpPhi[icolor].srcRatio()  = Point::Ones(2);
    m_relaxOpPhi[icolor].destShift() = s_colors[icolor];

    m_relaxOpRhs[icolor].destRatio() = Point::Ones(2);
    m_relaxOpRhs[icolor].srcRatio()  = Point::Ones(2);
    m_relaxOpRhs[icolor].destShift() = s_colors[icolor];

    m_prolong[icolor]  =  (1.0)*Shift(Point::Zeros());
    m_prolong[icolor].destRatio() = Point::Ones(2);
    m_prolong[icolor].destShift() = s_colors[icolor];

    m_restrict += (1.0/numpts)*Shift(s_colors[icolor]);
  }    
  m_restrict.srcRatio() = Point::Ones(2);
//  m_restrict.print();
}
/***/
//cheerfully stolen from the euler example
void
SGMultigridLevel::
enforceBoundaryConditions(BoxData<double, 1>& a_phi)
{
  PR_TIME("sgmglevel::enforcebc");
  for(int idir = 0; idir < DIM; idir++)
  {
    protocommon::enforceSGBoundaryConditions<double, 1>(a_phi, 1, idir);
  }
}
/****/
void
SGMultigridLevel::
residual(BoxData<double, 1>       & a_res,
         const BoxData<double, 1> & a_phi,
         const BoxData<double, 1> & a_rhs)
                    
{
  PR_TIME("sgmglevel::resid");
  BoxData<double,1>& phi = const_cast<BoxData<double, 1>& >(a_phi);
  enforceBoundaryConditions(phi);
  a_res |= m_negoperator(phi, m_domain);
  a_res += a_rhs;
}
/****/
void
SGMultigridLevel::
relax(BoxData<double, 1>       & a_phi,
      const BoxData<double, 1> & a_rhs)
{
  PR_TIME("sgmglevel::relax");
  if(SGMultigrid::s_usePointJacoby)
  {

    BoxData<double, 1> res(m_domain);
    residual(res, a_phi, a_rhs);
    res *= 0.5*m_lambda;
    a_phi += res;
  }
  else
  {
    Bx coarDom = m_domain.coarsen(2);
    int irel = 0;
    for(int icolor = 0; icolor < MG_NUM_COLORS; icolor++)
    {
      enforceBoundaryConditions(a_phi);
      //needs to be coarse domain because of the whole gsrb thing
      BoxData<double, 1> phisrc(a_phi.box());
      a_phi.copyTo(phisrc);
      a_phi += m_relaxOpPhi[icolor](phisrc, coarDom);
      a_phi += m_relaxOpRhs[icolor](a_rhs, coarDom);

      irel++;
    }
  }
}
/****/
void
SGMultigridLevel::
restrictResidual(BoxData<double, 1>       & a_resc,
                 const BoxData<double, 1> & a_res)
{
  PR_TIME("sgmglevel::restrict");
  //called by the coarser mg level
  a_resc |= m_restrict(a_res, m_domain);
}
/****/
void
SGMultigridLevel::
prolongIncrement(BoxData<double, 1>      & a_phi,
                 const BoxData<double, 1>& a_delta)
{
  PR_TIME("sgmglevel::prolong");
  //std::cout << "just inside prolong color phi = ";
  //a_phi.print();
  //std::cout << "just inside prolong color delta = ";
  //a_delta.print();
  //called by the coarser mg level
  for(int icolor = 0; icolor < MG_NUM_COLORS; icolor++)
  {
    a_phi += m_prolong[icolor](a_delta,m_domain);

//    iprong++;
//    std::cout << "prolong color phi = "  << icolor << std::endl;
//    a_phi.print();
//    std::cout << "prolong color delta = "  << icolor << std::endl;
//    a_phi.print();
  }
//  std::cout << "leaving prolong " << std::endl;
}
/****/
void 
SGMultigridLevel::
vCycle(BoxData<double, 1>         & a_phi,
       const BoxData<double, 1>   & a_rhs)
{

  PR_TIME("sgmglevel::vcycle");
  for(int irelax = 0; irelax < SGMultigrid::s_numSmoothDown; irelax++)
  {
    relax(a_phi,a_rhs); //don't do it
  }

  if (m_hasCoarser)
  {
    residual(m_resid,a_phi,a_rhs);                     // 
    m_coarser->restrictResidual(m_residC,m_resid);
    m_deltaC.setVal(0.);
    m_coarser->vCycle(m_deltaC,m_residC);
    m_coarser->prolongIncrement(a_phi,m_deltaC);
  }

  for(int irelax = 0; irelax < SGMultigrid::s_numSmoothUp; irelax++)
  {
    relax(a_phi,a_rhs);
  }

}
/****/
