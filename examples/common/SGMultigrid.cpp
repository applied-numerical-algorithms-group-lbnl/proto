#include "SGMultigrid.H"
#include "Proto.H"
#include "Proto_Timer.H"
#include "CommonTemplates.H"


using namespace std;
using namespace Proto;


int  SGMultigrid::s_numSmoothDown = 2;
int  SGMultigrid::s_numSmoothUp   = 2;

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
            const Box     & a_domain)
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
// Color offsets are grouped into "red"=even number of nonzeros (first 2^(DIM-1)) 
// and "black= odd number of nonzeros (the rest).
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
                 const Box     & a_domain)
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
    Box coardom = m_domain.coarsen(2);
    Box growdom = coardom.grow(1);
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
  m_negoperator = (-m_alpha)*Shift(Point::Zeros()) + (-m_beta/(m_dx*m_dx))*(Stencil<double>::Laplacian());

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
    m_updateOpPhi[icolor] = (1.0)*Shift(Point::Zeros());

    m_relaxOpPhi[icolor].destRatio() = Point::Ones(1);
    m_relaxOpPhi[icolor].srcRatio()  = Point::Ones(2);

    m_updateOpPhi[icolor].destRatio() = Point::Ones(2);
    m_updateOpPhi[icolor].srcRatio()  = Point::Ones(1);
    m_updateOpPhi[icolor].destShift() = s_colors[icolor];

    m_relaxOpRhs[icolor].destRatio() = Point::Ones(1);
    m_relaxOpRhs[icolor].srcRatio()  = Point::Ones(2);

    m_prolong[icolor]  =  (1.0)*Shift(Point::Zeros());
    m_prolong[icolor].destRatio() = Point::Ones(2);
    m_prolong[icolor].destShift() = s_colors[icolor];

    m_restrict += (1.0/numpts)*Shift(s_colors[icolor]);
  }    
  m_restrict.srcRatio() = Point::Ones(2);
//  m_restrict.print();
}
/****/
void
SGMultigridLevel::
enforceBoundaryConditions(BoxData<double, 1>& a_phi,int a_ghost)
{
  PR_TIMERS("EnforceBCs");
  Box ghostDir(Point::Ones()*(-1),Point::Ones());
  for (auto it = ghostDir.begin(); it != ghostDir.end(); ++it)
    {
      Point ghRegion = *it;
      if (ghRegion != Point::Zeros())
        {
          Point ptLow,ptHigh,ptShift;
          for (int dir = 0; dir < DIM; dir++)
            {
              if (ghRegion[dir] == -1)
                {
                  ptLow[dir] = m_domain.low()[dir] - a_ghost;
                  ptHigh[dir] = m_domain.low()[dir] -1;
                  ptShift[dir] = m_domain.size(dir);
                }
              else if (ghRegion[dir] == 1)
                {
                  ptLow[dir] = m_domain.high()[dir] + 1;
                  ptHigh[dir] = m_domain.high()[dir] + a_ghost;
                  ptShift[dir] = -m_domain.size(dir);
                }
              else
                {
                  ptLow[dir] = m_domain.low()[dir];
                  ptHigh[dir] = m_domain.high()[dir];
                  ptShift[dir] = 0;
                }
            }
          Box ghBox(ptLow,ptHigh);
          BoxData<double,1> ghostVals(ghBox);
          a_phi.copyTo(ghostVals,ghBox.shift(ptShift),ptShift*(-1));
          ghostVals.copyTo(a_phi,ghBox);
        }
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
  // GSRB. As implemented here, only correct for second-order 5 / 7 point operators. 
  // To go to higher order, need to use full multicolor algorithm. 
  Box coarDom = m_domain.coarsen(2);
  int irel = 0;
  BoxData<double, 1> phisrc(coarDom);
  // Loop over red, black colors.
  for (int evenOdd=0;evenOdd < 2 ; evenOdd++)
  {
    enforceBoundaryConditions(a_phi);
    // loop over 2^(DIM-1) coarsened domains in each color. 
    for(int icolor = evenOdd*MG_NUM_COLORS/2;
        icolor < evenOdd*MG_NUM_COLORS/2 + MG_NUM_COLORS/2; icolor++)
    {
      PR_TIMERS("GSMC");
      //needs to be coarse domain because of the whole gsrb thing            
      {
        PR_TIMERS("applyInRelax");
        phisrc |= m_relaxOpPhi[icolor](a_phi, coarDom);
      }
      {
        PR_TIMERS("rhsInRelax");
        phisrc += m_relaxOpRhs[icolor](a_rhs, coarDom);
      }
      {
        PR_TIMERS("updatePhi");
        a_phi += m_updateOpPhi[icolor](phisrc,coarDom);
      }
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
  //called by the coarser mg level
  for(int icolor = 0; icolor < MG_NUM_COLORS; icolor++)
  {
    a_phi += m_prolong[icolor](a_delta,m_domain);
  }
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
    residual(m_resid,a_phi,a_rhs);                      
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
