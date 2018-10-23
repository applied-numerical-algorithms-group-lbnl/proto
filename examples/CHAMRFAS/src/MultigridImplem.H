#include "AMRUtils.H" 

template<class OP, class DATA>
void Multigrid<OP, DATA>::define(const DisjointBoxLayout& a_layout,
                                 double     a_dx,
                                 int        a_level)  
{
    m_level = a_level;

    m_op.define(a_layout, a_dx);

    m_res.define(a_layout, OP::numcomps()); 
    cout << "dx, level = " << a_dx << " , " << m_level << endl;
    if (a_level > 0)
    {
        DisjointBoxLayout coarseLayout;
        // coarsened, possibly with a different patch layout
        coarsenGrid(coarseLayout, a_layout);
        
        m_phiC0.define(coarseLayout, OP::numcomps(), m_op.ghost());
        m_phiC.define(coarseLayout, OP::numcomps(), m_op.ghost());
        m_RC.define(coarseLayout, OP::numcomps());
        
        m_coarser = make_shared<Multigrid<OP,DATA>>(coarseLayout, a_dx*MG_REFRATIO, a_level - 1);
    }
    m_defined = true;    
}

template<class OP, class DATA>
void Multigrid<OP, DATA>::vcycle(LD& a_phi, const LD& a_R)
{
    //double resnorm = 0.0;
    char fname[100];
    CH_TIME("Multigrid::vcycle");
    if (m_level == 0)
    {
      m_op.relax(a_phi, m_res, a_R, BOTTOM_RELAX);

    } else {
      m_op.relax(a_phi, m_res, a_R, PRE_RELAX);
      m_op.residual(m_res,a_phi,a_R);
      m_op.coarsen(m_phiC, a_phi);
      m_phiC.copyTo(m_phiC0);
      m_op.coarseRhs(m_RC, m_phiC0, m_res);
      m_coarser->vcycle(m_phiC, m_RC);
      m_op.fineCorrection(a_phi, m_phiC, m_phiC0);
      m_op.relax(a_phi, m_res, a_R, POST_RELAX);
    }

    // return resnorm;
}

template<class OP, class DATA>
void Multigrid<OP, DATA>::vcycle(LD& a_phi, const LD& a_phiC, const LD& a_R)
{
  m_op.interpBoundary(a_phi, a_phiC);
  vcycle(a_phi, a_R);
}
