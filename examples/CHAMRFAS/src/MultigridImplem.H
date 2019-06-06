
template<class OP, class DATA>
void Multigrid<OP, DATA>::define(const DisjointBoxLayout& a_layout,
        double     a_dx,
        int        a_level,
        bool       a_amr,
        int        a_amrLevel)  
{
    m_level = a_level;

    m_op.define(a_layout, a_dx, false);
    if (a_amr && (a_amrLevel > 0))
    { 
        //how many MG levels are there going to be
        unsigned int num_mg = AMR_REFRATIO/MG_REFRATIO;
        //how far to coarsen to get to next AMR lvl
        unsigned int interpRatio = AMR_REFRATIO/ipow(MG_REFRATIO,num_mg - m_level - 1); 

        DisjointBoxLayout coarseAMRLayout;
        coarsen_dbl(coarseAMRLayout, a_layout, interpRatio);
        m_phiCAMR.define(coarseAMRLayout, OP::numcomps(), m_op.ghost() + m_op.interpGhost());
    }
    if (a_level > 0)
    {
        DisjointBoxLayout coarseLayout;
        coarsenGrid(coarseLayout, a_layout, MG_REFRATIO);

        m_phiC0.define(coarseLayout, OP::numcomps(), OP::ghost());
        m_phiC.define(coarseLayout, OP::numcomps(), OP::ghost());
        m_RC.define(coarseLayout, OP::numcomps(), Proto::Point::Zeros());

        m_coarser = make_shared<Multigrid<OP,DATA>>(coarseLayout, a_dx*MG_REFRATIO, a_level, a_amr, a_amrLevel);
    }
    m_defined = true;    
}

template<class OP, class DATA>
void Multigrid<OP, DATA>::vcycle(LD& a_phi, const LD& a_R)
{
    if (m_level == 0)
    {
        m_op.relax(a_phi, a_R, BOTTOM_RELAX);
    } else {
        m_op.relax(a_phi, a_R, PRE_RELAX);
        m_op.coarsen(m_phiC, a_phi);
        m_phiC.copyTo(m_phiC0);
        m_op.coarseRhs(m_RC, m_phiC, a_phi, a_R);
        m_coarser->vcycle(m_phiC, m_RC);
        m_op.fineCorrection(a_phi, m_phiC, m_phiC0);
        m_op.relax(a_phi, a_R, POST_RELAX);
    }
}

template<class OP, class DATA>
void Multigrid<OP, DATA>::vcycle(LD& a_phi, LD& a_phiC, const LD& a_R)
{
    a_phiC.copyTo(m_phiCAMR);
    m_op.interpBoundary(a_phi, m_phiCAMR, m_level);
    if (m_level == 0)
    {
        m_op.relax(a_phi, a_R, BOTTOM_RELAX);
    } else {
        m_op.relax(a_phi, a_R, PRE_RELAX);
        m_op.coarsen(m_phiC, a_phi);
        a_phiC.copyTo(m_coarser->m_phiCAMR);
        m_coarser->m_op.interpBoundary(m_phiC, m_coarser->m_phiCAMR, m_level-1);
        m_phiC.copyTo(m_phiC0);
        m_op.coarseRhs(m_RC, m_phiC, a_phi, a_R);
        
        m_coarser->vcycle(m_phiC, a_phiC, m_RC);
        
        m_op.fineCorrection(a_phi, m_phiC, m_phiC0);
        m_op.relax(a_phi, a_R, POST_RELAX);
    }
}
