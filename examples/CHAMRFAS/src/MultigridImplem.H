#include "AMRUtils.H" 

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

        //m_amrInterp = Proto::InterpStencil<Real>::Quadratic(interpRatio);
        m_amrInterp = Proto::InterpStencil<Real>::Build(
            3,
            Proto::Box::Cube(5).shift(Proto::Point::Ones(-2)),
            4,
            interpRatio);
        DisjointBoxLayout coarseAMRLayout;
        coarsen_dbl(coarseAMRLayout, a_layout, interpRatio);
        m_phiCAMR.define(coarseAMRLayout, OP::numcomps(), m_op.ghost() + Proto::Point::Ones());
    }
    if (a_level > 0)
    {
        DisjointBoxLayout coarseLayout;
        coarsenGrid(coarseLayout, a_layout, MG_REFRATIO);

        m_phiC0.define(coarseLayout, OP::numcomps(), OP::ghost());
        m_phiC.define(coarseLayout, OP::numcomps(), OP::ghost());
        m_RC.define(coarseLayout, OP::numcomps(), Proto::Point::Zeros());

        m_coarser = make_shared<Multigrid<OP,DATA>>(coarseLayout, a_dx*MG_REFRATIO, a_level - 1, a_amr, a_amrLevel);
    }
    m_defined = true;    
}

template<class OP, class DATA>
void Multigrid<OP, DATA>::vcycle(LD& a_phi, const LD& a_R)
{
    MG_LEVEL = m_level;
    CH_TIME("Multigrid::vcycle");
    //writeLevel (a_phi, "MG_Phi.%i.0.hdf5",m_level);
    //writeLevel (a_phi, "MG_R.%i.0.hdf5",m_level);
    if (m_level == 0)
    {
        m_op.relax(a_phi, a_R, BOTTOM_RELAX);
        //writeLevel (a_phi, "MG_Phi.%i.1.hdf5",m_level);

    } else {
        m_op.relax(a_phi, a_R, PRE_RELAX);
        //writeLevel (a_phi, "MG_Phi.%i.1.hdf5",m_level);
        m_op.coarsen(m_phiC, a_phi);
        m_phiC.copyTo(m_phiC0);
        m_op.coarseRhs(m_RC, m_phiC, a_phi, a_R);
        m_coarser->vcycle(m_phiC, m_RC);
        m_op.fineCorrection(a_phi, m_phiC, m_phiC0);
        //writeLevel (a_phi, "MG_Phi.%i.2.hdf5",m_level);
        m_op.relax(a_phi, a_R, POST_RELAX);
        //writeLevel (a_phi, "MG_Phi.%i.3.hdf5",m_level);
    }
}

template<class OP, class DATA>
void Multigrid<OP, DATA>::vcycle(LD& a_phi, LD& a_phiC, const LD& a_R)
{
    a_phiC.copyTo(m_phiCAMR);
    m_op.interpBoundary(a_phi, m_phiCAMR, m_amrInterp);
    //writeLevel(a_phi, "MG_phi.%i.0.hdf5",m_level);
    //writeLevel(a_R, "MG_R.%i.0.hdf5",m_level);
    //std::cout << "Integral of R at level " << m_level << ": "  << integrate(a_R, m_op.dx()) << std::endl;
    //writeLevel(a_phiC, "MG_a_phiC.%i.0.hdf5",m_level);
    LD res;
    res.define(a_phi.disjointBoxLayout(), 1, Proto::Point::Zeros());
    m_op.residual(res, a_phi, a_R);
    //writeLevel(res, "MG_res.%i.0.hdf5",m_level);
    //std::cout << "Integral of res at start of level " << m_level << ": " << integrate(res, m_op.dx()) << std::endl;
    MG_LEVEL = m_level;
    CH_TIME("Multigrid::vcycle");
    if (m_level == 0)
    {
        m_op.relax(a_phi, a_R, BOTTOM_RELAX);
        //writeLevel(a_phi, "MG_phi.%i.1.hdf5", m_level);
        m_op.residual(res, a_phi, a_R);
        //writeLevel(res, "MG_res.%i.1.hdf5",m_level);
        //std::cout << "Integral of res after bottom relax: " << integrate(res, m_op.dx()) << std::endl;
    } else {
        m_op.relax(a_phi, a_R, PRE_RELAX);
        m_op.residual(res, a_phi, a_R);
        //writeLevel(res, "MG_res.%i.1.hdf5",m_level);
        //writeLevel(a_phi, "MG_phi.%i.1.hdf5", m_level);
        //std::cout << "Integral of res after downward relax : " << integrate(res, m_op.dx()) << std::endl;
        m_op.coarsen(m_phiC, a_phi);
        a_phiC.copyTo(m_coarser->m_phiCAMR);
        m_coarser->m_op.interpBoundary(m_phiC, m_coarser->m_phiCAMR, m_coarser->m_amrInterp);
        m_phiC.copyTo(m_phiC0);
        m_op.coarseRhs(m_RC, m_phiC, a_phi, a_R);
        m_coarser->vcycle(m_phiC, a_phiC, m_RC);
        m_op.fineCorrection(a_phi, m_phiC, m_phiC0);
        m_op.residual(res, a_phi, a_R);
        //writeLevel(res, "MG_res.%i.2.hdf5",m_level);
        //writeLevel(a_phi, "MG_phi.%i.2.hdf5", m_level);
        //std::cout << "Integral of res after recursive call + fineCorrect: " << integrate(res, m_op.dx()) << std::endl;
        m_op.relax(a_phi, a_R, POST_RELAX);
        m_op.residual(res, a_phi, a_R);
        //writeLevel(res, "MG_res.%i.3.hdf5",m_level);
        //writeLevel(a_phi, "MG_phi.%i.3.hdf5", m_level);
        //std::cout << "Integral of res after upward relax: " << integrate(res, m_op.dx()) << std::endl;
    }
}
