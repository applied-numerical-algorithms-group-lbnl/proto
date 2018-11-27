template<class AMR_OP, class DATA>
inline void 
AMRFAS<AMR_OP, DATA>::define(const std::vector<DisjointBoxLayout>& a_grids,
        double a_dx,
        int    a_level,
        int    a_bottomDepth)
{
    m_level = a_level;
    if (a_level > 0)
    {
        m_mg.define(a_grids[a_level], a_dx, AMR_REFRATIO/MG_REFRATIO - 1);
        m_phiC0.define(a_grids[a_level-1], AMR_OP::numcomps(), AMR_OP::ghost());
        m_RC.define(a_grids[a_level-1], AMR_OP::numcomps(), Proto::Point::Zeros());
        m_coarser->define(a_grids, a_dx*AMR_REFRATIO, a_level -1, a_bottomDepth);
        m_reflux.define(a_grids[a_level], a_grids[a_level-1],
                a_grids[a_level].physDomain(), AMR_REFRATIO, 
                AMR_OP::numcomps());
    }
    else
    {
        m_mg.define(a_grids[a_level], a_dx, a_bottomDepth);
    }
    m_defined = true;
}

template<class AMR_OP, class DATA>
inline void
AMRFAS<AMR_OP, DATA>::vcycle(std::vector<LevelData<DATA>>& a_phi,
        const std::vector<LevelData<DATA>>& a_rho,
        const LevelData<DATA>& a_R)
{
    // alias values
    LevelData<DATA>& phi = a_phi[m_level];    
    LevelData<DATA>& rho = a_rho[m_level];

    if (m_level == 0)
    {
        // no coarser data exists; no need to interp BCs
        m_mg.vcycle(phi,a_R); //bottom solve
    }
    else
    {
        LevelData<DATA>& phiC = a_phi[m_level-1];
        LevelData<DATA>& rhoC = a_rho[m_level-1];
        m_mg.vcycle(phi, phiC, a_R);
        m_mg.op().coarsen(phiC, phi);
        if (m_level > 1) //at least 2 coarser levels exist
        {
            LevelData<DATA>& phiCC = a_phi[m_level-2];    
            computeRhs(m_RC, phiC, phi, a_R, rhoC, phiCC);//TODO
        } 
        else
        {
            computeRhs(m_RC, phiC, phi, a_R, rhoC);//TODO
        }
        phiC.copyTo(m_phiC0);

        m_coarser.vcycle(a_phi, a_rho, m_RC);

        m_mg.op().fineCorrection(phi, phiC, m_phiC0);
        m_mg.vcycle(phi, phiC, a_R); //TODO, phiC needed for BC

    }
}

template<class AMR_OP, class DATA>
inline void
AMRFAS<AMR_OP, DATA>::computeRhs(LevelData<DATA>& a_RC, 
        LevelData<DATA>& a_phiC,
        const LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_R,
        const LevelData<DATA>& a_rhoC,
        const LevelData<DATA>& a_phiCC)
{
    m_mg.op().interpBoundary(a_phiC, a_phiCC);
    m_mg.op().coarseRhs(a_RC, a_phiC, a_phi, a_R);


}

template<class AMR_OP, class DATA>
inline void
AMRFAS<AMR_OP, DATA>::computeRhs(LevelData<DATA>& a_RC, 
        LevelData<DATA>& a_phiC,
        const LevelData<DATA>& a_phi,
        const LevelData<DATA>& a_R,
        const LevelData<DATA>& a_rhoC)
{
}










