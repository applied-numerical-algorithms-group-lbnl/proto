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
        m_RC.define(a_grids[a_level-1], AMR_OP::numcomps(), Point::Zeros());
        m_coarser->define(a_grids, a_dx*AMR_REFRATIO, a_level -1, a_bottomDepth);
        m_reflux.define(a_grids[a_level], a_grids[a_level-1],
                        a_grids[a_level].physDomain(), AMR_REFRATIO, 
                        AMR_OP::numcomps());

        // used for interpolating BCs from coarser level
        coarsen_dbl(m_tempLayout, a_grids[a_level], AMR_REFRATIO);
        m_phiTemp.define(m_tempLayout, AMR_OP::numcomps(), AMR_OP::ghost());
    }
    else
    {
        m_mg.define(a_grids[a_level], a_dx, a_bottomDepth);
    }
    //m_op = m_mg.op();
    m_defined = true;
}

template<class AMR_OP, class DATA>
inline void
AMRFAS<AMR_OP, DATA>::vcycle(std::vector<LevelData<DATA>>& a_phi,
                             const std::vector<LevelData<DATA>>& a_rho,
                             const LevelData<DATA>& a_R)
{
    // alias values
    /*
    LevelData<DATA>& phi = a_phi[m_level];    
    LevelData<DATA>& rho = a_rho[m_level];

    if (m_level == 0)
    {
        m_mg.vcycle(phi,a_R); //bottom solve
    }
    else
    {
        LevelData<DATA>& phiC = a_phi[m_level-1];
        phiC.copyTo(m_phiTemp);    
        LevelData<DATA>& rhoC = a_rho[m_level-1];
        m_mg.vcycle(phi, m_phiTemp, a_R); //TODO, phiC needed for BC
        m_op.coarsen(phiC, phi);
        phiC.copyTo(m_phiC0);
        if (m_level > 1) //at least 2 coarser levels exist
        {
            LevelData<DATA>& phiCC = a_phi[m_level-2];    
            computeRhs(m_RC, phi, phiC, phiCC, rhoC);//TODO
        } 
        else
        {
            computeRhs(m_RC, phi, phiC, rhoC);//TODO
        }

        m_coarser.vcycle(a_phi, a_rho, m_RC);

        //TODO: CHANGE fineCorrection to NOT annihilate phiC!!!
        m_op.fineCorrection(phi, phiC, m_phiC0);
        m_mg.vcycle(phi, phiC, a_R); //TODO, phiC needed for BC

    }
    */
}

template<class AMR_OP, class DATA>
inline void
AMRFAS<AMR_OP, DATA>::computeRhs(LevelData<DATA>& a_RC, 
                                 const LevelData<DATA>& a_phi,
                                 const LevelData<DATA>& a_phiC,
                                 const LevelData<DATA>& a_phiCC,
                                 const LevelData<DATA>& a_rhoC)
{

}
