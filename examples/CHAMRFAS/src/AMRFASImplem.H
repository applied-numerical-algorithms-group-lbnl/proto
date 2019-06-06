template<class AMR_OP, class DATA>
inline void 
AMRFAS<AMR_OP, DATA>::define(const AMRLayout& a_layout,
        double a_dx,
        int    a_level,
        int    a_bottomDepth)

{
    m_level = a_level;
    m_op.define(a_layout[a_level], a_dx);
    if (a_level > 0)
    {
        m_mg.define(a_layout[a_level], a_dx, AMR_REFRATIO/MG_REFRATIO - 1, true, a_level);
        m_phiC0.define(a_layout[a_level-1], AMR_OP::numcomps(), m_op.ghost());
        m_RC.define(a_layout[a_level-1], AMR_OP::numcomps(), Proto::Point::Zeros());
        m_coarser = make_shared<AMRFAS<AMR_OP, DATA>>(a_layout, a_dx*AMR_REFRATIO, a_level - 1, a_bottomDepth);
        m_reflux.define(a_layout[a_level], a_layout[a_level-1],
                a_layout[a_level].physDomain(), AMR_REFRATIO, 
                AMR_OP::numcomps());
    }
    else
    {
        m_mg.define(a_layout[a_level], a_dx, a_bottomDepth, false);
    }
    m_defined = true;
}

template<class AMR_OP, class DATA>
inline void
AMRFAS<AMR_OP, DATA>::vcycle(AMR_DATA& a_phi, const AMR_DATA& a_rho, AMR_DATA& a_res, int a_n)
{
    this->vcycle(a_phi, a_rho, a_res, a_rho[a_rho.size()-1], a_n);
}

template<class AMR_OP, class DATA>
inline void
AMRFAS<AMR_OP, DATA>::vcycle(AMR_DATA& a_phi, const AMR_DATA& a_rho, AMR_DATA& a_res,
        const LevelData<DATA>& a_R, int a_n)
{
    auto& phi = a_phi[m_level];    
    auto& res = a_res[m_level];
    if (m_level == 0)
    {
        m_mg.vcycle(phi,a_R); //bottom solve
    }
    else
    {
        auto& phiC = a_phi[m_level-1];
        auto& rhoC = a_rho[m_level-1];
        auto& resC = a_res[m_level-1];
         
        m_mg.vcycle(phi, phiC, a_R);
        
        m_op.coarsen(phiC, phi);
        phiC.copyTo(m_phiC0);
      
        if (m_level > 1) //at least 2 coarser levels exist, need to interp BC of phiC
        {
            LevelData<DATA>& phiCC = a_phi[m_level-2];
            m_coarser->m_op.coarsen(phiCC, phiC);
            m_coarser->m_op.interpBoundary(phiC, phiCC);
        } 
        m_op.coarseResidual(resC, rhoC, phiC, a_R, phi, m_reflux);
        resC.copyTo(m_RC);
        m_op.plusApply(m_RC, phiC, m_op.dx()*AMR_REFRATIO);
        
        m_coarser->vcycle(a_phi, a_rho, a_res, m_RC,a_n);
        
        m_op.fineCorrection(phi, phiC, m_phiC0);
        
        m_op.interpBoundary(phi, phiC); 
        
        m_mg.vcycle(phi, phiC, a_R);
        
        m_op.coarseResidual(resC, rhoC, phiC, a_R, phi, m_reflux);
        m_op.residual(res, phi, a_R); //this line only matters at the finest level
    }
}

template<class AMR_OP, class DATA>
inline void
AMRFAS<AMR_OP, DATA>::residual(AMR_DATA& a_res, AMR_DATA& a_phi, const AMR_DATA& a_rho)
/*
AMRFAS<AMR_OP, DATA>::residual(std::vector<std::shared_ptr<LevelData<DATA>>> a_res,
        const std::vector<std::shared_ptr<LevelData<DATA>>> a_phi,
        const std::vector<std::shared_ptr<LevelData<DATA>>> a_rho)
*/
{
    // alias values
    auto& phi = a_phi[m_level];    
    auto& res = a_res[m_level];
    auto& rho = a_rho[m_level];
    if (m_level == 0)
    {
        m_op.residual(res, phi, rho);
    }
    else
    {
        auto& phiC = a_phi[m_level-1];
        auto& rhoC = a_rho[m_level-1];
        auto& resC = a_res[m_level-1];
         
        m_op.coarseResidual(resC, rhoC, phiC, rho, phi, m_reflux);
        m_op.residual(res, phi, rho); //this line only matters at the finest level
    }
}

// Mostly for testing; given phi -> compute RHS
template<class AMR_OP, class DATA>
inline void
AMRFAS<AMR_OP, DATA>::operator()(AMR_DATA& a_rhs, AMR_DATA& a_phi)
/*
AMRFAS<AMR_OP, DATA>::operator()(std::vector<std::shared_ptr<LevelData<DATA>>> a_rhs,
        std::vector<std::shared_ptr<LevelData<DATA>>> a_phi)
*/
{
    LevelData<DATA>& phi = a_phi[m_level];    
    LevelData<DATA>& rhs = a_rhs[m_level];
    if (m_level > 0)
    {   
        LevelData<DATA>& phiC = a_phi[m_level-1];    
        LevelData<DATA>& rhsC = a_rhs[m_level-1];    
        (*m_coarser)(a_rhs, a_phi);
        m_op.apply(rhs, phi, m_op.dx());
        m_op.reflux(rhsC, phiC, phi, m_reflux, -1);
        m_op.coarsen(rhsC, rhs); //to fix garbage from reflux in invalid regions
    } else {
        m_op.apply(rhs, phi, m_op.dx());
    }
}

/*
template<class AMR_OP, class DATA>
template<typename... Args>
inline void
AMRFAS<AMR_OP, DATA>::write(std::vector<std::shared_ptr<LevelData<FArrayBox>>> a_data,
    const char* a_fname, Args... a_params)
{
    char fname[100];
    sprintf(fname, a_fname, a_params...);
    Vector<DisjointBoxLayout> grids(m_level+1);
    Vector<LevelData<FArrayBox>*> data(m_level+1);
    for (int ii = 0; ii < grids.size(); ii++)
    {
        grids[ii] = (*(a_data[ii])).disjointBoxLayout();
        auto& temp = *(a_data[ii]);
        data[ii] = &temp;
    }

    Vector<string> varNames(AMR_OP::numcomps());
    for (int ii = 0; ii < varNames.size(); ii++)
    {
        char name[100];
        sprintf(name, "var_%i", ii);
        std::string s(name);
        varNames[ii] = s;
    }

    Box domain = (*(a_data[0])).disjointBoxLayout().physDomain().domainBox();
    Real dx = m_mg.op().dx()*ipow(AMR_REFRATIO, m_level);
    Real dt = 1.0;
    Real t = 0.0;
    Vector<int> ratios(m_level, AMR_REFRATIO);
    int numLevels = m_level + 1;
    WriteAMRHierarchyHDF5(fname, grids, data, varNames, domain, dx, dt, t, ratios, numLevels);
}
*/
